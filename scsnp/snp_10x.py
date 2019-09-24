import collections
import toolz
import numpy as np
import pysam
"""
SNPs in single cell 10x data
"""

def parse_chromium_bamread_metadata(alignment:pysam.AlignedSegment):
    """
    return the readname, error-corrected cellbarvode and error corrected UMI"
    according to:
    'https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam'
    """
    cellbarcode = alignment.get_tag('CB')
    umi = alignment.get_tag('UB')
    readname = alignment.query_name
    return readname, cellbarcode, umi

def get_bases_at_genomic_position(chrom:str, start:int, bamfile:pysam.AlignmentFile):
    """
    for each read covering the given position,
    return readname, cell-barcode, umi and base
    """
    n_reads = 0
    reads_covering_position = []
    for alignment in bamfile.fetch(chrom, start,start+1):
        n_reads+=1
        if start in alignment.get_reference_positions():
            if alignment.has_tag('CB') and alignment.has_tag('UB'):
                readname, cellbarcode, umi = parse_chromium_bamread_metadata(alignment)

                base_index = alignment.get_reference_positions().index(start) # this is the position of the base in the query sequence (due to .index())
                base = alignment.query_alignment_sequence[base_index]
#                 reference_base = alignment.reference_alignment_sequence[base_index]
    #             print(base)
                reads_covering_position.append((readname, cellbarcode,umi, base))
    return reads_covering_position


def quantify_SNP_in_position(chrom:str, pos:int, bamfile:pysam.AlignmentFile):
    """
    returns a dict of [cells] -> list[bases]
    the list has one element per molecule the covers this base in the cell
    i.e. if there's 5RNAs with a consistent 'A'nucleotide in cell `x`  cell_dict[x] == ['A', 'A', 'A', 'A', 'A']
    """
    reads_covering_snp = get_bases_at_genomic_position(chrom=chrom, start=pos, bamfile=bamfile)

    # the above are just the reads covering the base of interest, now aggregate UMI/Cells
    q = collections.defaultdict(list)
    for readname, cb, umi, base in reads_covering_snp:
        q[(cb,umi)].append(base)
    # up to here, each CB/UMI might have multiple reads, contradicting bases
#     toolz.dicttoolz.valfilter(lambda x: len(set(x))>1, q)

    # now collapse the reads from each molecule onto the most frequent base
    q = toolz.dicttoolz.valmap(lambda baselist: collections.Counter(baselist).most_common(1)[0][0], q)

    # now q contains the base for every cell/ every molecule
    # we still have to aggregate the molecules per cell: if its not a PCR error, all molecules in the cell should have the same base
    cell_dict = collections.defaultdict(list)
    for (cb, umi), base in q.items():
        cell_dict[cb].append(base)

    return cell_dict


def annotate_snp_2_adata(adata, chrom:str, pos:int, bam:pysam.AlignmentFile, refbase:str):
    """
    adds the SNP-status to the .obs dataframe of each cell

    adata: sc.AnnData containing the 10x data
    chrom: chromosome of the SNP
    pos: pos of the SNP
    bam: 10x bamfile
    refbase: the reference base to expect

    adds three fields to the .obs_names
    - _GT : the genotype:
        - '.' if no reads
        - '0/0' if only reference bases
        - '0/1' if at least one alt base
    - _nALT: number of molecules supporting the alt allele in the cell
    - _coverage: number of molecules covering the locus in the cell
    """
    GT_00 = 'GT 0/0'
    GT_01 = 'GT 0/1'
    GT_unknown = '.'
    the_dict = quantify_SNP_in_position(chrom, pos, bam)

    snp_field = f'SNP_{chrom}_{pos}_GT'
    snp_field_ALT = f'SNP_{chrom}_{pos}_nALT'
    snp_field_coverage = f'SNP_{chrom}_{pos}_coverage'

    adata.obs[snp_field] = GT_unknown
    adata.obs[snp_field_ALT] = 0
    adata.obs[snp_field_coverage] = 0


    unknown_cells = [] # keep track of cells where we called a mutation but its not part of the adata
    # iterate over all cells that have a read at the locus
    for CB, base_list in the_dict.items():
        if CB not in adata.obs_names:
            unknown_cells.append(CB)
            continue
        genotype = GT_00 if all([b == refbase for b in base_list]) else GT_01

        n_alt_alleles = np.sum([b != refbase for b in base_list])

        adata.obs.loc[CB, snp_field] = genotype
        adata.obs.loc[CB, snp_field_ALT] = n_alt_alleles
        adata.obs.loc[CB, snp_field_coverage] = len(base_list)
    
    if unknown_cells:
        print('Cells not in the adata:')
        print(' '.join(unknown_cells))

"""
do to speed considerations, we cant really go base by base (for each base, we
would again query the alginment file, loading the same reads over and over).
instead, we look at bigger chunks of the genome (e.g. exons), load all reads
covering that region and then (with all reads in mem), walk along the region
base by base
"""

class _EfficientRead():
    """
    just a helper to make the access of reads very fast
    """
    def __init__(self, alignment):
        readname, cellbarcode, umi = parse_chromium_bamread_metadata(alignment)

        self.readname = readname
        self.cellbarcode = cellbarcode
        self.umi = umi

        self.reference_positions = alignment.get_reference_positions()
        self.query_alignment_sequence = alignment.query_alignment_sequence
        self.reference_positions_set = set(self.reference_positions) # set for faster lookup
        self.start = alignment.reference_start
        self.end = alignment.reference_end


    def covers_base(self, i):
        if self.start <= i <= self.end and i in self.reference_positions_set:
            return True
        else:
            return False

def quantify_SNP_in_region(chrom, start, end, bamfile, MAPQ_threshold=255):
    """
    instead of going base-by-base as in `quantify_SNP_in_position()`
    we quantify all SNPs in a region. This should be much faster then base-by-base as we
    only have to load the reads once into memory

    MAPQ_threshold: only consider reads with mapping quality >= threshold;
           default=255 means uniquely mapped reads only (STAR alginer)
    """

    # only take reads having good quality and 10x metadata!
    reads = [alignment for alignment in bamfile.fetch(chrom, start, end) if alignment.mapping_quality>=MAPQ_threshold and alignment.has_tag('CB') and alignment.has_tag('UB')]

    # we have to do some more preprocessing to make the down there efficient
    # 1. get CB/UMI
    # 2. precalc reference_positions
    # 3. the reads base at a specific genomic position
    eff_reads = [_EfficientRead(r) for r in reads]

#     print(len(reads))
    for i in range(start, end):
        reads_covering_position = []

        """
        the old way:
        """
        """
        for r in reads:
            "i in r.get_reference_positions() is expensive. instead check if start,stop overlap with i, then check if the read actually covers!"
            if r.reference_start <= i <= r.reference_end and i in r.get_reference_positions():
#             if i in r.get_reference_positions():
                if r.has_tag('CB') and r.has_tag('UB'):
                    readname, cellbarcode, umi = parse_chromium_bamread_metadata(r)

                    base_index = r.get_reference_positions().index(i) # this is the position of the base in the query sequence (due to .index())
                    base = r.query_alignment_sequence[base_index]
        #                 reference_base = alignment.reference_alignment_sequence[base_index]
                    reads_covering_position.append((readname, cellbarcode,umi, base))
        """

        """
        faster with preprocess
        """
        for r in eff_reads:
            if r.covers_base(i):
                base_index = r.reference_positions.index(i) # this is the position of the base in the query sequence (due to .index())
                base = r.query_alignment_sequence[base_index]
                reads_covering_position.append((r.readname, r.cellbarcode, r.umi, base))


        # the above are just the reads covering the base of interest, now aggregate UMI/Cells
        q = collections.defaultdict(list)
        for readname, cb, umi, base in reads_covering_position:
            q[(cb,umi)].append(base)
        # up to here, each CB/UMI might have multiple reads, contradicting bases
        #     toolz.dicttoolz.valfilter(lambda x: len(set(x))>1, q)

        # now collapse the reads from each molecule onto the most frequent base
        q = toolz.dicttoolz.valmap(lambda baselist: collections.Counter(baselist).most_common(1)[0][0], q)

        # now q contains the base for every cell/ every molecule
        # we still have to aggregate the molecules per cell: if its not a PCR error, all molecules in the cell should have the same base
        cell_dict = collections.defaultdict(list)
        for (cb, umi), base in q.items():
            cell_dict[cb].append(base)
        yield i, cell_dict
