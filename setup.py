from setuptools import setup
import scsnp

setup(name='scsnp',
      version=0.1,
      description='Functions to handle time-lapse microscopy data',
      url='http://github.com/redst4r/scsnp/',
      author='redst4r',
      maintainer='redst4r',
      maintainer_email='redst4r@web.de',
      license='GNU GPL 3',
      keywords='scRNAseq SNP',
      packages=['scsnp'],
      install_requires=[
        'toolz',
        'numpy',
        'pysam'],
      zip_safe=False)
