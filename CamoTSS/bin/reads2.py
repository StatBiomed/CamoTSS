import numpy as np
import pandas as pd
from brie.utils import fetch_reads
import os
import anndata as ad
import multiprocessing 
import pysam
import time
import warnings
import sys
import matplotlib.pyplot as plt


global CACHE_CHROM
global CACHE_SAMFILE
CACHE_CHROM = None
CACHE_SAMFILE = None



def check_pysam_chrom(samFile,chrom=None):
    """
    check samFile is a file name or a pysam object, and if chrom format.
    """
    global CACHE_CHROM
    global CACHE_SAMFILE

    if CACHE_CHROM is not None:
        if (samFile==CACHE_SAMFILE) and (chrom==CACHE_CHROM):
            return CACHE_SAMFILE,CACHE_CHROM

    
    if type(samFile)==str or type(samFile)==np.str_:
        ftype=samFile.split(".")[-1]
        if ftype != "bam" and ftype!="sam" and ftype!="cram":
            print("Error: file type need suffix of bam, sam or cram")
            sys.exit(1)
        if ftype == "cram":
            samFile=pysam.AlignmentFile(samFile,'rc')
        elif ftype=="bam":
            samFile=pysam.AlignmentFile(samFile,'rb')
        else:
            samFile=pysam.AlignmentFile(samFile,'r')

    if chrom is not None:
        if chrom not in samFile.references:
            if chrom.startswith("chr"):
                chrom=chrom.split('chr')[1]
            else:
                chrom="chr"+chrom
        if chrom not in samFile.references:
            print("Can't find reference %s in samFile"%chrom)
            return samFile,None

    CACHE_CHROM=chrom
    CACHE_SAMFILE=samFile
    return samFile,chrom



def _getreads_2(bamfilePath,geneid,generefdf):
    #fetch reads2 in gene 
    samFile, _chrom = check_pysam_chrom(bamfilePath,'chr'+str(generefdf.loc[geneid]['Chromosome']))
    reads = fetch_reads(samFile, _chrom,  generefdf.loc[geneid]['Start'] , generefdf.loc[geneid]['End'],  trimLen_max=100)
    reads2_umi = reads["reads1u"]
#     print(geneid)

    #select according to GX tag and CB (filter according to user owned cell)
    reads2_umi=[r for r in reads2_umi if r.get_tag('GX')==geneid]
    reads2_umi=[r for r in reads2_umi if r.get_tag('CB') in self.cellBarcode]
    
#     print(reads2_umi)
    reads2_position = [r.positions[0] for r in reads2_umi]
    
    return reads2_position

    