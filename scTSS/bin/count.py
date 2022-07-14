from optparse import OptionParser,OptionGroup
from ..version import __version__
import sys
from ..utils.build_ref import get_TSSref,get_generef
from ..utils.get_counts import get_TSS_count
import pyranges as pr
import os
import pandas as pd


def main():
    parser = OptionParser()
    parser.add_option('--gtf','-g',dest='gtf_file',default=None,help='The annotation gtf file for your analysing species.')
    parser.add_option('--bam','-b',dest='bam_file',default=None,help='The bam file of aligned from Cellranger or other single cell aligned software.')
    parser.add_option('--outdir','-o',dest='out_dir',default=None,help='The directory for output [default : $bam_file]') #what should be after $
   
   
    group0=OptionGroup(parser,"Optional arguments")
    group0.add_option("--minCount",type="int",dest="minCount",default=30,
    help="Minimum counts for each transcript in all cells [default: 30]")
    # group0.add_option("--isoformNumber",type="int",dest="isoformNumber",default=2,
    # help="No. of isoform keeping in for each gene [default: 2]")
    group0.add_option('--nproc','-p',dest='nproc',default=4,
    help='Number of subprocesses [default: 4]')
    group0.add_option('--maxReadCount',dest='maxReadCount',default=50000,
    help='For each gene, the maxmium read count kept for clustering[default: 50000]')


    parser.add_option_group(group0)
    (options, args) = parser.parse_args()


    
    #this means that if users do not input any argument, then direct produce help. then end.
    if len(sys.argv[1:]) == 0:
        print('Welcome to scTSS-count v%s!\n'%(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)


    #output file 
    if options.out_dir is None:
        print("Warning: no outDir provided, we use $bamfilePath/scTSS\n")
        out_dir = os.path.dirname(os.path.abspath(options.bam_file)) + "/scTSS"
    elif os.path.dirname(options.out_dir) == "":
        out_dir= "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    #print(out_dir)
        
        
    #gtf file
    if options.gtf_file is None:
        print("Error: Need --gtf for annotation file.")
        sys.exit(1)
    else:
        gr = pr.read_gtf(options.gtf_file)
        grdf = gr.df
        ref_out_dir=str(out_dir)+'/ref_file/'
        if not os.path.exists(ref_out_dir):
            os.mkdir(ref_out_dir)
        tssrefpath=get_TSSref(grdf,ref_out_dir)
        tssdf=pd.read_csv(tssrefpath,delimiter='\t')
        generefpath=get_generef(grdf,tssdf,ref_out_dir)

    bam_file=options.bam_file
    minCount=options.minCount
    isoformNumber=options.isoformNumber



           
    #bam file
    if options.bam_file is None:
        print("Error: Need --bam for aligned file.")
        sys.exit(1)
    else:
        getTSScount=get_TSS_count(generefpath,tssrefpath,options.bam_file,out_dir,n_proc,minCount,maxReadCount)
        scadata=getTSScount.produce_sclevel()

        
    