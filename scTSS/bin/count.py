from optparse import OptionParser,OptionGroup
from ..version import __version__
import sys
from ..utils.build_ref import get_TSSref,get_generef,get_filter_TSS
from ..utils.get_counts import get_TSS_count
from ..utils.get_count_no_novel import get_old_TSS_count
import pyranges as pr
import os
import pandas as pd
import time 


START_TIME = time.time()


def main():
    parser = OptionParser()
    parser.add_option('--gtf','-g',dest='gtf_file',default=None,help='The annotation gtf file for your analysing species.')
    parser.add_option('--cdrFile','-c',dest='cdrFile',default=None,help='The file include cell barcode which users want to keep in the downstream analysis.Actually, it can be the same file input in the scTSS-quant')
    parser.add_option('--bam','-b',dest='bam_file',default=None,help='The bam file of aligned from Cellranger or other single cell aligned software.')
    parser.add_option('--outdir','-o',dest='out_dir',default=None,help='The directory for output [default : $bam_file]') #what should be after $
    parser.add_option('--refFastq','-r',dest='refFastq',default=None,help='The directory for reference fastq file') #what should be after $
    parser.add_option('--mode','-m',dest='mode',default=None,help='You can select run by finding novel TSS mode [New] or just based on gtf annotation file [Old]')

   
   
    group0=OptionGroup(parser,"Optional arguments")

    group0.add_option("--minCount",type="int",dest="minCount",default=30,
    help="Minimum counts for each transcript in all cells [default: 30]")
    # group0.add_option("--isoformNumber",type="int",dest="isoformNumber",default=2,
    # help="No. of isoform keeping in for each gene [default: 2]")
    group0.add_option('--nproc','-p',type="int",dest='nproc',default=4,
    help='Number of subprocesses [default: 4]')

    group0.add_option('--maxReadCount',type="int",dest='maxReadCount',default=10000,
    help='For each gene, the maxmium read count kept for clustering[default: 10000]')
    
    group0.add_option('--clusterDistance',type="float",dest='clusterDistance',default=300,
    help="The minimum distance between two cluster transcription start site [default: 300]")




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
    # elif os.path.dirname(options.out_dir) == "":
    #     out_dir= "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    #print(out_dir)

    if options.cdrFile is None:
        print("Error: Need --cdrFile for cell barcode file.")
        sys.exit(1)

    if options.refFastq is None:
        print("Error: Need --refFastq for reference fastq file.")
        sys.exit(1)


    #bam file
    if options.bam_file is None:
        print("Error: Need --bam for aligned file.")
        sys.exit(1)
    
   
        
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
    #isoformNumber=options.isoformNumber
    cellBarcodePath=options.cdrFile
    n_proc=options.nproc
    maxReadCount=options.maxReadCount
    clusterDistance=options.clusterDistance
    fastqFilePath=options.refFastq




        
    if options.mode == "Old":
        filterTssPath=get_filter_TSS(tssdf,ref_out_dir)


        getTSScount=get_old_TSS_count(generefpath,filterTssPath,bam_file,fastqFilePath,out_dir,cellBarcodePath,n_proc,minCount,maxReadCount,clusterDistance)
        scadata=getTSScount.produce_sclevel()

    else:
        getTSScount=get_TSS_count(generefpath,tssrefpath,bam_file,fastqFilePath,out_dir,cellBarcodePath,n_proc,minCount,maxReadCount,clusterDistance)
        scadata=getTSScount.produce_sclevel()

        run_time = time.time() - START_TIME
        print("[scTSS-count] All done: %d min %.1f sec" %(int(run_time / 60), 
                                                  run_time % 60))
