from optparse import OptionParser,OptionGroup
from ..version import __version__
import sys
from ..utils.build_ref import get_TSSref,get_generef,get_filter_TSS
from ..utils.get_counts import get_TSS_count
import pyranges as pr
import os
import pandas as pd
import time 


START_TIME = time.time()


def main():
    parser = OptionParser()
    parser.add_option('--gtf','-g',dest='gtf_file',default=None,help='The annotation gtf file for your analysing species.')
    parser.add_option('--cellbarcodeFile','-c',dest='cdrFile',default=None,help='The file include cell barcode which users want to keep in the downstream analysis.')
    parser.add_option('--bam','-b',dest='bam_file',default=None,help='The bam file of aligned from Cellranger or other single cell aligned software.')
    parser.add_option('--outdir','-o',dest='out_dir',default=None,help='The directory for output [default : $bam_file]') #what should be after $
    parser.add_option('--refFasta','-r',dest='refFasta',default=None,help='The directory for reference genome fasta file') #what should be after $
    parser.add_option('--mode','-m',dest='mode',default=None,help='You can select run by finding novel TSS cluster mode [TC]. If you also want to detect CTSS within one cluster, you can use [CTSS] mode')

   
   
    group0=OptionGroup(parser,"Optional arguments")

    group0.add_option("--minCount",type="int",dest="minCount",default=50,
    help="Minimum UMI counts for TC in all cells [default: 50]")

    group0.add_option('--nproc','-p',type="int",dest='nproc',default=4,
    help='Number of subprocesses [default: 4]')

    group0.add_option('--maxReadCount',type="int",dest='maxReadCount',default=10000,
    help='For each gene, the maxmium read count kept for clustering [default: 10000]')
    
    group0.add_option('--clusterDistance',type="float",dest='clusterDistance',default=300,
    help="The minimum distance between two cluster transcription start site [default: 300]")

    group0.add_option('--InnerDistance',type="float",dest='InnerDistance',default=100,
    help="The resolution of each cluster [default: 100]")

    group0.add_option('--windowSize',type="float",dest='windowSize',default=15,
    help="The width of sliding window [default: 15]")

    group0.add_option('--minCTSSCount',type="float",dest='minCTSSCount',default=100,
    help="The minimum UMI counts for each CTSS [default: 100]")

    group0.add_option('--minFC',type="float",dest='minFC',default=6,
    help="The minimum fold change for filtering CTSS [default: 6]")






    parser.add_option_group(group0)
    (options, args) = parser.parse_args()

    
    #this means that if users do not input any argument, then direct produce help. then end.
    if len(sys.argv[1:]) == 0:
        print('Welcome to CamoTSS v%s!\n'%(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)


    #output file 
    if options.out_dir is None:
        print("Warning: no outDir provided, we use $bamfilePath/CamoTSS\n")
        out_dir = os.path.dirname(os.path.abspath(options.bam_file)) + "/CamoTSS"
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

    if options.refFasta is None:
        print("Error: Need --refFasta for reference fasta file.")
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
    cellBarcodePath=options.cdrFile
    n_proc=options.nproc
    maxReadCount=options.maxReadCount
    clusterDistance=options.clusterDistance
    InnerDistance=options.InnerDistance
    fastqFilePath=options.refFasta
    windowSize=options.windowSize
    minCTSSCount=options.minCTSSCount
    minFC=options.minFC




        
    if options.mode == "TC":
        getTSScount=get_TSS_count(generefpath,tssrefpath,bam_file,fastqFilePath,out_dir,cellBarcodePath,n_proc,minCount,maxReadCount,clusterDistance,InnerDistance,windowSize,minCTSSCount,minFC)
        scadata=getTSScount.produce_sclevel()

    elif options.mode=="CTSS":
        getTSScount=get_TSS_count(generefpath,tssrefpath,bam_file,fastqFilePath,out_dir,cellBarcodePath,n_proc,minCount,maxReadCount,clusterDistance,InnerDistance,windowSize,minCTSSCount,minFC)
        scadata=getTSScount.produce_sclevel()
        twoctssadata=getTSScount.produce_CTSS_adata()

    else:
        print('Do not have this mode. Please check your spell!')
        run_time = time.time() - START_TIME
        print("[CamoTSS] All done: %d min %.1f sec" %(int(run_time / 60), 
                                                  run_time % 60))


