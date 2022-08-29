import os
import sys
import pandas as pd
import numpy as np 
import scanpy as sc
from optparse import OptionParser, OptionGroup
from ..utils.get_inputfile_toBrie import get_brie_input
import subprocess
from ..version import __version__



def main():
    parser=OptionParser()
    parser.add_option("--geneExpression","-g",dest="gene_expression",default=None,
    help="Gene expression file from cellranger count (The default directory name is "
     "filtered_feature_bc_matrix in cellranger count) ")

    parser.add_option('--cellInfo','-c',dest="cell_info",default=None,
    help="Seven columns tsv file. The first column is cell barcode. " 
    "The second column is corresponding cell clusters or condition info, such as normal and disease. "
    "This column info should be provided according to your mode selection "
    "The 3td to 7th column is the first 5 PCs for your gene expression data")

    parser.add_option('--countOut',dest="countOutPath",default=None,
    help="The path of scTSS-count output")


    parser.add_option('--mode','-m',dest='mode',default=None,
    help="Two options: cluster or disease. If you want to detect TSS markers in cell clusters, "
    "you can input cluster. At the same time, you should add cluster information in the --cellInfo file. "
    "If you want to detect differential TSS expression between two conditions, such as disease and normal, "
    "you can input disease.")

    parser.add_option('--outdir','-o',dest='out_dir',default=None,help='The directory for output [default : $alternativeTSSPath]')

    group0=OptionGroup(parser,"Optional arguments")

    group0.add_option("--mincellnum",type="int",dest="mincellnum",default=50,
    help="Minimum cell number for each cells [default: 50]")


    parser.add_option_group(group0)
    (options, args) = parser.parse_args()




    #this means that if users do not input any argument, then direct produce help. then end.
    if len(sys.argv[1:]) == 0:
        print('Welcome to scTSS-quant v%s!\n'%(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)


    #output file 
    if options.out_dir is None:
        print("Warning: no outDir provided, we use $cell_info_path/scTSS\n")
        out_dir = os.path.dirname(os.path.abspath(options.alternativeTSS)) + "/scTSS_quant"
    # elif os.path.dirname(options.out_dir) == "":
    #     out_dir= "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    # print(out_dir)


    
    #raw data 
    if options.gene_expression is None:
        print("Error: Need --geneExpression for 10x Cellranger output expression file")
        sys.exit(1)
    
    if options.cell_info is None:
        print("Error: Need --cellInfo for cell cluster or disease information")
        sys.exit(1)

    if options.mode is None:
        print("Error: Need --mode for cluster or disease")
        sys.exit(1)

    if options.countOutPath is None:
        print("Error: Need --countOut for path of scTSS-count output")
        sys.exit(1)


    

    
    rawExpFilePath=options.gene_expression
    cellInfoPath=options.cell_info
    mode=options.mode
    splicingFilePath=options.countOutPath
    cellnumThreshold=options.mincellnum
 


    quant_out_dir=str(out_dir)+'/quant/'
    if not os.path.exists(quant_out_dir):
        os.mkdir(quant_out_dir)
    # print(quant_out_dir)



    # get file to brie2
    getFile=get_brie_input(rawExpFilePath,splicingFilePath,cellInfoPath,quant_out_dir,cellnumThreshold)
    brie_h5ad_input,originadata=getFile.get_h5adFile()
    brie_cdr_input,numls=getFile.get_cluster_cdrFile(mode,originadata)
    outputFile=str(quant_out_dir)+'Diff_TSS_Gene.h5ad'

    #print(originadata)
    #print(brie_cdr_input)

    numls=[str(i) for i in numls]
    inputnum=','.join(numls)
    #print(inputnum)
    

    #run BRIE2
    bashCommond="brie-quant -i %s -c %s -o %s --batchSize 1000000 --minCell 250  --interceptMode gene --testBase null --LRTindex=All" %(brie_h5ad_input,brie_cdr_input,outputFile)
    print(bashCommond)
    process = subprocess.Popen(bashCommond.split(), stdout=subprocess.PIPE)
    #print(process)
    _output = process.communicate()[0]




    
        

    
