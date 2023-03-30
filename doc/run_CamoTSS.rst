==============
Run CamoTSS
==============

CamoTSS includes two kind of modes : TC mode and CTSS mode. 

The input files include:

* alignment file (bam file)
* annotation file (gtf file)
* cell list file and reference genome file (fasta file)
* cell barcode list file (csv file)

The output files include:

* cell by all TSSs matrix (h5ad)
* cell by two TSSs matrix (h5ad) 
* cell by CTSS matrix (h5ad)
* cell by CTSS matrix (h5ad) 
* ref file (reference gene and TSS csv)

Here is a quick test file. You can check it.
  
Download test file
===================

You can download test file from onedrive_.

.. _onedrive: https://connecthkuhk-my.sharepoint.com/:f:/g/personal/ruiyan_connect_hku_hk/Eqp1gYR5dlVIoWgH0udyJ5YB_9eVQ1e5WAxx3muAIeYdjw?e=SQ7fgb

Here, you can download some large file include genome.fa, possorted_genome_bam_filtered.bam.

Alternatively, you can also download the reference genome fasta file from Ensembl or Genecode or website of 10x Genomics.

Run CamoTSS
=============

Here are two modes in CamoTSS : **TC** and **CTSS**.

You can run CamoTSS by using test file according to the following code.

.. code-block:: bash

   #!/bin/bash
   gtfFile= $CamoTSS/test/Homo_sapiens.GRCh38.105.chr_test.gtf
   fastaFile = $download/genome.fa
   bamFile= $download/possorted_genome_bam_filtered.bam
   cellbarcodeFile=$CamoTSS/test/cellbarcode_to_CamoTSS

   CamoTSS --gtf gtfFile --refFasta fastaFile --bam bamFile -c cellbarcodeFile -o CamoTSS_out --mode TC


Options
========


There are more parameters for setting (``CamoTSS -h`` always give the version
you are using):


.. code-block:: html

   Usage: CamoTSS [options]

   Options:
        -h, --help            show this help message and exit
        -g GTF_FILE, --gtf=GTF_FILE
                        The annotation gtf file for your analysing species.
        -c CDRFILE, --cellbarcodeFile=CDRFILE
                        The file include cell barcode which users want to keep
                        in the downstream analysis.
        -b BAM_FILE, --bam=BAM_FILE
                        The bam file of aligned from Cellranger or other
                        single cell aligned software.
        -o OUT_DIR, --outdir=OUT_DIR
                        The directory for output [default : $bam_file]
        -r REFFASTA, --refFasta=REFFASTA
                        The directory for reference genome fasta file
        -m MODE, --mode=MODE  You can select run by finding novel TSS cluster mode
                        [TC]. If you also want to detect CTSS within one
                        cluster, you can use [CTSS] mode

   Optional arguments:
        --minCount=MINCOUNT
                        Minimum UMI counts for TC in all cells [default: 50]
        -p NPROC, --nproc=NPROC
                        Number of subprocesses [default: 4]
        --maxReadCount=MAXREADCOUNT
                        For each gene, the maxmium read count kept for
                        clustering [default: 10000]
        --clusterDistance=CLUSTERDISTANCE
                        The minimum distance between two cluster transcription
                        start site [default: 300]
        --InnerDistance=INNERDISTANCE
                        The resolution of each cluster [default: 100]
        --windowSize=WINDOWSIZE
                        The width of sliding window [default: 15]
        --minCTSSCount=MINCTSSCOUNT
                        The minimum UMI counts for each CTSS [default: 100]
        --minFC=MINFC       The minimum fold change for filtering CTSS [default:
                        6]


