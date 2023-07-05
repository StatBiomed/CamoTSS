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

* fetch_reads.pkl :A dictionary whose key is the gene id, the value is the reads information of this gene, including the position 0 of reads, cellbarcodes and cigar string, such as (805735, ‘CATGACATCTCAACTT-1’, ‘14S49M’).
* before_cluster_peak.pkl :A dictionary file whose key is gene id, the value is reads information of this gene including three numpy array. The first is  all the coordination of position 0 of reads. The second are all cell-barcodes of reads. The third are all cigar string of reads.  
* fourFeature.csv : A dataframe includes 8 columns: cluster_id; UMI_count of cluster; SD of cluster; summit count of cluster;the percentage of Unencoded G percentage; the order of TSS; gene id; summit position. 
* afterfiltered.csv : A dataframe have the same columns name as the fourFeature.csv. This file includes peaks filtered by classifier.
* keepdict.pkl :  A dictionary which includes details of peaks in afterfiltered.csv. 
* scTSS_count_all.h5ad: An anndata whose X is cell by TSS matrix. This file contained all TSS detected by CamoTSS. 
* scTSS_count_two.h5ad: An anndata whose X is cell by TSS matrix. This file exclusively includes genes that possess two or more TSSs. 
* CTSS_foldchange.pkl : A dictionary whose keys are peaks obtained at the first step and values are all CTSSs within this cluster and the related count and fold change. 
* all_ctss.h5ad : An anndata whose X is cell by CTSS matrix. This file contained all CTSS detected by CamoTSS. 
* all_ctss_two.h5ad : An anndata whose X is cell by CTSS matrix. This file contained TSS which have two or more CTSS.


Here is a quick test file. You can check it.
  
Download test file
===================

You can download test file from figshare_.

.. _figshare: https://figshare.com/articles/dataset/CamoTSS_test_data/22641031

Here, you can download some large file include genome.fa, possorted_genome_bam_filtered.bam.

Alternatively, you can also download the reference genome fasta file from Ensembl or Genecode or website of 10x Genomics.

Run CamoTSS
=============

Here are three modes in CamoTSS : **TC** , **CTSS** and **TC+CTSS**.

* TC : Just detect TSS cluster.  
* CTSS : Just detect CTSS within one cluster. But you should have the output from TC as the input to CTSS. The aim to add this mode is to prevent to rerun CamoTSS when user want to analysis CTSS.
* TC+CTSS : Directly to detect TSS cluster and CTSS within one TSS cluster. 



You can run CamoTSS by using test file according to the following code to run TC+CTSS mode.

.. code-block:: bash

   #!/bin/bash
   gtfFile= $download/Homo_sapiens.GRCh38.105.chr_test.gtf
   fastaFile = $download/genome.fa
   bamFile= $download/possorted_genome_bam_filtered.bam
   cellbarcodeFile=$download/cellbarcode_to_CamoTSS

   CamoTSS --gtf gtfFile --refFasta fastaFile --bam bamFile -c cellbarcodeFile -o CamoTSS_out --mode TC+CTSS


You can run CamoTSS by using test file according to the following code to run TC mode. 

.. code-block:: bash

   #!/bin/bash
   gtfFile= $download/Homo_sapiens.GRCh38.105.chr_test.gtf
   fastaFile = $download/genome.fa
   bamFile= $download/possorted_genome_bam_filtered.bam
   cellbarcodeFile=$download/cellbarcode_to_CamoTSS

   CamoTSS --gtf gtfFile --refFasta fastaFile --bam bamFile -c cellbarcodeFile -o CamoTSS_out --mode TC


You can run CamoTSS by using test file according to the following code to run CTSS mode. 

.. code-block:: bash

   #!/bin/bash 
   #note: the output file path should be the parent path of CamoTSS.
   outputfile=CamoTSS_out 

   CamoTSS -m CTSS -o $outputfile





Options
========


There are more parameters for setting (``CamoTSS -h`` always give the version
you are using):


.. code-block:: html

   Usage: CamoTSS [options]

   Options:
        -h, --help      show this help message and exit
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
        -m MODE, --mode=MODE  You can select run by finding novel TSS cluster and
                        CTSS within one cluster [TC+CTSS].
                        If you just want to detect TSS cluster, you can use
                        [TC] mode. If you just want to detect CTSS, you can
                        use [CTSS] mode which is based on the output of [TC
                        mode]

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

  Optional arguments:
        --minCTSSCount=MINCTSSCOUNT
                        The minimum UMI counts for each CTSS [default: 100]
        --minFC=MINFC       The minimum fold change for filtering CTSS [default:
                        6]



