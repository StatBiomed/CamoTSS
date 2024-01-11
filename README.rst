============================================================
CamoTSS for alternative TSS analysis in single cells
============================================================
|pypi| 

.. |pypi| image:: https://badge.fury.io/py/CamoTSS.svg
       :target: https://pypi.org/project/CamoTSS/

.. image:: https://zenodo.org/badge/497821671.svg
      :target: https://zenodo.org/badge/latestdoi/497821671


Note
============
Hi there, my github account did not notify me when there are issue. 
So if you are in a hurry, you can email me. ruiyan@connect.hku.hk.
I check email every day.  



Installation
============

You can install from this GitHub repository for latest (often development) 
version by following command line

.. code-block:: bash

  pip install -U git+https://github.com/StatBiomed/CamoTSS

In either case, add ``--user`` if you don't have the write permission for your 
Python environment.


Quick start
===========

Download test file
===================

You can download test file from figshare_.

.. _figshare: https://figshare.com/projects/CamoTSS/184603

Here, you can download some large file include genome.fa, possorted_genome_bam_filtered.bam.

Alternatively, you can also download the reference genome fasta file from Ensembl or Genecode or website of 10x Genomics. 
 
Run CamoTSS 
=============

Here are three modes in CamoTSS : **TC+CTSS** , **TC** and **CTSS**.

When you run **TC+CTSS** mode, you will get TC result and then get the CTSS result based on the TC.

When you run **TC** mode, you will only get the TC result.

The **TC+CTSS** and **TC** mode have the same required files.

The --outdir is the only required parameter for **CTSS** mode. But the outdir should include output of TC.  

If you want to run **CTSS** mode, you must based on the output of TC.

You can run CamoTSS **TC+CTSS** mode by using test file according to the following code.

For the remaining modes, you can check this document_.

.. _document: https://camotss.readthedocs.io/en/latest/run_CamoTSS.html

.. code-block:: bash

   #!/bin/bash 
   gtfFile= $download/Homo_sapiens.GRCh38.105.chr_test.gtf
   fastaFile = $download/genome.fa
   bamFile= $download/possorted_genome_bam_filtered.bam
   cellbarcodeFile=$download/cellbarcode_to_CamoTSS

   CamoTSS --gtf gtfFile --refFasta fastaFile --bam bamFile -c cellbarcodeFile -o CamoTSS_out --mode TC+CTSS


Alternative TSS or CTSS detecting
=================================

To identify alternative TSS usage or alternative CTSS usage, Brie2 (Huang & Sanguinetti, 2021) is recommend to be used. 

Here, we provide an example exploiting BRIE2 to detect alterntive TSS/CTSS usage. 

You can check it in our manual_.

.. _manual: https://camotss.readthedocs.io/en/latest/runBRIE.html  


Detailed Manual
================

The full manual is here_, including:

`Preprocess`_

`Run CamoTSS`_

`Detect alternative TSS/CTSS`_

.. _here: https://camotss.readthedocs.io/en/latest/index.html

.. _Preprocess: https://camotss.readthedocs.io/en/latest/preprocess.html

.. _Run CamoTSS: https://camotss.readthedocs.io/en/latest/run_CamoTSS.html

.. _Detect alternative TSS/CTSS: https://camotss.readthedocs.io/en/latest/runBRIE.html



Reference
===========

Hou, R., Hon, CC. & Huang, Y. CamoTSS: analysis of alternative transcription start sites for cellular phenotypes and regulatory patterns from 5' scRNA-seq data. Nat Commun 14, 7240 (2023). https://doi.org/10.1038/s41467-023-42636-1












