============================================================
CamoTSS for alternative TSS analysis in single cells
============================================================
|pypi| 

.. |pypi| image:: https://badge.fury.io/py/CamoTSS.svg
       :target: https://pypi.org/project/CamoTSS/




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

.. _figshare: https://figshare.com/articles/dataset/CamoTSS_test_data/22641031

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

For the remaining modes, you can check our manual_.

.. _manual: https://camotss.readthedocs.io/en/latest/run_CamoTSS.html

.. code-block:: bash

   #!/bin/bash 
   gtfFile= $CamoTSS/test/Homo_sapiens.GRCh38.105.chr_test.gtf
   fastaFile = $download/genome.fa
   bamFile= $download/possorted_genome_bam_filtered.bam
   cellbarcodeFile=$CamoTSS/test/cellbarcode_to_CamoTSS

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

Hou, R., Hon, C. C., & Huang, Y. (2023). CamoTSS: analysis of alternative transcription start sites for cellular phenotypes and regulatory patterns from 5'scRNA-seq data. bioRxiv, 2023-04.













