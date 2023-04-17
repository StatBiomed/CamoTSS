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

Here are two modes in CamoTSS : **TC** and **CTSS**. 

You can run CamoTSS by using test file according to the following code.

.. code-block:: bash

   #!/bin/bash 
   gtfFile= $CamoTSS/test/Homo_sapiens.GRCh38.105.chr_test.gtf
   fastaFile = $download/genome.fa
   bamFile= $download/possorted_genome_bam_filtered.bam
   cellbarcodeFile=$CamoTSS/test/cellbarcode_to_CamoTSS

   CamoTSS --gtf gtfFile --refFasta fastaFile --bam bamFile -c cellbarcodeFile -o CamoTSS_out --mode TC


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

Coming soon













