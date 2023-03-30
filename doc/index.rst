|PyPI| |Docs| |Build Status|

.. |PyPI| image:: https://badge.fury.io/py/CamoTSS.svg
       :target: https://pypi.org/project/CamoTSS/
.. |Docs| image:: https://readthedocs.org/projects/CamoTSS/badge/?version=latest
      :target: https://CamoTSS.readthedocs.io
.. |Build Status| image:: https://travis-ci.org/huangyh09/brie.svg?branch=master
      :target: https://travis-ci.org/huangyh09/brie



====
Home
====



About CamoTSS
==================

CamoTSS can  precisely identify TSS and quantify its expression by leveraging the cDNA on read 1, which enables effective detection of alternative TSS usage.

CamoTSS supports the analysis of alternative TSS at two levels, TSS cluster (TC mode) and CTSS (CTSS mode). 
The input files of CamoTSS include alignment file (bam file), annotation file (gtf file), cell list file and reference genome file (fasta file). 
The output files of CamoTSS include cell by all TSSs matrix (h5ad), cell by two TSSs matrix (h5ad), cell by CTSS matrix (h5ad) and cell by CTSS matrix (h5ad). 

CamoTSS identify TSS through the following steps.

.. image:: image/flow_chart.png
   :width: 600
   :alt: flow chart for CamoTSS
   :align: center


Specificlly, a convolutional neural network was applied to filter false positive peaks.

.. image:: image/classifier.png
   :width: 600
   :alt: classifier
   :align: center




It includes 3 steps to identify alternative TSS or CTSS usage: preprocessing, running CamoTSS and running Brie2.

Please refer to our tutorial for details.


* `Preprocess for one sample and multiple samples`_.

* `Run CamoTSS`_.

* `Run Brie2`_.

.. _Preprocess for one sample and multiple samples: melanoma.ipynb

.. _Run CamoTSS: differential_test_intestine.ipynb

.. _Run Brie2: differential_test_intestine.ipynb


    
