============================================================
scTSS: Detection and couting alternative TSS in single cells
============================================================

Installation
============

You can install from this GitHub repository for latest (often development) 
version by following command line

.. code-block:: bash

  pip install -U git+https://github.com/StatBiomed/scTSS

In either case, add ``--user`` if you don't have the write permission for your 
Python environment.


Quick start
===========

scTSS-count
===========

**STEP1:   Processing**


scTSS mainly deal with the output from cellranger (a common alignment tool for 10x data).

The preprocessing procedure based on the output file of cellranger. 

.. code-block:: bash

    1. cd /cellranger_out/outs
    2. samtools view  possorted_genome_bam.bam | LC_ALL=C grep "xf:i:25" > body_filtered_sam
    3. samtools view -H possorted_genome_bam.bam > header_filted_sam
    4. cat header_filted_sam body_filtered_sam > possorted_genome_bam_filterd.sam
    5. samtools view -b possorted_genome_bam_filterd.sam > possorted_genome_bam_filterd.bam
    6. samtools index possorted_genome_bam_filterd.bam possorted_genome_bam_filterd.bam.bai
 
**STEP2:   Run scTSS-count**

.. code-block:: bash

        scTSS-count --gtf $gtfFile --refFastq $fastFile --bam $possorted_genome_bam_filterd.bam -c $cluster_toscTSS.tsv  -o $output_fileFold --mode Unannotation

Want to learn about more parameter, you can use ``scTSS-count --help`` to check. 

You can find out the example file in the test folder. Please make sure you also have the same column name.

Here, you can select one of the mode from annotation and unannotation. 

Unannotation means that you can detect novel TSS. The distance between different TSS may be wide. 

Annotation means that you can detect TSS based on the annotation. The distance between different TSS may be narrow.

You can check our paper to learn more detail. 


scTSS-quant
===========

.. code-block:: bash

        scTSS-quant -g $filtered_feature_bc_matrix -c $cluster_toscTSS.tsv --countOut $scTSS_count_folder -m cluster -o $scTSS_quant_folder

Please use ``scTSS-quant --help`` to check more parameter.

Here, you can select one of the mode from cluster and disease. 

Cluster means that you can detect cell type-specific TSS. Any multiple groups detection can use this mode.

Disease mode help you select disease-specific TSS. Any two groups detection can use this mode. 





