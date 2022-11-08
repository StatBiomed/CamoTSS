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

STEP1:   Processing
===========

scTSS mainly deal with the output from cellranger (a common alignment tool for 10x data).

The preprocessing procedure based on the output file of cellranger. 

.. code-block:: bash

    1. cd /cellranger_out/outs
    2. samtools view  possorted_genome_bam.bam | LC_ALL=C grep "xf:i:25" > body_filtered_sam
    3. samtools view -H possorted_genome_bam.bam > header_filted_sam
    4. cat header_filted_sam body_filtered_sam > possorted_genome_bam_filterd.sam
    5. samtools view -b possorted_genome_bam_filterd.sam > possorted_genome_bam_filterd.bam
    6. samtools index possorted_genome_bam_filterd.bam possorted_genome_bam_filterd.bam.bai
 
STEP2:   Run scTSS-count
===========

scTSS-count --gtf $gtfFile --refFastq $fastFile --bam $possorted_genome_bam_filterd.bam -c $cell_idFile -o $output_fileFold --mode New



