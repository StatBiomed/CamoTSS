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

run cellranger count
====================
CamoTSS require 5'scRNA-seq data (10x Genomics) with the length of reads 1 more than 100bp.

That means the reads 1 should contain extra cDNA information except UMI and cell barcode. 

To align the extra cDNA, '--chemistry SC5P-PE' should be set during running cellranger count. 

 

CamoTSS
===========

**STEP1:   Processing**


CamoTSS mainly deal with the output from cellranger (a common alignment tool for 10x data).

The preprocessing procedure based on the output file of cellranger. 

.. code-block:: bash

    1. cd /cellranger_out/outs
    2. samtools view  possorted_genome_bam.bam | LC_ALL=C grep "xf:i:25" > body_filtered_sam
    3. samtools view -H possorted_genome_bam.bam > header_filted_sam
    4. cat header_filted_sam body_filtered_sam > possorted_genome_bam_filterd.sam
    5. samtools view -b possorted_genome_bam_filterd.sam > possorted_genome_bam_filterd.bam
    6. samtools index possorted_genome_bam_filterd.bam possorted_genome_bam_filterd.bam.bai
 
**STEP2:   Run CamoTSS**

.. code-block:: bash

        CamoTSS --gtf $gtfFile --refFastq $fastFile --bam $possorted_genome_bam_filterd.bam -c $cluster_toscTSS.tsv  -o $output_fileFold --mode Unannotation

Want to learn about more parameter, you can use ``CamoTSS --help`` to check. 

You can find out the example file in the test folder. Please make sure you also have the same column name.

Here, you can select one of the mode from "Unannotation" and "Unannotation_addCTSS". 

Unannotation means that you can detect novel TSS cluster. 

Unannotation_addCTSS  means that you can detect CTSS within one cluster. 

You can check our paper to learn more detail.



Multiple samples preprocessing
==============================

For most public single cell data, we can obtain the whole annotation of cell type from different samples. 

The sample ID information always show at the cell barcode for each cell.

In order to fully use the annotation described above, we can run cellranger count for each sample independently. 

Then manually add sample information to the cell barcode. We can implement it by using following script.

.. code-block:: python

        import pysam
        inputbamfile=$home+'/cellranger_out/outs/manual_filter/possorted_genome_bam_filterd.bam'
        outputbamfile=$home+'/cellranger_out/outs/manual_filter/possorted_genome_bam_filterd_add_suffix.bam'
        inputbam=pysam.Samfile(inputbamfile,'rb')
        outputbam=pysam.Samfile(outputbamfile,'wb',template=inputbam)
        for read in inputbam.fetch():
                cb=read.get_tag('CB')
                assert cb is not None
                cbfix=cb.replace('-1',"")
                cbfix=cbfix+'-sampleID'
                read.set_tag('CB',cbfix)
                outputbam.write(read)
        inputbam.close()
        outputbam.close()
        

Then the bam file with changed cellbarcode can be merged with samtools merge

.. code-block:: bash

        samtools merge $merged_bam -b $bamlist.fofn --write-index



Alternative TSS or CTSS detecting
=================================
In CamoTSS, one of output files is Tobrie.h5ad  which can be as input to Brie. 

To identify alternative TSS usage or alternative CTSS usage, Brie2 (Huang & Sanguinetti, 2021) is recommend to be used. 

For more information, please check https://brie.readthedocs.io/en/latest/ 






