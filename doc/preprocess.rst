===========
Preprocess
===========

run cellranger count
====================
CamoTSS require 5'scRNA-seq data (10x Genomics) with the length of reads 1 more than 100bp.

That means the reads 1 should contain extra cDNA information except UMI and cell barcode.

To align the extra cDNA, '--chemistry SC5P-PE' should be set during running cellranger count.


Filter
=======
CamoTSS usually take alignment file from cellranger count as input. 
Here, we try to filter 10x Genomic BAM file (i.e. possorted_genome_bam.bam) to
make it cooperate with filtered_feature_bc_matrix. 

Filtering UMI according to the specific criteria:

.. code-block:: html

* Have a MAPQ score of 255
* Maps to exactly one gene
* Overlaps an exon by at least 50% in a way consistent with annotated splice junctions and strand annotation. Records that align to exons will have an RE:A:E tag.
* Remove any records with matching UMI and Barcode values that map to different genes.


For more information, you can check this material_ .

.. _material: https://www.10xgenomics.com/resources/analysis-guides/tutorial-navigating-10x-barcoded-bam-files 

To get cleaner bam, you can process possorted_genome_bam.bam according to the following steps.

.. code-block:: linux

   1. cd /cellranger_out/outs
   2. samtools view  possorted_genome_bam.bam | LC_ALL=C grep "xf:i:25" > body_filtered_sam
   3. samtools view -H possorted_genome_bam.bam > header_filted_sam
   4. cat header_filted_sam body_filtered_sam > possorted_genome_bam_filterd.sam
   5. samtools view -b possorted_genome_bam_filterd.sam > possorted_genome_bam_filterd.bam
   6. samtools index possorted_genome_bam_filterd.bam possorted_genome_bam_filterd.bam.bai


For multiple samples
====================

In most instances, we want to detect alternative TSS usage among various samples.
We should run CamoTSS at the same batch for all of these samples. 
To accomplish this goal, it is necessary to merge bam file of these samples to one bam file. 
We should add suffix to cell barcode of each sample in order to distinguish origin of samples.
The following script can help you manually add sample information to cell barcode.

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

.. code-block:: linux 

   samtools merge $merged_bam -b $bamlist.fofn --write-index

