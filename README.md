# scATACseq_analysis_pipeline
This is a pipeline for analyzing single cell ATACseq.

The analysis takes the bam files from Cellranger, the mapped reads are quality filtered, duplicate removed, and barcode fixed. 

After that, MACS2 is used to call peaks and the blacklist (a bed file showing the known repeat regions) is removed from the called accessible chromatin regions (ACRs).