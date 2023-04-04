# **scATACseq_analysis_pipeline**
This is a pipeline for analyzing single cell ATACseq.

The analysis takes the bam files from Cellranger, the mapped reads are quality filtered, duplicate removed, and barcode fixed. 

After that, MACS2 is used to call peaks and the blacklist (a bed file showing the known repeat regions) is removed from the called accessible chromatin regions (ACRs).

## Prerequisite
Make sure the following softwares are installed/loaded in the PATH.


`module load SAMtools`

`module load picard`

`module load BEDTools`

`module load MACS2`

## Proccess the bam files 

`python bam2bed_processing_033023.py <cellranger_dir> <sample_name> <output_dir> <cpu#> <ref_index> <black_list>`
## Coverting processed bed file (peaks) to bw file for browser view

` python bed2bw_040423.py <script_dir> <ref_idx> <bedGraphToBigWig_dir> <sample_name> > `


