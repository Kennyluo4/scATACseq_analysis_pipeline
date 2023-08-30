# **scATACseq_analysis_pipeline**
Ziliang Luo


This is a pipeline for analyzing single cell ATACseq data.

It takes the mapping bam files (e.g. from Cellranger), filters the low quality reads by SAMtools, removes PCR duplicates by Picard, modifies the `CB` barcode, and makes Tn5 insertion bed files. 

## Prerequisite
Make sure the these software (samtools, picard, bedtools, and macs2) are installed/loaded in the PATH. Python module `pysam` is also required.

```
module load SAMtools
module load picard
module load BEDTools
module load MACS2
```


## Proccess the bam files 

This step calls `FixingBarcodeName.py` and `makeTn5bed.py` for the ananlysis, make sure their path is correct. 

```
usage: bam2bed_scatac_8.23.23.py [-h] -b <bam> -s <sample_name> -o <output_dir> [-x <cpu#>] -r <ref_index> [-l <black_list>]

This script processes raw mapping bam file (e.g. from Cellranger), and make Tn5 insertion bed file. 

options:
  -h, --help            show this help message and exit
  -b <bam>, --bam <bam>     the input bam file
  -s <sample_name>, --sample <sample_name>      the sanmple name
  -o <output_dir>, --out <output_dir>     the output directory
  -x <cpu#>, --cpu <cpu#>     the number of cores used
  -r <ref_index>, --ref <ref_index>     the genome reference index. Use samtools faidx to build the index
  -l <black_list>, --bklst <black_list>     a bed file black list that contains the low-complexity and homopolymeric regions, organelle sequence regions, Tn5 cutting bias regions and potential collapsed regions in the reference genome. (see methods in [Marand et al., 2021](https://doi.org/10.1016/j.cell.2021.04.014))
```

## Coverting processed bed file (peaks) to bw file for browser view

```
python bed2bw_scatac_8.23.23.py <script_dir> <ref_idx> <bedGraphToBigWig_dir> <sample_name>
```


