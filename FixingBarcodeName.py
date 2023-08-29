import pysam
import argparse
import ast
import os


## From pablo # Sohyun edited for cellranger v2 and do not consider non nuclear configs
## Usage: python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/4_BarcodeArrange/4-1_FixingBarcodeName.py \
# -BAM "${List[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam -exp_name Ex | samtools view -@ 12 -h - > ../4.Bam_FixingBarcode/"${List[SLURM_ARRAY_TASK_ID]}"_BarcodeFixed.sam

def read_bam_file(bam_file,exp_name):
    """
    First alter the CB tags as well as the other tags
    Next - count the tag and add tag to dictionary if not presenat and start at
    one. Add to total column, add to nuclear/nonnuclear column depending on the
    scaffold name and the list given above
    """
    #Read the File
    #save = pysam.set_verbosity(0)
    #read_bam_file = pysam.AlignmentFile(bam_file,"rb", ignore_truncation=True)
    #pysam.set_verbosity(save)
    read_bam_file = pysam.AlignmentFile(bam_file,"rb")

    outfile = pysam.AlignmentFile("-", "w", template=read_bam_file)
    for read in read_bam_file :
        #For each read alter CB tag
        try:
            original_tag = (read.get_tag("CB"))
            #print("here")
            #print(original_tag)
            exp_name_tag = "-" + exp_name
            new_tag = original_tag.replace("-1", exp_name_tag)
            read.set_tag("CB", new_tag, replace=True)
            outfile.write(read)
        except KeyError:
            pass


def get_parser():
    parser = argparse.ArgumentParser(description='Pull our reads aligning to a\
        region from multiple list of BAM files, puts them into a BAM file\
        for later assembly.')
    parser.add_argument('-BAM', "--bam_file", help="Bam file to \
        pull reads from.", required=True, dest='bam_f')
    parser.add_argument('-exp_name', "--experiment_name", help="10x config file to \
        pull scaffold names of nulcear and non-nuclear scaffolds",
        required=True, dest='exp')

    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    #Load all Bed files
    gathered_read_dict = read_bam_file(args.bam_f, args.exp)
