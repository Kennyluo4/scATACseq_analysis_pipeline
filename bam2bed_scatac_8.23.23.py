#!/usr/bin/env python
''' 
This script is to process bam files from the cell ranger output
    edited by Ziliang 2023.03 
        change: make the output directory if it doesn't exist;
                merge samtools processing steps; 
                merge remove blacklist step.
  '''
## 
import sys
import subprocess
import os
import argparse
  
def get_parser():
  ''' read the arguments'''
  parser = argparse.ArgumentParser(description='This script processes raw mapping bam file (e.g. from Cellranger), and make Tn5 insertion bed file.\n \
                                                ... make sure picard and SAMtools are loaded.\n \
                                                ... make sure the makeTn5bed.py and FixingBarcodeName.py are in the same folder')
  parser.add_argument('-b', '--bam', help=' the input bam file', dest='bam')
  parser.add_argument('-s', '--sample', help=' the sanmple name', dest='sample')
  parser.add_argument('-o', '--out', help=' the output directory', dest='out')
  parser.add_argument('-x', '--cpu', help=' the number of cores used', dest='cpu')
  parser.add_argument('-r', '--ref', help=' the reference index file', dest='ref')
  parser.add_argument('-l', '--bklst', help=' the bed file of blacklist regions in the genome', dest='bklst')
  return parser
  
def process_bam2bed(input_bam, output_dir, core_num, sample_name, qual):
    #make output sample folder in the output directory 
    output_dir = output_dir + '/' + sample_name
    isExist = os.path.exists(output_dir)
    if not isExist:
      os.mkdir(output_dir)

    ##############################################
    ## retain only mapped reads (q>1) and sort
    ##############################################
    cmd = 'samtools view -@ ' + str(core_num) + ' -bhq 1 ' + input_bam + ' | ' + \
          'samtools' + ' sort' + \
          ' -@ ' + str(core_num) + \
          ' -o ' + output_dir + '/temp_mapped_sorted.bam'
    print('run:\n %s' % cmd)
    subprocess.call(cmd, shell=True)
    print('samtools filtering mapped reads and sorting finished\n')

    ##############################################
    ## remove PCR duplicates by Picard
    ##############################################
    cmd = 'java -jar $EBROOTPICARD/picard.jar' + \
          ' MarkDuplicates' + \
          ' MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000' + \
          ' REMOVE_DUPLICATES=true' + \
          ' I=' + output_dir + '/temp_mapped_sorted.bam' + \
          ' O=' + output_dir + '/temp_mapped_rmpcr.bam'+ \
          ' METRICS_FILE=' + output_dir + '/dups.txt' + \
          ' BARCODE_TAG=CB' + \
          ' ASSUME_SORT_ORDER=coordinate' + \
          ' USE_JDK_DEFLATER=true' + \
          ' USE_JDK_INFLATER=true'
    print('run:\n %s' % cmd)
    subprocess.call(cmd, shell=True)
    print('Picard finished\n')

    ##############################################
    ##  filter bam MQ > qual and properly paired
    ##############################################
    cmd = 'samtools view -@ ' + str(core_num) + ' -f 3 -bhq ' + qual + ' ' + output_dir + '/temp_mapped_rmpcr.bam > ' + output_dir + '/temp_mq' + qual + '_rmpcr.bam'
    print('run:\n %s' % cmd)
    subprocess.call(cmd, shell=True)
    print('Samtools filtering finished\n')
    ##############################################
    ##    Fix barcodes
    ##############################################
    ###   add sample name to CB barcode
    ##############################################
    cmd = 'python FixingBarcodeName.py -BAM ' + output_dir + '/temp_mq' + qual + '_rmpcr.bam' + ' -exp_name ' + sample_name + ' | ' \
          'samtools view -@ ' + str(core_num) + ' -h - > ' + output_dir + '/temp_fixBC_mq' + qual + '_rmpcr.sam'
    print('run:\n %s' % cmd)
    subprocess.call(cmd, shell=True)
    print('FixingBarcodeName.py barcode filtering finished\n')

    ##############################################
    ##    make Tn5 bed files
    ##############################################
    out_bed = output_dir + '/' + sample_name + 'tn5_mq' + qual + '_processed.bed'
    # cmd = 'perl makeTn5bed.pl ' + output_dir + '/temp_fixBC_mq' + qual + '_rmpcr.bam | sort -k1,1 -k2,2n - | uniq - > ' + out_bed
    # alternative: use python script to generate the bed file:
    cmd = 'python makeTn5bed.py -sam' + output_dir + '/temp_fixBC_mq' + qual + '_rmpcr.sam' + ' -output_file ' + out_bed
    # makeTn5bed.py <sam_file> <output_directory>
    print('run:\n %s' % cmd)
    subprocess.call(cmd, shell=True)
    print('Making bed file finished')
    print('Output bed file: %s \n' % out_bed)
    
    ##############################################
    ##remove temporary files
    ##############################################
    cmd = 'rm temp*'
    print('...deleting temporary files')
    subprocess.call(cmd, shell=True)
    return out_bed      

def remove_blacklist(tn5_bed,output_dir,black_list_fl,ref_fl): 
    '''remove black list regions'''
    ##############################################
    ##    intersect with the bed files
    ##############################################
    cmd = 'bedtools intersect -a ' + tn5_bed + ' -b ' + black_list_fl + \
            ' -wa -wb -sorted -g ' + ref_fl + ' > ' + output_dir + '/temp_intersect.txt'
    print('run:\n %s \n' % cmd)
    subprocess.call(cmd,shell=True)
    print("Finished bedtools intersect, start to write filtered bed file.\n")
    
    store_black_region_dic = {}
    print('intersect regions:')
    with open (output_dir + '/temp_intersect.txt','r') as handle:
        for line in handle:
            col = line.strip().split()
            black_loc = col[0] + '_' + col[1] + '_' + col[2] + '_' + col[3] + '_' + col[4]
            store_black_region_dic[black_loc] = 1
    overlapNum = len(store_black_region_dic)
    print('Number of overlapped regions: %s \n' % overlapNum)
    
    ## name the file output bed file
    bed_file = open(tn5_bed, 'r')
    rmblack_bed = tn5_bed.replace('.bed', '') + '_rmblack.bed' 
    ouput_file = open(rmblack_bed, 'w')

    for ln in bed_file:
        ln = ln.strip('\n')
        col = ln.strip().split()
        loc = col[0] + '_' + col[1] + '_' + col[2] + '_' + col[3] + '_' + col[4]
        if loc not in store_black_region_dic:
            ouput_file.write(ln + '\n')

    bed_file.close()
    ouput_file.close()
    print('Finished bed file filtering')
    print('Final processed bed file: %s' % rmblack_bed)

if __name__=='__main__':  
## read argvs
  argvs = get_parser().parse_args()
  quality = '10'         ## change this based on your data/species
  bam_file = argvs.bam.rstrip('/')  
  sample_name = argvs.sample
  output_dir = argvs.out.rstrip('/')
  # check if output directory exists
  isExist = os.path.exists(output_dir)
  if isExist:
    print("output dir '%s' exists, will overwrite files\n" % output_dir)
  else:
    os.makedirs(output_dir)
  cpu_nm = argvs.cpu
  ref = argvs.ref
  black_list = argvs.bklst  
  print(" \
        The input bam file is:\t %s;\n \
        The sample folder is: %s;\n \
        The output folder is: %s;\n \
        The alignment quality threshold is %s; \n \
        The CPUs called for the jobs is %s; \n \
        The reference index file is %s\n \
        The blacklist file is %s." % (bam_file, sample_name, output_dir, quality, cpu_nm, ref, black_list))

## run the analysis
  tn5_bed = process_bam2bed(bam_file,output_dir,cpu_nm,sample_name,quality)
  remove_blacklist(tn5_bed,output_dir,black_list,ref)