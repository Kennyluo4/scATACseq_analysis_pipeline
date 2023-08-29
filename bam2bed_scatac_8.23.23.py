#!/usr/bin/env python
''' edited by Ziliang 2023.03 
    change: make the output directory if it doesn't exist.
            merge remove blacklist step
  '''
## this script is to process bam files from the cell ranger output
import sys
import subprocess
import os

def help():
  print("Use: python processbam_fromcellranger_033023.py <cellranger_dir> <sample_name> <output_dir> <cpu#> <ref_index> <black_list>")
  " ... make sure the following scripts are the same folder: makeTn5bed.pl; fixBC.pl"
  " ... make sure picard and SAMtools are loaded"

def process_bam2bed(CRanger_dir, output_dir, core_num, sample_name, qual):
    raw_bam_fl = CRanger_dir + '/' + sample_name + '/outs/possorted_bam.bam'
    print('Located bamfile: %s' % raw_bam_fl)
    
    #make output sample folder in the output directory 
    output_dir = output_dir.rstrip('/') + '/' + sample_name
    os.mkdir(output_dir)

    ##retain only mapped reads (q>1)
    cmd = 'samtools view -@ ' + str(core_num) + ' -bhq 1 ' + raw_bam_fl + ' > ' + output_dir + '/temp_mapped.bam'
    print('run:\n %s' % cmd)
    subprocess.call(cmd, shell=True)
    print('samtools filtering mapped reads finished\n')

    ##sort the bam
    cmd = 'samtools' + \
          ' sort ' + output_dir + '/temp_mapped.bam'+ \
          ' -@ ' + core_num + \
          ' -o ' + output_dir + '/temp_mapped_sorted.bam'
    print('run:\n %s' % cmd)
    subprocess.call(cmd, shell=True)
    print('Samtools sorting finished\n')

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

    ##filter bam MQ > qual and properly paired
    cmd = 'samtools view -@ ' + core_num + ' -f 3 -bhq ' + qual + ' ' + output_dir + '/temp_mapped_rmpcr.bam > ' + output_dir + '/temp_mq' + qual + '_rmpcr.bam'
    print('run:\n %s' % cmd)
    subprocess.call(cmd, shell=True)
    print('Samtools filtering finished\n')

    ##filter barcodes
    cmd = 'perl fixBC.pl ' + output_dir + '/temp_mq' + qual + '_rmpcr.bam | samtools view -bhS - > ' + output_dir + '/temp_fixBC_mq' + qual + '_rmpcr.bam'
    print('run:\n %s' % cmd)
    subprocess.call(cmd, shell=True)
    print('fixBC.pl barcode filtering finished\n')

#    ###!(alternative) using the makeTn5bed.py to make bed file
#     # bam to sam (!!Only required if using makeTn5bed.py:!)(without header)
#     cmd = "samtools view -o " +  output_dir + '/temp_fixBC_mq' + qual + '_rmpcr.sam\t' + output_dir+ '/temp_fixBC_mq' + qual + '_rmpcr.bam'
#     print('run:\n %s\n' % cmd)
#     subprocess.call(cmd, shell=True)

    ##make Tn5 bed files
    out_bed = output_dir + '/' + sample_name + 'tn5_mq' + qual + '_processed.bed'
    cmd = 'perl makeTn5bed.pl ' + output_dir + '/temp_fixBC_mq' + qual + '_rmpcr.bam | sort -k1,1 -k2,2n - | uniq - > ' + out_bed
    # alternative: use python script to generate the bed file:
    # cmd = 'python makeTn5bed.py ' + output_dir + '/temp_fixBC_mq' + qual + '_rmpcr.sam\t' + output_dir
    # makeTn5bed.py <sam_file> <output_directory>
    print('run:\n %s' % cmd)
    subprocess.call(cmd, shell=True)
    print('makeTn5bed.pl making bed file finished')
    print('Output bed file: %s \n' % out_bed)
    
    ##remove temporary files
    cmd = 'rm temp*'
    subprocess.call(cmd, shell=True)
    return out_bed      

def remove_blacklist(tn5_bed,output_dir,black_list_fl,ref_fl):  # input tn5_bed file with path
    print('Removing blacklist regions for: %s'% tn5_bed)
    ##intersect with the bed
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
  try:               ## read argvs
    argvs = sys.argv
    quality = '10'         ## change this based on your data/species
    CRanger_dir = argvs[1].rstrip('/') ##cellranger output directory, includes the cellranger result folders for each sample 
    sample_name = argvs[2].rstrip('/')
    output_dir = argvs[3].rstrip('/')
    isExist = os.path.exists(output_dir)
    if isExist:
      print("output dir '%s' exists, will overwrite files\n" % output_dir)
    else:
      os.makedirs(output_dir)
    cpu_nm = argvs[4]
    ref = argvs[5]
    black_list = argvs[6]  
    print(" \
      The input project directory is:\t %s;\n \
      The sample folder is: %s;\n \
      The output folder is: %s;\n \
      The alignment quality threshold is %s; \n \
      The CPUs called for the jobs is %s; \n \
      The reference index file is %s\n \
      The blacklist file is %s." % (CRanger_dir, sample_name, output_dir, quality, cpu_nm, ref, black_list))
  except IndexError:
    help()
  
  tn5_bed = process_bam2bed(CRanger_dir,output_dir,cpu_nm,sample_name,quality)
  remove_blacklist(tn5_bed,output_dir,black_list,ref)