#!/usr/bin/env python

##updating 050321 allow the flag_list is not empty
##updating this script for the specific ct genome
##tbis script transfers perl script from Alex to python format
import re
import sys
import numpy as np

input_sam_fl = sys.argv[1]
output_dir = sys.argv[2]

def makeTn5bed (input_sam_fl, output_dir):

    ##transfer the
    ##prepare the flag_dic
    flag_dic = {'0':'read_paired',
                '1':'read_properly',
                '2':'read_umapped',
                '3':'mate_unmapped',
                '4':'read_reverse',
                '5':'mate_reverse',
                '6':'first_pair',
                '7':'second_pair',
                '8':'secondary',
                '9':'fail_qc',
                '10':'duplicate',
                '11':'sup_align'
    }

    sam_fl = open(input_sam_fl,'r')
    bed_fl = open(output_dir + '/opt_Tn5bed.txt','w')

    #store_final_line_list = []

    count = 0
    #with open (input_sam_fl, 'r') as ipt:
    for eachline in sam_fl:
        eachline = eachline.strip('\n')
        col = eachline.strip().split()

        count += 1
        if (count/1000000) == 0:
            print('- iterated over $counts reads ... \n')

        ##extracct barcode information
        first_col = col[0].split(':')
        bc = 'CB:Z:' + first_col[0]
        #for i in range(9,len(col)):
        #    if col[i].startswith('CB:Z:'):
        #        bc = col[i]

        ##store flag list information
        flag_list = [] ##[read_reverse,second_pair,duplicate,sup_align]

        ##transfer the flag to binary information
        bin = np.binary_repr(int(col[1]), width=12)

        bin_list = list(bin)
        ##generate flag information
        for i in range(len(bin_list)):
            if bin_list[i] == '1':
                flag_list.append(flag_dic[str(i)])

        ##get start and end pos of read
        chr = col[2]
        pos1 = col[3]
        cigar = col[5]
        netdif = 1

        ##in order to make act list we need to do transfer the CIGAR string to another format:
        ##dic_list is eg. [{'M': 76}, {'I': 15}, {'M': 57}, {'S': 3}]
        list_CIGAR = re.findall('\d+|\D+', cigar)
        n = int(len(list_CIGAR) / 2)
        dic_list = []
        for i in range(0, int(n)):
            dic = {list_CIGAR[(2 * i + 1)]: int(list_CIGAR[(2 * i)])}
            dic_list.append(dic)


        cig_dic = {} ####cig_dic stores key and value. key is the 1_M and value is number besides the key in the CIGAR eg {'1_M':50,'2_S':30}
        act_list = [] ##transfer the cigar to [1_M,2_S] or others [1_S,2_M,3_S] stored in the act_list
        ##generate act_list
        item_count = 0
        for eachdic in dic_list:
            item_count += 1
            item_str = str(item_count) + '_' + list(eachdic.keys())[0]
            act_list.append(item_str)
            cig_dic[item_str] = str(eachdic[list(eachdic.keys())[0]])


        for eachact in act_list:
            ##eachact is 1_M or 2_S or others
            ##eachact_list = [1,M] or [2,S]
            ##eachact_list[1] is 'M'
            eachact_list = eachact.split('_')

            ##save the time value to the cig dictionary
            ##time indicates the value in the CIGAR eg. 50M. time is '50'
            time = cig_dic[eachact]

            if eachact_list[1] == 'M' or eachact_list[1] == 'S':
                netdif += int(time)
            elif eachact_list[1] == 'D':
                netdif += int(time)

        pos2 = netdif + int(pos1)

        ##shift
        ##updating 050321
        if flag_list != []:
            
            if flag_list[0] == 'read_reverse':
                end = pos2 - 4
                start = end - 1
                #print (chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '-')
                final_line = chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '-'
                bed_fl.write(final_line + '\n')
            else:
                start = int(pos1) + 5
                end = start + 1
                #print (chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '+')
                final_line = chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '+'
                bed_fl.write(final_line + '\n')

        #with open (output_dir + '/opt_Tn5bed.txt','w+') as opt:
        #    for eachline in store_final_line_list:
        #        opt.write(eachline + '\n')

makeTn5bed (input_sam_fl, output_dir)
