##### INPUT: .bam file to be sorted and output directory to place split BC
import os
import pysam
import time
import pandas as pd 
### Input varibles to set
# file to split on
unsplit_file = "/workspace/luot/data/EV_data/data/cellranger_out/100K5M/outs/possorted_genome_bam.bam"
# where to place output files
out_dir = "/workspace/luot/data/EV_data/data/QC/100K5M/"
barcode_list_ev = pd.read_csv("/workspace/luot/EV/script/QC/barcode_list.txt")["x"].tolist()


barcode_list = []
# variable to hold barcode index
CB_hold = 'unset'
flag_hold = False
itr = 0
# read in upsplit file and loop reads by line
samfile = pysam.AlignmentFile(unsplit_file, "rb")
for read in samfile.fetch(until_eof=True):
    # barcode itr for current read
    try:
        CB_itr = read.get_tag('CB')
    except:
        continue
    if CB_itr not in barcode_list_ev:
        continue
    # if change in barcode or first line; open new file  
    if( CB_itr!=CB_hold or itr==0 ):
        if CB_itr in barcode_list:
            file_name = CB_itr + "_tmp"
            flag_merge = True
        else:
            file_name = CB_itr
            barcode_list.append(CB_itr)
            flag_merge = False
        # close previous split file, only if not first read in file
        if(itr!=0):
            split_file.close()
            if (flag_hold):
                tmp_name = CB_hold + "_new"
                rm_name = CB_hold + "_tmp"
                os.system("samtools merge " + out_dir + "{}.bam".format(tmp_name) + " " + out_dir + "{}.bam".format(CB_hold)+" "+ out_dir + "{}.bam".format(rm_name))
                os.system("rm -rf "+ out_dir + "{}.bam".format(CB_hold))
                os.system("rm -rf "+ out_dir + "{}.bam".format(rm_name))
                os.system("mv "+out_dir + "{}.bam".format(tmp_name) + " " + out_dir + "{}.bam".format(CB_hold))
        CB_hold = CB_itr
        flag_hold = flag_merge
        itr+=1
        split_file = pysam.AlignmentFile(out_dir + "{}.bam".format(file_name), "wb", template=samfile)
    # write read with same barcode to file
    split_file.write(read)
split_file.close()
samfile.close()

