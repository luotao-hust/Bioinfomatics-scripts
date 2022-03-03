### Python 3.6.13
import os
import pysam
import time
import pandas as pd 
### Input varibles to set
# file to split on
unsplit_file = "/workspace/luot/data/EV_data/data/cellranger_out/100K5M/outs/possorted_genome_bam.bam"
# where to place output files
out_dir = "/workspace/luot/data/EV_data/data/QC/100K5M_memory_cost/"
barcode_list_ev = pd.read_csv("/workspace/luot/EV/script/QC/barcode_list.txt")["x"].tolist()
barcode_dict = {}
# read in upsplit file and loop reads by line
samfile = pysam.AlignmentFile(unsplit_file, "rb")
for read in samfile.fetch(until_eof=True):
    try:
        CB_itr = read.get_tag('CB')
    except:
        continue
    if CB_itr not in barcode_list_ev:
        continue
    if CB_itr not in barcode_dict.keys():
        split_file = pysam.AlignmentFile(out_dir + "{}.bam".format(CB_itr), "wb", template=samfile)
        barcode_dict[CB_itr] = split_file
    barcode_dict[CB_itr].write(read)
for key in barcode_dict.keys():
    barcode_dict[key].close()
samfile.close()