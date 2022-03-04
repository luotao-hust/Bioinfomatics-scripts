### Python 3.6.13
import pysam
from pathlib import Path
import sys
import os

def main(input_bam,barcodes_file,output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if Path(barcodes_file).is_file():
        with open(barcodes_file, "r") as fh:
            barcode_list = [l.rstrip() for l in fh.readlines()]
    else:
        barcode_list = [barcodes_file]
    barcode_dict = {}
    # read in upsplit file and loop reads by line
    samfile = pysam.AlignmentFile(input_bam, "rb")
    for read in samfile.fetch(until_eof=True):
        try:
            CB_itr = read.get_tag('CB')
        except:
            continue
        if CB_itr not in barcode_list:
            continue
        if CB_itr not in barcode_dict.keys():
            split_file = pysam.AlignmentFile(os.path.join(output_dir + "/{}.bam".format(CB_itr)), "wb", template=samfile,threads=8)
            barcode_dict[CB_itr] = split_file
        barcode_dict[CB_itr].write(read)
    samfile.close()
    for key in barcode_dict.keys():
        barcode_dict[key].close()
        
if __name__ == "__main__":
    main(input_bam=sys.argv[1],barcodes_file=sys.argv[2],output_dir=sys.argv[3])
# python split_bc_bam_memory_cost.py /workspace/luot/data/EV_data/data/cellranger_out/100K5M/outs/possorted_genome_bam.bam /workspace/luot/EV/QC/mapping_rates/K562_100_EV_barcode_list.txt  /workspace/luot/EV/QC/mapping_rates/K562_100_EV
