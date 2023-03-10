### Python 3.6.13
import os
import sys
import pysam
import argparse
import pandas as pd
from tqdm import trange
from pathlib import Path

def get_readcount(bam_file):
	"""
	Parses a bam file idxstats to get the number of reads.
	The BAM file MUST have an index.
	"""
	num_reads = pysam.idxstats(bam_file).split('\n')
	nums = {}
	for num in num_reads:
		try:
			chrom, chrlen, mapped, unmapped = num.split('\t')
			nums[chrom] = int(mapped) + int(unmapped)
		except ValueError:
			print(num)
	return pd.DataFrame(nums, index=['num']).T.sum().values[0]

def main(input_bam,barcodes_file,output_dir):
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	if Path(barcodes_file).is_file():
		fh = open(barcodes_file, "r")
		barcode_list = [l.rstrip() for l in fh.readlines()]
		fh.close()
	else:
		barcode_list = [barcodes_file]
	barcode_dict = {}
	# read in upsplit file and loop reads by line
	# For progress bar.
	progress = trange(get_readcount(input_bam))  # the total number of reads in the bam file
	samfile = pysam.AlignmentFile(input_bam, "rb")
	for read in samfile:
		progress.update(1)
		if (not read.is_unmapped) and (not read.is_secondary) and (not read.is_duplicate):  # get only primary mapped reads.
			try:
				CB_itr = read.get_tag('CB')
			except:
				continue
			if CB_itr not in barcode_list:
				continue
			if CB_itr not in barcode_dict.keys():
				split_file = pysam.AlignmentFile(os.path.join(output_dir + "/{}.bam".format(CB_itr)), "wb", template=samfile)
				barcode_dict[CB_itr] = split_file
			barcode_dict[CB_itr].write(read)
	samfile.close()
	print("close samfile")
	for key in barcode_dict.keys():
		barcode_dict[key].close()
	print("job done")
        
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--bam",  required=True, help="BAM file")
	parser.add_argument("--list", required=True, help="Barcode list (Without header)")
	parser.add_argument("--out",  required=True, help="Output folder")
	args = parser.parse_args()
	main(input_bam=args.bam,barcodes_file=args.list,output_dir=args.out)

# ulimit -HSn 1024000 && nohup python split_bc_bam_memory_cost.py /workspace/luot/data/EV_data/data/cellranger_out/100K5M/outs/possorted_genome_bam.bam /workspace/luot/EV/QC/mapping_rates/K562_100_EV_barcode_list_cb2_20220413.txt /workspace/luot/EV/QC/mapping_rates/K562_100_EV > K562_100_bam_split.log 2>&1 &
