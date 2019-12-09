import argparse
from qc import *
#reimagining fastqc in python
#WH
import inspect
from sortedcontainers import SortedList
import pybloomfilter
import zlib
import numpy as np
import fileinput
import sys
import gzip
from fastqandfurious import fastqandfurious
import collections
# read_count = 0
#bufsize = 10000
#do_dedup = True

parser = argparse.ArgumentParser(description="qc.py - fastQC inspired, experimental, quality reporting tool")
parser.add_argument("-f", action = "store", dest = "fastq_file", help = "fastq.gz file to be processed")
parser.add_argument("-p", action="store_true", default=False,dest = "pipe", help = "pipe fastq.gz data to stdout")

args = parser.parse_args()


iterate_fastq(args.fastq_file, args.pipe)
print_report()

#custom functions
extern_functions =inspect.getmembers(extern, inspect.isfunction) 

data_pack = [qual_dict, tile_avg_dict, avg_qual_count_dict, a_pos_count, t_pos_count, g_pos_count, c_pos_count, length_dict, dedup_dict]
for (f_name, f) in extern_functions:
    if "report_" in f_name:
        f(output_file, data_pack)




