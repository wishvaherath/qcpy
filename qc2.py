from qc import *
#reimagining fastqc in python
#WH

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



iterate_fastq(sys.argv[1].strip())
print_report()



