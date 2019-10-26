#reimagining fastqc in python
#WH

import numpy as np
import fileinput
import sys
import gzip
from fastqandfurious import fastqandfurious

bufsize = 20000
import collections
qual_dict = collections.defaultdict(dict)
#d['dict1']['innerkey'] = 'value'

from collections import defaultdict
# qual_dict = defaultdict(list)
#setting up the qual dict

quality_char = [ord(i)-33 for i in """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""]
quality_value = range(len(quality_char))
quality_dict = dict(zip(quality_char, quality_value))

#pupulating the dict with zeros
#for i in quality_char:
#    base_qual_dict[i] = 0

#use dict.get(k, default)

"""
def add_base_qual_dict2(qual):
    pos = 0
    for q in qual:
        counter = dict_pos_qual_count.get(pos, 0)
        counter = counter + 1
        dict_pos_qual_count[pos][

        if dict_pos_qual_count[pos][q]:
            dict_pos_qual_count[pos][q] = dict_pos_qual_count[pos][q] + 1
        else:
            dict_pos_qual_count[pos][q] = 1
        pos = pos + 1

"""


def add_base_qual_dict2(qual):
    pos = 0
    for q in qual:
        # key = str(pos) + "_" + str(q)
        key = (pos*1000) + q
        #counter = qual_dict.get(key,  0)
        #counter = counter + 1
        qual_dict[key] = qual_dict[key] + 1
        pos = pos + 1

def add_base_qual_dict(qual):
    pos =0

    for q in qual:
        if pos in qual_dict:
            if q in qual_dict[pos]:
                qual_dict[pos][q] = qual_dict[pos][q] + 1
            else:
                qual_dict[pos][q] = 1
        else:
            qual_dict[pos][q] = 1

with gzip.open(sys.argv[1].strip()) as fh:
    it = fastqandfurious.readfastq_iter(fh, bufsize, fastqandfurious.entryfunc)
    for entry in it:
        header = entry[0]
        seq = [i for i in entry[1]]
        qual = [quality_dict[i] for i in entry[2]]
        add_base_qual_dict(qual)


#per base sequence quality
#Base   Mean    Median  Lower Quartile  Upp    er Quartile  10th Percentile 90th Percentil    e

"""
qual_freq_arr = list(base_qual_dict.values())
cs = np.cumsum(qual_freq_arr)
Qn = np.searchsorted(cs, np.percentile(cs, 75)
"""



