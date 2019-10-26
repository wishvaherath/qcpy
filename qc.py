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
tile_dict = collections.defaultdict(dict)

#d['dict1']['innerkey'] = 'value'

from collections import defaultdict
# qual_dict = defaultdict(list)
#setting up the qual dict

quality_char = [ord(i)-33 for i in """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""]
quality_value = range(len(quality_char))
quality_dict = dict(zip(quality_char, quality_value))

qd = {}
for i in quality_dict.values():
    qd[i] = 0

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

def add_base_qual_dict3(qual):
    pos =0

    for q in qual:
        if pos in qual_dict:
            if q in qual_dict[pos]:
                qual_dict[pos][q] = qual_dict[pos][q] + 1
            else:
                qual_dict[pos][q] = 1
        else:
            qual_dict[pos][q] = 1
        pos = pos + 1

def add_base_qual_dict(tile, qual):
    pos=0
    for q in qual:
        if pos in qual_dict:
            qual_dict[pos][q] = qual_dict[pos][q] + 1
        else:
            qual_dict[pos] = qd.copy()
            qual_dict[pos][q] = 1
        pos = pos + 1


with gzip.open(sys.argv[1].strip()) as fh:
    it = fastqandfurious.readfastq_iter(fh, bufsize, fastqandfurious.entryfunc)
    for entry in it:
        header = str(entry[0])
        #extract tile numer
        tile = header.strip().split(":")[4].strip()
        print(tile)
        seq = [i for i in entry[1]]
        qual = [quality_dict[i] for i in entry[2]]

        add_base_qual_dict(tile, qual)



print("per base sequence quality")
print("Base   Mean    Median  Lower Quartile  Upper Quartile  10th Percentile 90th Percentile")
"""
qual_freq_arr = list(base_qual_dict.values())
cs = np.cumsum(qual_freq_arr)
Qn = np.searchsorted(cs, np.percentile(cs, 75)
"""
for p in qual_dict.keys():
    # print(p)
    quality_array = np.array(list(qual_dict[p].values()))
    # print(quality_array)
    cs = np.cumsum(quality_array)
    # print(cs)
    q_list = []
    for i in [10, 25, 50, 75, 90]:
        q_list.append(np.searchsorted(cs, np.percentile(cs,i)))

    # print(p, "->", q_list)
    #calculate mean
    quality_key = np.array(list(qual_dict[p].keys()))
    total = np.sum(quality_key*quality_array)
    quality_mean = total / np.sum(quality_array)

    print(p, quality_mean, q_list[0], q_list[1], q_list[2], q_list[3])



