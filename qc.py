#reimagining fastqc in python
#WH
import time
import queue
from threading import Thread
import multiprocessing
from sortedcontainers import SortedList
import pybloomfilter
import zlib
import numpy as np
import fileinput
import sys
# import gzip
# from fastqandfurious import fastqandfurious
import dnaio

bufsize = 20000
import collections


qual_dict = collections.defaultdict(dict)
tile_sum_dict = collections.defaultdict(dict)
tile_count_dict = collections.defaultdict(dict)

#d['dict1']['innerkey'] = 'value'

from collections import defaultdict
# qual_dict = defaultdict(list)
#setting up the qual dict


# quality_char = [ord(i)-33 for i in """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""]
quality_char = [i for i in """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""]
quality_value = range(len(quality_char))
quality_dict = dict(zip(quality_char, quality_value))


avg_qual_count_dict = {}

dedup_list = []
dedup_dict = collections.defaultdict(int)
dedup_original_dict={}

dedup_sorted_list = SortedList()


a_pos_count = collections.defaultdict(int)
t_pos_count = collections.defaultdict(int)
g_pos_count = collections.defaultdict(int)
c_pos_count = collections.defaultdict(int)
n_pos_count = collections.defaultdict(int)

length_dict = collections.defaultdict(int)




qd = {}
for i in quality_dict.values():
    qd[i] = 0

#pupulating the dict with zeros
#for i in quality_char:
#    base_qual_dict[i] = 0

#use dict.get(k, default)

def do_qc(in_data):
    print("called")
    def add_base_qual_dict(tile_listz, qual_listz):
        qual_dict = collections.defaultdict(dict)
        tile_sum_dict = collections.defaultdict(dict)
        tile_count_dict = collections.defaultdict(dict)

        for tile,qual in zip(tile_listz, qual_listz):
            pos=0
            for q in qual:
                if pos in tile_count_dict[tile]:
                    #if pos is in the tile dict then IS in the qual_dict
                    tile_sum_dict[tile][pos] = tile_sum_dict[tile][pos] + q
                    tile_count_dict[tile][pos] = tile_count_dict[tile][pos] + 1
                    qual_dict[pos][q] = qual_dict[pos][q] + 1
                else:
                    tile_sum_dict[tile][pos] = q
                    tile_count_dict[tile][pos] = 1
                    qual_dict[pos] = qd.copy()
                    qual_dict[pos][q] = 1
                pos = pos + 1

        return (qual_dict, tile_sum_dict, tile_count_dict)


    def avg_qual_count(qual_listz):

        avg_qual_count_dict = {}
        for qual in qual_listz:
            q = np.array(qual)
            m = int(np.mean(q))
            if m in avg_qual_count_dict:
                avg_qual_count_dict[m] = avg_qual_count_dict[m] + 1
            else:
                avg_qual_count_dict[m] = 1

        return avg_qual_count_dict


    def dedup(bseq):
        dedup_list = []
        dedup_dict = collections.defaultdict(int)
        dedup_original_dict={}

        dedup_sorted_list = SortedList()


        # k = zlib.crc32(bseq)
        k = bseq[0:50]
        if do_dedup:
            #populate dedup_dict:
            if k not in dedup_list:
                dedup_bloom.add(k)
                # dedup_list.append(k)
                # dedup_original_dict[k] = bseq
        else:
            #check for dedup
            if k not in dedup_bloom:
                pass
            else:
                # dedup_dict[k] = dedup_dict[k] + 1
                dedup_sorted_list.add(k)



    def base_level(seq_listz):

        a_pos_count = collections.defaultdict(int)
        t_pos_count = collections.defaultdict(int)
        g_pos_count = collections.defaultdict(int)
        c_pos_count = collections.defaultdict(int)
        n_pos_count = collections.defaultdict(int)

        length_dict = collections.defaultdict(int)

        for seq in seq_listz:
            l = len(seq)
            length_dict[l] = length_dict[l] + 1
            pos = 0
            for b in seq:
                # b = chr(b)
                if b == 'A':
                    a_pos_count[pos] = a_pos_count[pos] + 1
                elif b == 'T':
                    t_pos_count[pos] = t_pos_count[pos] + 1
                elif b == 'G':
                    g_pos_count[pos] = g_pos_count[pos] + 1
                elif b == 'C':
                    c_pos_count[pos] = c_pos_count[pos] + 1
                elif b == 'N':
                    n_pos_count[pos] = n_pos_count[pos] + 1
                pos = pos + 1

        return (a_pos_count, t_pos_count, g_pos_count, c_pos_count, n_pos_count, length_dict)

    if in_data =="EMPTY":
        print('done')
        #sys.exit(0)
    (seq_listz, qual_listz, tile_listz) = in_data
    print("sent")
    return [add_base_qual_dict(tile_listz, qual_listz), avg_qual_count(qual_listz), base_level(seq_listz)]







def read_fastq(fn):
    tile_list = []
    seq_list = []
    qual_list = []
    count = 0   
    read_count = 0
    add_bloom = False


    with dnaio.open(fn) as fh:
        for entry in fh:
            header = entry.name
            #extract tile numer
            tile = header.strip().split(":")[4].strip()
            seq = entry.sequence
            #convert qual from char to number
            qual = [quality_dict[i] for i in entry.qualities]

            tile_list.append(tile)
            seq_list.append(seq)
            qual_list.append(qual)


            count = count + 1
            read_count = read_count + 1
            if read_count > 10000:
                #this is the sample used for the dedup
                add_bloom=False

            if add_bloom:
                k = seq[0:50]
                if k not in dedup_list:
                    dedup_bloom.add(k)
                    dedup_list.append(k)
            if count == 10000:
                print("dispatching")
                data_queue.put((seq_list, qual_list, tile_list))

                # add_base_qual_dict(tile_list, qual_list)
                # avg_qual_count(qual_list)
                # base_level(seq_list)
                # dedup(seq_list)

                count = 0
                tile_list = []
                seq_list = []
                qual_list = []

    print("sending empty")
    data_queue.put("EMPTY")

    #push the final dataset
    sys.exit(1)




if __name__ == "__main__":
    #todo for dict entries with default int they start with 0 not 1. So make sure you add one when processing them in the end. 

    dedup_bloom = pybloomfilter.BloomFilter(100010, 0.1, "dedup.bloom")
    read_count = 0
    add_bloom = True
    dedup_list = []
    count = 0

    # tile_list = []
    # seq_list = []
    # qual_list = []

    pool = multiprocessing.Pool(processes=8)
    data_queue = queue.Queue()


    print("initiating")
    reader = Thread(target=read_fastq, args = (sys.argv[1].strip(),))
    # reader.setDaemon(True) #change this to process if you want
    reader.start()
    try_to_stop = False
    data_buffer = []
    while True:
        # time.sleep(1)
        if data_queue.qsize() > 0:
            #there is stuff in the queue
            temp_data = data_queue.get()

            if "EMPTY" in temp_data:
                print("empty detected in main loop")
                if len(data_buffer) > 0:
                    results = poolmap(do_qc, data_buffer)
                break
            else:
                data_buffer.append(temp_data)

                if len(data_buffer) == 4:
                    results = pool.map(do_qc, data_buffer)
                    data_buffer = []



        #if data_queue.qsize() >= 4:

        #    # add_base_qual_dict(tile_list, qual_list)
        #    # avg_qual_count(qual_list)
        #    # base_level(seq_list)
        #    # dedup(seq_list)
        #    temp_data = (data_queue.get(),data_queue.get(),data_queue.get(),data_queue.get())
        #    if "EMPTY" in temp_data:
        #        try_to_stop = True
        #        print("empty detected1")
        #        data_queue.join()
        #    results = pool.map(do_qc, temp_data)
        #    #reader.join()
        #elif try_to_stop:
        #        print("empty detected2222")
        #        data_queue.join()
        #        sys.exit(1)
        #        # results = pool.map(do_qc, temp_data)
        #        #reader.join()









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


print("per tile sequence quality")
print("Tile     Base    Mean")

for tile in tile_count_dict:
    for pos in tile_count_dict[tile]:
        m = tile_sum_dict[tile][pos] / tile_count_dict[tile][pos]
        print(tile, pos, m)

print("Per sequence quality scores")
print("Quality  Count")
for q in sorted(avg_qual_count_dict):
    print(q, avg_qual_count_dict[q])

print("Per base sequence content")
print("Base G   A   T   C")

for p in sorted(a_pos_count):
    a = a_pos_count[p]
    t = t_pos_count[p]
    g = g_pos_count[p]
    c = c_pos_count[p]
    n = n_pos_count[p]

    ff = (a+t+g+c+n) / 100

    print(p, g/ff, a/ff, t/ff, c/ff)

print("Per sequence GC content")
print("GC Content   count")
for p in sorted(a_pos_count):
    a = a_pos_count[p]
    t = t_pos_count[p]
    g = g_pos_count[p]
    c = c_pos_count[p]
    n = n_pos_count[p]

    gc = (g+c)/(a+t) * 100
    print(p, gc)

print("Per base N content")
print("Base N-Count")
for p in sorted(n_pos_count):
    print(p, n_pos_count[p])

print("Sequence Length Distribution")
print("Length   Count")

for l in sorted(length_dict):
    print(l, length_dict[l])

