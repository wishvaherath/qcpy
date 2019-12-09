# from extern import *
import inspect
import extern

# ereimagining fastqc in python
#WH
cimport cython
import operator
# from sortedcontainers import SortedList
import pybloomfilter
import zlib
import numpy as np
import fileinput
import sys
import gzip
from fastqandfurious import fastqandfurious
from collections import defaultdict
import collections
from collections import Counter


output_file = open("test.txt","w")

bufsize = 200000


qual_dict = collections.defaultdict(dict)
tile_count_dict = collections.defaultdict(dict)

tile_avg_dict = {}
tile_sum_dict = collections.defaultdict(int)
cdef int old_tile = 0

cdef int tile_read_count = 0

avg_qual_count_dict = {}
a_pos_count = collections.defaultdict(int)
t_pos_count = collections.defaultdict(int)
g_pos_count = collections.defaultdict(int)
c_pos_count = collections.defaultdict(int)
n_pos_count = collections.defaultdict(int)
length_dict = collections.defaultdict(int)
cdef list dedup_list = []
dedup_dict = collections.defaultdict(int)
# cdef dict dedup_dict
# dedup_sorted_list = SortedList()
dedup_bloom = pybloomfilter.BloomFilter(100010, 0.1, "dedup.bloom")
# read_count = 0
cdef int GC_count = 0
cdef int AT_count = 0 
seq_gc = collections.defaultdict(int)

cdef int read_count = 0
do_dedup = True


#setting up the qual dict
# quality_char = [ord(i)-33 for i in """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcpcpdefghijklmnopqrstuvwxyz{|}~"""]
# quality_value = range(len(quality_char))
# quality_dict = dict(zip(quality_char, quality_value))
# qd = {}
# for i in quality_dict.values():
#     qd[i] = 0

Q_symbol_ascii = [ord(i) for i in """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcpcpdefghijklmnopqrstuvwxyz{|}~"""]
Q_score = [i-33 for i in Q_symbol_ascii]


quality_dict = dict(zip(Q_symbol_ascii, Q_score))
qd = {}
for i in quality_dict.values():
     qd[i] = 0



@cython.profile(True)
cpdef iterate_fastq(fn, pipe):


    global read_count
    global do_dedup

    cdef str seq
    cdef list qual
    with gzip.open(fn) as fh:
        it = fastqandfurious.readfastq_iter(fh, bufsize, fastqandfurious.entryfunc)
        for entry in it:
            if read_count % 100000 == 0:
                print("Read count : ", read_count)
            read_count = read_count + 1
            if read_count >100000:
                do_dedup=False


            # print([i for i in entry[2]])
            header = str(entry[0])
            #extract tile numer
            tile = header.strip().split(":")[4].strip()
            # print(tile)
            seq = "".join([chr(i) for i in entry[1]])
            qual = [quality_dict[i] for i in entry[2]]

            add_base_qual_dict(tile, qual)
            avg_qual_count(qual)
            base_level(seq)
            dedup(entry[1],do_dedup)

            #processing the functions in extern
            extern_functions = inspect.getmembers(extern, inspect.isfunction) 
            for (f_name, f) in extern_functions:
                if "report_" not in f_name:
                    f(entry[1], qual, header)
            if pipe:
                print(str(header))
                print(str(seq))
                print(entry[2])





@cython.profile(True)
cpdef add_base_qual_dict(tile, qual):
    """
    popula tes three dicts
    tile_sum_dict[tile][pos]= sum quality
    tile_count_dict[tile][pos] = counts
    qual_dict[pos][quality] = count
    """
    tile = int(tile)
    # print(int(tile))
    global old_tile
    global tile_read_count
    global tile_sum_dict
    global qual_dict
    global tile_avg_dict
    if old_tile == 0:
        tile_read_count = 0
        #this is the fist iteration
        old_tile = tile


    tile_read_count = tile_read_count + 1
    #check the tile has changed
    if old_tile != tile:

        #convert sums to mean
        for p in tile_sum_dict:
            tile_sum_dict[p] = tile_sum_dict[p] / tile_read_count
            #add the avg quals to the qual_dict_large
        tile_avg_dict[old_tile] = tile_sum_dict.copy()
        tile_sum_dict = collections.defaultdict(int)
        #the tile has changed
        old_tile = int(tile)

    else:
        #same tile as old, bau
        pos = 0
        for q in qual:
            qual_dict[pos][q] = qual_dict[pos].get(q,0) + 1

            tile_sum_dict[pos] = tile_sum_dict[pos] + q
            pos = pos + 1





@cython.profile(True)
cpdef avg_qual_count(qual):
    """
    populates avg_qual_count_dict[mean_quality_of_read] = count

    """
    global avg_qual_count_dict
    # q = np.array(qual)
    # m = int(np.average(q, axis=0))
    cdef int m = int(sum(qual)/len(qual))
    if m in avg_qual_count_dict:
        avg_qual_count_dict[m] = avg_qual_count_dict[m] + 1
    else:
        avg_qual_count_dict[m] = 1


@cython.profile(True)
cpdef dedup(bseq, do_dedup):
    """
    populates dedup_dict[sequence_substring] = count
    """

    # k = zlib.crc32(bseq)
    global dedup_dict
    global dedup_bloom

    cdef str k = str(zlib.crc32(bseq[0:50]))
    # cdef str k = str(zlib.crc32(bseq[0:50]))
    if do_dedup:
        #populate dedup_dict:

        dedup_bloom.add(k)
        #dedup_list.append(k)
        # dedup_original_dict[k] = bseq
    else:
        #check for dedup
        if k not in dedup_bloom:
            pass
        else:
            # dedup_dict[k] = dedup_dict.get(k,0) + 1
            dedup_dict[k] = dedup_dict[k] + 1




@cython.profile(True)
cpdef dedup2(bseq, do_dedup):
    """
    populates dedup_dict[sequence_substring] = count
    """
    # k = zlib.crc32(bseq)
    

    #cdef str k = str(bseq[0:50])

    cdef str k = str(zlib.crc32(bseq[0:50]))
    if do_dedup:
        #populate dedup_dict:
        if k not in dedup_list:

            dedup_bloom.add(k)
            dedup_list.append(k)
            # dedup_original_dict[k] = bseq
    else:
        #check for dedup
        if k not in dedup_bloom:
            pass
        else:
            # dedup_dict[k] = dedup_dict.get(k,0) + 1
            dedup_dict[k] = dedup_dict[k] + 1


@cython.profile(True)
cpdef base_level(seq):
    global GC_count
    global AT_count
    cdef int l = len(seq)
    length_dict[l] = length_dict[l] + 1
    cdef int pos = 0
    global a_pos_count
    global t_pos_count
    global g_pos_count
    global c_pos_count
    old_AT_count = AT_count
    old_GC_count = GC_count
    for b in seq:
        if b == 'A':
            a_pos_count[pos] = a_pos_count[pos] + 1
            AT_count = AT_count + 1
        elif b == 'T':
            t_pos_count[pos] = t_pos_count[pos] + 1
            AT_count = AT_count + 1
        elif b == 'G':
            g_pos_count[pos] = g_pos_count[pos] + 1
            GC_count = GC_count + 1
        elif b == 'C':
            c_pos_count[pos] = c_pos_count[pos] + 1
            GC_count = GC_count + 1
        elif b == 'N':
            n_pos_count[pos] = n_pos_count[pos] + 1
        pos = pos + 1

    #calculating the GC% for this seq
    gc = GC_count - old_GC_count
    at = AT_count - old_AT_count
    this_seq_gc = int((gc / (at + gc)) * 100 )
    seq_gc[this_seq_gc] = seq_gc[this_seq_gc] + 1

#todo for dict entries with cpdefault int they start with 0 not 1. So make sure you add one when processing them in the end. 


@cython.profile(True)
cpdef write(txt):
    s = ""
    if type(txt) == list:
        for t in txt:
            s = s + str(t) + "\t"
    else:
        s = str(txt)
    output_file.write(s.strip()+ "\n")


cpdef print_basic_stats(read_count):
    write(["##qc.py", "0.1beta"])
    write([">>Basic Statistics", "pass"])
    write(["#Measure","Value"])
    write(["Filename", sys.argv[1].strip()])
    write(["File type", "?"])
    write(["Encoding", "?"])
    write(["Total Sequences",read_count])
    write(["Sequences flagged as poor quality", 0]) #[TODO]
    min_length = str(min(length_dict.keys()))
    max_length = str(max(length_dict.keys()))
    write(["Sequence length", min_length+"-"+max_length])
    gc_percentage = str((GC_count / (AT_count + GC_count))*100)
    write(["%GC", gc_percentage ])
    write(">>END_MODULE")



cpdef print_per_base_sequence_quality():
    write([">>Per base sequence quality", "pass"])
    write("Base   Mean    Median  Lower Quartile  Upper Quartile  10th Percentile 90th Percentile")
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

        write([p, quality_mean, q_list[0], q_list[1], q_list[2], q_list[3]])
    write(">>END_MODULE")


cpdef print_per_tile_sequence_quality():
    write([">>Per tile sequence quality", "pass"])
    write("Tile     Base    Mean")

    # print(tile_avg_dict)

    for tile in tile_avg_dict:
        for pos in tile_avg_dict[tile]:
            m = tile_avg_dict[tile][pos]
            # m = tile_sum_dict[tile][pos] / tile_count_dict[tile][pos]
            write([tile, pos, m])
    write(">>END_MODULE")


cpdef print_per_sequence_quality_scores():
    write(">>Per sequence quality scores")
    write(["#Quality","Count"])
    for q in sorted(avg_qual_count_dict):
        write([q, avg_qual_count_dict[q]])
    write(">>END_MODULE")



cpdef print_per_base_sequence_content():
    write([">>Per base sequence content", "pass"])
    write("#Base G   A   T   C")

    for p in sorted(a_pos_count):
        a = a_pos_count[p]
        t = t_pos_count[p]
        g = g_pos_count[p]
        c = c_pos_count[p]
        n = n_pos_count[p]

        ff = (a+t+g+c+n) / 100

        write([p, g/ff, a/ff, t/ff, c/ff])
    write(">>END_MODULE")

cpdef print_per_sequence_gc_content():
    write([">>Per sequence GC content", "pass"])
    write("#GC Content   count")
    gc_list = list(seq_gc.keys())
    gc_list.sort()
    for g in gc_list:
        write([g, seq_gc[g]])
    write(">>END_MODULE")


cpdef print_per_base_n_content(read_count):
    write([">>Per base N content", "pass"])
    write("#Base N-Count")
    for p in sorted(n_pos_count):
        write([p, (n_pos_count[p]/read_count)])
    write(">>END_MODULE")


cpdef print_sequence_length_distribution():
    write([">>Sequence Length Distribution", "pass"])
    write("#Length   Count")
    for l in sorted(length_dict):
        write([l, length_dict[l]])
    write(">>END_MODULE")

cpdef print_sequence_duplication_levels(read_count):
    dedup_counts = Counter(dedup_dict.values())
    total_dedup = sum(dedup_counts.values())

    print("dedup", dedup_counts, total_dedup)

    write([">>Sequence Duplication Levels", "pass"])
    write(["Total Deduplicated Percentage", "0"])
    write(["#Duplication Level", "Percentage of deduplicated", "Percentage of total"])
    c = list(dedup_counts.keys())
    c.sort()
    print(total_dedup, read_count)
    for i in c:
        write([i, dedup_counts[i]/total_dedup * 100, dedup_counts[i]/read_count*100])

    write(">>END_MODULE")




cpdef print_report():
    global read_count
    print_basic_stats(read_count)
    print_per_base_sequence_quality()
    print_per_tile_sequence_quality()
    print_per_sequence_quality_scores()
    print_per_base_sequence_content()
    print_per_sequence_gc_content()
    print_per_base_n_content(read_count)
    print_sequence_length_distribution()
    print_sequence_duplication_levels(read_count)






