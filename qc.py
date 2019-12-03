#reimagining fastqc in python
#WH
import pandas
import time
import collections
import queue
from threading import Thread
import multiprocessing
#from sortedcontainers import SortedList
import pybloomfilter
import zlib
import numpy as np
import fileinput
import sys

from collections import defaultdict
# import gzip
# from fastqandfurious import fastqandfurious
import dnaio

config = collections.defaultdict(dict)

def parse_config(config_folder):
    """
    Parses the files in the Configuration folder
    """
    for line in open(config_folder+"/limits.txt"):
        line = line.strip()
        if "#" not in line and len(line) != 0:
            a,b,c=line.strip().split()
            try:
                c = int(c)
            except:
                c = float(c)
            config[a][b] = c

    for line in open(config_folder+"/adapter_list.txt"):
        line = line.strip()
        if "#" not in line and len(line) != 0:
            temp = line.split("\t")
            a = temp[0].strip()
            b = temp[-1].strip()
            config["adapter_list"][a] = b.strip()

    for line in open(config_folder+"/contaminant_list.txt"):
        line = line.strip()
        if "#" not in line and len(line) != 0:
            temp = line.split("\t")
            a = temp[0].strip()
            b = temp[-1].strip()
            config["contaminant_list"][a] = b.strip()
    return config





def do_qc(in_data):
    #print(config)

    #building qual_dict[qual_char] = qual_value to convert qual values to scores
    qd = {}
    quality_char = [i for i in """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""]
    quality_value = range(len(quality_char))
    quality_dict = dict(zip(quality_char, quality_value))

    for i in quality_dict.values():
        qd[i] = 0


    def add_base_qual_dict(header_list, qual_list):
        """
        Propogates the following
        1. tile_sum_dict[tile][pos] = sum quality
        2. tile_count_dict[tile][pos] = counts
        3. qual_dict[pos][quality] = count

        returns a list of [1.,2.,3.]

        """
        qual_dict = collections.defaultdict(dict)
        tile_sum_dict = collections.defaultdict(dict)
        tile_count_dict = collections.defaultdict(dict)

        for tile,qual in zip(tile_list, qual_list):
            pos = 0
            for q in qual:
                if pos in tile_count_dict[tile]:
                    #we have seen this pos before
                    tile_sum_dict[tile][pos] = tile_sum_dict[tile][pos] + q
                    tile_count_dict[tile][pos] = tile_count_dict[tile][pos] + 1
                    qual_dict[pos][q] = qual_dict[pos][q] + 1
                else:
                    #this pos is new
                    tile_sum_dict[tile][pos] = q
                    tile_count_dict[tile][pos] = 1
                    qual_dict[pos] = qd.copy()
                    qual_dict[pos][q] = 1
                pos = pos + 1
        return [qual_dict, tile_sum_dict, tile_count_dict]


    def avg_qual_count(qual_list):
        """
        Propogates avg_qual_count_dict[mean_quality_of_read] = count
        """

        avg_qual_count_dict = {}
        for qual in qual_list:
            q = np.array(qual)
            m = int(np.mean(q))
            if m in avg_qual_count_dict:
                avg_qual_count_dict[m] = avg_qual_count_dict[m] + 1
            else:
                avg_qual_count_dict[m] = 1

        return [avg_qual_count_dict]


    def dedup(seq_list, bloom_start=0, bloom_stop=50):
        """
        Propogates dedup_dict[sequence_key] = count for duplicated sequences.
        sequence_key = seq[bloom_start:bloom_stop]

        """
        dedup_dict = collections.defaultdict(int)

        # dedup_sorted_list = SortedList()

        for seq in seq_list:
            seq_key = seq[bloom_start:bloom_stop]
            if seq_key not in dedup_bloom:
                # the sequence is not a duplicate
                dedup_dict[seq_key] = dedup_dict[seq_key] + 1

        return [dedup_dict]

        # k = zlib.crc32(bseq)
        #k = bseq[0:50]
        #if do_dedup:
        #    #populate dedup_dict:
        #    if k not in dedup_list:
        #        dedup_bloom.add(k)
        #        # dedup_list.append(k)
                ## dedup_original_dict[k] = bseq
        #else:
            ##check for dedup
            #if k not in dedup_bloom:
                #pass
            #else:
                # dedup_dict[k] = dedup_dict[k] + 1
                #dedup_sorted_list.add(k)



    def base_level(seq_listz):
        """
        propogates the <base>_pos_count[pos] = count dicts. 
        and 
        length_dict[length] = count
        """

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

        return [a_pos_count, t_pos_count, g_pos_count, c_pos_count, n_pos_count, length_dict]

    def dummy(seq_list):
        i = 0
        d = {}
        for j in seq_list:
            i = i + 1

        d["count"] = i
        return [d]

    if in_data =="EMPTY":
        print('done')
        #sys.exit(0)
    (seq_list, qual_list, header_list, tile_list) = in_data
    # print("sent")

    return [[add_base_qual_dict(tile_list, qual_list) + avg_qual_count(qual_list) +  base_level(seq_list)+  dedup(seq_list)],[dummy(seq_list)]]



def read_fastq(fn, chunk_size = 10000 ,bloom_limit = 100000, bloom_start = 0, bloom_stop = 50):
    """
    read given fastq.gz file (fn) and put chunks of reads into data_queue as (seq_list, qual_list, header_list)

    fn : file name
    chunk_size : max sequences per chunk in dna_queue
    bloom_limit: sequences needs to be added into the bloom filter

    """

    quality_char = [i for i in """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""]
    quality_value = range(len(quality_char))
    quality_dict = dict(zip(quality_char, quality_value))

    tile_list = []
    seq_list = []
    qual_list = []
    header_list = []
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
            header_list.append(entry.name)
            seq_list.append(seq)
            qual_list.append(qual)


            count = count + 1
            read_count = read_count + 1

            #propogating the bloom filter
            if read_count > bloom_limit:
                #this is the sample used for the dedup
                add_bloom=False
                #resetting dedup_list to free memory
                dedup_list = []


            if add_bloom:
                k = seq[bloom_start:bloom_length]
                if k not in dedup_list:
                    dedup_bloom.add(k)
                    dedup_list.append(k)
            if count == chunk_size:
                # chunk size was hit. Adding the chunk to data_queue
                if seq_list == []:
                    pass
                else:
                    data_queue.put((seq_list, qual_list, header_list, tile_list))
                    #reset chunk count and lists

                    count = 0
                    tile_list = []
                    seq_list = []
                    qual_list = []

    print("File read complete")
    #The EMPTY flag denotes the end of the read.
    data_queue.put("EMPTY")
    #exiting thread
    sys.exit(1)



def join_dict2(one,two, level=0):
    """
    joins two dictionaries with one key
    """
    if level == 0 :
        #autodetermine the type
        if type(one) == dict:
            #standard dict
            level = 1
        elif type(one) == collections.defaultdict:
            if one.default_factory == dict:
                level =2
            else:
                level=1

    if level == 1:
        return pandas.Series(one).add(pandas.Series(two), fill_value=0)
    elif level == 2:
        return collections.defaultdict(dict, pandas.DataFrame(one).add(pandas.DataFrame(two), fill_value=0).to_dict())

def join_dict(one,two, level=0):
    """
    joins two dictionaries with one key
    """
    if level == 0 :
        #autodetermine the type
        if type(one) == dict:
            #standard dict
            level = 1
        elif type(one) == collections.defaultdict:
            if one.default_factory == dict:
                level =2
            else:
                level=1

    if level == 1:
        #standard dict
        for (k,v) in one.items():
            if k in two:
                two[k] = two[k] + v
            else:
                two[k] = v
        return two
    elif level == 2:
        print("length", len(one), len(two))
        for (k1,v1) in one.items():
            for (k2,v2) in v1.items():
                if k2 in two[k1]:
                    two[k1][k2] = two[k1][k2] + v2
                else:
                    two[k1][k2] = v2


                # try:
                #     two[k1][k2] = two[k1][k2] + v2
                # except:
                #     print("KEY error", k1, k2,v2, type(one), one.default_factory, type(two), two.default_factory, v2)
                #     print(one[k1].keys())
                #     print(two[k1].keys())
                #     print(v1)
                #     two[k1][k2] = v2

        #return collections.defaultdict(dict, pandas.DataFrame(one).add(pandas.DataFrame(two), fill_value=0).to_dict())
        print("length2", len(two))
        return two



def merge_results2(m_results, results):


    [m_qual_dict, m_tile_sum_dict, m_tile_count_dict, m_avg_qual_count_dict, m_a_pos_count, m_t_pos_count, m_g_pos_count, m_c_pos_count, m_n_pos_count, m_length_dict] = m_results
    [qual_dict, tile_sum_dict, tile_count_dict, avg_qual_count_dict, a_pos_count, t_pos_count, g_pos_count, c_pos_count, n_pos_count, length_dict] = results


    #[(m_qual_dict, m_tile_sum_dict, m_tile_count_dict), m_avg_qual_count_dict, (m_a_pos_count, m_t_pos_count, m_g_pos_count, m_c_pos_count, m_n_pos_count, m_length_dict)] = m_results
    #[(qual_dict, tile_sum_dict, tile_count_dict), avg_qual_count_dict, (a_pos_count, t_pos_count, g_pos_count, c_pos_count, n_pos_count, length_dict)] = results


    m_qual_dict = join_dict(m_qual_dict, qual_dict)
    m_tile_sum_dict = join_dict(m_tile_sum_dict, tile_sum_dict)

    m_tile_count_dict=join_dict(m_tile_count_dict, tile_count_dict)
    m_avg_qual_count_dict=join_dict(m_avg_qual_count_dict, avg_qual_count_dict)
    m_a_pos_count = join_dict(m_a_pos_count, a_pos_count)
    m_t_pos_count = join_dict(m_t_pos_count, t_pos_count)
    m_g_pos_count = join_dict(m_g_pos_count, g_pos_count)
    m_c_pos_count = join_dict(m_c_pos_count, c_pos_count)
    m_n_pos_count = join_dict(m_n_pos_count, n_pos_count)
    m_length_dict = join_dict(m_length_dict, length_dict)

    return  [m_qual_dict, m_tile_sum_dict, m_tile_count_dict, m_avg_qual_count_dict, m_a_pos_count, m_t_pos_count, m_g_pos_count, m_c_pos_count, m_n_pos_count, m_length_dict]


def merge_results(results_queue):
    """
    collapese the results_queue by combining the dicts.
    """
    global merging_complete
    global final_results
    while True:
        time.sleep(1)
        if results_queue.qsize() > 2:
            #There are mroe than two entries in the results queue.
            rr1 = results_queue.get()
            rr2 = results_queue.get()

            print("merging", len(rr1), len(rr2))
            if rr1 == "EMPTY":
                #you cant have this
                print("rr1 is emptry. ERROR")
                #print(rr2)

            if rr2 == "EMPTY":
                print("empty recieved by merge_results")
                #the do_qc process is complete. 
                #time to wrap up
                final_results = rr1
                print("bye from merge_results")
                merging_complete = True
                sys.exit(1)
            #print(len(rr1), len(rr2))
            [[s_r1], [c_r1]] = rr1
            [[s_r2], [c_r2]] = rr2
            #print(len(rr1), len(rr2), len(s_r1), len(s_r2))
            # [m1_1,m2_1, m3_1, m4_1, c_r1] = rr1
            # [m1_2,m2_2, m3_2, m4_2, c_r2] = rr2





            s_r3 = []
            for (d1,d2) in zip(s_r1, s_r2):
                joined = join_dict(d1,d2)
                if joined == None:
                    print ("ERROR", len(s_r1), len(s_r2), len(d1), len(d2))
                s_r3.append(joined)
            c_r3 = []
            for (d1,d2) in zip(c_r1, c_r2):
                c_r3.append(join_dict(d1,d2))

            results_queue.put([[s_r3],[c_r3]])


def print_stats():

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



if __name__ == "__main__":
    final_results = []
    merging_complete = False
    pool_size = 4
    config = parse_config("Configuration")

    m_qual_dict = collections.defaultdict(dict)
    m_tile_sum_dict = collections.defaultdict(dict)
    m_tile_count_dict = collections.defaultdict(dict)
    m_avg_qual_count_dict = {}
    m_a_pos_count = collections.defaultdict(int)
    m_t_pos_count = collections.defaultdict(int)
    m_g_pos_count = collections.defaultdict(int)
    m_c_pos_count = collections.defaultdict(int)
    m_n_pos_count = collections.defaultdict(int)
    m_length_dict = collections.defaultdict(int)
    dedup_dict = collections.defaultdict(int)
    m_results = [m_qual_dict, m_tile_sum_dict, m_tile_count_dict, m_avg_qual_count_dict, m_a_pos_count, m_t_pos_count, m_g_pos_count, m_c_pos_count,m_n_pos_count, m_length_dict, dedup_dict]


    #todo for dict entries with default int they start with 0 not 1. So make sure you add one when processing them in the end. 

    #Setting up bloom filter
    dedup_bloom = pybloomfilter.BloomFilter(100010, 0.1, "dedup.bloom")

    read_count = 0
    add_bloom = True
    dedup_list = []
    count = 0
    results = []

    #setting up the thread to read fastq
    pool = multiprocessing.Pool(processes=pool_size)

    #data_queue contains the contents of the def read_fastq output.
    data_queue = queue.Queue(maxsize=pool_size)

    #results_que contains the contents of the def do_qc output. 
    results_queue = queue.Queue()

    print("initiating")
    reader = Thread(target=read_fastq, args = (sys.argv[1].strip(),))
    # reader.setDaemon(True) #change this to process if you want
    reader.start()

    merger = Thread(target=merge_results, args = (results_queue,))
    merger.start()

    try_to_stop = False
    data_buffer = []

    while True:
        time.sleep(1)
        if merging_complete:
            #print_stats()
            print(results_queue.get())
            print(results_queue.qsize())

            sys.exit()
        if data_queue.qsize() > 0:
            #there is stuff in the queue
            temp_data = data_queue.get()

            if "EMPTY" in temp_data:
                #this is the last entry of the do_qc thread.
                print("File read is complete")
                #push the rest of the entries to do_qc

                results_n = pool.map(do_qc, data_buffer)
                for r in results_n:
                    results_queue.put(r)
                data_buffer = []

                #push "EMPTY to results_queue
                results_queue.put("EMPTY")
                print("push empty to merge_results")




            else:
                data_buffer.append(temp_data)
                if len(data_buffer) == pool_size:
                    #there are entries equal to the size of the multiprocessing pool.
                    #so push it to pool
                    print("dispatchin do_qc")
                    results_n = pool.map(do_qc, data_buffer)

                    print("do_qc results recieved", len(results_n))
                    #clear data buffer
                    data_buffer = []
                    for r in results_n:
                        results_queue.put(r)







#                if len(data_buffer) > 0:

#                    # results = poolmap(do_qc, data_buffer)
#                    results.append(poolmap(do_qc, data_buffer))

#                    #m_results = merge_results(m_results, results)
#                    sys.exit(0)
#                    break
#                else:
#                    sys.exit(0)
#            else:
#                data_buffer.append(temp_data)
#                if len(data_buffer) == 4:
#                    results_4 = pool.map(do_qc, data_buffer)
#                    for r in results_4:
#                        results_queue.put(r)
#                        #m_results = merge_results(m_results, r)
#                        #results.append(r)
#                    data_buffer = []



    #while True:
    #    if data_queue.qsize() > 0:
    #        #there is stuff in the queue
    #        temp_data = data_queue.get()

    #        if "EMPTY" in temp_data:
    #            print("empty detected in main loop")
    #            if len(data_buffer) > 0:
    #                results = pool.map(do_qc, data_buffer)
    #            break
    #        else:
                #data_buffer.append(temp_data)
                ##wait until there are 4 chunks in the data_buffer
                #if len(data_buffer) == 4:
                #    results = pool.map(do_qc, data_buffer)
                #    data_buffer = []



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






print(len(results))
#sys.exit(1)


