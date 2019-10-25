#reimagining fastqc in python
#WH



import fileinput
import sys
import gzip
from fastqandfurious import fastqandfurious

# from collections import defaultdict
# base_qual_dict = defaultdict(list)
#setting up the qual dict

quality_char = [ord(i) for i in """ !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""]
quality_value = range(len(quality_char))
quality_dict = dict(zip(quality_char, quality_value))


print(quality_dict)








bufsize = 20000
base_qual_dict = {}
def add_base_qual_dict(qual):
    for q in qual:
        if q in base_qual_dict:
            base_qual_dict[q] = base_qual_dict[q]+ 1
        else:
            base_qual_dict[q] = 1


with gzip.open(sys.argv[1].strip()) as fh:
    it = fastqandfurious.readfastq_iter(fh, bufsize, fastqandfurious.entryfunc)
    for entry in it:
        header = entry[0]
        seq = [i for i in entry[1]]
        qual = [quality_dict[i] for i in entry[2]]
        add_base_qual_dict(qual)







