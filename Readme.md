# qc.py - a fastqc inspired, **experimental**, quality control tool written in python. 

### Dependencies
qc.py was written in Python 3.7 and converted using cython 0.29. 

It requires the following modules. 
* pybloomfilter
* fastqandfurious
* numpy
* cython


### Running
To start a run clone this repository.
Then do the following 
`python setup.py build_ext --inplace`

Then run 
`python qc2.py -f your_file.fastq.gz`

This will create the **QCreport.txt** file. 

### Adding custom functions

qc.py will run all functions in extern.py. For each feature you need, create two functions. 

~~~python
def *any_name*(seq, qual, header):
	<code here>

~~~

To summerize the output and save it into the QC report use

~~~python
report_any_name(report, data_pack):
	<code here>

~~~

*report* is a file handle for the report. Use it to write to it. 
*data_pack* is a list defined as following.

~~~python
 32 data_pack = [qual_dict, tile_avg_dict, avg_qual_count_dict, a_pos_count, t_pos_count, g_pos_count, c_pos_count, length_dict, dedup_dict]



~~~




