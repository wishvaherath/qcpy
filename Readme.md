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


