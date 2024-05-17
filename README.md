# KerTAD

KerTAD is an algorithm for identifying topologically associating domains in Hi-C maps. There are two versions of the code: the 'KerTAD' folder is written in Matlab and 'KerTAD python' is a python version. 

KerTAD takes as input a symmetric NxN matrix of positive values (not necessarily integers) and returns a list of TAD locations by their corner points in the upper triangular matrix. The Matlab version indices begin at 1 and the Python version indices begin at 0. 

### Matlab
The Matlab code requires the 'Image Processing Toolbox' and the 'Statistics and Machine Learning Toolbox' to run.
Download the 'KerTAD' folder and change current folder to the location of the 'KerTAD' folder. Then, given an NxN matrix (input_map)
run:
```
[output_map, TAD_list] =callTADs(input_map);
```
The outputs are:

-output_map: a binary NxN matrix where 1s represent corner points of TADs in the upper triangular matrix

-TAD_list: a two column list of indices for each corner TAD point in the upper triangular matrix


### Python
The Python code requires the numpy and scipy packages. Download the 'KerTAD python' folder and in terminal (assuming you've changed to the appropriate directory) run: 
```
python3 callTADs.py input output
```
input refers to the input matrix (ideally in tab-delimited .txt format) and output refers to output file. For example:

```
python3 callTADs.py  'chr5.txt' 'output.txt'
```

By default, the output is a two column list of indices for each corner TAD point in the upper triangular matrix. 


For the Hi-C maps used in the manuscript see:
https://figshare.com/s/c92bd17f5bd0882fa3e0
