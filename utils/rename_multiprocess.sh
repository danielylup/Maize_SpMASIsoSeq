#!/bin/bash

# Shell script that enable multi-thread processing 
# To run : ./rename_multiprocess.sh <input.file> <rename_dict.file> <path/to/output.file> <rename.awk> <number of thread>

inputfile=$1
dictfile=$2
outputfile=$3
script=$4
count=$5
split -n l/$count $inputfile /tmp/_pawk$$
for file in /tmp/_pawk$$*; do
            awk -f script.awk $dictfile $file > ${file}.out &
    done
    wait
cat /tmp/_pawk$$*.out > $outputfile
rm /tmp/_pawk$$*
