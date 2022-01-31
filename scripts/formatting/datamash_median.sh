#!/usr/bin/bash

module load datamash

while read x
do    tr -s '\t' '\n' <<< "$x" | \
      datamash  median 1
done <input.txt &> output.txt
