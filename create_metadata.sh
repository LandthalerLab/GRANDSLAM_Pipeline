#!/bin/bash

path=$(echo $1)
samples=$(ls -1 "$path"/*gz | xargs -n 1 basename | awk -F '_R(1|2).' '{ print $1 }' | uniq)

echo "sample read1 read2" >> metadata.tsv
for i in $(echo $samples)
do
         # uncomment below
	 echo $i
	 reads=$(ls -1 "$path/$i"* | xargs -n 1 basename)
	 echo $i $reads >> metadata.tsv
done
