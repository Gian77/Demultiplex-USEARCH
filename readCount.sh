#!/bin/bash

# Copyright © 2023 Gian M.N. Benucci, Ph.D.
# email: benucci@msu.edu URL: https://github.com/Gian77

# This script will loop throug a series of fastq files and count the number of sequences present
# in the fiels, then sum all the counts into a total count value for all files.

# Usage example
#sh readCount.sh <input-dir-path> R1

if test -f Sample_R1.counts | test -f Sample_R2.counts ; then 
	rm Sample_R1.counts; rm Sample_R2.counts;
	echo -e "\noverwriting files!\n"
fi

input_dir=$1
read=$2
total_counts=0

# Loop through each FASTQ file
for file in "$input_dir"/*_"$read".fastq; do
	line_count=$(wc -l < "$file")
	sequence_count=$((line_count / 4))

	echo "$file:" "$sequence_count" >> "$input_dir"/Sample_"$read".counts

	# Generate the total counts
	total_counts=$((total_counts + sequence_count))
done

echo "Total:" "$total_counts" >> "$input_dir"/Sample_"$read".counts
