# Copyright © 2023 Gian M.N. Benucci, Ph.D.
# email: benucci@msu.edu URL: https://github.com/Gian77

# Header example: @M03127:552:000000000-C9R4D:1:1101:19227:1880;sample=R18_T1R2FBR6Root;

# Usage example
# python splitFastqReads.py <pathg-to-dir>/file.fastq 

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os

input_file = sys.argv[1]
output_directory = sys.argv[2]
read = sys.argv[3]

def split_fastq_by_sample(input_file, output_directory, read):
	records = SeqIO.parse(input_file, 'fastq')
	sample_files = {}
	sample_counts = {}

	# Create the output directory if it doesn't exist
	os.makedirs(output_directory, exist_ok=True)

	for record in records:
		header_parts = record.description.split(';')
		for part in header_parts:
			if part.strip().startswith('sample='):
				sample_name = part.strip().split('=')[1]
				if sample_name not in sample_files:
					sample_filename = os.path.join(output_directory, f'{sample_name}_{read}.fastq')
					sample_files[sample_name] = open(sample_filename, 'w')
					sample_counts[sample_name] = 1		

				# Handle different quality encoding formats
				if 'solexa_quality' in record.letter_annotations:
					qualities = record.letter_annotations['solexa_quality']
				elif 'phred_quality' in record.letter_annotations:
					qualities = record.letter_annotations['phred_quality']
				else:
					qualities = None

				# Create a new SeqRecord with the updated sample name and quality scores
				sequence_number = sample_counts[sample_name]
				renamed_record = SeqRecord(record.seq, id=f'{sample_name}.{sequence_number}', name='', description='')
				if qualities:
					renamed_record.letter_annotations['phred_quality'] = qualities

				SeqIO.write(renamed_record, sample_files[sample_name], 'fastq')

				# Increment the sequence count for the sample
				sample_counts[sample_name] += 1

				break

	# Close all the sample files
	for file in sample_files.values():
		file.close()

# Split the FASTQ file by sample
split_fastq_by_sample(input_file, output_directory, read)
