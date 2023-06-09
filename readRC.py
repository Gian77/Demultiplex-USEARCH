# Copyright © 2023 Gian M.N. Benucci, Ph.D.
# email: benucci@msu.edu URL: https://github.com/Gian77

# For each record, it reverses the sequence using the reverse_complement method 
# of Seq from Biopython. It also reverses the quality scores by using 
# slicing ([::-1]) on the quality list.
#
# Use: pyhton readRC.py input.fastq output.fastq


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

def reverse_complement_fastq(input_file, output_file):
	with open(output_file, 'w') as output_handle:
		for record in SeqIO.parse(input_file, "fastq"):
			seq = record.seq.reverse_complement()
			qual = record.letter_annotations["phred_quality"][::-1]
			reversed_record = SeqRecord(seq, id=record.id, name=record.name, description=record.description)
			reversed_record.letter_annotations["phred_quality"] = qual
			SeqIO.write(reversed_record, output_handle, "fastq")

	print(f"Reversed complement and quality of {input_file} written to {output_file}")


reverse_complement_fastq(input_file, output_file)
