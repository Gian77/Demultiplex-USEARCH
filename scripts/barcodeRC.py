# Copyright © 2023 Gian M.N. Benucci, Ph.D.
# email: benucci@msu.edu URL: https://github.com/Gian77

# For each record, it reverses the sequence using the reverse_complement method 
# of Seq from Biopython. 
# The os.path.split shoudl make it so it writes the output in the same directory as
# the imput file.

# Use: pyhton barcodeRC.py input.fasta

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os

input_file = sys.argv[1]
output_file = os.path.splitext(input_file)[0] + "_rc.fasta"

def reverse_complement_fasta(input_file, output_file):
	with open(input_file, 'r') as input_handle, open(output_file, 'w') as output_handle:
		for record in SeqIO.parse(input_handle, 'fasta'):
			
			reverse_complement_seq = record.seq.reverse_complement()
			
			reverse_complement_record = SeqRecord(reverse_complement_seq, id=record.id,
				name=record.name, description=record.description)

			SeqIO.write(reverse_complement_record, output_handle, 'fasta')

reverse_complement_fasta(input_file, output_file)
