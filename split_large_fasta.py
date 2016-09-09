from Bio import SeqIO
import sys

program, fasta_filename, size_lim = sys.argv
size_lim = int(size_lim)

part = 1

record_iterator = SeqIO.parse(fasta_filename, "fasta")
while True:
	outfile = open(fasta_filename + ".part" + str(part), 'a')
	for counter in range(size_lim):
		try:
			record = next(record_iterator)
		except StopIteration:
			outfile.close()
			sys.exit()
		SeqIO.write(record, outfile, "fasta")
	outfile.close()
	part += 1
