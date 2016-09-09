import csv
import os

def CombineTsv(virus, numberOfParts):
	os.chdir("/Users/cdisuchicago/Documents/" + virus + "_database")

	out_filename = "all_" + virus + "_nucleotide.fasta.compiled.xml.summary.tsv"
	outfile = open(out_filename, 'w')
	outWriter = csv.writer(outfile, delimiter='\t')
	outWriter.writerow(["Query_ID", "Query_Start", "Query_Finish", 
						"Subject_Start", "Subject_Finish"])
	for part in range(1, numberOfParts+1):
		part_filename = "all_" + virus + "_nucleotide.fasta.part" + str(part) + ".xml.summary.tsv"
		partfile = open(part_filename, 'r')
		partReader = csv.reader(partfile, delimiter = '\t')
		for row in partReader:
			outWriter.writerow(row)
		partfile.close()
	outfile.close()

CombineTsv("west_nile", 21)
