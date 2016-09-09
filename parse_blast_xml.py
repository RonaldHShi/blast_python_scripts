from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from operator import itemgetter
from pdb import set_trace
import sys
import time
import csv
import urllib2

E_THRESH = 10**(-10)

def WriteARow(writerObj, blast_record, filteredHspStartEnds):
	for hspStartEnd in filteredHspStartEnds:
		hspStartEnd.insert(0, blast_record.query)
		writerObj.writerow(hspStartEnd)

def FilterBlastRecord(blast_record):
	if blast_record.alignments:
		hspStartEnds = []
		# data stored in form [[query_start, query_finish, sbjct_start, sbjct_finish]...]
		hspList = blast_record.alignments[0].hsps
		for hsp in hspList:
			if hsp.frame[1] == -1:
				break
			if hsp.expect < E_THRESH:
				query_algn_len = hsp.align_length - hsp.query.count('-')
				if (hsp.frame[1] == 1):
					sbjct_finish = hsp.sbjct_start + (hsp.align_length - hsp.sbjct.count('-')) - 1
				else:
					sbjct_finish = hsp.sbjct_start - (hsp.align_length - hsp.sbjct.count('-')) + 1
				query_finish = hsp.query_start + query_algn_len - 1
				hspStartEnds.append([hsp.query_start, query_finish, hsp.sbjct_start, sbjct_finish])
		sortedHspList = sorted(hspStartEnds, key = itemgetter(0))
		if len(sortedHspList) > 1:
			filteredHspList = []
			for index in range(len(sortedHspList) - 1):
				currentHsp = sortedHspList[index]
				try:
					nextHsp = sortedHspList[index+1]
				except IndexError:
					break
				if currentHsp[1]-100 < nextHsp[0]:
					filteredHspList.append(currentHsp)
				else:
					currentHsp[1] = nextHsp[1]
					currentHsp[3] = nextHsp[3]
					sortedHspList.remove(nextHsp)
					index -= 1
			try:
				filteredHspList.append(sortedHspList[index+1])
			except IndexError:
				pass
			return filteredHspList
		else:
			return False
	else:
		return False

def CheckPossibleRecomb(sortedHspStartEnds, range_start = None, range_end = None):
	if sortedHspStartEnds:
		for indx in range(len(sortedHspStartEnds)-1):
			hspdata = sortedHspStartEnds[indx]
			nexthspdata = sortedHspStartEnds[indx+1]
			query_gap = nexthspdata[0] - hspdata[1]
			sbjct_gap = nexthspdata[2] - hspdata[3]
			if abs(query_gap - sbjct_gap) < 100:
				return False
			else:
				return True

def BlastFastaXmlIndv(fasta_filename = None, xml_filename = None):
	if fasta_filename:
		record_iterator = SeqIO.parse(fasta_filename, "fasta")
		output_table = open(fasta_filename + ".summary.tsv", 'w')
		outputWriter = csv.writer(output_table, delimiter = "\t")
		for seq_record in record_iterator:
			wait_time = 1
			while True:
				print seq_record.id
				try:
					result_handle = NCBIWWW.qblast("blastn", "nr", seq_record.seq, entrez_query="KM204118.1")
					break
				except ValueError:
					print "Error encountered"
					print "Trying again in " + str(wait_time) + " seconds"
					if wait_time > 100:
						sys.exit()
					time.sleep(wait_time)
					wait_time *= 2

			blast_record = NCBIXML.read(result_handle)
			filteredHspStartEnds = FilterBlastRecord(blast_record)
			if filteredHspStartEnds and CheckPossibleRecomb(filteredHspStartEnds):
				WriteARow(outputWriter, blast_record, filteredHspStartEnds)

			result_handle.close()

	elif xml_filename:
		output_table = open(xml_filename + ".summary.tsv", 'w')
		outputWriter = csv.writer(output_table, delimiter = "\t")
		result_handle = open(xml_filename)
		blast_records = NCBIXML.parse(result_handle)
		for blast_record in blast_records:
			filteredHspStartEnds = FilterBlastRecord(blast_record)
			if filteredHspStartEnds and CheckPossibleRecomb(filteredHspStartEnds):
				WriteARow(outputWriter, blast_record, filteredHspStartEnds)

		result_handle.close()

	output_table.close()

# Searching downloaded data from NCBI
###### DENGUE ######
'''
for i in range(1, 21):
	print i
	BlastFastaXmlIndv(xml_filename = "/Users/cdisuchicago/Documents/dengue_database/all_dengue_nucleotide.fasta.part" + str(i) + ".xml")
	'''

###### WEST NILE ######
for i in range(1,22):
	print i
	BlastFastaXmlIndv(xml_filename = "/Users/cdisuchicago/Documents/west_nile_database/all_west_nile_nucleotide.fasta.part" + str(i) + ".xml")
#######################

###### TICK-BORNE ENCEPHALITIS ######
#BlastFastaXmlIndv(xml_filename = "/Users/cdisuchicago/Documents/tick-borne_encephalitis_database/all_tick-borne_encephalitis_nucleotide.fasta.xml")
#####################################

###### YELLOW FEVER ######
#BlastFastaXmlIndv(xml_filename = "/Users/cdisuchicago/Documents/yellow_fever_database/all_yellow_fever_nucleotide.fasta.xml")
##########################
