from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import csv
import os
import pdb
import urllib2
import time

E_THRESH = 10**(-20)

def WriteARow(WriterObj, blastRec, gap = False):
	for hsp in blastRec.alignments[0].hsps:
		query_algn_len = hsp.align_length - hsp.query.count('-')
		coverage = float(query_algn_len) / blastRec.query_length
		if (hsp.frame[1] == 1):
			sbjct_finish = hsp.sbjct_start + (hsp.align_length - hsp.sbjct.count('-')) - 1
		else:
			sbjct_finish = hsp.sbjct_start - (hsp.align_length - hsp.sbjct.count('-')) + 1
		query_finish = hsp.query_start + query_algn_len - 1
		if hsp.expect < E_THRESH:
			if gap:
				WriterObj.writerow([
					blastRec.query, 
					hsp.query_start, query_finish, 
					hsp.sbjct_start, sbjct_finish, 
					coverage, gap])
			else:
				WriterObj.writerow([
					blastRec.query, 
					hsp.query_start, query_finish, 
					hsp.sbjct_start, sbjct_finish, 
					coverage])

def CompareQryAndSbjctGap(blastRec):
	hspList = blastRec.alignments[0].hsps
	hspStartEnds = []
	# data stored in form [[query_start, query_finish, sbjct_start, sbjct_finish]...]
	for hsp in hspList:
		query_algn_len = hsp.align_length - hsp.query.count('-')
		if (hsp.frame[1] == 1):
			sbjct_finish = hsp.sbjct_start + (hsp.align_length - hsp.sbjct.count('-')) - 1
		else:
			sbjct_finish = hsp.sbjct_start - (hsp.align_length - hsp.sbjct.count('-')) + 1
		query_finish = hsp.query_start + query_algn_len - 1
		hspStartEnds.append([hsp.query_start, query_finish, hsp.sbjct_start, sbjct_finish])
	for indx in range(len(hspStartEnds)-1):
		hspdata = hspStartEnds[indx]
		nexthspdata = hspStartEnds[indx+1]
		if hspdata[0] < nexthspdata[0]:
			query_gap = nexthspdata[0] - hspdata[1]
		else:
			query_gap = hspdata[0] - nexthspdata[1]
		if nexthspdata[2] > max(hspdata[2:]) and nexthspdata[3] > max(hspdata[2:]):
			sbjct_gap = min(nexthspdata[2:]) - max(hspdata[2:])
		elif nexthspdata[2] < min(hspdata[2:]) and nexthspdata[3] < min(hspdata[2:]):
			sbjct_gap = min(hspdata[2:]) - max(nexthspdata[2:])
		else:
			return True
		if abs(sbjct_gap - query_gap) > 250:
			return abs(sbjct_gap - query_gap)
		else:
			return False

def FilterLargeXml(xml_filename, output_table_filename):
	result_handle = open(xml_filename)
	blast_records = NCBIXML.parse(result_handle)

	output_table = open(output_table_filename, 'w')
	csvWriter = csv.writer(output_table, delimiter = "\t")
	csvWriter.writerow(["Query_ID", 
						"Query_Start", "Query_Finish", 
						"Subject_Start", "Subject_Finish", 
						"Coverage", "Difference in Query Gap and Subject Gap"])

	total_reads = 0
	dengue_aligned = 0

	try:
		for blast_record in blast_records:
			total_reads += 1
			row_written = False
			if blast_record.alignments:
				dengue_aligned += 1
				if len(blast_record.alignments[0].hsps) > 1:
					original_hsps = list(blast_record.alignments[0].hsps)
					for hsp in original_hsps:
						if hsp.expect > E_THRESH:
							blast_record.alignments[0].hsps.remove(hsp)
					if len(blast_record.alignments[0].hsps) <= 1:
						continue
					frame = blast_record.alignments[0].hsps[0].frame[1]
					for hsp in blast_record.alignments[0].hsps[1:]:
						# check if all aren't +/+ or +/-
						if hsp.frame[1] != frame:
							row_written = True
							break
					if row_written == False:
						gapOrNot = CompareQryAndSbjctGap(blast_record)
						if gapOrNot == True:
							WriteARow(csvWriter, blast_record)
						elif gapOrNot != False:
							WriteARow(csvWriter, blast_record, gapOrNot)
				else:
					query_algn_len = blast_record.alignments[0].hsps[0].align_length - blast_record.alignments[0].hsps[0].query.count('-')
					coverage = float(query_algn_len) / blast_record.query_length
					if coverage < 0.85 and blast_record.descriptions[0].e < E_THRESH:
						WriteARow(csvWriter, blast_record)
		csvWriter.writerow(["Total Reads: " + str(total_reads), "Dengue Reads: " + str(dengue_aligned)])
	except:
		pass
	result_handle.close()
	output_table.close()

#FilterLargeXml("/Users/cdisuchicago/Downloads/WXWNZ5KE01R-Alignment.xml", "/Users/cdisuchicago/Documents/DENV_Minion_fasta/testcase.filtered_summary")
#FilterLargeXml("/Users/cdisuchicago/Downloads/WXXV4Y9E01R-Alignment.xml", "/Users/cdisuchicago/Documents/DENV_Minion_fasta/testcase2.filtered_summary")
FilterLargeXml("/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_blastn_DENV2NGS.xml", "/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_blastn_DENV2NGS.xml.filtered_summary")

def BlastFastaXmlIndv(fasta_filename, xml_folder, starting_qname = False, local = False):
	output_table_filename = fasta_filename + ".tsv"

	fasta_string = open(fasta_filename)
	if starting_qname:
		startline = 1
		while fasta_string.readline()[1:-2] != starting_qname:
			startline += 1
		fasta_string.close()
		fasta_string = open(fasta_filename)
		for i in range(1, startline):
			fasta_string.readline()
	else:
		output_table = open(output_table_filename, 'w')
		csvWriter = csv.writer(output_table, delimiter = "\t")
		csvWriter.writerow(["Query_ID", 
							"Query_Start", "Query_Finish", 
							"Subject_Start", "Subject_Finish", 
							"Coverage"])
		output_table.close()
	if local:
		xml_filename = fasta_filename + ".local.xml"
		output_table = open(output_table_filename, 'a')
		blastn_cline = NcbiblastnCommandline(query = fasta_filename, db = local, outfmt = 5, out = xml_filename)
		stdout, stderr = blastn_cline()
		FilterLargeXml(xml_filename, output_table_filename)
	else:
		while True:
			query_name = fasta_string.readline()
			if query_name == '':
				break
			query_seq = fasta_string.readline()
			print query_name
			wait_time = 1
			while True:
				try:
					result_handle = NCBIWWW.qblast("blastn", "nr", query_name + query_seq, entrez_query="\"Dengue virus\"[porgn:__txid12637]")
					break
			
				except(ValueError, urllib.error.URLError, urllib.error.HTTPError):
					print "Error encountered"
					print "Trying again in " + str(wait_time) + " seconds"
					if wait_time > 500:
						break
					time.sleep(wait_time)
					wait_time *= 2
					
			xml_filename = os.path.join(xml_folder, query_name[1:-2] + ".xml")
			save_file = open(xml_filename, 'w')
			save_file.write(result_handle.read())
			save_file.close()
			result_handle.close()

			FilterIndvXml(xml_filename, output_table_filename)

