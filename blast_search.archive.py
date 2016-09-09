from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import csv
import os
import pdb
import urllib2
import time

E_THRESH = 10**(-20)

#fasta_filename1 = "/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta.part1"
#fasta_filename2 = "/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta.part2"
#fasta_filename3 = "/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta.part3"
#fasta_filename4 = "/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta.part4"
#fasta_filename5 = "/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta.part5"
#fasta_filename6 = "/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta.part6"
#fasta_filename7 = "/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta.part7"
#fasta_filename8 = "/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta.part8"
#fasta_filename9 = "/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta.part9"
#xml_filename1 = fasta_filename1 + ".xml"
#xml_filename2 = fasta_filename2 + ".xml"
#xml_filename3 = fasta_filename3 + ".xml"
#xml_filename4 = fasta_filename4 + ".xml"
#xml_filename5 = fasta_filename5 + ".xml"
#xml_filename6 = fasta_filename6 + ".xml"
#xml_filename7 = fasta_filename7 + ".xml"
#xml_filename8 = fasta_filename8 + ".xml"
#xml_filename9 = fasta_filename9 + ".xml"

def SplitLargeFasta(fasta_filename, size_lim):
	in_file = open(fasta_filename, 'r')
	line_num = 1
	part_num = 1
	line_lim = size_lim * 2 + 1
	out_file = open(fasta_filename + ".part" + str(part_num), 'w')
	for line in in_file:
		if line_num == line_lim:
			out_file.close()
			line_num = 1
			part_num += 1
			out_file = open(fasta_filename + ".part" + str(part_num), 'w')
		out_file.write(line)
		line_num += 1
	out_file.close()

#SplitLargeFasta("/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta", 1000)

def BlastFastaXml(fasta_filename):
	fasta_string = open(fasta_filename).read()
	result_handle = NCBIWWW.qblast("blastn", "nr", fasta_string)

	xml_filename = fasta_filename + ".xml"
	save_file = open(xml_filename, 'w')
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()

def WriteARow(WriterObj, blastRec):
	for hsp in blastRec.alignments[0].hsps:
		query_algn_len = hsp.align_length - hsp.query.count('-')
		coverage = float(query_algn_len) / blastRec.query_length
		if (hsp.frame[1] == 1):
			sbjct_finish = hsp.sbjct_start + (hsp.align_length - hsp.sbjct.count('-')) - 1
		else:
			sbjct_finish = hsp.sbjct_start - (hsp.align_length - hsp.sbjct.count('-')) + 1
		query_finish = hsp.query_start + query_algn_len - 1
		if hsp.expect < E_THRESH:
			WriterObj.writerow([
				blastRec.query, 
				hsp.query_start, query_finish, 
				hsp.sbjct_start, sbjct_finish, 
				coverage])

def CompareQryAndSbjctGap(blastRec):
	hspList = blastRec.alignments[0].hsps
	hspStartEnds = []
	# data stored in form [[query_start, query_finish, sbjct_start, sbjct_finish]...]
	#pdb.set_trace()
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
		if abs(sbjct_gap - query_gap) > 100:
			return True
		else:
			return False
		
def ParseBlastXmlV2(xml_filename):
	result_handle = open(xml_filename)
	blast_records = NCBIXML.parse(result_handle)
	output_table_filename = xml_filename + ".v2.tsv"

	output_table = open(output_table_filename, 'w')
	csvWriter = csv.writer(output_table, delimiter = "\t")
	csvWriter.writerow(["Query_ID", 
						"Query_Start", "Query_Finish", 
						"Subject_Start", "Subject_Finish", 
						"Coverage", "Sequence #"])

	i = 0
	for blast_record in blast_records:
		i += 1
		if blast_record.alignments and ("Dengue" in blast_record.alignments[0].title):
			if len(blast_record.alignments[0].hsps) > 1:
				hspStartEndList = []
				row_written = False
				frame = blast_record.alignments[0].hsps[0].frame[1]
				for hsp in blast_record.alignments[0].hsps[1:]:
					# check if all aren't +/+ or +/-
					if hsp.frame[1] != frame:
						WriteARow(csvWriter, blast_record, i)
						row_written = True
						break
				if row_written == False:
					if CompareQryAndSbjctGap(blast_record):
						WriteARow(csvWriter, blast_record, i)
			else:
				coverage = float(query_algn_len) / blast_record.query_length
				if coverage < 0.95 and blast_record.descriptions[0].e < .001:
					WriteARow(csvWriter, blast_record, i)
	output_table.close()

def ParseBlastXmlV1(xml_filename):
	result_handle = open(xml_filename)
	blast_records = NCBIXML.parse(result_handle)
	output_table_filename = xml_filename + ".v1.tsv"

	output_table = open(output_table_filename, 'w')
	csvWriter = csv.writer(output_table, delimiter = "\t")
	csvWriter.writerow(["Query_ID", 
						"Query_Start", "Query_Finish", 
						"Subject_Start", "Subject_Finish", 
						"Coverage", "Sequence #"])

	i = 0
	for blast_record in blast_records:
		i += 1
		if blast_record.alignments and ("Dengue" in blast_record.alignments[0].title):
			for hsp in blast_record.alignments[0].hsps:
				query_algn_len = hsp.align_length - hsp.query.count('-')
				coverage = float(query_algn_len) / blast_record.query_length
				if blast_record.descriptions[0].e < .001 and coverage < 0.95:
					if (hsp.frame[1] == 1):
						sbjct_finish = hsp.sbjct_start + (hsp.align_length - hsp.sbjct.count('-')) - 1
					else:
						sbjct_finish = hsp.sbjct_start - (hsp.align_length - hsp.sbjct.count('-')) + 1
					query_finish = hsp.query_start + query_algn_len - 1
					csvWriter.writerow([
						blast_record.query, 
						hsp.query_start, query_finish, 
						hsp.sbjct_start, sbjct_finish, 
						coverage, i])

	output_table.close()

def FindInsertions(xml_filename):
	result_handle = open(xml_filename)
	blast_records = NCBIXML.parse(result_handle)

	i = 0
	for blast_record in blast_records:
		#print float(hsp.align_length - hsp.query.count('-')) / blast_record.query_length
		#for alignment in blast_record.alignments:
		i += 1
		if blast_record.alignments and ("Dengue" in blast_record.alignments[0].title):
			if len(blast_record.alignments[0].hsps) > 1:
				print blast_record.query
				print i

def FilterIndvXml(xml_filename, output_table_filename):
	result_handle = open(xml_filename)
	blast_record = NCBIXML.read(result_handle)
	result_handle.close()

	output_table = open(output_table_filename, 'a')
	csvWriter = csv.writer(output_table, delimiter = "\t")
	#pdb.set_trace()

	row_written = False

	if blast_record.alignments and ("Dengue" in blast_record.alignments[0].title):
		if len(blast_record.alignments[0].hsps) > 1:
			frame = blast_record.alignments[0].hsps[0].frame[1]
			for hsp in blast_record.alignments[0].hsps[1:]:
				# check if all aren't +/+ or +/-
				if hsp.frame[1] != frame:
					WriteARow(csvWriter, blast_record)
					row_written = True
					break
			if row_written == False:
				if CompareQryAndSbjctGap(blast_record):
					WriteARow(csvWriter, blast_record)
					row_written = True
		else:
			query_algn_len = blast_record.alignments[0].hsps[0].align_length - blast_record.alignments[0].hsps[0].query.count('-')
			coverage = float(query_algn_len) / blast_record.query_length
			if coverage < 0.85 and blast_record.descriptions[0].e < .0001:
				WriteARow(csvWriter, blast_record)
				row_written = True
		if not row_written:
			os.remove(xml_filename)
	else:
		os.remove(xml_filename)
	output_table.close()

def FilterLargeXml(xml_filename, output_table_filename):
	result_handle = open(xml_filename)
	blast_records = NCBIXML.parse(result_handle)

	output_table = open(output_table_filename, 'w')
	csvWriter = csv.writer(output_table, delimiter = "\t")

	try:
		for blast_record in blast_records:
			#pdb.set_trace()
			row_written = False
			if blast_record.alignments:
				if len(blast_record.alignments[0].hsps) > 1:
					for hsp in blast_record.alignments[0].hsps:
						if hsp.expect > E_THRESH:
							blast_record.alignments[0].hsps.remove(hsp)
					#pdb.set_trace()
					if len(blast_record.alignments[0].hsps) <= 1:
						continue
					frame = blast_record.alignments[0].hsps[0].frame[1]
					for hsp in blast_record.alignments[0].hsps[1:]:
						# check if all aren't +/+ or +/-
						if hsp.frame[1] != frame:
							row_written = True
							break
					if row_written == False:
						if CompareQryAndSbjctGap(blast_record):
							WriteARow(csvWriter, blast_record)
				else:
					query_algn_len = blast_record.alignments[0].hsps[0].align_length - blast_record.alignments[0].hsps[0].query.count('-')
					coverage = float(query_algn_len) / blast_record.query_length
					if coverage < 0.85 and blast_record.descriptions[0].e < E_THRESH:
						WriteARow(csvWriter, blast_record)
	except:
		pass
	result_handle.close()
	output_table.close()
FilterLargeXml("/Users/cdisuchicago/Downloads/WXWNZ5KE01R-Alignment.xml", "/Users/cdisuchicago/Documents/DENV_Minion_fasta/testcase.filtered_summary")
FilterLargeXml("/Users/cdisuchicago/Downloads/WXXV4Y9E01R-Alignment.xml", "/Users/cdisuchicago/Documents/DENV_Minion_fasta/testcase2.filtered_summary")
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

#FilterIndvXml("gi|699980880|gb|KM204118.1| Dengue virus 2 strain New Guinea C, complete genom.v2.xml", "test.fasta.tsv")
#BlastFastaXmlIndv("test.fasta", '.')
#BlastFastaXmlIndv("/Users/cdisuchicago/Documents/DENV_Minion_fasta/DENV2_Ron_AllFail_LowQual_gt2000.fasta.part7", "/Users/cdisuchicago/Documents/DENV_Minion_fasta/xmls", "316f3c37-faa2-402c-a487-f711728a83a6_Basecall_2D_template")
'''
base_folder = "/Users/cdisuchicago/Documents/DENV_Minion_fasta"
for dirpath, dirnames, filenames in os.walk(base_folder):
	for f in filenames:
		if "part" in f:
			BlastFastaXml(os.path.join(base_folder, f))
			ParseBlastXmlV2(os.path.join(base_filder, f + ".xml"))
			'''
