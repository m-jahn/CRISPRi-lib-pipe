#!usr/bin/ncbi-blast-2.6.0+
#This script intends to map reads from deep sequencing to our sgRNA library. This script is run in sgRNA_Library/seqdata/scripts. KS

import sys
import datetime

# This is the main function. This function goes through all reads and map them back to the reference using the local_blast function. 
def map_reads(library,reads_file):

	perfect_counts = {} #100% match no embiguity		24 + 20 (max 25) + 27(out of 30) = 71  76 max seq data
	read_number = 0

	fh = open(reads_file)

	for line in fh.readlines():
		if line[0] == 'C':
			read_number +=1
			#print '** Read'

			line=line.strip()
			if 'N' in line: continue #No embiguities
			#print line[0:70]
			
			alignment_results = align(line,library) #We will modify the length of the line with the length of ths sgRNA are at in the loop

			for alignment in alignment_results:
				#print alignment
				if alignment not in perfect_counts.keys(): perfect_counts[alignment] =1
				else : perfect_counts[alignment] += 1

	fh.close()

	return perfect_counts

# Align to every entry in ths sgRNA library
def align(query_string,library):

	alignment_results=[]


	for entry in library.keys():
		if len(library[entry]) >= len(query_string):
			if query_string in library[entry]: alignment_results.append(str(entry))
		else: #Case where some bp non specific to library in sequencing query
			if library[entry] == query_string[0:(len(library[entry]))]: alignment_results.append(str(entry))

	return alignment_results

def print_results(alignment_results, file_to_write):
	fh = open(file_to_write,'w')
	now = datetime.datetime.now()
	fh.write('Library sequencing results - KS %s \n ******************\n Displaying perfect match reads \n ****************** \n'%(str(now)))

	for alignment in alignment_results:
		fh.write('>'+str(alignment)+'\t'+str(alignment_results[alignment])+'\n')

	fh.close()
	return 0

# Read library
def read_library(fasta_file):

    entry = {}
    fh = open(fasta_file)
    for line in fh.readlines():
        if line[0] == '>':
            ID = line.split('|')[0]+'_'+line.split('|')[1]
            entry[ID.split('>')[1]] = ''

        else:
            entry[ID.split('>')[1]] += line.strip()

    fh.close()
    return entry

def main():

	sample_name = str(sys.argv[1])
	library_file = str(sys.argv[2])


	seq_data_file = 'data/filtered/' + sample_name
	counts_file = 'data/counts/' +sample_name.split('.fastq')[0]+'.txt' 

	library=read_library(library_file)
	alignment_results=map_reads(library,seq_data_file)
	print_results(alignment_results,counts_file)

	return 0

main()
