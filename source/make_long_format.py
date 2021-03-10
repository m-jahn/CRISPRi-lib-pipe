# This script assembles read counts into a single table in a csv format.
# Author: Kiyan Shabestary
# Date: 2021-03-10

import os
import sys
import re

def main():

	outfile = sys.argv[1]+sys.argv[2]

	fh=open(outfile,'w')
	fh.write('sgRNA\tsampleID\tcounts\n')

	#Loop through all sample files and append all counts from all samples to a master count file
	dirname=sys.argv[1]
	counts_files = [f for f in os.listdir(dirname) if re.match('.*counts.txt$', f)]

	for count_file in counts_files:
		if count_file[0]=='.': continue 
		ifh=open(os.path.join(dirname,count_file))

		ifh.readline()
		ifh.readline()
		ifh.readline()
		ifh.readline()

		for line in ifh.readlines():
			if line[0] == '>': fh.write(line[1:].strip()+'\t')
			else: fh.write(count_file.split('_filtered_counts.txt')[0]+'\t'+line[1:])

		ifh.close()
	fh.close()
	return 0

main()
