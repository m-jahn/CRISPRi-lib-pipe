# This script assembles counts into a single table in a csv format. KS

import os
import sys

def main():

	outfile = sys.argv[1]+'.txt'

	fh=open(outfile,'w')
	fh.write('sgRNA\tsampleID\tcounts\n')

	#List through all samples counts within the counts folder and append all counts from all samples to a master count file
	#for sample in count_folder:
	dirname="data/counts"
	counts_files = os.listdir(dirname)

	for count_file in counts_files:
		if count_file[0]=='.': continue 
		ifh=open(os.path.join(dirname,count_file))

		ifh.readline()
		ifh.readline()
		ifh.readline()
		ifh.readline()

		for line in ifh.readlines():
			if line[0] == '>': fh.write(line[1:].strip()+'\t')
			else: fh.write(count_file.split('.txt')[0]+'\t'+line[1:]+'\n')

		ifh.close()



	fh.close()
	return 0


main()