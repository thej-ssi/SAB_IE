#!/usr/bin/env python	
"""
Goes through all fasta files that match the pattern, converts them to a dictionary and counts the occurances of each kmer for each pattern
input:
	./kmer_states.py -i <dir> -p <pattern1>,<pattern2>,...
output:
	dict file with kmer, pattern1 count, pattern2, count
	file with files in pattern1, files in pattern2, unused files

Citation:
Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3 http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878
"""

__author__ = "Kim Lee Ng"
__copyright__ = "Copyright 2017, Statens Serum Institut"
__credits__ = ["Kim Lee Ng"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Kim Lee Ng"
__email__ = "kimn@ssi.dk"
__status__ = "Prototype"

import os
import sys
import argparse
import re
from Bio.Seq import Seq

def program_initialization():	
	#Handle arguements
	parser = argparse.ArgumentParser(description='Kmer states')		
	#will not require any inputs other than config for batch version
	parser.add_argument("-i","--input_directory",
						help="Give a working directory",
						required=True)
	parser.add_argument("-o","--output_directory",
						help="Give a output directory",
						required=True)
	parser.add_argument("-p","--patterns", 
						help="State patterns (pattern1,pattern2,...)",
						type=str,
						required=True)
	parser.add_argument("-k","--kmer_size", 
						help="kmer_size",
						type=int,
	                    default=33)
	args = parser.parse_args()
	args.input_directory = args.input_directory.strip()
	args.output_directory = args.output_directory.strip()

	if(args.patterns != ""):
		args.patterns = args.patterns.split(",")

	return args

def add_to_kmer_dict(kmer_dict,k,fasta_text,number_of_states,pattern_number,sample_dict):	
	for j in range(len(fasta_text)-k+1):
		reverse_complement = str(Seq(fasta_text[j:j+k]).reverse_complement())

		if fasta_text[j:j+k] in kmer_dict or reverse_complement in kmer_dict:			
			if not (fasta_text[j:j+k] in sample_dict or reverse_complement in sample_dict):
				sample_dict[fasta_text[j:j+k]] = True
				if reverse_complement in kmer_dict:
					kmer_dict[reverse_complement][pattern_number]+=1
				else:
					kmer_dict[fasta_text[j:j+k]][pattern_number] +=1
		else:
			sample_dict[fasta_text[j:j+k]] = True
			kmer_dict[fasta_text[j:j+k]] = [0]*number_of_states
			kmer_dict[fasta_text[j:j+k]][pattern_number] +=1
	return ""

def main(argv):
	args = program_initialization()	

	pattern_list = []
	pattern_file_list = []
	for i in range(len(args.patterns)):
		pattern_file_list.append([])
	unused_file_list = []	
	kmer_dict = {}	

	for pattern in args.patterns:
		pattern_list.append(re.compile(pattern))

	for read_file in sorted(os.listdir(args.input_directory)):
		any_match = False
		for i, pattern in enumerate(pattern_list):			
			result = re.match(pattern,read_file)
			if(result):				
				print("{}\t{}\t{}").format(i, read_file, pattern.pattern) 
				pattern_file_list[i].append(read_file)				
				any_match = True
				with open(os.path.join(args.input_directory,read_file),"r") as fasta_file:
					fasta_text=""
					sample_dict={}
					for line in fasta_file:
						if line[0]==">":
							if fasta_text != "":
								fasta_text=add_to_kmer_dict(kmer_dict,args.kmer_size,fasta_text,len(args.patterns),i,sample_dict)							
						else:
							fasta_text+=line.strip()

					if fasta_text !="":
						fasta_text=add_to_kmer_dict(kmer_dict,args.kmer_size,fasta_text,len(args.patterns),i,sample_dict)

		if any_match == False:
			unused_file_list.append(read_file)
	
	with open(os.path.join(args.output_directory,"kmer_state_samples_list.txt"),"w+") as kmer_state_samples_list_file:
		for i in range(len(pattern_file_list)):
			kmer_state_samples_list_file.write("Group: {}\n".format(i))
			for j in pattern_file_list[i]:
				kmer_state_samples_list_file.write("\t{}\n".format(j))
		kmer_state_samples_list_file.write("Group: unused\n".format())
		for j in unused_file_list:			
			kmer_state_samples_list_file.write("\t{}\n".format(j))

	with open(os.path.join(args.output_directory,"kmer_state_dict.txt"),"w+") as kmer_state_dict:
		for key in kmer_dict:
			kmer_state_dict.write("{}\t{}\n".format(key,"\t".join(str(x) for x in kmer_dict[key])))
	print("Done!")

if __name__ == "__main__":   
	main(sys.argv)   