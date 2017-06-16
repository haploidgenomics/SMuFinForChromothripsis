#! /usr/bin/env python

import sys

#######################################################################
##
## Han Tan, University of Maine, 2017
##
## This script is used to generate bin files for de novo PRICE assembly
## for breakpoint junctions as a result of chromothripsis
## from SMuFin's somatic_large_SVs output files
##
## PRICE assembly will typically result in a larger junction sequence to
## aid in the design of PCR amplicons
##
## USAGE: junction-bins-fromSMuFin.py somatic_large_SVs_$file output.txt
##
#########################################################################

largeSV = open(sys.argv[1])
o = open(sys.argv[2], 'w')

window = 500 ## Adjust this to required bin size around breakpoint

o.write("Chrom"+'\t'+"BinStart"+'\t'+"BinEnd"+'\n')

for line in largeSV:
	j = line.split('\t')
	if j[0] != "Mut_ID":
		chr1 = j[2]
		pos1 = j[3]
		chr2 = j[4]
		pos2 = j[5]
		o.write(chr1+'\t'+str(int(pos1)-int(window))+'\t'+str(int(pos1)+int(window))+'\n')
	
o.close()