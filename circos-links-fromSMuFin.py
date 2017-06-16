#! /usr/bin/env python

import sys

#######################################################################
##
## Han Tan, University of Maine, 2017
##
## This script is used to generate the link file for circos plots
## that represents the junctions as a result of chromothripsis
## from SMuFin's somatic_large_SVs output files
##
## USAGE: circos-links-fromSMuFin.py somatic_large_SVs_$file output.txt
##
#########################################################################

largeSV = open(sys.argv[1])
o = open(sys.argv[2], 'w')

for line in largeSV:
	j = line.split('\t')
	if j[0] != "Mut_ID":
		chr1 = j[2]
		pos1 = j[3]
		chr2 = j[4]
		pos2 = j[5]
		o.write(chr1.lower()+'\t'+pos1+'\t'+str(int(pos1)+100)+'\t'+chr2.lower()+'\t'+pos2+'\t'+str(int(pos2)+100)+'\n')
	
o.close()