#!/usr/bin/env python

# Han Tan
# UC Davis
# November 2015

# Uninterleaves interleaved [ ].fq files from batch-sam-junction script
# These are reads that are within specific edges of duplicated/triplicated blocks
# The script then assembles each of the junction files using Price to generate a series of .fa files

# Usage programname.py

import os, sys, math

li = os.listdir(os.getcwd())
todo = filter(lambda x: x[10:15] == 'junct' and x[-3:] == '.fq', li)
todo.sort()
print todo

os.system("mkdir Assemblies")

for file in todo:
    print file
    name = file.split('.')[0]
    u = open(file)
    f = open(name+"-1.fq",'w')
    r = open(name+"-2.fq",'w')

    while True:
        name1 = u.readline()
        if name1 == '':
            break
        seq1 = u.readline()
        p1 = u.readline()
        qual1 = u.readline()
        name2 = u.readline()
        seq2 = u.readline()
        p2 = u.readline()
        qual2 = u.readline()
        f.write(name1+seq1+p1+qual1)
        r.write(name2+seq2+p2+qual2)

    f.close()
    r.close()
    u.close()
    os.system("PriceTI -fpp "+name+"-1.fq"" "+name+"-2.fq"" 400 95 -icf "+name+"-1.fq"" 1 1 5 -nc 20 -a 2 -lenf 200 15 -o "+name+".fa"" ")
    os.system("mv "+name+"*.cycle20.fa"" Assemblies/ ")
    os.system("rm "+name+"*.fa"" ")
    os.system("rm "+name+"-1.fq"" ")
    os.system("rm "+name+"-2.fq"" ")
