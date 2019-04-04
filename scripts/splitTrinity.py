#!/usr/bin/env python
# This code uses splits Trinity fasta output according to the number of genes indicated by a parameter 
# Trinity fasta contains transcripts indicated as COMPONENT_ID, GENE_ID, ISOFORM_ID.  
# The unique identifier for a gene is given by combination of COMPONENT_ID and GENE_ID
# Example:
#>c0_g1_i4
#>c0_g1_i5
#>c0_g1_i6
#...
#>c1_g13_i1
#>c1_g12_i1
#>c1_g11_i1
#>c1_g10_i1
#>c1_g17_i1


__author__ = 'luca.cozzuto@crg.eu'
# -*- coding utf-8 -*-

#MODULES
import sys
import re
import optparse
import gzip
import pprint
import uuid
import os
from Bio import SeqIO
from collections import defaultdict

# Define options
def options_arg():
    usage = "usage: %prog "
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-n', '--seqname', help='Input sequence name', dest="seqname", default="Trinity.fasta")
    parser.add_option('-s', '--size', help='Size of chunks', dest="size", default=3000)
    (opts,args) = parser.parse_args()
    if opts.seqname and opts.size:pass
    else: parser.print_help()
    return (opts)

		
opts = options_arg()
seqname = opts.seqname
size = opts.size



num = 0
old_gene_id = ""

# make fasta index
records = SeqIO.index(seqname, "fasta")

# fasta_sequences = SeqIO.parse(open(seqname),'fasta')
# for fasta in fasta_sequences:
# 	name, sequence = fasta.id, str(fasta.seq)
# 	pieces = name.split("_")
# 	gene_id = pieces[0] + "_" + pieces[1]
# 	if (gene_id != old_gene_id):
# 		print gene_id + " " + str(num)
# 		num=+1
# 		old_gene_id=gene_id


		


