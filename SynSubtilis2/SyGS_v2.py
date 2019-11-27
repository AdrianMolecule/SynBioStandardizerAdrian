#!/usr/bin/env python
#####
#
# Synthetic Biology Gene Standardizer
# Copyright (c) 2015, Tyson R. Shepherd, PhD
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies,
# either expressed or implied, of Uppsala University.
#
#####
import copy
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import genestand2
#
# Open file to read in fasta sequences for modified and original records
#
fIn=open('parsed.txt','r')
geneName=[]
geneSeqIn=[]
geneSeqOut=[]
genes=[]
for line in fIn:
	genes.append(line);
for i in range(0,int(len(genes)/2)):
	geneName.append(genes[int(i*2)].split('>')[1].rstrip());
	geneSeqIn.append(genes[int(i*2+1)].rstrip());
fIn.close()
genes=[];
#
# Start your engines
#
stnds = ['N','BioB','BglB','MoClo','GB','Chi']
SynthRecs = []
q = 0
for p in geneSeqIn:
#
# Look for: non-ATG start codons, non-TAA stop codons,
#           NdeI:        NdeI
#           BioBrick:    EcoRI, SpeI, XbaI, PstI, mfeI, avrII, NheI, NsiI, SbfI, NotI, ApoI
#           BglBrick:    EcoRI, XhoI, BglII, BamHI,
#           MoClo:       BbsI, BsaI, MlyI
#           GoldenBraid: BsmI, BtgZI
#           Chi sites
# Then makes non-conflicting point mutations to highest allowed codon usage
#
	tmpGeneSeqOut=genestand2.refactor(geneName[q], p, stnds, 1)
	geneSeqOut.append(genestand2.mutatePromoters(geneName[q], tmpGeneSeqOut))
	q=q+1
#	print(str(q)+'/'+str(len(geneName)));
#
# Print results
#
print('Success!')
#
genestand2.statistics();
#
outFile = "SyntheticSequence.txt"
output_handle = open(outFile, "w")
for i in range(len(geneName)):
	output_handle.write('>'+geneName[i]+"\n"+geneSeqOut[i]+"\n")
