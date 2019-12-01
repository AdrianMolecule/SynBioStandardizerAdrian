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
import genestand
#
# Open file to read in fasta sequences for modified and original records
#
handle = open("sequence.txt", "rU")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()
handle = open("sequence.txt", "rU")
recsO = list(SeqIO.parse(handle, "fasta"))
handle.close()
muts_handle = open("Mutations.txt", "a")
#
# Start your engines
#
stnds = ['N','BioB','BglB','MoClo','GB','Chi']
SynthRecs = []
q = 0
for p in records:
	outSeq = copy.deepcopy(records[q])
	recSeq = copy.deepcopy(recsO[q])
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
	changes = []
	genestand.refactor(outSeq, recSeq, changes, stnds)
	SynthRecs.append(copy.deepcopy(outSeq))
	muts_handle.write(str(changes)+'\n')
	q=q+1
#
# Print results
#
print 'Success!'
muts_handle.close()
outFile = "SyntheticSequence.txt"
output_handle = open(outFile, "w")
SeqIO.write(SynthRecs, output_handle, "fasta")
output_handle.close()
