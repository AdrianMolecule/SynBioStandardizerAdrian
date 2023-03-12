#!/usr/bin/env python
import copy
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import zerorpc
#
# Open file to read in fasta sequences for modified and original records
#
fName=sys.argv[1]
fIn=open(fName,'r')
geneName=[]
geneSeqIn=[]
geneSeqOut={}
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
c = zerorpc.Client()
c.connect("tcp://192.168.1.180:4242")
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
	print("process gene:",p)
	if (((p[:3]=='GTG') or (p[:3]=='TTG') or (p[:3]=='ATT')) and ((p[-3:]=='TAA') or (p[-3:]=='TGA') or (p[-3:]=='TAG'))):
		p='ATG'+p[3:]
	tmpGeneSeqOut=c.MinimizeCodonUsage(p)
	tmpGeneSeqOut=c.refactor(geneName[q], p, stnds)
	tmpGeneSeqOut=c.mutatePromoters(tmpGeneSeqOut)
	if (c.tlCheck(p,tmpGeneSeqOut)):
		geneSeqOut[geneName[q]]=tmpGeneSeqOut
	else:
		"There was an error in the translation check"
		geneSeqOut[geneName[q]]="Error in translation check"
	q=q+1
	print(str(q)+'/'+str(len(geneName)));
#
# Print results
#
print('Success!')
#
# genestand2.statistics();
#
outFile = "SyntheticSequence.txt"
output_handle = open(outFile, "w")
for i in geneSeqOut:
	output_handle.write('>'+i+"\n"+geneSeqOut[i]+"\n\n")
