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
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Restriction import *
from math import *
#
# Defining codon tables for codon minimization
#
RLD = {
	'TGC': 'C1', 'TGT': 'C2', 'GAT': 'D1', 'GAC': 'D2',
	'AGC': 'S1', 'TCT': 'S2', 'AGT': 'S3', 'TCC': 'S4', 'TCA': 'S5', 'TCG': 'S6',
	'CAG': 'Q1', 'CAA': 'Q2', 'ATG': 'M1', 'AAC': 'N1', 'AAT': 'N2',
	'CCG': 'P1', 'CCA': 'P2', 'CCT': 'P3', 'CCC': 'P4', 'AAA': 'K1', 'AAG': 'K2',
	'ACC': 'T1', 'ACG': 'T2', 'ACT': 'T3', 'ACA': 'T4', 'TTT': 'F1', 'TTC': 'F2',
	'GCG': 'A1', 'GCC': 'A2', 'GCA': 'A3', 'GCT': 'A4',
	'GGC': 'G1', 'GGT': 'G2', 'GGG': 'G3', 'GGA': 'G4',
	'ATT': 'I1', 'ATC': 'I2', 'ATA': 'I3',
	'CTG': 'L1', 'TTA': 'L2', 'TTG': 'L3', 'CTT': 'L4', 'CTC': 'L5', 'CTA': 'L6',
	'CAT': 'H1', 'CAC': 'H2',
	'CGT': 'R1', 'CGC': 'R2', 'CGG': 'R3', 'CGA': 'R4', 'AGA': 'R5', 'AGG': 'R6',
	'TGG': 'W1', 'GTG': 'V1', 'GTT': 'V2', 'GTC': 'V3', 'GTA': 'V4',
	'GAA': 'E1', 'GAG': 'E2', 'TAT': 'Y1', 'TAC': 'Y2',
	'TAA': 'STOP', 'TGA': 'STOP', 'TAG': 'STOP'
};
ReducedCodonTable = {
	'C': ['TGC', 'TGT'],
	'D': ['GAT', 'GAC'],
	'S': ['AGC'],
	'Q': ['CAG', 'CAA'],
	'M': ['ATG'],
	'N': ['AAC', 'AAT'],
	'P': ['CCG', 'CCA'],
	'K': ['AAA', 'AAG'],
	'STOP': ['TAA'],
	'T': ['ACC', 'ACG'],
	'F': ['TTT', 'TTC'],
	'A': ['GCG', 'GCC', 'GCA'],
	'G': ['GGC', 'GGT'],
	'I': ['ATT', 'ATC'],
	'L': ['CTG'],
	'H': ['CAT', 'CAC'],
	'R': ['CGT', 'CGC'],
	'W': ['TGG'],
	'V': ['GTG', 'GTT', 'GTC'],
	'E': ['GAA', 'GAG'],
	'Y': ['TAT', 'TAC']
};
KilledCodons = {
	'C': [],
	'D': [],
	'S': ['TCT', 'AGT', 'TCC', 'TCA', 'TCG'],
	'Q': [],
	'M': [],
	'N': [],
	'P': ['CCT', 'CCC'],
	'K': ['AAA', 'AAG'],
	'STOP': ['TAA', 'TGA', 'TAG'],
	'T': ['ACT', 'ACA'],
	'F': [],
	'A': ['GCT'],
	'G': ['GGG', 'GGA'],
	'I': ['ATA'],
	'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA'],
	'H': [],
	'R': ['CGG', 'CGA', 'AGA', 'AGG'],
	'W': [],
	'V': ['GTA'],
	'E': [],
	'Y': []
};
ReducedCodonStat = {
	'C': [54, 46],
	'D': [63, 37],
	'S': [100],
	'Q': [66, 34],
	'M': [100],
	'N': [51, 49],
	'P': [71, 29],
	'K': [74, 26],
	'STOP': [100],
	'T': [62, 38],
	'F': [58, 42],
	'A': [40, 31, 29],
	'G': [51, 49],
	'I': [56, 44],
	'L': [100],
	'H': [57, 43],
	'R': [50, 50],
	'W': [100],
	'V': [42, 34, 24],
	'E': [68, 32],
	'Y': [59, 41]
};
#
# Meat and potatoes
#
def refactor( seqIn, stnds, minCodon ):
	"""Refactor sequences: outSeq is output sequence, record is input, changes is mutations"""
	changes=[]
	if minCodon==1:
		newSeqIn=''
		for i in range(0,int(len(seqIn)/3)):
			tmpCodon=seqIn[i]+seqIn[i+1]+seqIn[i+2]
			newSeqIn=newSeqIn+MinimizeCodonUsage(tmpCodon)
		seqIn=newSeqIn
	m = 1
	while m == 1:
		stpCdn = str(seqIn[len(seqIn)-3:len(seqIn)])
		m = 0
		if str(seqIn[0:3]) != 'ATG':
			A = 'Start Codon: '+str(seqIn[0:3])+'1'+'ATG'
			changes.append(A)
			seqIn='ATG'+seqIn[3:]
			# recSeq.seq='ATG'+recSeq.seq[3:]
		if 'N' in stnds:
			ndeInst = NdeI.search(seqIn)
			for i in ndeInst:
				m = 1
				j = (i - 3) % 3
				k = i - 2
				if j == 0:
					k = k + 2
					A = "NdeI: T"+str(k)+"C"
					seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK == "CT" or stK == "CC"):
						A = "NdeI: C"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
						changes.append(A)
					A = "NdeI: A"+str(k)+"T"
					seqIn = seqIn[:k-1]+"T"+seqIn[k:]
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if (stK == "C" or stK == "A" or stK == "G"):
						A = "NdeI: A"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					elif stK == "T":
						A = "NdeI: TCA"+str(k-5)+"AGC"
						seqIn = seqIn[:k-6]+"AGC"+seqIn[k-3:]
					else:
						A = "NdeI: T"+str(k)+"C"
						seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				changes.append(A)
		if 'BglB' in stnds:
			xhoInst = XhoI.search(seqIn)
			for i in xhoInst:
				m = 1
				j = (i - 2) % 3
				k = i - 1
				if j == 0:
					k = k + 2
					A = "XhoI: C"+str(k)+"G"
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					stJ = str(seqIn[k+2])
					if (stJ == "A" or stJ == "G"):
						A = "XhoI: AG"+stJ+str(k+1)+"CGC"
						seqIn = seqIn[:k]+"CGC"+seqIn[k+3:]
					elif (stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC"):
						A = "XhoI: C"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					else:
						A = "XhoI: TCG"+str(k-2)+"AGC"
						seqIn = seqIn[:k-3]+"AGC"+seqIn[k:]
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if stK == "C":
						A = "XhoI: T"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
						changes.append(A)
					A = "XhoI: A"+str(k)+"C"
					seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				changes.append(A)
		if (('BioB' in stnds) or ('BglB' in stnds)):
			ecoInst = EcoRI.search(seqIn)
			for i in ecoInst:
				m = 1
				j = (i - 2) % 3
				k = i - 1
				if j == 0:
					k = k + 2
					A = "EcoRI: A"+str(k)+"G"
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if stK == "AA":
						A = "EcoRI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GA":
						A = "EcoRI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GG":
						A = "EcoRI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "AG":
						A = "EcoRI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
					elif stK == "CG":
						A = "EcoRI: CGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					else:
						A = "EcoRI: T"+str(k)+"C"
						seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if stK == "C" or stK == "A":
						A = "EcoRI: "+stK+"GA"+str(k-5)+"CGT"
						seqIn = seqIn[:k-6]+"CGT"+seqIn[k-3:]
					elif stK == "G":
						A = "EcoRI: A"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "T":
						A = "EcoRI: G"+str(k-4)+"A"
						seqIn = seqIn[:k-6]+"TAA"+seqIn[k-3:]
					else:
						A = "EcoRI: T"+str(k)+"C"
						seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				changes.append(A)
		if 'BioB' in stnds:
			xbaInst = XbaI.search(seqIn)
			for i in xbaInst:
				m = 1
				j = (i - 2) % 3
				k = i - 1
				if j == 0:
					k = k + 3
					A = "XbaI: AGA"+str(k)+"CGT"
					changes.append(A)
					seqIn = seqIn[:k-1]+"CGT"+seqIn[k+2:]
				elif j == 2:
					k = k + 3
					A = "XbaI: A"+str(k)+"G"
					changes.append(A)
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 1:
					k = k + 4
					A = "XbaI: G"+str(k)+"A"
					changes.append(A)
					seqIn = seqIn[:k-1]+"A"+seqIn[k:]
			speInst = SpeI.search(seqIn)
			for i in speInst:
				m = 1
				j = (i - 2) % 3
				k = i - 1
				if j == 0:
					k = k + 2
					A = "SpeI: T"+str(k)+"C"
					changes.append(A)
					seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				elif j == 2:
					k = k + 3
					A = "SpeI: A"+str(k)+"G"
					changes.append(A)
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 1:
					k = k + 4
					A = "SpeI: G"+str(k)+"A"
					changes.append(A)
					seqIn = seqIn[:k-1]+"A"+seqIn[k:]
			pstInst = PstI.search(seqIn)
			for i in pstInst:
				m = 1
				j = (i - 6) % 3
				k = i - 5
				if j == 0:
					k = k + 5
					stK = str(seqIn[k-1:k+5])
					if (stK == "GGATCT" or stK == "GTGCAT"):
						A = "PstI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					else:
						A = "PstI: G"+str(k)+"A"
						seqIn = seqIn[:k-1]+"A"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					stJ = str(seqIn[k-8:k-4])
					if stJ == "GAAT" or stJ == "AAAT":
						A = "PstI: C"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC":
						A = "PstI: C"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					elif stK == "AG" or stK=="AC" or stK=="CA":
						A = "PstI: C"+str(k)+"T"
						seqIn = seqIn[:k-1]+"T"+seqIn[k:]
					else:
						A = "PstI: C"+str(k-3)+"T"
						seqIn = seqIn[:k-4]+"T"+seqIn[k-3:]
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if stK == "C" or stK=="G":
						A = "PstI: T"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					else:
						A = "PstI: A"+str(k)+"G"
						seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				changes.append(A)
			mfeInst = MfeI.search(seqIn)
			for i in mfeInst:
				m = 1
				j = (i - 2) % 3
				k = i - 1
				if j == 0:
					k = k + 3
					A = "MfeI: T"+str(k)+"C"
					seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC"):
						A = "MfeI: C"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					elif (stK == "AG" or stK == "AC" or stK == "CA"):
						A = "MfeI: T"+str(k)+"C"
						seqIn = seqIn[:k-1]+"C"+seqIn[k:]
					else:
						A = "MfeI: C"+str(k-3)+"T"
						seqIn = seqIn[:k-4]+"T"+seqIn[k-3:]
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if (stK == "C" or stK == "A" or stK == "G"):
						A = "MfeI: A"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					elif stK == "T":
						A = "MfeI: TCA"+str(k-5)+"AGC"
						seqIn = seqIn[:k-6]+"AGC"+seqIn[k-3:]
					else:
						A = "MfeI: T"+str(k)+"C"
						seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				changes.append(A)
			avrInst = AvrII.search(seqIn)
			for i in avrInst:
				m = 1
				j = (i - 2) % 3
				k = i - 1
				if j == 0:
					k = k + 3
					A = "AvrII: AGG"+str(k)+"CGT"
					changes.append(A)
					seqIn = seqIn[:k-1]+"CGT"+seqIn[k+2:]
				elif j == 2:
					k = k + 3
					A = "AvrII: A"+str(k)+"G"
					changes.append(A)
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 1:
					k = k + 4
					A = "AvrII: G"+str(k)+"A"
					changes.append(A)
					seqIn = seqIn[:k-1]+"A"+seqIn[k:]
			nheInst = NheI.search(seqIn)
			for i in nheInst:
				m = 1
				j = (i - 2) % 3
				k = i - 1
				if j == 0:
					k = k + 2
					A = "NheI: T"+str(k)+"G"
					changes.append(A)
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					A = "NheI: A"+str(k)+"G"
					changes.append(A)
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 1:
					k = k + 4
					A = "NheI: G"+str(k)+"A"
					changes.append(A)
					seqIn = seqIn[:k-1]+"A"+seqIn[k:]
			nsiInst = NsiI.search(seqIn)
			for i in nsiInst:
				m = 1
				j = (i - 6) % 3
				k = i - 5
				if j == 0:
					k = k + 5
					A = "NsiI: T"+str(k)+"C"
					changes.append(A)
					seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK=="CA" or stK=="GC" or stK=="CT" or stK=="CC" or stK=="AC"):
						A = "NsiI: A"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					elif stK == "GT":
						A = "NsiI: A"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					elif (stK == "AT" or stK == "CG"):
						A = "NsiI: A"+str(k-3)+"T"
						seqIn = seqIn[:k-4]+"T"+seqIn[k-3:]
					elif (stK == "GG"):
						A = "NsiI: A"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "AG":
						A = "NsiI: AGA"+str(k-3)+"CGT"
						seqIn = seqIn[:k-6]+"CGT"+seqIn[k-3:]
					elif stK == "TT":
						A = "NsiI: TTA"+str(k-3)+"CTG"
						seqIn = seqIn[:k-6]+"CTG"+seqIn[k-3:]
					elif stK == "TC":
						A = "NsiI: TCA"+str(k-3)+"AGC"
						seqIn = seqIn[:k-6]+"AGC"+seqIn[k-3:]
					else:
						A = "NsiI: C"+str(k)+"T"
						seqIn = seqIn[:k-1]+"T"+seqIn[k:]
					changes.append(A)
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k+1:k+3])
					if (stK == "TA" or stK == "TG"):
						A = "NsiI: T"+stK+str(k+1)+"CTG"
						seqIn = seqIn[:k]+"CTG"+seqIn[k+3:]
					else:
						A = "NsiI: A"+str(k)+"G"
						seqIn = seqIn[:k-1]+"G"+seqIn[k:]
					changes.append(A)
		if 'BglB' in stnds:
			bglInst = BglII.search(seqIn)
			for i in bglInst:
				m = 1
				j = (i - 2) % 3
				k = i - 1
				if j == 0:
					k = k + 0
					A = "BglII: AGA"+str(k)+"CGT"
					changes.append(A)
					seqIn = seqIn[:k-1]+"CGT"+seqIn[k+2:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK=="CA" or stK=="GC" or stK=="CT" or stK=="CC" or stK=="AC"):
						A = "BglII: A"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					elif stK == "GT":
						A = "BglII: A"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					elif (stK == "AT" or stK == "CG"):
						A = "BglII: A"+str(k-3)+"T"
						seqIn = seqIn[:k-4]+"T"+seqIn[k-3:]
					elif (stK == "GG"):
						A = "BglII: A"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "AG":
						A = "BglII: AGA"+str(k-3)+"CGT"
						seqIn = seqIn[:k-6]+"CGT"+seqIn[k-3:]
					elif stK == "TT":
						A = "BglII: TTA"+str(k-3)+"CTG"
						seqIn = seqIn[:k-6]+"CTG"+seqIn[k-3:]
					elif stK == "TC":
						A = "BglII: TCA"+str(k-3)+"AGC"
						seqIn = seqIn[:k-6]+"AGC"+seqIn[k-3:]
					else:
						A = "BglII: T"+str(k)+"C"
						seqIn = seqIn[:k-1]+"C"+seqIn[k:]
					changes.append(A)
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k+1:k+3])
					if (stK == "TA" or stK == "TG"):
						A = "BglII: T"+stK+str(k+1)+"CTG"
						seqIn = seqIn[:k]+"CTG"+seqIn[k+3:]
					else:
						A = "BglII: C"+str(k)+"T"
						seqIn = seqIn[:k-1]+"T"+seqIn[k:]
					changes.append(A)
			bamInst = BamHI.search(seqIn)
			for i in bamInst:
				m = 1
				j = (i - 2) % 3
				k = i - 1
				if j == 0:
					k = k + 2
					A = "BamHI: A"+str(k)+"C"
					changes.append(A)
					seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK == "TT" or stK == "TA" or stK == "AA" or stK == "GA"):
						A = "BamHI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif (stK == "CG" or stK == "AC" or stK == "GC" or stK == "GG"):
						A = "BamHI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "GT":
						A = "BamHI: G"+str(k-3)+"T"
						seqIn = seqIn[:k-4]+"T"+seqIn[k-3:]
					elif stK == "TC":
						A = "BamHI: TCG"+str(k-5)+"AGC"
						seqIn = seqIn[:k-6]+"AGC"+seqIn[k-3:]
					elif stK == "AG":
						A = "BamHI: AGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
					else:
						A = "BamHI: T"+str(k)+"C"
						seqIn = seqIn[:k-1]+"C"+seqIn[k:]
					changes.append(A)
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if (stK == "C" or stK == "A"):
						A = "BamHI: "+stK+"GG"+str(k-5)+"CGT"
						seqIn = seqIn[:k-6]+"CGT"+seqIn[k-3:]
					elif stK=="G":
						A = "BamHI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					else:
						A = "BamHI: C"+str(k)+"T"
						seqIn = seqIn[:k-1]+"T"+seqIn[k:]
					changes.append(A)
		if 'BioB' in stnds:
			notInst = NotI.search(seqIn)
			for i in notInst:
				m = 1
				j = (i - 3) % 3
				k = i - 2
				if j == 0:
					k = k + 5
					A = "NotI: C"+str(k)+"G"
					changes.append(A)
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK == "AA" or stK == "GA"):
						A = "NotI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GG":
						A = "NotI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "AG":
						A = "NotI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
					elif stK == "CG":
						A = "NotI: CGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					else:
						A = "NotI: G"+str(k)+"C"
						seqIn = seqIn[:k-1]+"C"+seqIn[k:]
					changes.append(A)
				elif j == 1:
					k = k + 4
					A = "NotI: C"+str(k-3)+"T"
					changes.append(A)
					seqIn = seqIn[:k-4]+"T"+seqIn[k-3:]
			apoInst = ApoI.search(seqIn)
			for i in apoInst:
				m = 1
				j = (i - 2) % 3
				k = i - 1
				if j == 0:
					k = k + 2
					A = "ApoI: A"+str(k)+"G"
					changes.append(A)
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					A = "ApoI: T"+str(k)+"C"
					changes.append(A)
					seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				elif j == 1:
					k = k + 4
					A = "ApoI: T"+str(k)+"C"
					changes.append(A)
					seqIn = seqIn[:k-1]+"C"+seqIn[k:]
		if 'MoClo' in stnds:
			bbsInst = seqIn.find("GAAGAC")
			if bbsInst > 0:
				m = 1
				j = bbsInst % 3
				k = bbsInst + 1
				if j == 0:
					k = k + 2
					A = "BbsI: A"+str(k)+"G"
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK == "AG"):
						A = "BbsI: AGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
						changes.append(A)
					A = "BbsI: G"+str(k)+"A"
					seqIn = seqIn[:k-1]+"A"+seqIn[k:]
				elif j == 1:
					k = k + 2
					A = "BbsI: AGA"+str(k)+"CGT"
					seqIn = seqIn[:k-1]+"CGT"+seqIn[k+2:]
				changes.append(A)
				bbsInst = 0
			bbsInst = seqIn.find("GTCTTC")
			if bbsInst > 0:
				m = 1
				j = bbsInst % 3
				k = bbsInst + 1
				if j == 0:
					k = k + 2
					A = "BbsI: C"+str(k)+"G"
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if stK == "AA":
						A = "BbsI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GA":
						A = "BbsI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GG":
						A = "BbsI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "AG":
						A = "BbsI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
					elif stK == "CG":
						A = "BbsI: CGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					else:
						A = "BbsI: TCT"+str(k-2)+"AGC"
						seqIn = seqIn[:k-3]+"AGC"+seqIn[k:]
				elif j == 1:
					k = k + 4
					A = "BbsI: T"+str(k)+"G"
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				changes.append(A)
				bbsInst = 0
			bsaInst = seqIn.find("GGTCTC")
			if bsaInst > 0:
				m = 1
				j = bsaInst % 3
				k = bsaInst + 1
				if j == 0:
					k = k + 5
					A = "BsaI: C"+str(k)+"G"
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if stK == "AA":
						A = "BsaI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GA":
						A = "BsaI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GG":
						A = "BsaI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "AG":
						A = "BsaI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
					elif stK == "CG":
						A = "BsaI: CGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					else:
						A = "BsaI: C"+str(k)+"G"
						seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if (stK == "C" or stK == "G"):
						A = "BsaI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "A":
						A = "BsaI: AGG"+str(k-5)+"CGT"
						seqIn = seqIn[:k-6]+"CGT"+seqIn[k-3:]
					else:
						A = "BsaI: TCT"+str(k-2)+"AGC"
						seqIn = seqIn[:k-3]+"AGC"+seqIn[k:]
				changes.append(A)
				bsaInst = 0
			bsaInst = seqIn.find("GAGACC")
			if bsaInst > 0:
				m = 1
				j = bsaInst % 3
				k = bsaInst + 1
				if j == 0:
					k = k + 2
					A = "BsaI: G"+str(k)+"A"
					seqIn = seqIn[:k-1]+"A"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if stK == "AA":
						A = "BsaI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
						changes.append(A)
					elif stK == "GA":
						A = "BsaI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
						changes.append(A)
					elif stK == "AG":
						A = "BsaI: AGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
						changes.append(A)
					A = "BsaI: AGA"+str(k-2)+"CGT"
					seqIn = seqIn[:k-3]+"CGT"+seqIn[k:]
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if (stK == "C" or stK == "G"):
						A = "BsaI: A"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "A":
						A = "BsaI: AGA"+str(k-5)+"CGT"
						seqIn = seqIn[:k-6]+"CGT"+seqIn[k-3:]
					else:
						A = "BsaI: C"+str(k)+"T"
						seqIn = seqIn[:k-1]+"T"+seqIn[k:]
				changes.append(A)
				bsaInst = 0
			mlyInst = seqIn.find("GAGTC")
			if mlyInst > 0:
				m = 1
				j = mlyInst % 3
				k = mlyInst + 1
				if j == 0:
					k = k + 2
					A = "MlyI: G"+str(k)+"A"
					seqIn = seqIn[:k-1]+"A"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if stK == "AA":
						A = "MlyI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GA":
						A = "MlyI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GG":
						A = "MlyI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "AG":
						A = "MlyI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
					elif stK == "CG":
						A = "MlyI: CGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					else:
						A = "MlyI: T"+str(k)+"C"
						seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if stK == "C" or stK == "A":
						A = "MlyI: "+stK+"GA"+str(k-5)+"CGT"
						seqIn = seqIn[:k-6]+"CGT"+seqIn[k-3:]
					elif stK == "G":
						A = "MlyI: A"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "T":
						A = "MlyI: G"+str(k-4)+"A"
						seqIn = seqIn[:k-6]+"TAA"+seqIn[k-3:]
					else:
						A = "MlyI: C"+str(k)+"G"
						seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				changes.append(A)
				mlyInst = 0
			mlyInst = seqIn.find("GACTC")
			if mlyInst > 0:
				m = 1
				j = mlyInst % 3
				k = mlyInst + 1
				if j == 0:
					k = k + 2
					A = "MlyI: C"+str(k)+"T"
					seqIn = seqIn[:k-1]+"T"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if stK == "AA":
						A = "MlyI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GA":
						A = "MlyI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif stK == "GG":
						A = "MlyI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "AG":
						A = "MlyI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
					elif stK == "CG":
						A = "MlyI: CGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					else:
						A = "MlyI: T"+str(k)+"G"
						seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if stK == "C" or stK == "A":
						A = "MlyI: "+stK+"GA"+str(k-5)+"CGT"
						seqIn = seqIn[:k-6]+"CGT"+seqIn[k-3:]
						changes.append(A)
					A = "MlyI: C"+str(k)+"G"
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				changes.append(A)
				mlyInst = 0
		if 'GB' in stnds:
			bsmInst = seqIn.find("CGTCTC")
			if bsmInst > 0:
				m = 1
				j = bsmInst % 3
				k = bsmInst + 1
				if j == 0:
					k = k + 5
					A = "BsmBI: C"+str(k)+"G"
					changes.append(A)
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC"):
						A = "BsmBI: C"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					else:
						A = "BsmBI: C"+str(k)+"G"
						seqIn = seqIn[:k-1]+"G"+seqIn[k:]
					changes.append(A)
				elif j == 1:
					k = k + 4
					A = "BsmBI: TCT"+str(k)+"AGC"
					changes.append(A)
					seqIn = seqIn[:k-3]+"AGC"+seqIn[k:]
				bsmInst = 0
			bsmInst = seqIn.find("GCAGAG")
			if bsmInst > 0:
				m = 1
				j = bsmInst % 3
				k = bsmInst + 1
				if j == 0:
					k = k + 2
					A = "BsmBI: G"+str(k)+"G"
					changes.append(A)
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					stJ = str(seqIn[k+2])
					if (stJ == "A" or stJ == "G"):
						A = "BsmBI: AG"+stJ+str(k+1)+"CGT"
						seqIn = seqIn[:k]+"CGT"+seqIn[k+3:]
					elif (stJ == "T"):
						A = "BsmBI: AGT"+str(k+1)+"TCT"
						seqIn = seqIn[:k]+"TCT"+seqIn[k+3:]
					elif (stK == "TT" or stK == "TA" or stK == "AA" or stK == "GA"):
						A = "BsmBI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif (stK == "CG" or stK == "AC" or stK == "GC" or stK == "GG"):
						A = "BsmBI: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "GT":
						A = "BsmBI: G"+str(k-3)+"T"
						seqIn = seqIn[:k-4]+"T"+seqIn[k-3:]
					elif stK == "TC":
						A = "BsmBI: TCG"+str(k-5)+"AGC"
						seqIn = seqIn[:k-6]+"AGC"+seqIn[k-3:]
					elif stK == "AG":
						A = "BsmBI: AGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
					else:
						A = "BsmBI: G"+str(k)+"A"
						seqIn = seqIn[:k-1]+"A"+seqIn[k:]
					changes.append(A)
				elif j == 1:
					k = k + 2
					A = "BsmBI: AGA"+str(k)+"CGT"
					changes.append(A)
					seqIn = seqIn[:k-1]+"CGT"+seqIn[k+2:]
				bsmInst = 0
			btgzInst = seqIn.find("GCGATG")
			if btgzInst > 0:
				m = 1
				j = btgzInst % 3
				k = btgzInst + 1
				if j == 0:
					k = k + 2
					A = "BtgZI: G"+str(k)+"C"
					changes.append(A)
					seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK == "AA" or stK == "GA"):
						A = "BtgZI: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
						changes.append(A)
					elif (stK == "AG"):
						A = "BtgZI: AGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
						changes.append(A)
					A = "BtgZI: A"+str(k)+"T"
					changes.append(A)
					seqIn = seqIn[:k-1]+"T"+seqIn[k:]
				elif j == 1:
					k = k + 4
					A = "BtgZI: C"+str(k-3)+"T"
					changes.append(A)
					seqIn = seqIn[:k-4]+"T"+seqIn[k-3:]
				btgzInst = 0
			btgzInst = seqIn.find("CATCGC")
			if btgzInst > 0:
				m = 1
				j = btgzInst % 3
				k = btgzInst + 1
				if j == 0:
					k = k + 2
					A = "BtgZI: T"+str(k)+"C"
					changes.append(A)
					seqIn = seqIn[:k-1]+"C"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC"):
						A = "BtgZI: C"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					else:
						A = "BtgZI: C"+str(k)+"T"
						seqIn = seqIn[:k-1]+"T"+seqIn[k:]
					changes.append(A)
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if stK == "C":
						A = "BtgZI: A"+str(k-3)+"G"
						seqIn = seqIn[:k-6]+"CCG"+seqIn[k-3:]
						changes.append(A)
					A = "BtgZI: TCG"+str(k-2)+"AGC"
					changes.append(A)
					seqIn = seqIn[:k-3]+"AGC"+seqIn[k:]
				btgzInst = 0
		if 'chi' in stnds:
			chiInst = seqIn.find("GCTGGTGG")
			if chiInst > 0:
				m = 1
				j = chiInst % 3
				k = chiInst + 1
				if j == 0:
					k = k + 2
					A = "Chi site: T"+str(k)+"G"
					changes.append(A)
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK == "TT" or stK == "TA" or stK == "AA" or stK == "GA"):
						A = "Chi site: G"+str(k-3)+"A"
						seqIn = seqIn[:k-4]+"A"+seqIn[k-3:]
					elif (stK == "CG" or stK == "AC" or stK == "GC" or stK == "GG"):
						A = "Chi site: G"+str(k-3)+"C"
						seqIn = seqIn[:k-4]+"C"+seqIn[k-3:]
					elif stK == "GT":
						A = "Chi site: G"+str(k-3)+"T"
						seqIn = seqIn[:k-4]+"T"+seqIn[k-3:]
					elif stK == "TC":
						A = "Chi site: TCG"+str(k-5)+"AGC"
						seqIn = seqIn[:k-6]+"AGC"+seqIn[k-3:]
					elif stK == "AG":
						A = "Chi site: AGG"+str(k-5)+"CGC"
						seqIn = seqIn[:k-6]+"CGC"+seqIn[k-3:]
					else:
						A = "Chi site: G"+str(k+3)+"T"
						seqIn = seqIn[:k+2]+"T"+seqIn[k+3:]
					changes.append(A)
				elif j == 1:
					k = k + 1
					A = "Chi site: C"+str(k)+"T"
					changes.append(A)
					seqIn = seqIn[:k-1]+"T"+seqIn[k:]
				chiInst = 0
			chiInst = seqIn.find("CCACCAGC")
			if chiInst > 0:
				m = 1
				j = chiInst % 3
				k = chiInst + 1
				if j == 0:
					k = k + 2
					A = "Chi site: A"+str(k)+"G"
					seqIn = seqIn[:k-1]+"G"+seqIn[k:]
				elif j == 2:
					k = k + 3
					stK = str(seqIn[k-6:k-4])
					if (stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC"):
						A = "Chi site: C"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					else:
						A = "Chi site: C"+str(k)+"T"
						seqIn = seqIn[:k-1]+"T"+seqIn[k:]
				elif j == 1:
					k = k + 4
					stK = str(seqIn[k-6])
					if (stK == "C" or stK == "A" or stK == "G"):
						A = "Chi site: C"+str(k-3)+"G"
						seqIn = seqIn[:k-4]+"G"+seqIn[k-3:]
					else:
						A = "Chi site: TCC"+str(k-3)+"AGC"
						seqIn = seqIn[:k-6]+"AGC"+seqIn[k-3:]
				changes.append(A)
				chiInst = 0
		if (stpCdn == 'TGA' or stpCdn == 'TAG'):
			A = 'Stop Codon: '+str(seqIn[len(seqIn)-3:len(seqIn)])+str(len(seqIn)-3)+'TAA'
			changes.append(A)
			seqIn=seqIn[:len(seqIn)-3]+'TAA'
		if (stpCdn != 'TAA' and stpCdn != 'TAG' and stpCdn != 'TGA'):
			A = 'Add Stop Codon: '+str(len(seqIn))+'TAA'
			changes.append(A)
			seqIn=seqIn[:len(seqIn)]+'TAA'
		if len(changes) > (len(seqIn)/5):
			m = 0
			print("Error: Conflicting mutations")
			for i in changes:
				print(i)
			# print(recSeq.seq)
			sys.exit("Error: Mutation limit exceded")
#
# Test the output Protein sequence vs. the input Protein Sequence
#
	# outProt=seqIn.translate()
	# inProt=recSeq.seq.translate()
	# if str(outProt) != str(inProt):
	# 	print(recSeq.seq)
	# 	print(seqIn)
	# 	print(inProt)
	# 	print(outProt)
	# 	print("Error in silent mutation")
	# 	sys.exit("Error in silent mutation")
	fOut=open('Mutations.txt','a')
	for i in changes:
		fOut.write(i+"\n")
	fOut.close()
	outSeq=seqIn
	seqIn=newSeqIn
	return outSeq

def MinimizeCodonUsage( inCodon ):
	inAA = RLD[str(inCodon)][:1];
	inNum = int(RLD[str(inCodon)][1:]);
	outCodon = inCodon
	tmpCList = []
	if ((inAA != 'M') and (inAA != 'W') and (inAA != 'STOP') and (inAA != 'S') and (inAA != 'L')):
		for i in range(0, len(ReducedCodonTable[inAA])):
			for j in range (0, ReducedCodonStat[inAA][i]):
				tmpCList.append(str(ReducedCodonTable[inAA][i]));
		while outCodon in KilledCodons[inAA]:
			outCodon = random.choice(tmpCList);
	if inAA == 'STOP':
		outCodon = 'TAA'
	return outCodon;
