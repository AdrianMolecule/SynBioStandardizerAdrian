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
import random
SynCodons = { 
	'C': ['TGC', 'TGT'], 
	'D': ['GAT', 'GAC'], 
	'S': ['AGC', 'TCT', 'AGT', 'TCC', 'TCA', 'TCG'], 
	'Q': ['CAG', 'CAA'], 
	'M': ['ATG'], 
	'N': ['AAC', 'AAT'], 
	'P': ['CCG', 'CCA', 'CCT', 'CCC'], 
	'K': ['AAA', 'AAG'], 
	'STOP': ['TAA', 'TGA', 'TAG'], 
	'T': ['ACC', 'ACG', 'ACT', 'ACA'], 
	'F': ['TTT', 'TTC'], 
	'A': ['GCG', 'GCC', 'GCA', 'GCT'], 
	'G': ['GGC', 'GGT', 'GGG', 'GGA'], 
	'I': ['ATT', 'ATC', 'ATA'], 
	'L': ['CTG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA'], 
	'H': ['CAT', 'CAC'], 
	'R': ['CGT', 'CGC', 'CGG', 'CGA', 'AGA', 'AGG'], 
	'W': ['TGG'], 
	'V': ['GTG', 'GTT', 'GTC', 'GTA'], 
	'E': ['GAA', 'GAG'], 
	'Y': ['TAT', 'TAC'] 
};

CodonStat = { 
	'C': [54, 46], 
	'D': [63, 37], 
	'S': [25, 17, 16, 15, 14, 14], 
	'Q': [66, 34], 
	'M': [100], 
	'N': [51, 49], 
	'P': [49, 20, 18, 13], 
	'K': [74, 26], 
	'STOP': [61, 30, 9], 
	'T': [40, 25, 19, 17], 
	'F': [58, 42], 
	'A': [33, 26, 23, 18], 
	'G': [37, 35, 15, 13], 
	'I': [49, 39, 11], 
	'L': [47, 14, 13, 12, 10, 4], 
	'H': [57, 43], 
	'R': [36, 36, 11, 7, 7, 4], 
	'W': [100], 
	'V': [35, 28, 20, 17], 
	'E': [68, 32], 
	'Y': [59, 41] 
};

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

def RanReCodon( inCodon ):
	inAA = RLD[str(inCodon)][:1];
	outCodon = inCodon
	tmpCList = []
	if ((inAA != 'M') and (inAA != 'W') and (inAA != 'STOP')):
		for i in range(0, len(SynCodons[inAA])):
			for j in range (0, CodonStat[inAA][i]):
				tmpCList.append(str(SynCodons[inAA][i]));
		while outCodon == inCodon:
			outCodon = random.choice(tmpCList);
	if inAA == 'STOP':
		outCodon = 'TAA'
	return outCodon;
	
def MaxReCodon( inCodon ):
	inAA = RLD[str(inCodon)][:1];
	inNum = int(RLD[str(inCodon)][1:]);
	outCodon = inCodon
	if ((inAA != 'M') and (inAA != 'W') and (inAA != 'STOP')):
		if inNum > 1:
			outCodon = SynCodons[inAA][0];
		else:
			outCodon = SynCodons[inAA][1];
	if inAA == 'STOP':
		outCodon = 'TAA'
	return outCodon;
