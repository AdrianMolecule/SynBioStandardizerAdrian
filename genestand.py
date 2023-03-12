import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Restriction
from math import *
import random
from ast import literal_eval
import json
import zerorpc
#
# Defining codon tables for codon minimization
#

#
# Meat and potatoes
#
class RefactorRPC(object):
	TxUnits = [['glyW','cysT','leuZ'],['serV','argV'],['metU','glnV'],['glnU'],
		['asnT'],['proK'],['hisR','leuT','proM'],['valT','lysW'],['thrU','tyrU','glyT','thrT'],
		['pheU'],['ileT','alaT'],['alaW'],['rrsC','gltU','rrlC','rrfC','aspT','trpT'],['valW'],['metV']]
	muts={}

	promMuts=[]
	print("Starting SERVER")
	def refactor( self, gName, inSeq, stnds ):
		"""Refactor sequences: outSeq is output sequence, record is input, changes is mutations"""
		RNASeqs = {
			'rrlA': '',
			'rrfA': '',
			'rrsA': '',
			'rrlB': '',
			'rrfB': '',
			'rrsB': '',
			'rrlC': '',
			'rrfC': '',
			'rrsC': '',
			'rrlD': '',
			'rrfD': '',
			'rrsD': '',
			'rrlE': '',
			'rrfE': '',
			'rrsE': '',
			'rrlG': '',
			'rrfG': '',
			'rrsG': '',
			'rrlH': '',
			'rrfH': '',
			'rrsH': '',
			'rnpB': 'GAAGCTGACCAGACAGTCGCCGCTTCGTCGTCGTCCTCTTCGGGGGAGACGGGCGGAGGGGAGGAAAGTCCGGGCTCCATAGGGCAGGGTGCCAGGTAACGCCTGGGGGGGAAACCCACGACCAGTGCAACAGAGAGCAAACCGCCGATGGCCCGCGCAAGCGGGATCAGGTAAGGGTGAAAGGGTGCGGTAAGAGCGCACCGCGCGGCTGGTAACAGTCCGTGGCACGGTAAACTCCACCCGGAGCAAGGCCAAATAGGGGTTCATAAGGTACGGCCCGTACTGAACCCGGGTAGGCTGCTTGAGCCAGTGAGCGATTGCTGGCCTAGATGAATGACTGTCCACGACAGAACCCGGCTTATCGGTCAGTTTCACCT',
			'cysT': 'GGCGCGTTAACAAAGCGGTTATGTAGCGGATTGCAAATCCGTCTAGTCCGGTTCGACTCCGGAACGCGCCTCCA',
			'aspT': 'GGAGCGGTAGTTCAGTCGGTTAGAATACCTGCCTGTCACGCAGGGGGTCGCGGGTTCGAGTCCCGTCCGTTCCGCCA',
			'aspU': 'GGAGCGGTAGTTCAGTCGGTTAGAATACCTGCCTGTCACGCAGGGGGTCGCGGGTTCGAGTCCCGTCCGTTCCGCCA',
			'aspV': 'GGAGCGGTAGTTCAGTCGGTTAGAATACCTGCCTGTCACGCAGGGGGTCGCGGGTTCGAGTCCCGTCCGTTCCGCCA',
			'serV': 'GGTGAGGTGGCCGAGAGGCTGAAGGCGCTCCCCTGCTAAGGGAGTATGCGGTCAAAAGCTGCATCCGGGGTTCGAATCCCCGCCTCACCGCCA',
			'serW': 'GGTGAGGTGTCCGAGTGGCTGAAGGAGCACGCCTGGAAAGTGTGTATACGGCAACGTATCGGGGGTTCGAATCCCCCCCTCACCGCCA',
			'serX': 'GGTGAGGTGTCCGAGTGGCTGAAGGAGCACGCCTGGAAAGTGTGTATACGGCAACGTATCGGGGGTTCGAATCCCCCCCTCACCGCCA',
			'serT': 'GGAAGTGTGGCCGAGCGGTTGAAGGCACCGGTCTTGAAAACCGGCGACCCGAAAGGGTTCCAGAGTTCGAATCTCTGCGCTTCCGCCA',
			'serU': 'GGAGAGATGCCGGAGCGGCTGAACGGACCGGTCTCGAAAACCGGAGTAGGGGCAACTCTACCGGGGGTTCAAATCCCCCTCTCTCCGCCA',
			'glnV': 'TGGGGTATCGCCAAGCGGTAAGGCACCGGATTCTGATTCCGGCATTCCGAGGTTCGAATCCTCGTACCCCAGCCA',
			'glnX': 'TGGGGTATCGCCAAGCGGTAAGGCACCGGATTCTGATTCCGGCATTCCGAGGTTCGAATCCTCGTACCCCAGCCA',
			'glnU': 'TGGGGTATCGCCAAGCGGTAAGGCACCGGTTTTTGATACCGGCATTCCCTGGTTCGAATCCAGGTACCCCAGCCA',
			'glnW': 'TGGGGTATCGCCAAGCGGTAAGGCACCGGTTTTTGATACCGGCATTCCCTGGTTCGAATCCAGGTACCCCAGCCA',
			'metV': 'CGCGGGGTGGAGCAGCCTGGTAGCTCGTCGGGCTCATAACCCGAAGGTCGTCGGTTCAAATCCGGCCCCCGCAACCA',
			'metW': 'CGCGGGGTGGAGCAGCCTGGTAGCTCGTCGGGCTCATAACCCGAAGGTCGTCGGTTCAAATCCGGCCCCCGCAACCA',
			'metY': 'CGCGGGGTGGAGCAGCCTGGTAGCTCGTCGGGCTCATAACCCGAAGATCGTCGGTTCAAATCCGGCCCCCGCAACCA',
			'metZ': 'CGCGGGGTGGAGCAGCCTGGTAGCTCGTCGGGCTCATAACCCGAAGGTCGTCGGTTCAAATCCGGCCCCCGCAACCA',
			'metT': 'GGCTACGTAGCTCAGTTGGTTAGAGCACATCACTCATAATGATGGGGTCACAGGTTCGAATCCCGTCGTAGCCACCA',
			'metU': 'GGCTACGTAGCTCAGTTGGTTAGAGCACATCACTCATAATGATGGGGTCACAGGTTCGAATCCCGTCGTAGCCACCA',
			'asnT': 'TCCTCTGTAGTTCAGTCGGTAGAACGGCGGACTGTTAATCCGTATGTCACTGGTTCGAGTCCAGTCAGAGGAGCCA',
			'asnU': 'TCCTCTGTAGTTCAGTCGGTAGAACGGCGGACTGTTAATCCGTATGTCACTGGTTCGAGTCCAGTCAGAGGAGCCA',
			'asnV': 'TCCTCTGTAGTTCAGTCGGTAGAACGGCGGACTGTTAATCCGTATGTCACTGGTTCGAGTCCAGTCAGAGGAGCCA',
			'asnW': 'TCCTCTGTAGTTCAGTCGGTAGAACGGCGGACTGTTAATCCGTATGTCACTGGTTCGAGTCCAGTCAGAGGAGCCA',
			'proK': 'CGGTGATTGGCGCAGCCTGGTAGCGCACTTCGTTCGGGACGAAGGGGTCGGAGGTTCGAATCCTCTATCACCGACCA',
			'proM': 'CGGCGAGTAGCGCAGCTTGGTAGCGCAACTGGTTTGGGACCAGTGGGTCGGAGGTTCGAATCCTCTCTCGCCGACCA',
			'proL': 'CGGCACGTAGCGCAGCCTGGTAGCGCACCGTCATGGGGTGTCGGGGGTCGGAGGTTCAAATCCTCTCGTGCCGACCA',
			'lysQ': 'GGGTCGTTAGCTCAGTTGGTAGAGCAGTTGACTTTTAATCAATTGGTCGCAGGTTCGAATCCTGCACGACCCACCA',
			'lysT': 'GGGTCGTTAGCTCAGTTGGTAGAGCAGTTGACTTTTAATCAATTGGTCGCAGGTTCGAATCCTGCACGACCCACCA',
			'lysV': 'GGGTCGTTAGCTCAGTTGGTAGAGCAGTTGACTTTTAATCAATTGGTCGCAGGTTCGAATCCTGCACGACCCACCA',
			'lysW': 'GGGTCGTTAGCTCAGTTGGTAGAGCAGTTGACTTTTAATCAATTGGTCGCAGGTTCGAATCCTGCACGACCCACCA',
			'lysY': 'GGGTCGTTAGCTCAGTTGGTAGAGCAGTTGACTTTTAATCAATTGGTCGCAGGTTCGAATCCTGCACGACCCACCA',
			'lysZ': 'GGGTCGTTAGCTCAGTTGGTAGAGCAGTTGACTTTTAATCAATTGGTCGCAGGTTCGAATCCTGCACGACCCACCA',
			'thrT': 'GCTGATATAGCTCAGTTGGTAGAGCGCACCCTTGGTAAGGGTGAGGTCGGCAGTTCGAATCTGCCTATCAGCACCA',
			'thrV': 'GCTGATATGGCTCAGTTGGTAGAGCGCACCCTTGGTAAGGGTGAGGTCCCCAGTTCGACTCTGGGTATCAGCACCA',
			'thrW': 'GCCGATATAGCTCAGTTGGTAGAGCAGCGCATTCGTAATGCGAAGGTCGTAGGTTCGACTCCTATTATCGGCACCA',
			'thrU': 'GCCGACTTAGCTCAGTAGGTAGAGCAACTGACTTGTAATCAGTAGGTCACCAGTTCGATTCCGGTAGTCGGCACCA',
			'pheU': 'GCCCGGATAGCTCAGTCGGTAGAGCAGGGGATTGAAAATCCCCGTGTCCTTGGTTCGATTCCGAGTCCGGGCACCA',
			'pheV': 'GCCCGGATAGCTCAGTCGGTAGAGCAGGGGATTGAAAATCCCCGTGTCCTTGGTTCGATTCCGAGTCCGGGCACCA',
			'alaT': 'GGGGCTATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCATAGCTCCACCA',
			'alaU': 'GGGGCTATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCATAGCTCCACCA',
			'alaV': 'GGGGCTATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCATAGCTCCACCA',
			'alaW': 'GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA',
			'alaX': 'GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA',
			'glyV': 'GCGGGAATAGCTCAGTTGGTAGAGCACGACCTTGCCAAGGTCGGGGTCGCGAGTTCGAGTCTCGTTTCCCGCTCCA',
			'glyW': 'GCGGGAATAGCTCAGTTGGTAGAGCACGACCTTGCCAAGGTCGGGGTCGCGAGTTCGAGTCTCGTTTCCCGCTCCA',
			'glyX': 'GCGGGAATAGCTCAGTTGGTAGAGCACGACCTTGCCAAGGTCGGGGTCGCGAGTTCGAGTCTCGTTTCCCGCTCCA',
			'glyY': 'GCGGGAATAGCTCAGTTGGTAGAGCACGACCTTGCCAAGGTCGGGGTCGCGAGTTCGAGTCTCGTTTCCCGCTCCA',
			'glyU': 'GCGGGCGTAGTTCAATGGTAGAACGAGAGCTTCCCAAGCTCTATACGAGGGTTCGATTCCCTTCGCCCGCTCCA',
			'glyT': 'GCGGGCATCGTATAATGGCTATTACCTCAGCCTTCCAAGCTGATGATGCGGGTTCGATTCCCGCTGCCCGCTCCA',
			'ileT': 'AGGCTTGTAGCTCAGGTGGTTAGAGCGCACCCCTGATAAGGGTGAGGTCGGTGGTTCAAGTCCACTCAGGCCTACCA',
			'ileU': 'AGGCTTGTAGCTCAGGTGGTTAGAGCGCACCCCTGATAAGGGTGAGGTCGGTGGTTCAAGTCCACTCAGGCCTACCA',
			'ileV': 'AGGCTTGTAGCTCAGGTGGTTAGAGCGCACCCCTGATAAGGGTGAGGTCGGTGGTTCAAGTCCACTCAGGCCTACCA',
			'leuP': 'GCGAAGGTGGCGGAATTGGTAGACGCGCTAGCTTCAGGTGTTAGTGTTCTTACGGACGTGGGGGTTCAAGTCCCCCCCCTCGCACCA',
			'leuQ': 'GCGAAGGTGGCGGAATTGGTAGACGCGCTAGCTTCAGGTGTTAGTGTCCTTACGGACGTGGGGGTTCAAGTCCCCCCCCTCGCACCA',
			'leuT': 'GCGAAGGTGGCGGAATTGGTAGACGCGCTAGCTTCAGGTGTTAGTGTCCTTACGGACGTGGGGGTTCAAGTCCCCCCCCTCGCACCA',
			'leuV': 'GCGAAGGTGGCGGAATTGGTAGACGCGCTAGCTTCAGGTGTTAGTGTCCTTACGGACGTGGGGGTTCAAGTCCCCCCCCTCGCACCA',
			'leuZ': 'GCCCGGATGGTGGAATCGGTAGACACAAGGGATTTAAAATCCCTCGGCGTTCGCGCTGTGCGGGTTCAAGTCCCGCTCCGGGTACCA',
			'leuX': 'GCCGAAGTGGCGAAATCGGTAGACGCAGTTGATTCAAAATCAACCGTAGAAATACGTGCCGGTTCGAGTCCGGCCTTCGGCACCA',
			'leuU': 'GCCGAGGTGGTGGAATTGGTAGACACGCTACCTTGAGGTGGTAGTGCCCAATAGGGCTTACGGGTTCAAGTCCCGTCCTCGGTACCA',
			'leuW': 'GCGGGAGTGGCGAAATTGGTAGACGCACCAGATTTAGGTTCTGGCGCCGCAAGGTGTGCGAGTTCAAGTCTCGCCTCCCGCACCA',
			'hisR': 'GGTGGCTATAGCTCAGTTGGTAGAGCCCTGGATTGTGATTCCAGTTGTCGTGGGTTCGAATCCCATTAGCCACCCCA',
			'argQ': 'GCATCCGTAGCTCAGCTGGATAGAGTACTCGGCTACGAACCGAGCGGTCGGAGGTTCGAATCCTCCCGGATGCACCA',
			'argV': 'GCATCCGTAGCTCAGCTGGATAGAGTACTCGGCTACGAACCGAGCGGTCGGAGGTTCGAATCCTCCCGGATGCACCA',
			'argY': 'GCATCCGTAGCTCAGCTGGATAGAGTACTCGGCTACGAACCGAGCGGTCGGAGGTTCGAATCCTCCCGGATGCACCA',
			'argZ': 'GCATCCGTAGCTCAGCTGGATAGAGTACTCGGCTACGAACCGAGCGGTCGGAGGTTCGAATCCTCCCGGATGCACCA',
			'argX': 'GCGCCCGTAGCTCAGCTGGATAGAGCGCTGCCCTCCGGAGGCAGAGGTCTCAGGTTCGAATCCTGTCGGGCGCGCCA',
			'argU': 'GCGCCCTTAGCTCAGTTGGATAGAGCAACGACCTTCTAAGTCGTGGGCCGCAGGTTCGAATCCTGCAGGGCGCGCCA',
			'argW': 'GTCCTCTTAGTTAAATGGATATAACGAGCCCCTCCTAAGGGCTAATTGCAGGTTCGATTCCTGCAGGGGACACCA',
			'trpT': 'AGGGGCGTAGTTCAATTGGTAGAGCACCGGTCTCCAAAACCGGGTGTTGGGAGTTCGAGTCTCTCCGCCCCTGCCA',
			'valT': 'GGGTGATTAGCTCAGCTGGGAGAGCACCTCCCTTACAAGGAGGGGGTCGGCGGTTCGATCCCGTCATCACCCACCA',
			'valU': 'GGGTGATTAGCTCAGCTGGGAGAGCACCTCCCTTACAAGGAGGGGGTCGGCGGTTCGATCCCGTCATCACCCACCA',
			'valX': 'GGGTGATTAGCTCAGCTGGGAGAGCACCTCCCTTACAAGGAGGGGGTCGGCGGTTCGATCCCGTCATCACCCACCA',
			'valY': 'GGGTGATTAGCTCAGCTGGGAGAGCACCTCCCTTACAAGGAGGGGGTCGGCGGTTCGATCCCGTCATCACCCACCA',
			'valZ': 'GGGTGATTAGCTCAGCTGGGAGAGCACCTCCCTTACAAGGAGGGGGTCGGCGGTTCGATCCCGTCATCACCCACCA',
			'valW': 'GCGTCCGTAGCTCAGTTGGTTAGAGCACCACCTTGACATGGTGGGGGTCGGTGGTTCGAGTCCACTCGGACGCACCA',
			'valV': 'GCGTTCATAGCTCAGTTGGTTAGAGCACCACCTTGACATGGTGGGGGTCGTTGGTTCGAGTCCAATTGAACGCACCA',
			'gltT': 'GTCCCCTTCGTCTAGAGGCCCAGGACACCGCCCTTTCACGGCGGTAACAGGGGTTCGAATCCCCTAGGGGACGCCA',
			'gltU': 'GTCCCCTTCGTCTAGAGGCCCAGGACACCGCCCTTTCACGGCGGTAACAGGGGTTCGAATCCCCTAGGGGACGCCA',
			'gltV': 'GTCCCCTTCGTCTAGAGGCCCAGGACACCGCCCTTTCACGGCGGTAACAGGGGTTCGAATCCCCTAGGGGACGCCA',
			'gltW': 'GTCCCCTTCGTCTAGAGGCCCAGGACACCGCCCTTTCACGGCGGTAACAGGGGTTCGAATCCCCTAGGGGACGCCA',
			'tyrT': 'GGTGGGGTTCCCGAGCGGCCAAAGGGAGCAGACTGTAAATCTGCCGTCATCGACTTCGAAGGTTCGAATCCTTCCCCCACCACCA',
			'tyrU': 'GGTGGGGTTCCCGAGCGGCCAAAGGGAGCAGACTGTAAATCTGCCGTCACAGACTTCGAAGGTTCGAATCCTTCCCCCACCACCA',
			'tyrV': 'GGTGGGGTTCCCGAGCGGCCAAAGGGAGCAGACTGTAAATCTGCCGTCATCGACTTCGAAGGTTCGAATCCTTCCCCCACCACCA'
		}
		if gName in RNASeqs:
			changes=[]
			m = 0
			mySeq=SeqRecord(Seq(inSeq))
		else:
			changes=[]

			# if ('ACCAATTGCAG' in inSeq):
			# 	inSeq=inSeq[:inSeq.find('ACCAATTGCAG')]+'ACGAATTGCAG'+inSeq[inSeq.find('ACCAATTGCAG')+len('ACCAATTGCAG'):]
			# elif ('ACCAACTGCAG' in inSeq):
			# 	inSeq=inSeq[:inSeq.find('ACCAACTGCAG')]+'ACGAATTGCAG'+inSeq[inSeq.find('ACCAACTGCAG')+len('ACCAACTGCAG'):]
			mySeq=SeqRecord(Seq(inSeq))
			m = 1

		while m == 1:
			stpCdn = str(mySeq.seq[len(mySeq.seq)-3:len(mySeq.seq)])
			m = 0
			if str(mySeq.seq[0:3]) != 'ATG':
				A = 'Start Codon: '+str(mySeq.seq[0:3])+'1'+'ATG'
				changes.append(A)
				mySeq.seq='ATG'+mySeq.seq[3:]
				# recSeq.seq='ATG'+recSeq.seq[3:]
			if 'N' in stnds:
				ndeInst = Restriction.NdeI.search(mySeq.seq)
				for i in ndeInst:
					m = 1
					j = (i - 3) % 3
					k = i - 2
					if j == 0:
						k = k + 2
						A = "NdeI: T"+str(k)+"C"
						mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK == "CT" or stK == "CC"):
							A = "NdeI: C"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
							changes.append(A)
						A = "NdeI: A"+str(k)+"T"
						mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if (stK == "C" or stK == "A" or stK == "G"):
							A = "NdeI: A"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						elif stK == "T":
							A = "NdeI: TCA"+str(k-5)+"AGC"
							mySeq.seq = mySeq.seq[:k-6]+"AGC"+mySeq.seq[k-3:]
						else:
							A = "NdeI: T"+str(k)+"C"
							mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					changes.append(A)
			if 'BglB' in stnds:
				xhoInst = Restriction.XhoI.search(mySeq.seq)
				for i in xhoInst:
					m = 1
					j = (i - 2) % 3
					k = i - 1
					if j == 0:
						k = k + 2
						A = "XhoI: C"+str(k)+"G"
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						stJ = str(mySeq.seq[k+2])
						if (stJ == "A" or stJ == "G"):
							A = "XhoI: AG"+stJ+str(k+1)+"CGC"
							mySeq.seq = mySeq.seq[:k]+"CGC"+mySeq.seq[k+3:]
						elif (stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC"):
							A = "XhoI: C"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						else:
							A = "XhoI: TCG"+str(k-2)+"AGC"
							mySeq.seq = mySeq.seq[:k-3]+"AGC"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if stK == "C":
							A = "XhoI: T"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
							changes.append(A)
						A = "XhoI: A"+str(k)+"C"
						mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					changes.append(A)
			if (('BioB' in stnds) or ('BglB' in stnds)):
				ecoInst = Restriction.EcoRI.search(mySeq.seq)
				for i in ecoInst:
					m = 1
					j = (i - 2) % 3
					k = i - 1
					if j == 0:
						k = k + 2
						A = "EcoRI: A"+str(k)+"G"
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if stK == "AA":
							A = "EcoRI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GA":
							A = "EcoRI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GG":
							A = "EcoRI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "EcoRI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
						elif stK == "CG":
							A = "EcoRI: CGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						else:
							A = "EcoRI: T"+str(k)+"C"
							mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if stK == "C" or stK == "A":
							A = "EcoRI: "+stK+"GA"+str(k-5)+"CGT"
							mySeq.seq = mySeq.seq[:k-6]+"CGT"+mySeq.seq[k-3:]
						elif stK == "G":
							A = "EcoRI: A"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "T":
							A = "EcoRI: G"+str(k-4)+"A"
							mySeq.seq = mySeq.seq[:k-6]+"TAA"+mySeq.seq[k-3:]
						else:
							A = "EcoRI: T"+str(k)+"C"
							mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					changes.append(A)
			if 'BioB' in stnds:
				xbaInst = Restriction.XbaI.search(mySeq.seq)
				for i in xbaInst:
					m = 1
					j = (i - 2) % 3
					k = i - 1
					if j == 0:
						k = k + 3
						A = "XbaI: AGA"+str(k)+"CGT"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"CGT"+mySeq.seq[k+2:]
					elif j == 2:
						k = k + 3
						A = "XbaI: A"+str(k)+"G"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						A = "XbaI: G"+str(k)+"A"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"A"+mySeq.seq[k:]
				speInst = Restriction.SpeI.search(mySeq.seq)
				for i in speInst:
					m = 1
					j = (i - 2) % 3
					k = i - 1
					if j == 0:
						k = k + 2
						A = "SpeI: T"+str(k)+"C"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						A = "SpeI: A"+str(k)+"G"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						A = "SpeI: G"+str(k)+"A"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"A"+mySeq.seq[k:]
				pstInst = Restriction.PstI.search(mySeq.seq)
				for i in pstInst:
					m = 1
					j = (i - 6) % 3
					k = i - 5
					if j == 0:
						k = k + 5
						stK = str(mySeq.seq[k-1:k+5])
						if (stK == "GGATCT" or stK == "GTGCAT"):
							A = "PstI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						else:
							A = "PstI: G"+str(k)+"A"
							mySeq.seq = mySeq.seq[:k-1]+"A"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						stJ = str(mySeq.seq[k-8:k-4])
						if stJ == "GAAT" or stJ == "AAAT":
							A = "PstI: C"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC":
							A = "PstI: C"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						elif stK == "AG" or stK=="AC" or stK=="CA":
							A = "PstI: C"+str(k)+"T"
							mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
						else:
							A = "PstI: C"+str(k-3)+"T"
							mySeq.seq = mySeq.seq[:k-4]+"T"+mySeq.seq[k-3:]
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if stK == "C" or stK=="G":
							A = "PstI: T"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						else:
							A = "PstI: A"+str(k)+"G"
							mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					changes.append(A)
				# mfeInst = MfeI.search(mySeq.seq)
				# for i in mfeInst:
				# 	m = 1
				# 	j = (i - 2) % 3
				# 	k = i - 1
				# 	if j == 0:
				# 		k = k + 3
				# 		A = "MfeI: T"+str(k)+"C"
				# 		mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
				# 	elif j == 2:
				# 		k = k + 3
				# 		stK = str(mySeq.seq[k-6:k-4])
				# 		if (stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC"):
				# 			A = "MfeI: C"+str(k-3)+"G"
				# 			mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
				# 		elif (stK == "AG" or stK == "AC" or stK == "CA"):
				# 			A = "MfeI: T"+str(k)+"C"
				# 			mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
				# 		else:
				# 			A = "MfeI: C"+str(k-3)+"T"
				# 			mySeq.seq = mySeq.seq[:k-4]+"T"+mySeq.seq[k-3:]
				# 	elif j == 1:
				# 		k = k + 4
				# 		stK = str(mySeq.seq[k-6])
				# 		if (stK == "C" or stK == "A" or stK == "G"):
				# 			A = "MfeI: A"+str(k-3)+"G"
				# 			mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
				# 		elif stK == "T":
				# 			A = "MfeI: TCA"+str(k-5)+"AGC"
				# 			mySeq.seq = mySeq.seq[:k-6]+"AGC"+mySeq.seq[k-3:]
				# 		else:
				# 			A = "MfeI: T"+str(k)+"C"
				# 			mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
				# 	changes.append(A)
				avrInst = Restriction.AvrII.search(mySeq.seq)
				for i in avrInst:
					m = 1
					j = (i - 2) % 3
					k = i - 1
					if j == 0:
						k = k + 3
						A = "AvrII: AGG"+str(k)+"CGT"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"CGT"+mySeq.seq[k+2:]
					elif j == 2:
						k = k + 3
						A = "AvrII: A"+str(k)+"G"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						A = "AvrII: G"+str(k)+"A"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"A"+mySeq.seq[k:]
				nheInst = Restriction.NheI.search(mySeq.seq)
				for i in nheInst:
					m = 1
					j = (i - 2) % 3
					k = i - 1
					if j == 0:
						k = k + 2
						A = "NheI: T"+str(k)+"G"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						A = "NheI: A"+str(k)+"G"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						A = "NheI: G"+str(k)+"A"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"A"+mySeq.seq[k:]
				nsiInst = Restriction.NsiI.search(mySeq.seq)
				for i in nsiInst:
					m = 1
					j = (i - 6) % 3
					k = i - 5
					if j == 0:
						k = k + 5
						A = "NsiI: T"+str(k)+"C"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK=="CA" or stK=="GC" or stK=="CT" or stK=="CC" or stK=="AC"):
							A = "NsiI: A"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						elif stK == "GT":
							A = "NsiI: A"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						elif (stK == "AT" or stK == "CG"):
							A = "NsiI: A"+str(k-3)+"T"
							mySeq.seq = mySeq.seq[:k-4]+"T"+mySeq.seq[k-3:]
						elif (stK == "GG"):
							A = "NsiI: A"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "NsiI: AGA"+str(k-3)+"CGT"
							mySeq.seq = mySeq.seq[:k-6]+"CGT"+mySeq.seq[k-3:]
						elif stK == "TT":
							A = "NsiI: TTA"+str(k-3)+"CTG"
							mySeq.seq = mySeq.seq[:k-6]+"CTG"+mySeq.seq[k-3:]
						elif stK == "TC":
							A = "NsiI: TCA"+str(k-3)+"AGC"
							mySeq.seq = mySeq.seq[:k-6]+"AGC"+mySeq.seq[k-3:]
						else:
							A = "NsiI: C"+str(k)+"T"
							mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
						changes.append(A)
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k+1:k+3])
						if (stK == "TA" or stK == "TG"):
							A = "NsiI: T"+stK+str(k+1)+"CTG"
							mySeq.seq = mySeq.seq[:k]+"CTG"+mySeq.seq[k+3:]
						else:
							A = "NsiI: A"+str(k)+"G"
							mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
						changes.append(A)
			if 'BglB' in stnds:
				bglInst = Restriction.BglII.search(mySeq.seq)
				for i in bglInst:
					m = 1
					j = (i - 2) % 3
					k = i - 1
					if j == 0:
						k = k + 0
						A = "BglII: AGA"+str(k)+"CGT"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"CGT"+mySeq.seq[k+2:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK=="CA" or stK=="GC" or stK=="CT" or stK=="CC" or stK=="AC"):
							A = "BglII: A"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						elif stK == "GT":
							A = "BglII: A"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						elif (stK == "AT" or stK == "CG"):
							A = "BglII: A"+str(k-3)+"T"
							mySeq.seq = mySeq.seq[:k-4]+"T"+mySeq.seq[k-3:]
						elif (stK == "GG"):
							A = "BglII: A"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "BglII: AGA"+str(k-3)+"CGT"
							mySeq.seq = mySeq.seq[:k-6]+"CGT"+mySeq.seq[k-3:]
						elif stK == "TT":
							A = "BglII: TTA"+str(k-3)+"CTG"
							mySeq.seq = mySeq.seq[:k-6]+"CTG"+mySeq.seq[k-3:]
						elif stK == "TC":
							A = "BglII: TCA"+str(k-3)+"AGC"
							mySeq.seq = mySeq.seq[:k-6]+"AGC"+mySeq.seq[k-3:]
						else:
							A = "BglII: T"+str(k)+"C"
							mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
						changes.append(A)
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k+1:k+3])
						if (stK == "TA" or stK == "TG"):
							A = "BglII: T"+stK+str(k+1)+"CTG"
							mySeq.seq = mySeq.seq[:k]+"CTG"+mySeq.seq[k+3:]
						else:
							A = "BglII: C"+str(k)+"T"
							mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
						changes.append(A)
				bamInst = Restriction.BamHI.search(mySeq.seq)
				for i in bamInst:
					m = 1
					j = (i - 2) % 3
					k = i - 1
					if j == 0:
						k = k + 2
						A = "BamHI: A"+str(k)+"C"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK == "TT" or stK == "TA" or stK == "AA" or stK == "GA"):
							A = "BamHI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif (stK == "CG" or stK == "AC" or stK == "GC" or stK == "GG"):
							A = "BamHI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "GT":
							A = "BamHI: G"+str(k-3)+"T"
							mySeq.seq = mySeq.seq[:k-4]+"T"+mySeq.seq[k-3:]
						elif stK == "TC":
							A = "BamHI: TCG"+str(k-5)+"AGC"
							mySeq.seq = mySeq.seq[:k-6]+"AGC"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "BamHI: AGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
						else:
							A = "BamHI: T"+str(k)+"C"
							mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
						changes.append(A)
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if (stK == "C" or stK == "A"):
							A = "BamHI: "+stK+"GG"+str(k-5)+"CGT"
							mySeq.seq = mySeq.seq[:k-6]+"CGT"+mySeq.seq[k-3:]
						elif stK=="G":
							A = "BamHI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						else:
							A = "BamHI: C"+str(k)+"T"
							mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
						changes.append(A)
			if 'BioB' in stnds:
				notInst = Restriction.NotI.search(mySeq.seq)
				for i in notInst:
					m = 1
					j = (i - 3) % 3
					k = i - 2
					if j == 0:
						k = k + 5
						A = "NotI: C"+str(k)+"G"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK == "AA" or stK == "GA"):
							A = "NotI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GG":
							A = "NotI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "NotI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
						elif stK == "CG":
							A = "NotI: CGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						else:
							A = "NotI: G"+str(k)+"C"
							mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
						changes.append(A)
					elif j == 1:
						k = k + 4
						A = "NotI: C"+str(k-3)+"T"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-4]+"T"+mySeq.seq[k-3:]
				apoInst = Restriction.ApoI.search(mySeq.seq)
				for i in apoInst:
					m = 1
					j = (i - 2) % 3
					k = i - 1
					if j == 0:
						k = k + 2
						A = "ApoI: A"+str(k)+"G"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						A = "ApoI: T"+str(k)+"C"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						A = "ApoI: T"+str(k)+"C"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
			if 'MoClo' in stnds:
				bbsInst = mySeq.seq.find("GAAGAC")
				if bbsInst > 0:
					m = 1
					j = bbsInst % 3
					k = bbsInst + 1
					if j == 0:
						k = k + 2
						A = "BbsI: A"+str(k)+"G"
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK == "AG"):
							A = "BbsI: AGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
							changes.append(A)
						A = "BbsI: G"+str(k)+"A"
						mySeq.seq = mySeq.seq[:k-1]+"A"+mySeq.seq[k:]
					elif j == 1:
						k = k + 2
						A = "BbsI: AGA"+str(k)+"CGT"
						mySeq.seq = mySeq.seq[:k-1]+"CGT"+mySeq.seq[k+2:]
					changes.append(A)
					bbsInst = 0
				bbsInst = mySeq.seq.find("GTCTTC")
				if bbsInst > 0:
					m = 1
					j = bbsInst % 3
					k = bbsInst + 1
					if j == 0:
						k = k + 2
						A = "BbsI: C"+str(k)+"G"
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if stK == "AA":
							A = "BbsI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GA":
							A = "BbsI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GG":
							A = "BbsI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "BbsI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
						elif stK == "CG":
							A = "BbsI: CGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						else:
							A = "BbsI: TCT"+str(k-2)+"AGC"
							mySeq.seq = mySeq.seq[:k-3]+"AGC"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						A = "BbsI: T"+str(k)+"G"
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					changes.append(A)
					bbsInst = 0
				bsaInst = mySeq.seq.find("GGTCTC")
				if bsaInst > 0:
					m = 1
					j = bsaInst % 3
					k = bsaInst + 1
					if j == 0:
						k = k + 5
						A = "BsaI: C"+str(k)+"G"
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if stK == "AA":
							A = "BsaI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GA":
							A = "BsaI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GG":
							A = "BsaI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "BsaI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
						elif stK == "CG":
							A = "BsaI: CGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						else:
							A = "BsaI: C"+str(k)+"G"
							mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if (stK == "C" or stK == "G"):
							A = "BsaI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "A":
							A = "BsaI: AGG"+str(k-5)+"CGT"
							mySeq.seq = mySeq.seq[:k-6]+"CGT"+mySeq.seq[k-3:]
						else:
							A = "BsaI: TCT"+str(k-2)+"AGC"
							mySeq.seq = mySeq.seq[:k-3]+"AGC"+mySeq.seq[k:]
					changes.append(A)
					bsaInst = 0
				bsaInst = mySeq.seq.find("GAGACC")
				if bsaInst > 0:
					m = 1
					j = bsaInst % 3
					k = bsaInst + 1
					if j == 0:
						k = k + 2
						A = "BsaI: G"+str(k)+"A"
						mySeq.seq = mySeq.seq[:k-1]+"A"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if stK == "AA":
							A = "BsaI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
							changes.append(A)
						elif stK == "GA":
							A = "BsaI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
							changes.append(A)
						elif stK == "AG":
							A = "BsaI: AGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
							changes.append(A)
						A = "BsaI: AGA"+str(k-2)+"CGT"
						mySeq.seq = mySeq.seq[:k-3]+"CGT"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if (stK == "C" or stK == "G"):
							A = "BsaI: A"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "A":
							A = "BsaI: AGA"+str(k-5)+"CGT"
							mySeq.seq = mySeq.seq[:k-6]+"CGT"+mySeq.seq[k-3:]
						else:
							A = "BsaI: C"+str(k)+"T"
							mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
					changes.append(A)
					bsaInst = 0
				mlyInst = mySeq.seq.find("GAGTC")
				if mlyInst > 0:
					m = 1
					j = mlyInst % 3
					k = mlyInst + 1
					if j == 0:
						k = k + 2
						A = "MlyI: G"+str(k)+"A"
						mySeq.seq = mySeq.seq[:k-1]+"A"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if stK == "AA":
							A = "MlyI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GA":
							A = "MlyI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GG":
							A = "MlyI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "MlyI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
						elif stK == "CG":
							A = "MlyI: CGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						else:
							A = "MlyI: T"+str(k)+"C"
							mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if stK == "C" or stK == "A":
							A = "MlyI: "+stK+"GA"+str(k-5)+"CGT"
							mySeq.seq = mySeq.seq[:k-6]+"CGT"+mySeq.seq[k-3:]
						elif stK == "G":
							A = "MlyI: A"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "T":
							A = "MlyI: G"+str(k-4)+"A"
							mySeq.seq = mySeq.seq[:k-6]+"TAA"+mySeq.seq[k-3:]
						else:
							A = "MlyI: C"+str(k)+"G"
							mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					changes.append(A)
					mlyInst = 0
				mlyInst = mySeq.seq.find("GACTC")
				if mlyInst > 0:
					m = 1
					j = mlyInst % 3
					k = mlyInst + 1
					if j == 0:
						k = k + 2
						A = "MlyI: C"+str(k)+"T"
						mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if stK == "AA":
							A = "MlyI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GA":
							A = "MlyI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif stK == "GG":
							A = "MlyI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "MlyI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
						elif stK == "CG":
							A = "MlyI: CGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						else:
							A = "MlyI: T"+str(k)+"G"
							mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if stK == "C" or stK == "A":
							A = "MlyI: "+stK+"GA"+str(k-5)+"CGT"
							mySeq.seq = mySeq.seq[:k-6]+"CGT"+mySeq.seq[k-3:]
							changes.append(A)
						A = "MlyI: C"+str(k)+"G"
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					changes.append(A)
					mlyInst = 0
			if 'GB' in stnds:
				bsmInst = mySeq.seq.find("CGTCTC")
				if bsmInst > 0:
					m = 1
					j = bsmInst % 3
					k = bsmInst + 1
					if j == 0:
						k = k + 5
						A = "BsmBI: C"+str(k)+"G"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC"):
							A = "BsmBI: C"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						else:
							A = "BsmBI: C"+str(k)+"G"
							mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
						changes.append(A)
					elif j == 1:
						k = k + 4
						A = "BsmBI: TCT"+str(k)+"AGC"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-3]+"AGC"+mySeq.seq[k:]
					bsmInst = 0
				bsmInst = mySeq.seq.find("GCAGAG")
				if bsmInst > 0:
					m = 1
					j = bsmInst % 3
					k = bsmInst + 1
					if j == 0:
						k = k + 2
						A = "BsmBI: A"+str(k)+"G"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						stJ = str(mySeq.seq[k+2])
						if (stJ == "A" or stJ == "G"):
							A = "BsmBI: AG"+stJ+str(k+1)+"CGT"
							mySeq.seq = mySeq.seq[:k]+"CGT"+mySeq.seq[k+3:]
						elif (stJ == "T"):
							A = "BsmBI: AGT"+str(k+1)+"TCT"
							mySeq.seq = mySeq.seq[:k]+"TCT"+mySeq.seq[k+3:]
						elif (stK == "TT" or stK == "TA" or stK == "AA" or stK == "GA"):
							A = "BsmBI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif (stK == "CG" or stK == "AC" or stK == "GC" or stK == "GG"):
							A = "BsmBI: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "GT":
							A = "BsmBI: G"+str(k-3)+"T"
							mySeq.seq = mySeq.seq[:k-4]+"T"+mySeq.seq[k-3:]
						elif stK == "TC":
							A = "BsmBI: TCG"+str(k-5)+"AGC"
							mySeq.seq = mySeq.seq[:k-6]+"AGC"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "BsmBI: AGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
						else:
							A = "BsmBI: G"+str(k)+"A"
							mySeq.seq = mySeq.seq[:k-1]+"A"+mySeq.seq[k:]
						changes.append(A)
					elif j == 1:
						k = k + 2
						A = "BsmBI: AGA"+str(k)+"CGT"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"CGT"+mySeq.seq[k+2:]
					bsmInst = 0
				btgzInst = mySeq.seq.find("GCGATG")
				if btgzInst > 0:
					m = 1
					j = btgzInst % 3
					k = btgzInst + 1
					if j == 0:
						k = k + 2
						A = "BtgZI: G"+str(k)+"C"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK == "AA" or stK == "GA"):
							A = "BtgZI: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
							changes.append(A)
						elif (stK == "AG"):
							A = "BtgZI: AGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
							changes.append(A)
						A = "BtgZI: A"+str(k)+"T"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						A = "BtgZI: C"+str(k-3)+"T"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-4]+"T"+mySeq.seq[k-3:]
					btgzInst = 0
				btgzInst = mySeq.seq.find("CATCGC")
				if btgzInst > 0:
					m = 1
					j = btgzInst % 3
					k = btgzInst + 1
					if j == 0:
						k = k + 2
						A = "BtgZI: T"+str(k)+"C"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"C"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC"):
							A = "BtgZI: C"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						else:
							A = "BtgZI: C"+str(k)+"T"
							mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
						changes.append(A)
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if stK == "C":
							A = "BtgZI: A"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-6]+"CCG"+mySeq.seq[k-3:]
							changes.append(A)
						A = "BtgZI: TCG"+str(k-2)+"AGC"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-3]+"AGC"+mySeq.seq[k:]
					btgzInst = 0
			if 'chi' in stnds:
				chiInst = mySeq.seq.find("GCTGGTGG")
				if chiInst > 0:
					m = 1
					j = chiInst % 3
					k = chiInst + 1
					if j == 0:
						k = k + 2
						A = "Chi site: T"+str(k)+"G"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK == "TT" or stK == "TA" or stK == "AA" or stK == "GA"):
							A = "Chi site: G"+str(k-3)+"A"
							mySeq.seq = mySeq.seq[:k-4]+"A"+mySeq.seq[k-3:]
						elif (stK == "CG" or stK == "AC" or stK == "GC" or stK == "GG"):
							A = "Chi site: G"+str(k-3)+"C"
							mySeq.seq = mySeq.seq[:k-4]+"C"+mySeq.seq[k-3:]
						elif stK == "GT":
							A = "Chi site: G"+str(k-3)+"T"
							mySeq.seq = mySeq.seq[:k-4]+"T"+mySeq.seq[k-3:]
						elif stK == "TC":
							A = "Chi site: TCG"+str(k-5)+"AGC"
							mySeq.seq = mySeq.seq[:k-6]+"AGC"+mySeq.seq[k-3:]
						elif stK == "AG":
							A = "Chi site: AGG"+str(k-5)+"CGC"
							mySeq.seq = mySeq.seq[:k-6]+"CGC"+mySeq.seq[k-3:]
						else:
							A = "Chi site: G"+str(k+3)+"T"
							mySeq.seq = mySeq.seq[:k+2]+"T"+mySeq.seq[k+3:]
						changes.append(A)
					elif j == 1:
						k = k + 1
						A = "Chi site: C"+str(k)+"T"
						changes.append(A)
						mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
					chiInst = 0
				chiInst = mySeq.seq.find("CCACCAGC")
				if chiInst > 0:
					m = 1
					j = chiInst % 3
					k = chiInst + 1
					if j == 0:
						k = k + 2
						A = "Chi site: A"+str(k)+"G"
						mySeq.seq = mySeq.seq[:k-1]+"G"+mySeq.seq[k:]
					elif j == 2:
						k = k + 3
						stK = str(mySeq.seq[k-6:k-4])
						if (stK == "CT" or stK == "CC" or stK == "GT" or stK == "GC"):
							A = "Chi site: C"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						else:
							A = "Chi site: C"+str(k)+"T"
							mySeq.seq = mySeq.seq[:k-1]+"T"+mySeq.seq[k:]
					elif j == 1:
						k = k + 4
						stK = str(mySeq.seq[k-6])
						if (stK == "C" or stK == "A" or stK == "G"):
							A = "Chi site: C"+str(k-3)+"G"
							mySeq.seq = mySeq.seq[:k-4]+"G"+mySeq.seq[k-3:]
						else:
							A = "Chi site: TCC"+str(k-3)+"AGC"
							mySeq.seq = mySeq.seq[:k-6]+"AGC"+mySeq.seq[k-3:]
					changes.append(A)
					chiInst = 0
			if (stpCdn == 'TGA' or stpCdn == 'TAG'):
				A = 'Stop Codon: '+str(mySeq.seq[len(mySeq.seq)-3:len(mySeq.seq)])+str(len(mySeq.seq)-3)+'TAA'
				changes.append(A)
				mySeq.seq=mySeq.seq[:len(mySeq.seq)-3]+'TAA'
			if (stpCdn != 'TAA' and stpCdn != 'TAG' and stpCdn != 'TGA'):
				A = 'Add Stop Codon: '+str(len(mySeq.seq))+'TAA'
				changes.append(A)
				mySeq.seq=mySeq.seq[:len(mySeq.seq)]+'TAA'
			if len(changes) > (len(mySeq.seq)/5):
				m = 0
				print("Error: Conflicting mutations")
				prntSeq=''
				for l in range(len(mySeq.seq)):
					prntSeq=prntSeq+mySeq.seq[l];
				print(inSeq)
				print(prntSeq)
				for i in changes:
					print(i)
				# print(recSeq.seq)
				sys.exit("Error: Mutation limit exceded")
#		self.muts.append(changes)
		outSeq=''
		for l in range(len(mySeq.seq)):
			outSeq=outSeq+mySeq.seq[l];
		return outSeq
#
	def MinimizeCodonUsage( self, inSeq ):
		"""Minimizes sequences to 2 codons per amino acid, serviced by single tRNAs when possible"""
		RLD = {
			'TGC': 'C', 'TGT': 'C', 'GAT': 'D', 'GAC': 'D',
			'AGC': 'S', 'TCT': 'S', 'AGT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
			'CAG': 'Q', 'CAA': 'Q', 'ATG': 'M', 'AAC': 'N', 'AAT': 'N',
			'CCG': 'P', 'CCA': 'P', 'CCT': 'P', 'CCC': 'P', 'AAA': 'K', 'AAG': 'K',
			'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'ACA': 'T', 'TTT': 'F', 'TTC': 'F',
			'GCG': 'A', 'GCC': 'A', 'GCA': 'A', 'GCT': 'A',
			'GGC': 'G', 'GGT': 'G', 'GGG': 'G', 'GGA': 'G',
			'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
			'CTG': 'L', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
			'CAT': 'H', 'CAC': 'H',
			'CGT': 'R', 'CGC': 'R', 'CGG': 'R', 'CGA': 'R', 'AGA': 'R', 'AGG': 'R',
			'TGG': 'W', 'GTG': 'V', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
			'GAA': 'E', 'GAG': 'E', 'TAT': 'Y', 'TAC': 'Y',
			'TAA': 'STOP', 'TGA': 'STOP', 'TAG': 'STOP'
		};
		ReducedCodonTable = {
			'C': ['TGC', 'TGT'],
			'D': ['GAT', 'GAC'],
			'S': ['AGC', 'AGT'],
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
			'L': ['CTG', 'TTA'],
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
			'S': ['TCT', 'TCC', 'TCA', 'TCG'],
			'Q': [],
			'M': [],
			'N': [],
			'P': ['CCT', 'CCC'],
			'K': [],
			'STOP': ['TGA', 'TAG'],
			'T': ['ACT', 'ACA'],
			'F': [],
			'A': ['GCT'],
			'G': ['GGG', 'GGA'],
			'I': ['ATA'],
			'L': ['TTG', 'CTT', 'CTC', 'CTA'],
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
			'S': [61, 39],
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
			'L': [77, 23],
			'H': [57, 43],
			'R': [50, 50],
			'W': [100],
			'V': [42, 34, 24],
			'E': [68, 32],
			'Y': [59, 41]
		};
		newInSeq=''
		for i in range(0,int(len(inSeq)/3)):
			inCodon=inSeq[i*3]+inSeq[i*3+1]+inSeq[i*3+2]
			inAA = RLD[str(inCodon)];
			outCodon = inCodon
			tmpCList = []
			if ((inAA != 'M') and (inAA != 'W') and (inAA != 'STOP')):
				for i in range(0, len(ReducedCodonTable[inAA])):
					for j in range (0, ReducedCodonStat[inAA][i]):
						tmpCList.append(str(ReducedCodonTable[inAA][i]));
				while outCodon in KilledCodons[inAA]:
					outCodon = random.choice(tmpCList);
			if inAA == 'STOP':
				outCodon = 'TAA'
			newInSeq=newInSeq+outCodon
		outSeq=newInSeq
		return outSeq;
#
	def tlCheck(self, inSeq, outSeq):
		"""Checks the translation of the engineered sequence against the wild-type sequence"""
		myInSeq=SeqRecord(Seq(inSeq))
		myOutSeq=SeqRecord(Seq(outSeq))
		if myInSeq.translate().seq==myOutSeq.translate().seq:
			successFlag=True
		else:
			successFlag=False
		return successFlag
#
	def mutatePromoters(self, inSeq):
		"""Mutates known E. coli promoters by swapping codons if needed"""
		mainlist=[]
		with open('SynColiProMuts.txt') as f:
			mainlist=[list(literal_eval(line)) for line in f]
		promoterDict={}
		for i in mainlist:
			for j in i:
				promoterDict[j.split(': ')[1].split(' to ')[0]]=j.split(' to ')[1]
		promMuts=[]
		outSeq=inSeq
		for i in promoterDict:
			if i in inSeq:
				outSeq=inSeq[:inSeq.find(i)]+promoterDict[i]+inSeq[(inSeq.find(i)+len(i)):]
				if len(inSeq) != len(outSeq):
					print('Error')
#				promMuts.append(gName+": "+i+str(inSeq.find(i))+str(promoterDict[i]))
		return outSeq
#
	def CodonUsage(self, inSeq):
		CodonCount = {
			'TGC': 0, 'TGT': 0, 'GAT': 0, 'GAC': 0,
			'AGC': 0, 'TCT': 0, 'AGT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0,
			'CAG': 0, 'CAA': 0, 'ATG': 0, 'AAC': 0, 'AAT': 0,
			'CCG': 0, 'CCA': 0, 'CCT': 0, 'CCC': 0, 'AAA': 0, 'AAG': 0,
			'ACC': 0, 'ACG': 0, 'ACT': 0, 'ACA': 0, 'TTT': 0, 'TTC': 0,
			'GCG': 0, 'GCC': 0, 'GCA': 0, 'GCT': 0,
			'GGC': 0, 'GGT': 0, 'GGG': 0, 'GGA': 0,
			'ATT': 0, 'ATC': 0, 'ATA': 0,
			'CTG': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 'CTC': 0, 'CTA': 0,
			'CAT': 0, 'CAC': 0,
			'CGT': 0, 'CGC': 0, 'CGG': 0, 'CGA': 0, 'AGA': 0, 'AGG': 0,
			'TGG': 0, 'GTG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
			'GAA': 0, 'GAG': 0, 'TAT': 0, 'TAC': 0,
			'TAA': 0, 'TGA': 0, 'TAG': 0
		};
		for i in range(0,int(len(inSeq)/3)):
			tmpCodon=inSeq[i*3]+inSeq[i*3+1]+inSeq[i*3+2]
			CodonCount[tmpCodon]=CodonCount[tmpCodon]+1
		return CodonCount

s = zerorpc.Server(RefactorRPC())
s.bind("tcp://0.0.0.0:4242")
s.run()


#def statistics():
	# """Runs edit statistics and codon usage"""
	# # fOut=open('CodonUsageTable.txt','w')
	# # fOut.write(json.dumps(CodonCountBefore))
	# # fOut.write("\n\n")
	# # fOut.write(json.dumps(CodonCountAfter))
	# # fOut.close();
	# import matplotlib.pyplot as plt; plt.rcdefaults()
	# import numpy as np
	# import matplotlib.pyplot as plt
	# codons=()
	# usageA=[]
	# usageB=[]
	# muts=[]
	# print(len(promMuts))
	# fIn=open('Mutations.txt','w')
	# for line in muts:
	# 	fIn.write(line+"\n");
	# fIn.close();
	# for i in CodonCountAfter:
	# 	codons=codons+(i,)
	# 	usageA.append(CodonCountAfter[i])
	# 	usageB.append(CodonCountBefore[i])
	# y_pos = np.arange(len(codons))
	# fig, ax = plt.subplots()
	# bar_width = 0.5
	# opacity = 0.8
	# rects1=plt.bar(y_pos, usageB, bar_width, alpha=0.5,label='Before')
	# rects2=plt.bar(y_pos+bar_width, usageA, bar_width, alpha=0.5,label='After')
	# plt.xticks(y_pos+bar_width, codons, rotation=45)
	# plt.ylabel('Number of Instances')
	# plt.title('Codon Usage')
	# plt.legend()
	# plt.tight_layout()
	# plt.show()
	# return
