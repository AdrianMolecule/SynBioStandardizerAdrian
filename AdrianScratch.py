
#from Bio.Restriction import NHel
from Bio import Restriction
from Bio.Seq import Seq
seq = Seq( 'GGTACCGGGCCCCCCCTCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTCCTG') 
#print(Restriction.EcoRI.search(seq))
#print(Restriction.CommOnly('EcoRI').search(seq))
print(Restriction.AllEnzymes.elements)
print(Restriction.AllEnzymes.get("EcoRI"))
print(Restriction.AllEnzymes.get("EcoRI").site)   # or .elucidate()
print(Restriction.AllEnzymes.get("NheI"))
print(Restriction.AllEnzymes.get("NheI").site)
print(Restriction.AllEnzymes.get("NheI"))
print(Restriction.AllEnzymes.get("NheI"))
print(Restriction.AllEnzymes.get("EcoRI"))
print(Restriction.AllEnzymes.get("EcoRI"))



# from Bio.Restriction import *
# from Bio.Restriction import EcoRI
# EcoRI.search("asdsafdsa")