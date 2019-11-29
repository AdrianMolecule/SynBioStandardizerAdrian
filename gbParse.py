from Bio import SeqIO
import sys
fOut=open(sys.argv[2],'w')
for record in SeqIO.parse(sys.argv[1], "genbank"):
    for feature in record.features:
        if feature.type == "gene":
            try:
                fOut.write(">%s\n%s\n" % (feature.qualifiers['gene'][0], feature.location.extract(record).seq))
            except KeyError as e:
                print("error")
fOut.close()
