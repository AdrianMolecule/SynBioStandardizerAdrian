#!/usr/bin/python
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
# of the authors. 
#
#####
import cgi
import cgitb; cgitb.enable()
import copy
import sys
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from Bio import SeqIO
from Bio.Restriction import *
from math import *
import genestand
#
#
#
print("Content-Type: text/html\n\n")
print("""
<html>
<head>
<style>
   body {
        font-family:courier;
   }
</style>
<title>Synthetic Biology Gene Standardizer</title>
</head>
<body >
<h3>Synthetic Biology Gene Standardizer</h3>
""")
form = cgi.FieldStorage()
stnds = []
i=0
j=0
if "gene" in form:
        organism = form.getvalue("organism")
        usergene = form.getvalue("gene")
        if organism == "coli":
                handle = open("SynColiSeqs.txt", "rU")
                records = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                for p in records:
                        tmpst=str(p.description)
                        tmpg1 = tmpst.split('[gene=', 1);
                        tmpg2 = str(tmpg1[1])
                        tmpg3 = tmpg2.split('] [')
                        p.description = tmpg3[0]
                        gdesc = str(p.description)
                        if gdesc.lower() == usergene.lower():
                                j=i
                                outSeq = copy.deepcopy(p)
                        i=i+1
                if j==0:
                        print("Check gene name - not in database</body></html>")
                        sys.exit("Gene not in Database")
                i=0
                with open("SynColiProMuts.txt") as f:
                        chngtmpP = f.readlines();
                        i=i+1;
                chP1 = chngtmpP[j].split('[\'');
                if len(chP1) == 1:
                        changesP = '0';
                if len(chP1) > 1:
                        chP2 = chP1[1].split('\']');
                        changesP = chP2[0].split('\', \'');
                with open("SynColiMuts.txt") as f:
                        chngtmp = f.readlines()
                        i=i+1
                ch1 = chngtmp[j].split('[\'')
                ch2 = ch1[1].split('\']')
                changes = ch2[0].split('\', \'')
        elif organism == "subtilis":
                handle = open("SynSubtilisSeqs.txt", "rU")
                records = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                for p in records:
                        tmpst=str(p.description)
                        tmpg1 = tmpst.split('[gene=', 1);
                        tmpg2 = str(tmpg1[1])
                        tmpg3 = tmpg2.split('] [')
                        p.description = tmpg3[0]
                        gdesc = str(p.description)
                        if gdesc.lower() == usergene.lower():
                                j=i
                                outSeq = copy.deepcopy(p)
                        i=i+1
                if j==0:
                        print("Check gene name - not in database</body></html>")
                        sys.exit("Gene not in Database")
                i=0
                with open("SynSubtilisMuts.txt") as f:
                        chngtmp = f.readlines()
                        i=i+1
                ch1 = chngtmp[j].split('[\'')
                ch2 = ch1[1].split('\']')
                changes = ch2[0].split('\', \'')
elif "seq" in form:
        organism = "";
        userseq = form.getvalue("seq")
        useq = ''.join(v for v in userseq if v.isalnum())
        xyy = re.search('[^ACGTacgt]', useq)
#               xyy = re.compile('[^ACGT]', re.IGNORECASE)
#               m = xyy.match(useq)
        if xyy:
                print("Only A, C, T, or G in the sequence.")
                sys.exit("ATCG error")
        sseq = Seq(useq)
        record = SeqRecord(sseq)
        usergene = "User"
        record.description = str(usergene)
        outSeq = copy.deepcopy(record)
        record.description = str(usergene)+" Original Sequence<br>\n"
        changes = []
        stnds.append(form.getvalue("N"))
        stnds.append(form.getvalue("BioB"))
        stnds.append(form.getvalue("BglB"))
        stnds.append(form.getvalue("MoClo"))
        stnds.append(form.getvalue("GB"))
        stnds.append(form.getvalue("chi"))
        genestand.refactor(outSeq, record, changes, stnds)
else:
        print("Please go back and fill in either a gene name or sequence</body></html>")
        sys.exit("Please fill in either gene or seq")
print('<br><hr>');
outSeq.description = str(usergene)+"'s Refactored Sequence<br>\n";
print("<h4>Synthetic gene</h4>");
print(outSeq.format("fasta"));
print("<br><hr><br>");
print("<h4>Mutations</h4>");
for xyz in changes:
        print(xyz+"<br>");
if (organism == 'coli'):
        print("<br><h4>E. coli Promoter Mutations:</h4>")
        for xyw in changesP:
                print(xyw+"<br>");
print("""
</body>
</html>
""")
