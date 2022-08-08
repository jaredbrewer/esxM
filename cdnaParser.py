from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os

script_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(script_path)

fasta = open("Mycobacterium_tuberculosis_H37Rv_genes_v4.fasta")
fastas = {}

for f in SimpleFastaParser(fasta):
	name = ">" + f[0].split("|")[1]
	seq = f[1]
	fastas[name] = seq

# This code block is good for writing out the intermediate, but is not strictly necessary.

# with open('/Users/jared/Downloads/Mtb_H37Rv.fa', 'w') as cdna:
# 	for gene, seq in fastas.items():
# 		cdna.write(gene + "\n" + seq + "\n")

### Everything above is very general in parsing these files for any use, the subsequent lines are to find and replace the esxMJKPW sequences with their corresponding 5' UTR sequences for the purposes of the subsequent alignment. ###

esx = open("esx_5.fa")

esxs = {}

for f in SimpleFastaParser(esx):
	name = ">" + f[0]
	seq = f[1]
	esxs[name] = seq

esx_genes = ["esxM", "esxJ", "esxP", "esxW", "esxK"]

[fastas.pop(key) for key in esx_genes]
new = fastas | esxs

with open('Mtb_H37Rv.fa', 'w') as cdna:
	for gene, seq in new.items():
		cdna.write(gene + "\n" + seq + "\n")


