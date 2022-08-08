from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

fasta = open("/Users/jared/Downloads/Mycobacterium_tuberculosis_H37Rv_genes_v4.fasta")
fastas = {}

for f in SimpleFastaParser(fasta):
	name = ">" + f[0].split("|")[1]
	seq = f[1]
	fastas[name] = seq

with open('/Users/jared/Downloads/Mtb_H37Rv.fa', 'w') as cdna:
	for gene, seq in fastas.items():
		cdna.write(gene + "\n" + seq + "\n")
