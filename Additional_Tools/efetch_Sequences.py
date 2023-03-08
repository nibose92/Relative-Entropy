#!/usr/bin/python

from Bio import SeqIO
import pandas as pd
import os, argparse

parser = argparse.ArgumentParser()

parser.add_argument('-in', '--input', dest='input', help='Enter .csv file containing genomic accession, start, and stop positions. example: output.csv from target.py output')
parser.add_argument('-out', '--output', dest='output', help='Enter output filename with .fasta extension, default = ./output.fasta', default='./output.fasta')

args = parser.parse_args()

inDF = pd.read_csv(args.input, delimiter = ',')	#Reading csv file
x=open(args.output, 'w')	#Created and opened an empty fasta file where fasta sequences will be dumped   

for index, row in inDF.iterrows():	#For each indexed row, I will access elements within each row, iterratively going through df, store full rows in  
	
	os.system("efetch -db 'nucleotide' -id "+str(row['genomic_accession'])+" -seq_start "+str(row['start'])+" -seq_stop "+str(row['end'])+" -format fasta > temp.fasta")
	handle = open("temp.fasta", 'r')
	record = SeqIO.read(handle, "fasta")
	handle.close()

	if row['strand'] == "+":
		string1= str(record.id)
		string2= str(record.seq)
		x.write(">"+string1+"\n"+string2+"\n")
		
	elif row['strand'] == "-":
		string1= str(record.id)
		string2= str(record.seq.reverse_complement())
		x.write(">"+string1+"\n"+string2+"\n")
	
	print(">"+string1+"\n"+string2+"\n")
	
x.close()

print("Complete!")

os.remove("./temp.fasta")
