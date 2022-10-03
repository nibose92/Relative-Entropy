#!/usr/bin/python3

########################
# This script will determine the relative entropy at each position of the multiple sequence by comparing the frequency of a residue in a strain vs. the frequency of that residue in the entire population. The multiple sequence alignment includes sequences of all allelic copies of the gene per strain. Allelic copies from the same strain are represented by the same accession number in the FASTA header, followed by ":" delimiter and an allele-specific designation. In the example multiple sequence alignments used in our study, this allele-specific designation is the allele locus in the strain genome.

# Usage with bash: $python RelativeEntropy.py -in /path/to/multiple_sequence_alignment.fasta -pssm /path/to/PSSM.csv -entropy /path/to/Relative_Entropy.csv

# The input file is a multiple sequence alignment of all sequences, stored as .fasta. The output files are: 1) PSSM.csv -- A csv file showing the position specific residue occurrences, and 2) The computed position-specific relative entropy

# Note that this method only considers nucleotides A, C, G, and T, as well as - (deletions). U and non-conventional nucleotide designations are ignored.

# Requirements for this script to function: 
# 1) python (2.0 or higher)
# 2) NumPy
# 3) pandas
# 4) biopython (1.7 or higher)  

#######################



from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
import pandas as pd
import numpy as np
import math, os, argparse

parser = argparse.ArgumentParser()

parser.add_argument('-in', '--input', dest='input', help='Enter the .fasta file with multiple sequence alignments that will be used to calculate relative entropy')
parser.add_argument('-pssm', '--pssm', dest='pssm', help='Enter the csv output file name for Position Specific Scoring (Frequency) Matrix-- default=./PSSM.csv', default='./PSSM.csv')
parser.add_argument('-entropy', '--entropy', dest='entropy', help='Enter the csv output file name for calculated relative entropy -- default=./Relative_Entropy.csv', default='./Relative_Entropy.csv')

args = parser.parse_args()



####
# Initialize files from argparse and initialize what constitutes nucleotides

inFile = open(args.input, "r")
pssm_csv = open(args.pssm, "w")
Entropy = open(args.entropy, "w")

nucleotides = ['A', 'C', 'G', 'T']




#####
# Obtain the position specific frequency matrix for residues per position in the full set multiple sequence alignment and save it as a csv file  

align_all = AlignIO.read(args.input, "fasta")                                                                                                               # Read alignment file
summary_align1 = AlignInfo.SummaryInfo(align_all)
consensus_all = summary_align1.dumb_consensus()                                                                                                             # Generate consensus for alignment
full_pssm = summary_align1.pos_specific_score_matrix(consensus_all, chars_to_ignore=["Y", "R", "W", "S", "K", "M", "D", "V", "H", "B", "N"])                # Generate position-specific frequency matrix for the full alignment
align_length = align_all.get_alignment_length()


pssm_csv.write(str(full_pssm))                                                                                                                               # Store position-specific frequency matrix for full aligment
pssm_csv.close()




#####
# Write the first row of the relative entropy output file

Entropy.write("ID_#,")
i = 1
while i <= align_length:
    Entropy.write(str(i)+",")
    i+=1
Entropy.write("\n")




#####
# Create a temporary csv file of all Accession numbers, input that file into a pandas dataframe, and deduplicate accessions (per strain). I have chosen to create a temporary file so that in the event of an error due to an unrecognizable character, you can identify which accession to correct 

csv1 = open("./temp1.csv", "w")
align_count_all = 0

for line in inFile:
    if ">" in line:
        csv1.write(line+"\n")
        align_count_all+=1                                                                                                                                  # Variable to obtain number of FASTA sequences
csv1.close()

df1 = pd.read_csv("./temp1.csv", delimiter = "[:,>]", engine='python')
df2 = df1.iloc[:,1]                                                                                                                                         # Extract accessions
df3 = df2.drop_duplicates(keep='first').reset_index(drop=True)                                                                                              # Removing duplicate accessions. reset_index resets that index number system.





#####
# Cycle through each accession, determine position specific frequency matrix for each strain, and then compute relative entropy by comparing frequency of residue across alleles in a strain vs. frequency across all strains  

frequency_ratio = 0.0

for row in df3:
    align_count_strain = 0
    csv2 = open("./temp2.fasta", "w")                                                                                                                       # Create a temporary FASTA file for multiple sequence alignment per strain 
    for seq_record1 in SeqIO.parse(args.input, "fasta"):
        if str(row) in seq_record1.id:
            csv2.write(">"+seq_record1.id+"\n"+str(seq_record1.seq)+"\n")
            align_count_strain+=1
    csv2.close()
    align_strain = AlignIO.read("./temp2.fasta", "fasta")                                                                                                   # Read alignment of sequences from the strain
    summary_align2 = AlignInfo.SummaryInfo(align_strain)
    consensus_strain = summary_align2.dumb_consensus()                                                                                                      # Generate consensus for alignment
    strain_pssm = summary_align2.pos_specific_score_matrix(consensus_strain, chars_to_ignore=["Y", "R", "W", "S", "K", "M", "D", "V", "H", "B", "N"])       # Determine position specific frequency matrix of each residue within the strain
    i = 0                                                                                                                                                   # Count variable to cycle through each position
    Entropy.write(str(row)+",")
    while i < align_length:
        RE = 0.0
        ntd_count_strain = 0.0
        ntd_count_all = 0.0
        for n in nucleotides:                                                                                                                               # Determine relative entropy for nucleotides
            if float(full_pssm[i][n])!=0.0 and float(strain_pssm[i][n])!=0.0:
                frequency_ratio = float(((strain_pssm[i][n])/align_count_strain)/((full_pssm[i][n])/align_count_all))
                RE = round((RE + (((strain_pssm[i][n])/align_count_strain)*(math.log(frequency_ratio, 10)))), 5)
                ntd_count_strain = ntd_count_strain + strain_pssm[i][n]
                ntd_count_all = ntd_count_all + full_pssm[i][n]


        # Determine relative entropy when deletions are encountered
        if ntd_count_strain < align_count_strain and ntd_count_all < align_count_all:
            frequency_ratio = float(((align_count_strain - ntd_count_strain)/align_count_strain)/((align_count_all - ntd_count_all)/align_count_all))
            RE = round((RE + (((align_count_strain - ntd_count_strain)/align_count_strain)*(math.log(frequency_ratio, 10)))), 5)

        Entropy.write(str(RE)+",")
        i+=1
    Entropy.write("\n")


Entropy.close()
print("Complete!")




#####
# Delete temporary files

os.remove("./temp1.csv")
os.remove("./temp2.fasta")

