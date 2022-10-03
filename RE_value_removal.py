#!/usr/bin/python

########################
# This script will truncate a Relative entropy (RE) file for a large population set by removing accessions and their computed RE of a subpopulation that has already been evaluated. This script requires that relative entropy of the larger population (Reference) and the population to be removed (Subpopulation) are separately determined. For example, in the example set, to get the RE for non-coli Escherichia, we removed RE corresponding to E. coli (already computed) from the RE file corresponding to Escherichia.  

# Usage with bash: $python RE_value_removal.py -ref /path/to/Reference_Relative_Entropy.csv -sub /path/to/Subpopulation_Relative_Entropy.csv -out /path/to/Truncated_RE.csv

# The output file generated is the RE for accession numbers in the large population that are not included in the subpopulation.

# Requirements for this script to function: 
# 1) python (2.0 or higher)
# 2) Relative entropy files for the Reference large population and subpopulation (using the RelativeEntropy.py script)
# 3) NumPy
# 4) pandas

#######################


import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-ref', '--reference', dest='reference', help='Enter the name of Reference .csv Relative Entropy file -- This is usually the overall population')
parser.add_argument('-sub', '--subpopulation', dest='subpopulation', help='Enter the name of the .csv file that contains relative entropy of the subpopulation group that will be removed from the overall population')
parser.add_argument('-out', '--output', dest='output', help='Enter the name of .csv output file, default = Truncated_RE.csv', default='./Truncated_RE.csv')

args = parser.parse_args()



#####
# Initialize argparse files and place Input and Reference files into pandas dataframes

SubFile = args.subpopulation
OutFile = args.output
RefFile = args.reference

df1 = pd.read_csv(RefFile)
df2 = pd.read_csv(SubFile)

df1.drop(df1.columns[len(df1.columns)-1], axis=1, inplace=True)
df2.drop(df2.columns[len(df2.columns)-1], axis=1, inplace=True)




#####
# Write the first row of the output file listing residue positions

out = open(args.output, 'w') 

size = df1.shape[1] -1

for col in df1.columns:
    out.write(str(col) + ",")

out.write("\n")




#####
# Removal of RE for already evaluated subpopulation

for index1, row1 in df1.iterrows():
    a = 0
    for index2, row2 in df2.iterrows():                                             # Remove already known REs
        if str(row1['Accession_#']) == str(row2['Accession_#']):
            a = 1

    if a == 0:                                                                      # Keep REs for unevaluated subpopulation
        out.write(str(row1['Accession_#'] + ','))
        for i in range(size):
            out.write(str(row1[str(i+1)]) + ',')
        out.write("\n")

out.close()

print("Complete!")
