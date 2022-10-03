Two codes are listed here:

1) RelativeEntropy.py

This script will determine the relative entropy at each position of the multiple sequence by comparing the frequency of a residue in a strain vs. the frequency of that residue in the entire population. The multiple sequence alignment includes sequences of all allelic copies of the gene per strain. Allelic copies from the same strain are represented by the same accession number in the FASTA header, followed by ":" delimiter and an allele-specific designation. In the example multiple sequence alignments used in our study, this allele-specific designation is the allele locus in the strain genome.


Usage with bash: $python RelativeEntropy.py -in /path/to/multiple_sequence_alignment.fasta -pssm /path/to/PSSM.csv -entropy /path/to/Relative_Entropy.csv

Display script help: $python RelativeEntropy.py -h

The input file is a multiple sequence alignment of all sequences, stored as .fasta. The output files are: 1) PSSM.csv -- A csv file showing the position specific residue occurrences, and 2) The computed position-specific relative entropy

Note that this method only considers nucleotides A, C, G, and T, as well as - (deletions). U and non-conventional nucleotide designations are ignored.



Requirements for this script to function: 
1) python (2.0 or higher)
2) NumPy
3) pandas
4) biopython (1.7 or higher)








2) RE_value_removal.py

This script will truncate a Relative entropy (RE) file for a large population set by removing accessions and their computed RE of a subpopulation that has already been evaluated. This script requires that relative entropy of the larger population (Reference) and the population to be removed (Subpopulation) are separately determined. For example, in the example set, to get the RE for non-coli Escherichia, we removed RE corresponding to E. coli (already computed) from the RE file corresponding to Escherichia.  

Usage with bash: $python RE_value_removal.py -ref /path/to/Reference_Relative_Entropy.csv -sub /path/to/Subpopulation_Relative_Entropy.csv -out /path/to/Truncated_RE.csv

Display script help: $python RE_value_removal.py -h

The output file generated is the RE for accession numbers in the large population that are not included in the subpopulation.



Requirements for this script to function:
1) python (2.0 or higher)
2) Relative entropy files for the Reference large population and subpopulation (using the RelativeEntropy.py script)
3) NumPy
4) pandas
