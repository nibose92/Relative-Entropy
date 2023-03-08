1) Retrieve_Regions.py


Explanation

This script mines through the NCBI Assembly feature table files that you have downloaded (see requirements below for details) and then extracts the gene information that you are looking for, based on keywords in a "search.csv" file, while avoiding keywords in a second "not_search.csv" file. The extracted information gets placed in an output.csv file. Also the script ultimately prints some statistics about the gene regions, such as the total number of genes it has extracted, the average length of the genes it has extracted and the standard deviation of the length.




Requirements for Retrieve_Regions.py to work:

a) NCBI Assembly derived Feature table files of the relevant genomes. To accomplish this, go to the NCBI Assembly website (https://www.ncbi.nlm.nih.gov/assembly), search the organism you are looking for (in my case it was the species "Escherichia coli", or the genus "Escherichia" or "Shigella"), apply necessary filters in the search such as "Complete Genome" and "Latest RefSeq", then select "Download Assemblies" within which select the "Feature table.txt" option.

b) Unpack all the feature tables into a single folder

c) Make sure you have the following tools: glob, and pandas, python (2.0 or higher)

d) Make sure you have two csv files, one that lists what keywords you are looking for, and one that lists keywords to avoid. I have found that having these keywords in lists are much more convenient than typing them into the command line. The default .csv file for search terms I have uploaded ("search.csv") lists the keywords "16S ribosomal RNA" and "ribosomal RNA-16S" because I was wanted to obtain the gene information from the feature table text files corresponding to the 16S rRNA genes. The default .csv file for terms to avoid I have uploaded ("not_search") lists "transferase", "pseudo", and "methyltransferase" because I do not want information about genes encoding enzymes that modify the 16S rRNA. You can download the "search.csv" and the "not_search.csv" files and change the keywords you wish to use, PROVIDED YOU DO NOT CHANGE THE HEADING (row 1) IN EACH FILE!



Example Usage with bash: $python Retrieve_Regions.py -path /path/to/featuretables -in *_feature_table.txt -out output.csv -search search.csv -not_search not_search.csv

Display script help: $python Retrieve_Regions.py -h


NOTE: The output file must have the following information under their respective columns to be able to use the next code efetch_Sequences.py: a) NCBI genomic accession (for example "NZ_..." or "NC..."), b) the gene start site, and c) the gene stop site. 

If you have used NCBI Assembly derived feature tables for complete genomes of your organism(s) of interest, this should have automatically occurred.











2) efetch_Sequences.py

Explanation

The script takes the output.csv file obtained from Retrieve_Regions.py and retrieves the "+" strand sequences from NCBI nucleotide and puts them in a FASTA format. Again, you must have NCBI genomic accessions, gene start site, and gene stop sites listed in order for this to work.
NOTE: This script can take a while and can sometimes disconnect if you are using too much of NCBI's server time. It may be a good idea to split your input files into multiple ones (making sure you retain the first row containing titles in each), obtaining multiple FASTA outputs, and then concatenating them.




Requirements for efetch_Sequences.py

a) conda or anaconda - using conda, install the entrez-direct tool (https://anaconda.org/bioconda/entrez-direct)

b) pandas

c) python (2.0 or higher)

d) biopython (1.7 or higher)

Example Usage with bash: $python efetch_Sequences.py -in /path/to/Retrieve_Regions_outputfile.csv -out /path/to/outputfile.fasta

Display script help: $python efetch_Sequences.py -h


