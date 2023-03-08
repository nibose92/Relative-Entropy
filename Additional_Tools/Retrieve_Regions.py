#!/usr/bin/python

#Make sure you have all input files stored in a single folder -- Like "infile/"


import os, glob, argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument('-path', '--path', dest='path', help='Enter path folder that contains NCBI Assembly feature table files; example: /folder1/folder2/folder3/. Default=./', default='./')
parser.add_argument('-in', '--input', dest='input', nargs='+', help='Enter filename(s). Note that wildcards may be used. Default: *_feature_table.txt', default='*_feature_table.txt')
parser.add_argument('-out', '--output', dest='output', help='Enter output filename with .csv extension. Default = ./output.csv', default='./output.csv')
parser.add_argument('-search', '--search', dest='search', help='Enter .csv file containing list of search phrases', default='./search.csv')
parser.add_argument('-not_search', '--not_search', dest='not_search', help='Enter .csv file containing list of phrases to avoid', default='./not_search.csv')

args = parser.parse_args()
FileList = []

for each in args.input:
    FileList+=glob.glob(args.path + each) #Generate a list of all input files with their paths

outFile = open(args.output, 'w')
searchDF = pd.read_csv(args.search)
NotSearchDF = pd.read_csv(args.not_search)

outFile.write("feature,class,assembly,assembly_unit,seq_type,chromosome,genomic_accession,start,end,strand,product_accession,non-redundant_refseq,related_accession,name,symbol,GeneID,locus_tag,feature_interval_length,product_length,attributes\n")

Count1 = 0  #Count number of targets
Count2 = 0  #Count number of bacteria with target

for each in FileList:
    print("\n\n" + each)
    a = open(each, "r")
    check1 = False  #Check if line contains search term and does not contain terms to avoid
    check2 = False
    for line in a:
        for index1, search in searchDF.iterrows(): 
            if str(search['search_terms']) in line:
                check1 = True
                for index2, not_search in NotSearchDF.iterrows():
                    if str(not_search['not_search_terms']) in line:
                        check1 = False
                if check1 == True:
                    Count1 +=1
                    outFile.write(line.replace('\t', ','))
                    print(line)
                    check2 = True
    if check2 == True:
        Count2 +=1

    a.close()

print("Complete!        Total number of bugs = " + str(Count2) + "      Total number of 16S = " + str(Count1))

outFile.close()

#Getting stats about target

fileDF = pd.read_csv(args.output, delimiter=',')
statsDF = fileDF['end']-fileDF['start']

print("Average target size = " + str(statsDF.mean(axis=0)) + "      Standard deviation of target size = " +str(statsDF.std(axis=0)))
