Files in the directory:

01_intermediateIndex.cpp
this file is to produce intermediate postings and page table
How to run:
first you need to put the dataset in the path of "collection/collection.tsv"
second you need to create a directory "output/"
g++ -o 01_intermediateIndex 01_intermediateIndex.cpp
./01_intermediateIndex


02_mergeIndex.cpp
this file merge the intermediate postings, produce index, lexicon and block metadata file
How to run:
g++ -o 02_mergeIndex 02_mergeIndex.cpp
./02_mergeIndex

03_queryProcessing.cpp
this file is to run query processing. It will prompt the user to input queries and show the result
How to run:
g++ -o 03_queryProcessing 03_queryProcessing.cpp
./03_queryProcessing