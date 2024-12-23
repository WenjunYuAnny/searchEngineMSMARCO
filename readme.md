# Search Engine for 8.8M MS MARCO Passage Dataset

This repository is a search engine built on the MS MARCO Passage Dataset using inverted index and BM25 ranking scores.  
The index is compressed with Varbyte.  
It supports both conjunctive (AND) and disjunctive (OR) queries.  
The query processing uses DAAT (Document-At-A-Time).  
It's written in C++ for IO and computation efficiency.  

## Sample Query Result:
### dog and cat and mouse
![image](https://github.com/user-attachments/assets/9240d4b1-ccf5-4d70-a701-fd3c93ae9afa)  
![image](https://github.com/user-attachments/assets/e37fcca2-9342-4d0f-b8fc-a63de38a5922)  

### Dog or cat or armadillo
![image](https://github.com/user-attachments/assets/0a845bbc-8cc4-48b5-8141-92a64f9fd335)  
![image](https://github.com/user-attachments/assets/a0d12707-5801-4be4-82fc-f260eeb97f8a)

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
