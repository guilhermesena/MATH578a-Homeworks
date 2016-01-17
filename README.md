# MATH578a-Homeworks
Solutions for the Computational Molecular Biology course

# About the homeworks

## Homework 1

This was a homework on multiple alignment
**Input**: A FASTA text file with several protein/DNA sequences
**Output**: The score and global alignment of all pairs of sequences using Match = 3, Mismatch = -1, Gap = -3

###How to Run
First build the CPP file:

> g++ -o hw1 hw1.cpp

Then run with the FASTA input as an argument. For example, if your fasta is called **input1.txt**, do

> ./hw1 input1.txt

If you want to save the output file:

> ./hw1 input1.txt >output1.txt
