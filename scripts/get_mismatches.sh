#!/bin/bash

seqmap 0 ../data/PNA_sequences/PNA_sequences.fasta \
       ../data/reference_sequences/oligos_cds.fasta \
       ../data/mismatches/mismatches_PNA.dat /output_all_matches\
       /available_memory:200




seqmap 1 ../data/PNA_sequences/PNA_sequences.fasta \
       ../data/reference_sequences/oligos_cds.fasta \
       ../data/mismatches/mismatches_1_PNA.dat /output_all_matches\
       /available_memory:200

awk 'NR == 1; NR > 1 {print $0 | "sort -k4,4 -k1,1"}' ../data/mismatches/mismatches_PNA.dat  > ../data/mismatches/mismatches_PNA_0.dat
awk 'NR == 1; NR > 1 {print $0 | "sort -k4,4 -k1,1"}' ../data/mismatches/mismatches_1_PNA.dat >  ../data/mismatches/mismatches_PNA_1.dat


rm -rf  ../data/mismatches/mismatches_PNA.dat  ../data/mismatches/mismatches_1_PNA.dat
