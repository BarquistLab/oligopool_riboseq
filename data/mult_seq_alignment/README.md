## multiple sequence alignments

First, I downloaded all genomes related to *E. coli* (*Salmonella, Klebsiella, Citrobacter*) and blasted the TIR of our possible alternative start sites against it:

```bash
# combine fastas:
cat FQ312003_wplasmids.fasta klebsiella.fasta citrobacter.fasta > all_seqs.fasta
# create BLAST DB:
makeblastdb -in all_seqs.fasta -dbtype nucl


```

Now we can run for each region this (after checking for the start site): 

```bash
GNAME=b3888_yiiD
echo ">Escherichia coli K12" > query_$GNAME.fasta
echo "AAACAAGAGAGAGTATCGCTATGTATCACC" >> query_$GNAME.fasta
# create result file:
cat query_$GNAME.fasta > result_$GNAME.fasta
# run blastn and get it into fasta format:
blastn -task blastn  -query query_$GNAME.fasta -db all_seqs.fasta -num_threads 4  -evalue 0.1  | grep -E ">|Sbjct" | sed -E "s/Sbjct[^ATGC]+([ATGC]+).*/\\1/" >> result_$GNAME.fasta 
```

