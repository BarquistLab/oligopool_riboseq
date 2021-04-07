#!/bin/bash

#~/bin/bbmap/bbduk.sh -Xmx1g in=SRR8156056.fastq \
#			     ref=~/bin/bbmap/resources/adapters.fa t=20\
#			     out=SRR8156056_trimmed.fastq ktrim=r k=23 mink=11\
#			     hdist=1 qtrim=r trimq=10

#~/bin/bbmap/bbmap.sh in=SRR8156056_trimmed.fastq trimreaddescription=t  t=20 \
#			     ref=../../reference_sequences/NC_000913.3.fasta \
#			     k=13 ambig=toss outm=SRR8156056.sam

~/bin/bbmap/bbmap.sh in=SRR8156056_trimmed.fastq trimreaddescription=t  t=20 \
			     ref=../../reference_sequences/ecolibw25113.fasta \
			     k=13 ambig=toss outm=SRR8156056.sam


# sort sam file, create BAM file:
samtools sort -O BAM -@ 40 SRR8156056.sam > SRR8156056.bam
# remove sam file: (not actually needed)
rm SRR8156056.sam
samtools index SRR8156056.bam
bamCoverage -b SRR8156056.bam \
	    -o SRR8156056_normal.bw -of bigwig -bs 1
bamCoverage -b SRR8156056.bam \
	    -o SRR8156056_3primeend.bw \
	    -of bigwig -bs 1 --Offset -1
