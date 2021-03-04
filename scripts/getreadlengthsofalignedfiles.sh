#!/bin/bash

for f in ../data/rna_align/*bam; do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    samtools view -F 4 "$f" |  cut -f 10 |\
	perl -ne 'chomp;print length($_) . "\n"' | \
	sort | uniq -c |\
	sed -E 's/[^0-9]*([0-9]*)[^0-9]([0-9]*).*/\1\t\2/g'\
	> ../analysis/readlengths/"$filename".txt
    done
