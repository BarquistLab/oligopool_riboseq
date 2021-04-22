#!/bin/bash

PROJECT=../data

DIR=$PROJECT/rna_align
for i in $(ls $DIR/*no_PNA.bam)
do
NAME=${i##*/}
echo "Starting coverage for sample: $NAME"
#bamCoverage -b $DIR/$NAME\
#	    -o ../analysis/bamcoverage/"$NAME"_normal_new.bw -of bigwig -bs 1\
#	    --normalizeUsing CPM
bamCoverage -b $DIR/$NAME \
		    -o ../analysis/bamcoverage/"$NAME"_3primeend_new.bw \
		    -of bigwig -bs 1 --Offset -1 --normalizeUsing CPM
done
