#!/bin/bash

PROJECT=../data

DIR=$PROJECT/rna_align
for i in $(ls $DIR/*.bam)
do
NAME=${i##*/}
echo "Starting mapping for sample: $NAME"
bamCoverage -b $DIR/$NAME\
	    -o ../analysis/bamcoverage/"$NAME"_normal.bw -of bigwig -bs 1\
	    --normalizeUsing CPM
bamCoverage -b $DIR/$NAME \
		    -o ../analysis/bamcoverage/"$NAME"_3primeend.bw \
		    -of bigwig -bs 1 --Offset -1 --normalizeUsing CPM
done
