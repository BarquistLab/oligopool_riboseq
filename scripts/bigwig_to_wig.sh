#!/bin/bash

main(){
    # go to project directory where the reads and reference
    # sequences are stored:
    PROJECT=../analysis/bamcoverage
    echo "Convert bigwig to wiggle files!"
    bigwig_to_wiggle
}

bigwig_to_wiggle(){
    for i in $(ls $PROJECT/*_3primeend.bw)
    do
        NAME=${i%.*}
        echo $NAME
	# convert:
	bigWigToWig $i "$NAME.wig"
    done
}

main
