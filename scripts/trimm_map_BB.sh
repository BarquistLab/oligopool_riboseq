#!/bin/bash

main(){
    # go to project directory where the reads and reference
    # sequences are stored:
    PROJECT=../data
    echo "Start trimming"
    rename_trim_rna_libs
    echo "Trimming done. Start maopping"
    align_rna_reads_genome
    echo "Finished mapping. Start connecting all tab files"
    # here I count reads for sRNA and 
    featureCounts -T 5 -t CDS,sRNA -g locus_tag \
		  -a $PROJECT/reference_sequences/oligos_cds.gff \
		  -o $PROJECT/rna_align/counttable.txt \
		  $PROJECT/rna_align/*.bam
    
    # Now find nr of rRNAs:
    featureCounts -T 5 -t rRNA -g gene \
		  -a $PROJECT/reference_sequences/oligos_cds.gff \
		  -o $PROJECT/rna_align/counttable_rRNAs.txt \
		  $PROJECT/rna_align/*.bam
    # Now find nr of tRNAs:
    featureCounts -T 5 -t tRNA -g gene \
		  -a $PROJECT/reference_sequences/oligos_cds.gff \
		  -o $PROJECT/rna_align/counttable_tRNAs.txt \
		  $PROJECT/rna_align/*.bam
    
}

rename_trim_rna_libs(){
    mkdir -pv $PROJECT/libs
    for NAME in $(ls $PROJECT/fastq/*.fq.gz)
    do
        echo "$NAME starts trimming nowwwwwww"
        NEWNAME=${NAME##*/}
        NEWNAME=${NEWNAME%.fq.gz}_trimmed.fq.gz
        echo $NEWNAME
       # bbduk trims low quality bases and removes adapters:
        ~/bin/bbmap/bbduk.sh -Xmx1g in=$NAME \
			     ref=~/bin/bbmap/resources/adapters.fa t=20\
			     out=$PROJECT/libs/${NEWNAME} ktrim=r k=23 mink=11\
			     hdist=1 qtrim=r trimq=10
	# add fastqc of reads:
	fastqc -o ../analysis/fastqc -f $PROJECT/libs/${NEWNAME}
    done
}


align_rna_reads_genome(){
    mkdir -pv $PROJECT/rna_align
    DIR=$PROJECT/rna_align
    for i in $(ls $PROJECT/libs/*.fq.gz)
    do
        NAME=${i##*/}
        NAME=${NAME%_trimmed.fq.gz}
        echo "Starting mapping for sample: $NAME"
        ~/bin/bbmap/bbmap.sh in=$i trimreaddescription=t  t=20 \
			     ref=$PROJECT/reference_sequences/oligos_cds.fasta \
			     k=13 ambig=toss outm=$DIR/$NAME.sam
        # sort sam file, create BAM file:
        samtools sort -O BAM -@ 40 $DIR/$NAME.sam > $DIR/$NAME.bam
        # remove sam file: (not actually needed)
        rm $DIR/$NAME.sam
	samtools index $DIR/$NAME.bam
	bamCoverage -b $DIR/$NAME.bam \
		    -o ../analysis/bamcoverage/"$NAME"_normal.bw -of bigwig -bs 1
	bamCoverage -b $DIR/$NAME.bam \
		    -o ../analysis/bamcoverage/"$NAME"_3primeend.bw \
		    -of bigwig -bs 1 --Offset -1
    done
}

main  
