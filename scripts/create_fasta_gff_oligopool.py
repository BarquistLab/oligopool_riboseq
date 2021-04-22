"""
Title:      Create Oligos
Author:     Jakob Jung
Date:       10-06-2020
Details:    functions which take in fasta & gff file and create fasta and gffs for our oligo sequences (only specific ones)
"""

from make_oligos import chromosomes_parsed, get_seqdict
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
import os
import re
import pandas as pd


def make_oligos(cds, feature,  record, output_fasta, length_prom, length_orf, locus_tags):
    """ appends an oligo to fasta file per cds

    :param length_orf: length of ORF (we use 120)
    :param length_prom: lengt of promoter (we use 30)
    :param feature: feature of sequence (usually gene)
    :param cds: coding sequence, SeqRecord
    :param record: record (of chromosome), SeqRecord
    :param file: fasta file, path
    """
    # only extract cds regions:

    if cds.type != "CDS":
        print(cds.type)
        raise KeyError

    # create sequence start and end depending on strandedness:
    try:
        gene = feature.qualifiers["gene"][0]
    except:
        gene = locus_tag
    if feature.strand == -1:
        end = cds.location.end + length_prom
        start = cds.location.end - 3 - length_orf
        # get locus tags and gene name (if name is available)
        locus_tag = feature.qualifiers["locus_tag"][0]

        locus_tag += ("_" + gene + "_" + str(cds.location.end) + "_{}".format(feature.strand))
    else:
        start = cds.location.start
        end = start + 3 + length_orf
        start = start - length_prom
        # get locus tags and gene name (if name is available)
        locus_tag = feature.qualifiers["locus_tag"][0]
        locus_tag += ("_" + gene + "_" + str(cds.location.start+1) + "_{}".format(feature.strand))

    # generate possible sequences, which have to include start codon (3 nt), orf (120 nt) and promoter (30 nt):
    sequence = record.seq
    f = open(output_fasta, "a")
    oligo_location = FeatureLocation(start, end)
    oligo_target = oligo_location.extract(sequence)
    if feature.strand == -1:
        oligo_target = oligo_target.reverse_complement()

    oligo_rec = SeqIO.SeqRecord(oligo_target, id="{}".format(locus_tag),
                              name="", description="")

    # add T7 promoter, stop codon and 3'UTR:
    t7 = "GTTTTTTTTAATACGACTCACTATAGGG"
    stop_utr_3 = "TAATCTTTTGTAGATTGCACTTGCTTAAAAT"
    oligo_final = t7 + oligo_rec + stop_utr_3

    SeqIO.write(oligo_final, f, "fasta")
    # check length (if it is as wanted):
    if len(oligo_final) != (len(t7 + stop_utr_3) + 3 + length_orf + length_prom):
        print("length of cds {} not normal: {}".format(locus_tag, len(oligo_final)))
    f.close()
    return locus_tag


def create_gff_fasta_oligos(in_gff, in_fasta, output_fasta, length_prom, length_orf):
    """creates specific oligos_2 from gff file and saves them in a fasta file

    :param length_orf: length of ORF (we use 120)
    :param length_prom: lengt of promoter (we use 30)
    :param output_file: file to save the foutput fasta file
    :param in_gff: gff annotation file of genes, path
    :param in_fasta: fasta files of chromosomes, path
    """
    locus_tags = []
    # delete output file if it already exists:
    if os.path.exists(output_fasta):
        os.remove(output_fasta)

    seqdict = get_seqdict(in_fasta)
    with open(in_gff) as in_handle:
        chroms = chromosomes_parsed(in_handle, seqdict)
        genes = 0
        # loop through all chromosome records within the gff file:
        for record in chroms:
            # loop through all features (genes/sRNAs/...) for each chromosome record (record):
            for feature in record.features:
                # now loop through all sub-features (cds) of the genes:
                for cds in feature.sub_features:
                    # I use try/except for the case that there's no locus tag available (-> KeyError):
                    try:
                        lt = make_oligos(cds, feature, record, output_fasta, length_prom, length_orf, locus_tags)
                        genes += 1
                        locus_tags += [lt]
                        print(cds)
                    except KeyError:
                        continue
    print("{} oligos_2 created and written to fasta file: {}".format(genes, output_fasta))
    return locus_tags


if __name__ == '__main__':
    # define paths for annotation GFF file and FASTA (should be input by user):
    ecoli_fa = "../data/reference_sequences/NC_000913.3.fasta"
    ecoli_gff = "../data/reference_sequences/NC_000913.3.gff"
    # create the fasta file with oligo sequences
    ltags = create_gff_fasta_oligos(ecoli_gff, ecoli_fa, "../data/reference_sequences/oligos_cds_new.fasta", 30, 150)


    # change gff:
    with open("../data/reference_sequences/NC_000913.3.gff", 'r') as gff:
        with open("../data/reference_sequences/oligos_cds_new.gff", 'w') as new_gff:
            for line in gff:
                line_splitted = line.split("\t")
                if len(line_splitted) != 9:
                    continue
                annotations = line_splitted[8]
                match = ".*;locus_tag=(b\\d\\d\\d\\d).*"
                match2 = ".*;gene=([^;]+).*"
                m = re.search(match, annotations)
                m2 = re.search(match2, annotations)  # gene name
                if m is None:
                    continue
                if m2 is None:
                    m2 = m
                lt = m.group(1)
                g = m2.group(1)
                if line_splitted[6] == "+":
                    id = lt + "_" + g + "_" + line_splitted[3] + "_" + "1"
                else:
                    id = lt + "_" + g + "_" + line_splitted[4] + "_" + "-1"
                line_splitted[3] = "1"
                line_splitted[4] = "242"
                line_splitted[0] = id
                newline = "\t".join(line_splitted)
                if line_splitted[0] in ltags:
                    new_gff.write(newline)

    gff_new = pd.read_csv("../data/reference_sequences/oligos_cds_new.gff", sep="\t",  header=None)
    all_gff_entries = gff_new[[0]].to_numpy()
    with open("../data/reference_sequences/oligos_cds_new.fasta") as fasta:
        for line in fasta:
            m = re.match(">(.*)\n", line)
            if m:
                if m.group(1) not in all_gff_entries:
                    print(m.group(1))
                    f = open("../data/reference_sequences/oligos_cds_new.gff", "a+")
                    if m.group(1).split("_")[3] == "-1":
                        f.write(m.group(1) +
                                "\tRefSeq\tgene\t1\t242\t.\t-\t.\tlocus_tag={};gene={}".format(m.group(1).split("_")[0],
                                                                                          m.group(1).split("_")[1]))
                    else:
                        f.write(m.group(1) +
                                "\tRefSeq\tgene\t1\t242\t.\t+\t.\tlocus_tag={};gene={}".format(m.group(1).split("_")[0],
                                                                                          m.group(1).split("_")[1]))
                    f.write("\n")
                    f.close()

    # change other gff (for trnas and rrnas):
    with open("../data/reference_sequences/trnas_rrnas.gff", 'r') as gff:
        with open("../data/reference_sequences/trnas_rrnas_mod.gff", 'w') as new_gff:
            for line in gff:
                line_splitted = line.split("\t")
                print(line_splitted)
                if len(line_splitted) != 9:
                    continue
                start = int(line_splitted[3])
                end = int(line_splitted[4])
                newstart = 1
                newend = end-start+1
                newid = line_splitted[0] + ":" + str(start-1) + "-" + str(end)

                line_splitted[3] = "1"
                line_splitted[4] = str(newend)
                line_splitted[0] = newid
                newline = "\t".join(line_splitted)
                print(newline)
                new_gff.write(newline)


os.system("cat ../data/reference_sequences/trnas_rrnas_mod.gff >> ../data/reference_sequences/oligos_cds_new.gff")
os.system("cat ../data/reference_sequences/trnas_rrnas.fasta >> ../data/reference_sequences/oligos_cds_new.fasta")
