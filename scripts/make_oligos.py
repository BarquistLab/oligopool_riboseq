"""
Title:      Create Oligos
Author:     Jakob Jung
Date:       10-06-2020
Details:    functions which take in fasta & gff file and create specific oligo sequences  in FASTA format
"""

from BCBio import GFF
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
import os
import re
import pandas as pd


def chromosomes_parsed(in_gff, ref_dict):
    """Parse gff output, generating SeqRecord and SeqFeatures for chromosome & plasmids of fasta
    """
    for rec in GFF.parse(in_gff,  base_dict=ref_dict):
        yield rec


def get_seqdict(fasta_file):
    """returns sequence dictionary of inputted fasta file for later use.

    :param fasta_file: input file of chromosomes in fasta format
    :return: returns sequence dictionary in dict format
    """
    with open(fasta_file) as fasta:
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    return seq_dict


def make_oligos(cds, feature,  record, file, length_prom, length_orf):
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
        raise KeyError

    # create sequence start and end depending on strandedness:
    if feature.strand == -1:
        end = cds.location.end + length_prom
        start = cds.location.end - 3 - length_orf
        start_codon = FeatureLocation(cds.location.end - 3, cds.location.end).extract(record.seq).reverse_complement()
    else:
        start = cds.location.start
        end = start + 3 + length_orf
        start = start - length_prom
        start_codon = FeatureLocation(cds.location.start, cds.location.start + 3).extract(record.seq)
    # do warning if start codon is not usual (ATG/GTG/TTG)
    if start_codon not in ["ATG", "GTG", "TTG", "ATT", "CTG"]:
        print("Start codon for {} not ATG, GTG, TTG, ATT and CTG!!!".format(cds.qualifiers['gene'][0]))
    # get locus tags and gene name (if name is available)
    locus_tag = feature.qualifiers["locus_tag"][0]
    try:
        gene = cds.qualifiers['gene'][0]
    except KeyError:
        gene = "unknown gene name"
    # generate possible sequences, which have to include start codon (3 nt), orf (120 nt) and promoter (30 nt):
    sequence = record.seq
    f = open(file, "a")
    oligo_location = FeatureLocation(start, end)
    oligo_target = oligo_location.extract(sequence)
    if feature.strand == -1:
        oligo_target = oligo_target.reverse_complement()
    oligo_rec = SeqIO.SeqRecord(oligo_target, id="{}; {}".format(locus_tag, gene),
                              name="{}".format(locus_tag), description="{}".format(feature.strand))
    # add T7 promoter, stop codon and 3'UTR:
    t7 = "GTTTTTTTTAATACGACTCACTATAGGG"
    stop_utr_3 = "TAATCTTTTGTAGATTGCACTTGCTTAAAAT"
    oligo_final = t7 + oligo_rec + stop_utr_3
    #SeqIO.write(oligo_final, f, "fasta")
    f.write("JVOpool-001" +"\t" + str(oligo_final.seq)+"\n")
    # check length (if it is as wanted):
    if len(oligo_final) != (len(t7 + stop_utr_3) + 3 + length_orf + length_prom):
        print("length of cds {} not normal: {}".format(locus_tag, len(oligo_final)))
    f.close()


def create_oligos(in_gff, in_fasta, output_file, length_prom, length_orf):
    """creates specific oligos_2 from gff file and saves them in a fasta file

    :param length_orf: length of ORF (we use 120)
    :param length_prom: lengt of promoter (we use 30)
    :param output_file: file to save the foutput fasta file
    :param in_gff: gff annotation file of genes, path
    :param in_fasta: fasta files of chromosomes, path
    """
    # delete output file if it already exists:
    if os.path.exists(output_file):
        os.remove(output_file)
    # create header for output tabfile:
    with open(output_file, "w") as f:
        f.write("Pool name" + "\t" + "Sequence" + "\n")

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
                        make_oligos(cds, feature, record, output_file, length_prom, length_orf)
                        genes += 1
                    except KeyError:
                        continue
    print("{} oligos_2 created and written to fasta file: {}".format(genes, output_file))


if __name__ == '__main__':
    # define paths for annotation GFF file and FASTA (should be input by user):
    ecoli_fa = "../data/e_coli_ref_seqs/NC_000913.3.fasta"
    ecoli_gff = "../data/e_coli_ref_seqs/NC_000913.3.gff"
    # create the fasta file with oligo sequences
    create_oligos(ecoli_gff, ecoli_fa, "../data/oligo_order/oligos_ecoli.tab", 30, 150)

    with open("../data/e_coli_ref_seqs/NC_000913.3.gff", 'r') as gff:
        with open("../data/oligos_cds.gff", 'w') as new_gff:
            for line in gff:
                line_splitted = line.split("\t")
                #print(line_splitted)
                if len(line_splitted) != 9:
                    continue
                line_splitted[3] = "1"
                line_splitted[4] = "183"
                annotations = line_splitted[8]
                match = ".*;locus_tag=(b\\d\\d\\d\\d).*"
                m = re.search(match, annotations)
                if m is None:
                    continue
                if not line_splitted[2] == "gene":
                    continue
                lt = m.group(1)
                print(lt)
                if line_splitted[6] == "+":
                    id = lt + " " + "1"
                else:
                    id = lt + " " + "-1"
                print(id)
                line_splitted[0] = id
                #print(annotations)
                #print(line_splitted)
                newline = "\t".join(line_splitted)
                new_gff.write(newline)













