"""
    Script to identify novel translation initiation sites using ribosome profiling of an oligo pool containing
    translation start regions of all coding genes.
"""

import pandas as pd
from Bio import SeqIO


def import_wiggles(wig_file):
    df = pd.read_csv(wig_file, sep="\t", names=["gene", "position", "pos2", "count"])
    df = df.iloc[:, [0, 1, 3]]
    return df


def find_alt_starts(wig_file):
    data_frame = import_wiggles(wig_file)
    alternate_sites = {}
    gene_list = data_frame.gene.unique()
    # loop through genes:
    for gene in gene_list:
        df = data_frame[data_frame["gene"] == gene]
        df_ass = df[(df["count"] > 20) & (df["position"] > 100)]
        if df_ass.__len__() > 0:
            df_ass = df_ass.reset_index(drop=True)
            maxind = df_ass["count"].idxmax()
            alternate_sites[gene] = [df_ass.iloc[maxind]["position"], df_ass.iloc[maxind]["count"]]
    return alternate_sites


def check_start_codons(fasta_file, wig_file):
    alt_start_sites = find_alt_starts(wig_file)
    as_w_sc = {}
    start_codons = ["ATG", "GTG", "TTG"]
    for record in SeqIO.parse(fasta_file, format="fasta"):
        if record.id in alt_start_sites.keys():
            pos_peak = alt_start_sites[record.id][0]
            area_before = record.seq[pos_peak-40:pos_peak]
            if any(sc in area_before for sc in start_codons):
                as_w_sc[record.id] = alt_start_sites[record.id]
    return as_w_sc


def compare_to_in_vivo(in_wig, fasta_file_in, other_wig_pos, other_wig_neg):
    list_alternate_starts = []
    alt_sites = check_start_codons(fasta_file_in, in_wig)
    print(alt_sites, len(alt_sites))
    for alt_site in alt_sites.keys():
        locus_tag, start_pos, strand = alt_site.split("_")
        start_pos = int(start_pos)
        pos_oligo = alt_sites[alt_site][0]
        if strand == "-1":
            tp = pd.read_csv(other_wig_neg, skiprows=2, names=["counts"])
            peak_pos = start_pos + 58 - pos_oligo
        else:
            tp = pd.read_csv(other_wig_pos, skiprows=2, names=["counts"])
            peak_pos = start_pos - 58 + pos_oligo
        sum_counts_other_wigs = tp["counts"][peak_pos-20:peak_pos+20].sum()
        if sum_counts_other_wigs > 10:
            print(locus_tag)
            list_alternate_starts += [locus_tag]
    return list_alternate_starts


alt_starts = compare_to_in_vivo("../data/wigglefiles/01_no_PNA.bam_3primeend_mod.wig",
                                "../data/reference_sequences/oligos_cds.fasta",
                                "../data/wigglefiles/GSM3509320_Onc112_plus.wig",
                                "../data/wigglefiles/GSM3509320_Onc112_minus.wig")



