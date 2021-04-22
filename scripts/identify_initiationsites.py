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
    '''
    identifies possible alternative start sites in our oligo-pool.
    :param wig_file: wiggle file containing ribosome toeprints
    :type wig_file: str
    :return: genes with position and height of alternative start site (start site not at start codon)
    :rtype: dict
    '''
    data_frame = import_wiggles(wig_file)
    alternate_sites = pd.DataFrame(columns=["gene", "position", "count"])
    gene_list = data_frame.gene.unique()
    # loop through genes:
    for gene in gene_list:
        df = data_frame[data_frame["gene"] == gene]  # make df for each oligo (gene)
        threshold = 5
        df_ass = df[(df["count"] > threshold) & (df["position"] > 80) & (df["position"] < 208)]
        df_ass = df_ass.reset_index(drop=True)  # drop row names
        df_new = df_ass.copy()
        # identify high cov. points:
        if df_ass.__len__() > 0:
            for i in df_ass.index[1:]:
                row = df_ass.iloc[i]
                prev_row = df_ass.iloc[i - 1]
                # add max if there are more high values next to each other (less than 5 NT away)
                if row["position"] - prev_row["position"] < 5:
                    index_del = df_ass.iloc[i - 1:i]["count"].idxmin()
                    df_new = df_new.drop(index_del)
            df_ass = df_new
        alternate_sites = alternate_sites.append(df_ass)
    alternate_sites = alternate_sites.reset_index(drop=True)
    return alternate_sites


def check_start_codons(fasta_file, wig_file):
    # get alt. start site df:
    alt_start_sites = find_alt_starts(wig_file)
    as_w_sc = alt_start_sites.copy()
    start_codons = ["ATG", "GTG", "TTG"]
    for record in SeqIO.parse(fasta_file, format="fasta"):
        if record.id in frozenset(alt_start_sites["gene"]):
            gene_peaks = alt_start_sites[alt_start_sites["gene"].str.match(record.id)]
            pos_peak = gene_peaks["position"]
            areas_before = [record.seq[p-20:p-10] for p in pos_peak]
            print(areas_before)
            for s in range(len(areas_before)):
                if not any(sc in areas_before[s] for sc in start_codons):
                    idx_del = gene_peaks.iloc[s].name
                    as_w_sc = as_w_sc.drop(idx_del)
                    print(idx_del)
    print(len(alt_start_sites))
    print(len(as_w_sc))
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


sc = check_start_codons("../data/reference_sequences/oligos_cds.fasta", "../data/wigglefiles/01_no_PNA.bam_3primeend_mod.wig")

#alt_starts = compare_to_in_vivo("../data/wigglefiles/01_no_PNA.bam_3primeend_mod.wig",
#                                "../data/reference_sequences/oligos_cds.fasta",
#                                "../data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_F_new.wig",
#                                "../data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_F_new.wig")



