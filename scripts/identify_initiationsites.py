"""
    Script to identify novel translation initiation sites using ribosome profiling of an oligo pool containing
    translation start regions of all coding genes.
"""

import pandas as pd
from Bio import SeqIO


def import_wiggles(wig_file_in_paths):
    wigfiles = []
    for w in wig_file_in_paths:
        df = pd.read_csv(w, sep="\t", names=["gene", "position", "pos2", "count"])
        df = df.iloc[:, [0, 1, 3]]
        wigfiles += [df]
    return wigfiles


def find_alt_starts(wig_file_in_paths):
    '''
    identifies possible alternative start sites in our oligo-pool.
    :param wig_file_in_paths: list of wiggle file paths containing ribosome toeprints
    :type wig_file_in_paths: list
    :return: genes with position and height of alternative start site (start site not at start codon)
    :rtype: dict
    '''
    data_frames = import_wiggles(wig_file_in_paths)
    alternate_sites = [pd.DataFrame(columns=["gene", "position", "count"]) for i in range(len(wig_file_in_paths))]
    gene_list = data_frames[0].gene.unique()
    # loop through genes:
    for gene in gene_list:
        dfs = [data_frame[data_frame["gene"] == gene] for data_frame in data_frames]  # make df for each oligo (gene)
        threshold = 5
        df_ass = [df[(df["count"] > threshold) & df.position.isin(list(range(30, 60)) + list(range(78, 208)))] for df in dfs]
        df_ass = [d.reset_index(drop=True) for d in df_ass]  # drop row names
        dfs_new = df_ass.copy()
        # identify high cov. points per sample:
        for ind in range(len(df_ass)):
            if df_ass[ind].__len__() > 0:
                for i in df_ass[ind].index[1:]:
                    row = df_ass[ind].iloc[i]
                    prev_row = df_ass[ind].iloc[i - 1]
                    # add max if there are more high values next to each other (less than 5 NT away)
                    if row["position"] - prev_row["position"] < 5:
                        index_del = df_ass[ind].iloc[i - 1:i]["count"].idxmin()
                        dfs_new[ind] = dfs_new[ind].drop(index_del)
                df_ass[ind] = dfs_new[ind]
            alternate_sites[ind] = alternate_sites[ind].append(df_ass[ind])
    alternate_sites = [alternate_sites[n].reset_index(drop=True) for n in range(len(alternate_sites))]
    return alternate_sites


def check_start_codons(fasta_file, wig_file_in_paths):
    # get alt. start site df:
    alt_start_sites = find_alt_starts(wig_file_in_paths)
    as_w_sc = alt_start_sites.copy()
    start_codons = ["ATG", "GTG", "TTG"]
    for record in SeqIO.parse(fasta_file, format="fasta"):
        for ind in range(len(alt_start_sites)):
            if record.id in frozenset(alt_start_sites[ind]["gene"]):
                gene_peaks = alt_start_sites[ind][alt_start_sites[ind]["gene"].str.match(record.id)]
                pos_peak = gene_peaks["position"]
                areas_before = [record.seq[p-20:p-10] for p in pos_peak]
                print(areas_before)
                for s in range(len(areas_before)):
                    if not any(sc in areas_before[s] for sc in start_codons):
                        idx_del = gene_peaks.iloc[s].name
                        as_w_sc[ind] = as_w_sc[ind].drop(idx_del)
                        print(idx_del)
    print(alt_start_sites)
    print(as_w_sc)
    return as_w_sc


def summarize_wigs(wiggle_list):
    all_samples = pd.DataFrame(columns=["gene", "position", "count"])
    for df in wiggle_list:
        all_samples = all_samples.append(df)
    all_samples = all_samples.sort_values(["gene", "position"])
    all_samples = all_samples.reset_index(drop=True)
    collapsed_samples = pd.DataFrame(columns=["gene", "position", "count"])
    for ind in range(1, len(all_samples)):
        row_before = all_samples.loc[ind-1]
        row = all_samples.loc[ind]
        if row.gene == row_before.gene and abs(row_before["position"] - row["position"]) < 5:
            new_row_idx = all_samples.loc[ind-1:ind]["count"].idxmax()
            new_row = all_samples.loc[new_row_idx]
            collapsed_samples = collapsed_samples.append(new_row)
            if len(collapsed_samples) > 1:
                if new_row["gene"] == collapsed_samples.iloc[-2]["gene"] and \
                        abs(collapsed_samples.iloc[-2]["position"] - new_row["position"]) < 5:
                    print(new_row, collapsed_samples.iloc[-2])
                    if new_row["count"] > collapsed_samples.iloc[-2]["count"]:
                        collapsed_samples = collapsed_samples.iloc[:-1]
                    else:
                        collapsed_samples = collapsed_samples.iloc[:-1]
        else:
            print(row)
    return collapsed_samples


def compare_to_in_vivo(wig_file_in_paths, fasta_file_in, other_wig_pos, other_wig_neg):
    list_alternate_starts = []
    alt_sites = check_start_codons(fasta_file_in, wig_file_in_paths)
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


wiggles_in = ["../data/wigglefiles/01_no_PNA_3primeend_new_mod.wig",
                "../data/wigglefiles/07_no_PNA_3primeend_new_mod.wig",
                "../data/wigglefiles/13_no_PNA_3primeend_new_mod.wig"]
data_frames = import_wiggles(wiggles_in)

ss = find_alt_starts(wiggles_in)

sco = check_start_codons("../data/reference_sequences/oligos_cds.fasta", wiggles_in)


#alt_starts = compare_to_in_vivo("../data/wigglefiles/01_no_PNA_3primeend_mod.wig",
#                                "../data/reference_sequences/oligos_cds.fasta",
#                                "../data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_F_new.wig",
#                                "../data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_F_new.wig")



