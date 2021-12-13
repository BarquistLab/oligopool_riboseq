
from BCBio import GFF
from make_oligos import chromosomes_parsed, get_seqdict
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
import pandas as pd
import re


def import_wiggles(wig_file_in_paths):
    wigfiles = []
    for w in wig_file_in_paths:
        df = pd.read_csv(w, sep="\t", names=["gene", "position", "pos2", "count"])
        df = df.iloc[:, [0, 1, 3]]
        wigfiles += [df]
    return wigfiles


in_gff = "../data/reference_sequences/NC_000913.3.gff"
in_fasta = "../data/reference_sequences/NC_000913.3.fasta"
in_oligos_fasta = "../data/reference_sequences/oligos_cds_new.fasta"
wiggles_in = ["../data/wigglefiles/01_no_PNA_3primeend_new_mod.wig",
              "../data/wigglefiles/07_no_PNA_3primeend_new_mod.wig",
              "../data/wigglefiles/13_no_PNA_3primeend_new_mod.wig",
              "../data/wigglefiles/79_no_PNA_3primeend_new_mod.wig",
              "../data/wigglefiles/85_no_PNA_3primeend_new_mod.wig"]


wig_oligos = import_wiggles(wiggles_in)
all_oligos = set(wig_oligos[0]["gene"])
wig_invivo_F = "../data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_F_new.wig"
wig_invivo_R = "../data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_R_new.wig"
threshold = 5

df_invivo_F = pd.read_table(wig_invivo_F, skiprows=2, index_col=None, names=["position", "count"])
df_invivo_R = pd.read_table(wig_invivo_R, skiprows=2, index_col=None, names=["position", "count"])
df_invivo_R[1] = df_invivo_R["count"] * (-1)


#out_tsv = open("../analysis/check_ann_ss.tsv", "w")
#out_tsv.write("\t".join(["locus_tag", "gene", "start_position", "strand", "peak_R1", "peak_R2", "peak_R3", "peak_R4",
#                         "peak_R5", "peak_found_Hör", "peak_wmeydan" "peak_found_meydan"]) + "\n")
out_df = pd.DataFrame(columns=["locus_tag", "gene", "start_position", "start_codon", "strand", "peak_R1", "peak_R2",
                               "peak_R3", "peak_R4", "peak_R5", "peak_found_Hör", "peak_meydan", "peak_found_meydan"])


for oligo in SeqIO.parse(in_oligos_fasta, "fasta"):
    if oligo.id.startswith("b"):
        start_codon = oligo.seq[58:61]
        id = oligo.id
        strand = id.split("_")[3]
        locus_tag = id.split("_")[0]
        gene = id.split("_")[1]
        start_pos = int(id.split("_")[2])

        newline = [locus_tag, gene, start_pos, str(start_codon), strand, 0, 0, 0, 0, 0, False, 0, False]
        dfs = [data_frame[data_frame["gene"] == id] for data_frame in wig_oligos]

        # check for known start sites:
        df_kss = [df[(df["count"] > threshold) & df.position.isin(list(range(70, 77)))] for df in dfs]
        rec_peaks = 0
        for df in range(len(df_kss)):
            if len(df_kss[df]) > 0:
                df_kss[df]["position"] = df_kss[df]["position"] - 58
                peak_int = df_kss[df].loc[df_kss[df]["count"].idxmax()]["count"]
                peak_pos = df_kss[df].loc[df_kss[df]["count"].idxmax()]["position"]
                newline[df+5] = peak_int
                rec_peaks += 1
        if rec_peaks > 1:
            newline[10] = True
            print(peak_int, gene, peak_pos)

        if strand == "1":
            start_seq = start_pos - 58
            end_seq = start_pos + 184
            d = df_invivo_F[df_invivo_F.iloc[:, 0].isin(list(range(start_seq, end_seq)))]
            newcol = d.iloc[:, 0] - start_pos
            d.iloc[:, 0] = newcol
        else:
            end_seq = start_pos + 58
            start_seq = start_pos - 184
            d = df_invivo_R[df_invivo_R.iloc[:, 0].isin(list(range(start_seq, end_seq)))]
            newcol = (d.iloc[:, 0] - start_pos) * (-1)
            d.iloc[:, 0] = newcol

        df_invivo_kss = d[(d["count"] > threshold) & d.position.isin(list(range(12, 19)))]
        if len(df_invivo_kss) > 0:
            peak_int = df_invivo_kss.loc[df_invivo_kss["count"].idxmax()]["count"]
            peak_pos = df_invivo_kss.loc[df_invivo_kss["count"].idxmax()]["position"]
            print(peak_int, gene, peak_pos)
            newline[11] = peak_int
            newline[12] = True




        out_df = out_df.append(pd.DataFrame(columns=out_df.columns, data=[newline]))
        #out_tsv.write("\t".join(newline)+"\n")

out_df.to_csv("../analysis/annotated_sites.csv")

print("done with first comp.")

tot_length = len(out_df)
found_in_our_ds = sum(out_df["peak_found_Hör"])
found_in_other_ds = sum(out_df["peak_found_meydan"])
found_both = sum(out_df["peak_found_meydan"] & out_df["peak_found_Hör"])
start_codons = ["ATG", "GTG", "TTG", "CTG", "ATC", "ATT"]
stop_codons = ["TGA", "TAA", "TAG"]
sc_re = r"(ATG)|(GTG)|(TTG)|(CTG)|(ATC)|(ATT)"

out_df_ss = pd.DataFrame(columns=["locus_tag", "gene", "start_position", "start_codon", "pos_from_annot_start",
                                  "in_frame", "strand", "peak", "nr_peaks_Hör", "peak_found_Hör", "peak_meydan",
                                  "peak_found_meydan"])

for oligo in SeqIO.parse(in_oligos_fasta, "fasta"):
    if oligo.id.startswith("b"):
        id = oligo.id
        print(id)
        strand = id.split("_")[3]
        locus_tag = id.split("_")[0]
        gene = id.split("_")[1]
        start_pos = int(id.split("_")[2])

        dfs = [data_frame[data_frame["gene"] == id] for data_frame in wig_oligos]
        # get total read counts for relative density:
        tot_counts = [sum(data_frame["count"]) for data_frame in dfs]
        for n in range(len(dfs)):
            dfs[n]["norm_count"] = dfs[n]["count"] / sum(dfs[n]["count"])

        # check for start sites which have norm_value (rel_density) more than 0.1:
        df_ss = [df[(df["count"] > threshold) & (df["norm_count"] > 0.1) &
                    df.position.isin(list(range(30, 70)) + list(range(79, 210)))] for df in dfs]
        df_ss = [df.reset_index(drop=True) for df in df_ss]
        rec_peaks = 0

        for df in range(len(df_ss)):
            df_ss[df]["norm_counts"] = df_ss[df]["count"] / tot_counts[df]
            if len(df_ss[df]) > 0:
                nr = 1
                # just get the peaks:
                while nr < len(df_ss[df]):
                    row = df_ss[df].iloc[nr]
                    prev_row = df_ss[df].iloc[nr - 1]
                    if row["position"] - prev_row["position"] < 5:
                        index_del = df_ss[df].iloc[nr - 1:nr + 1]["count"].idxmin()
                        df_ss[df] = df_ss[df].drop(index_del)
                        df_ss[df].reset_index(drop=True)
                    else:
                        nr += 1
                # remove if theres no start codon in region before:
                df_ss[df] = df_ss[df].reset_index(drop=True)

        for df in range(len(df_ss)):
            for peak in range(len(df_ss[df])):
                peakpos = df_ss[df].loc[peak]["position"]
                region_sc = oligo.seq[peakpos-18:peakpos-11]
                if not any(sc in region_sc for sc in start_codons):
                    df_ss[df] = df_ss[df].drop(peak)

            df_ss[df] = df_ss[df].reset_index(drop=True)
            df_ss[df]["position"] = df_ss[df]["position"] - 58

            for peak in range(len(df_ss[df])):
                peakpos = df_ss[df].loc[peak]["position"]
                region_sc = oligo.seq[peakpos + 58 - 18:peakpos + 58 - 11]
                print(region_sc, gene)
                start_c = re.search(sc_re, str(region_sc))[0]
                pos_from_peak = re.search(sc_re, str(region_sc)).start() + peakpos - 18

                newline = [locus_tag, gene, start_pos, start_c, 0, "in_frame", strand, 0, 1, False, 0, False]
                if strand == "1":
                    pos_pg = df_ss[df].loc[peak].position + start_pos
                    iv_d = df_invivo_F[df_invivo_F["position"].isin(list(range(pos_pg - 5, pos_pg + 5)))]
                else:
                    pos_pg = df_ss[df].loc[peak].position*(-1) + start_pos
                    iv_d = df_invivo_F[df_invivo_F["position"].isin(list(range(pos_pg - 5, pos_pg + 5)))]
                if any(iv_d["count"] > 5):
                    newline[11] = True
                    newline[10] = iv_d.loc[iv_d["count"].idxmax()]["count"]

                peak_int = df_ss[df].iloc[peak]["count"]
                peak_pos = df_ss[df].iloc[peak]["position"]
                newline[4] = pos_from_peak
                newline[7] = peak_int
                rec_peaks += 1
                # check in frame:
                if any(sc in oligo.seq[peak_pos+58:210] for sc in stop_codons):
                    newline[5] = "out_of_frame"
                out_df_ss = out_df_ss.append(pd.DataFrame(columns=out_df_ss.columns, data=[newline]))

out_df_ss = out_df_ss.sort_values(["locus_tag", "start_position", "pos_from_annot_start"])
out_df_ss = out_df_ss.reset_index(drop=True)


reduced_df_ss = out_df_ss.copy()
# just get the peaks:
nr = 1
# just get the peaks:
while nr < len(reduced_df_ss):
    row = reduced_df_ss.iloc[nr]
    prev_row = reduced_df_ss.iloc[nr - 1]
    if abs(row["pos_from_annot_start"] - prev_row["pos_from_annot_start"]) < 5 and row["locus_tag"] == prev_row["locus_tag"]:
        index_del = reduced_df_ss.iloc[nr - 1:nr + 1]["peak"].idxmin()
        index_max = reduced_df_ss.iloc[nr - 1:nr + 1]["peak"].idxmax()
        reduced_df_ss = reduced_df_ss.drop(index_del)
        reduced_df_ss.loc[index_max, "peak_found_Hör"] = True
        reduced_df_ss.loc[index_max, "nr_peaks_Hör"] = reduced_df_ss.loc[index_max, "nr_peaks_Hör"] + 1
        print(row["locus_tag"], prev_row["locus_tag"])
        reduced_df_ss.reset_index(drop=True)
    else:
        nr += 1

reduced_df_ss = reduced_df_ss.reset_index(drop=True)

reduced_df_ss.to_csv("../analysis/alt_tis_w_threshold.csv")
sum(reduced_df_ss["peak_found_Hör"])


meydan_genes = pd.read_csv("../data/wigglefiles/meydan_alt_tis.csv")














chrom = get_seqdict("../data/reference_sequences/NC_000913.3.fasta")
gffs = chromosomes_parsed("../data/reference_sequences/NC_000913.3.gff", chrom)

comparison_start_sites = open("../analysis/comparison_startsites.tab", "w")

with open(in_gff) as in_handle:
    chroms = chromosomes_parsed(in_handle, chrom)
    genes = 0
    for record in chroms:
        # loop through all features (genes/sRNAs/...) for each chromosome record (record):
        for feature in record.features:
            # now loop through all sub-features (cds) of the genes:
            if feature.type == "gene":
                for cds in feature.sub_features:
                    if cds.type == "CDS":

                        print(feature.qualifiers['locus_tag'])
                        print(cds.location.strand)
                        print(record.seq[cds.location.start:cds.location.start+3])
                        genes += 1
    print(genes)


comparison_start_sites.close()