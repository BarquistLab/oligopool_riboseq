import subprocess

from BCBio import GFF
from make_oligos import chromosomes_parsed, get_seqdict
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
import pandas as pd
import numpy as np
import re
from ushuffle import shuffle, Shuffler
import random


def import_wiggles(wig_file_in_paths):
    wigfiles = []
    for w in wig_file_in_paths:
        df = pd.read_csv(w, sep="\t", names=["gene", "position", "pos2", "count"])
        df = df.iloc[:, [0, 1, 3]]
        wigfiles += [df]
    return wigfiles


# import files:
in_gff = "../data/reference_sequences/NC_000913.3.gff"
in_fasta = "../data/reference_sequences/NC_000913.3.fasta"
in_oligos_fasta = "../data/reference_sequences/oligos_cds_new.fasta"
wiggles_in = ["../data/wigglefiles/01_no_PNA_3primeend_new_mod.wig",
              "../data/wigglefiles/07_no_PNA_3primeend_new_mod.wig",
              "../data/wigglefiles/13_no_PNA_3primeend_new_mod.wig",
              "../data/wigglefiles/79_no_PNA_3primeend_new_mod.wig",
              "../data/wigglefiles/85_no_PNA_3primeend_new_mod.wig"]


# import the wiggle files:
wig_oligos = import_wiggles(wiggles_in)
all_oligos = set(wig_oligos[0]["gene"])
wig_invivo_F = "../data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_F_new.wig"
wig_invivo_R = "../data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_R_new.wig"
threshold = 5   # threshold above which a peak is called

df_invivo_F = pd.read_table(wig_invivo_F, skiprows=2, index_col=None, names=["position", "count"])
df_invivo_R = pd.read_table(wig_invivo_R, skiprows=2, index_col=None, names=["position", "count"])
df_invivo_R[1] = df_invivo_R["count"] * (-1)


# out_tsv = open("../analysis/check_ann_ss.tsv", "w")
# out_tsv.write("\t".join(["locus_tag", "gene", "start_position", "strand", "peak_R1", "peak_R2", "peak_R3", "peak_R4",
#                         "peak_R5", "peak_found_Hör", "peak_wmeydan" "peak_found_meydan"]) + "\n")

# output dataframeis created, capturing all peaks
out_df = pd.DataFrame(columns=["locus_tag", "gene", "start_position", "start_codon", "strand", "peak_R1", "peak_R2",
                               "peak_R3", "peak_R4", "peak_R5", "peak_found_Hör", "peak_meydan", "peak_found_meydan",
                               "pred_RBS_region", "p_value", "sequence_start", "upstream_seq", "MFE"])


gibbs_df = pd.DataFrame(columns=["locus_tag", "gene", "SD_seq", "pos_from_stc", "shuffled", "delta_G", "p_value"])

# add a table showing enriched start codons/sd motifs:
motif_df = pd.DataFrame(np.zeros((35, 8)), columns=["ATG", "GTG", "TTG", "CTG", "ATC", "ATT", "AGG", "GGA"])
re_motiv = r"(ATG)|(GTG)|(TTG)|(CTG)|(ATC)|(ATT)|(AGG)|(GGA)"

# go through all oligos and check for annotated TIS
for oligo in SeqIO.parse(in_oligos_fasta, "fasta"):
    if oligo.id.startswith("b"):
        start_codon = oligo.seq[58:61]  # extract start codon
        id = oligo.id
        strand = id.split("_")[3]
        locus_tag = id.split("_")[0]
        gene = id.split("_")[1]
        start_pos = int(id.split("_")[2])

        newline = [locus_tag, gene, start_pos, str(start_codon), strand, 0, 0, 0, 0, 0, False, 0, False, False, False,
                   False, False, False]
        dfs = [data_frame[data_frame["gene"] == id] for data_frame in wig_oligos]


        # calculate sd-binding regions:
        region_sd = oligo.seq[58 - 25: 58]

        print(region_sd)
        # calculate sd-binding regions:
        fasta_pathname = "../data/sd_search/freescan_files_annot_ss/seq" + "_" + locus_tag
        with open(fasta_pathname + ".fasta", "w") as fasta:
            fasta.write(">Seq_SD\n")
            fasta.write(region_sd.__str__())
        sd_bits = str.encode(region_sd.__str__())

        # freescan determine delta gibbs etc:
        ASD = "auuccuccacuag"
        output_freescan = subprocess.run(["free_scan", "-e", ASD, fasta_pathname + ".fasta"], capture_output=True).stdout
        # get list with all free energies:
        list_values = output_freescan.decode("utf-8").split("\n")
        list_values = list_values[:-2]
        list_values = [float(i) for i in list_values]
        lowest_val = min(list_values)
        ind_lowest_val = list_values.index(lowest_val)
        best_binding_region = region_sd[ind_lowest_val:ind_lowest_val + len(ASD)]
        print(locus_tag, best_binding_region, lowest_val, ind_lowest_val, "original")
        with open(fasta_pathname + ".fasta", "w") as fasta:
            fasta.write(">Seq_SD\n")
            fasta.write(best_binding_region.__str__())


        shuffled_gibbs = np.array([])
        for i in range(5):
            bits_best_binding = str.encode(best_binding_region.__str__())
            shuffled = shuffle(bits_best_binding, 2).decode("utf-8")
            print(gene, shuffled, "shuffled", i)
            with open(fasta_pathname + ".fasta", "w") as fasta:
                fasta.write(">Seq_SD\n")
                fasta.write(shuffled + "\n")
            output_freescan = subprocess.run(["free_scan", "-e", "auuccuccacuag", fasta_pathname + ".fasta"], capture_output=True).stdout  # freescan
            print(output_freescan)
            # get list with all free energies:
            val = float(output_freescan.decode("utf-8")[:-2])
            newline_gibbs = [locus_tag, gene, shuffled, ind_lowest_val - 13, True, val]
            # gibbs_df = pd.concat([gibbs_df, pd.DataFrame(columns=gibbs_df.columns, data=[newline_gibbs])])
            shuffled_gibbs = np.append(shuffled_gibbs, val)
        # end freescan determine delta gibbs etc:

        less_eq = np.less_equal(shuffled_gibbs, lowest_val)
        print(lowest_val)
        pval = sum(less_eq)/len(shuffled_gibbs)
        print("pvalue for gene", gene, pval)

        newline_gibbs = [locus_tag, gene, best_binding_region.__str__(), ind_lowest_val - 13, False, lowest_val, pval]
        gibbs_df = pd.concat([gibbs_df, pd.DataFrame(columns=gibbs_df.columns, data=[newline_gibbs])])
        newline[13] = best_binding_region.__str__()
        newline[14] = pval
        newline[15] = region_sd.__str__()


        # check for known start sites (at position 15+-3):
        df_kss = [df[(df["count"] > threshold) & df.position.isin(list(range(70, 77)))] for df in dfs]
        rec_peaks = 0

        # loop through 5 dfs and check whether theres a peak per replicate:
        for df in range(len(df_kss)):
            if len(df_kss[df]) > 0:
                df_kss[df]["position"] = df_kss[df]["position"] - 58
                peak_int = df_kss[df].loc[df_kss[df]["count"].idxmax()]["count"]
                peak_pos = df_kss[df].loc[df_kss[df]["count"].idxmax()]["position"]
                newline[df+5] = peak_int
                rec_peaks += 1
        # check whether there are 2 or more peaks -> validate TIS:
        if rec_peaks > 1:
            newline[10] = True
            print(peak_int, gene, peak_pos, locus_tag)
        # check in vivo dataframes to get whether it was also detected in their dataset:
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
        # check for peaks above threshold:
        df_invivo_kss = d[(d["count"] > threshold) & d.position.isin(list(range(12, 19)))]
        if len(df_invivo_kss) > 0:
            peak_int = df_invivo_kss.loc[df_invivo_kss["count"].idxmax()]["count"]
            peak_pos = df_invivo_kss.loc[df_invivo_kss["count"].idxmax()]["position"]
            newline[11] = peak_int
            newline[12] = True

        if newline[10]:
            upstream_seq = oligo.seq[int(peak_pos) + 58 - 37:int(peak_pos) + 58].__str__()

            for i in range(len(upstream_seq)-2):
                codon = upstream_seq[i:i+3]
                found_motif = re.search(re_motiv, codon)
                if found_motif is not None:
                    fm = found_motif[0]
                    motif_df.loc[i, fm] = motif_df.loc[i, fm] + 1
        else:
            upstream_seq = oligo.seq[58 + 15 - 37: 58+15].__str__()
        newline[16] = upstream_seq
        # calculate MFE of upstream_seq:
        mfe_seq = oligo.seq[58 - 30: 58 + 15].__str__()
        cmd = " ".join(["echo", "'" + mfe_seq + "'", "|", "RNAfold", "|", "grep", "-Ev", "'A|U|G|C'", "|",
                        "sed", "-E", "'s/.* \\(([^\\)]+).*$/\\1/'"])
        ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        MFE = float(ps.communicate()[0].decode()[:-2].strip())
        newline[17] = MFE


        out_df = out_df.append(pd.DataFrame(columns=out_df.columns, data=[newline]))
        #out_tsv.write("\t".join(newline)+"\n")

motif_df_normalized = motif_df / motif_df.sum().sum()
motif_df.to_csv("../data/motif_df.tsv", sep="\t")
motif_df_normalized.to_csv("../data/motif_df_normalizd.tsv", sep="\t")


out_df.to_csv("../analysis/annotated_sites.csv")

print("done with first comp.")


# now we want to identify alt. TIS:
tot_length = len(out_df)
found_in_our_ds = sum(out_df["peak_found_Hör"])
found_in_other_ds = sum(out_df["peak_found_meydan"])
found_both = sum(out_df["peak_found_meydan"] & out_df["peak_found_Hör"])
start_codons = ["ATG", "GTG", "TTG", "CTG", "ATC", "ATT"]
stop_codons = ["TGA", "TAA", "TAG"]
sc_re = r"(ATG)|(GTG)|(TTG)|(CTG)|(ATC)|(ATT)"

# initiate DF for identification of novel TIS:
out_df_ss = pd.DataFrame(columns=["locus_tag", "gene", "start_position", "start_codon", "pos_from_annot_start",
                                  "in_frame", "strand", "peak", "nr_peaks_Hör", "peak_found_Hör", "peak_meydan",
                                  "peak_found_meydan", "pred_RBS_region", "p_value", "sequence_start", "upstream_seq",
                                  "MFE"])
# initiate DF for SD search:
sd_df = pd.DataFrame(columns=["location", "gene", "Gibbs"])


# add a table showing enriched start codons/sd motifs in alternative TIS:
motif_df_alt = pd.DataFrame(np.zeros((35, 8)), columns=["ATG", "GTG", "TTG", "CTG", "ATC", "ATT", "AGG", "GGA"])

# screen for TIS
for oligo in SeqIO.parse(in_oligos_fasta, "fasta"):
    if oligo.id.startswith("b"):
        id = oligo.id
        strand = id.split("_")[3]
        locus_tag = id.split("_")[0]
        gene = id.split("_")[1]
        start_pos = int(id.split("_")[2])


        dfs = [data_frame[data_frame["gene"] == id] for data_frame in wig_oligos]
        # get total read counts for relative density:
        tot_counts = [sum(data_frame["count"]) for data_frame in dfs]
        for n in range(len(dfs)):
            dfs[n]["norm_count"] = dfs[n]["count"] / sum(dfs[n]["count"])

        # check for start sites which have norm_value (rel_density) more than 0.1 and >5 cpm etc.:
        df_ss = [df[(df["count"] > threshold) & (df["norm_count"] > 0.1) &
                    df.position.isin(list(range(30, 70)) + list(range(79, 210)))] for df in dfs]
        df_ss = [df.reset_index(drop=True) for df in df_ss]
        rec_peaks = 0
        # loop through replicates, get peaks
        for df in range(len(df_ss)):
            df_ss[df]["norm_counts"] = df_ss[df]["count"] / tot_counts[df]
            if len(df_ss[df]) > 0:
                nr = 1
                # just get the peaks ,if within 5 nt there's another peak, the highest is taken:
                while nr < len(df_ss[df]):
                    row = df_ss[df].iloc[nr]
                    prev_row = df_ss[df].iloc[nr - 1]
                    if row["position"] - prev_row["position"] < 5:
                        index_del = df_ss[df].iloc[nr - 1:nr + 1]["count"].idxmin()
                        df_ss[df] = df_ss[df].drop(index_del)
                        df_ss[df].reset_index(drop=True)
                    else:
                        nr += 1
                df_ss[df] = df_ss[df].reset_index(drop=True)

        for df in range(len(df_ss)):
            for peak in range(len(df_ss[df])):
                peakpos = df_ss[df].loc[peak]["position"]
                region_sc = oligo.seq[peakpos-18:peakpos-11]
                # remove if there's no start codon in region before:
                if not any(sc in region_sc for sc in start_codons):
                    df_ss[df] = df_ss[df].drop(peak)

            df_ss[df] = df_ss[df].reset_index(drop=True)
            df_ss[df]["position"] = df_ss[df]["position"] - 58

            # write to resulting df (if theres start codon)!
            for peak in range(len(df_ss[df])):
                peakpos = df_ss[df].loc[peak]["position"]
                region_sc = oligo.seq[peakpos + 58 - 18:peakpos + 58 - 11]

                start_c = re.search(sc_re, str(region_sc))[0]
                pos_from_peak = re.search(sc_re, str(region_sc)).start() + peakpos - 18

                startsite_pos = peakpos + 58 - 18 + re.search(sc_re, str(region_sc)).start()
                region_sd = oligo.seq[startsite_pos-25: startsite_pos]
                sequence_start = oligo.seq[startsite_pos-20: startsite_pos + 10]
                upstream_seq = oligo.seq[int(peakpos) + 58 - 37:int(peakpos) + 58].__str__()
                # calculate MFE of upstream_seq:
                mfe_seq = oligo.seq[int(peakpos) + 58 - 45:int(peakpos) + 58].__str__()
                cmd = " ".join(["echo", "'" + mfe_seq + "'", "|", "RNAfold", "|", "grep", "-Ev", "'A|U|G|C'", "|",
                                "sed", "-E", "'s/.* \\(([^\\)]+).*$/\\1/'"])
                ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                MFE = float(ps.communicate()[0].decode()[:-2].strip())





                # calculate sd-binding regions:
                fasta_pathname = "../data/sd_search/freescan_files_alt_ss/seq" + "_" + locus_tag + "_" + str(peakpos)
                with open(fasta_pathname + ".fasta", "w") as fasta:
                    fasta.write(">Seq_SD\n")
                    fasta.write(region_sd.__str__())
                sd_bits = str.encode(region_sd.__str__())
      #          for i in range(10):
      #              with open("../data/sd_search/freescan_shuffled_alt/seq" + "_" + locus_tag + "_" + str(
      #                      i) + ".fasta",
       #                       "w") as fasta:
       #                 fasta.write(">Seq_SD\n")
       #                 fasta.write(shuffle(sd_bits, 2).decode("utf-8"))

                # freescan determine delta gibbs etc:
                ASD = "auuccuccacuag"
                output_freescan = subprocess.run(["free_scan", "-e", ASD, fasta_pathname + ".fasta"],
                                                 capture_output=True).stdout
                # get list with all free energies:
                list_values = output_freescan.decode("utf-8").split("\n")
                list_values = list_values[:-2]
                list_values = [float(i) for i in list_values]
                lowest_val = min(list_values)
                ind_lowest_val = list_values.index(lowest_val)
                best_binding_region = region_sd[ind_lowest_val:ind_lowest_val + len(ASD)]
                print(locus_tag)
                with open(fasta_pathname + ".fasta", "w") as fasta:
                    fasta.write(">Seq_SD\n")
                    fasta.write(best_binding_region.__str__())

                shuffled_gibbs = np.array([])
                for i in range(10):
                    bits_best_binding = str.encode(best_binding_region.__str__())
                    shuffled = shuffle(bits_best_binding, 2).decode("utf-8")
                    with open(fasta_pathname + ".fasta", "w") as fasta:
                        fasta.write(">Seq_SD\n")
                        fasta.write(shuffled + "\n")
                    output_freescan = subprocess.run(["free_scan", "-e", "auuccuccacuag", fasta_pathname + ".fasta"],
                                                     capture_output=True).stdout  # freescan
                    # get list with all free energies:
                    val = float(output_freescan.decode("utf-8")[:-2])
                    newline_gibbs = [locus_tag, gene, shuffled, ind_lowest_val - 13, True, val]
                    # gibbs_df = pd.concat([gibbs_df, pd.DataFrame(columns=gibbs_df.columns, data=[newline_gibbs])])
                    shuffled_gibbs = np.append(shuffled_gibbs, val)
                # end freescan determine delta gibbs etc:

                less_eq = np.less_equal(shuffled_gibbs, lowest_val)
                pval = sum(less_eq) / len(shuffled_gibbs)

                newline_gibbs = [id, gene, best_binding_region.__str__(), ind_lowest_val - 13, False, lowest_val,
                                 pval]
                gibbs_df = pd.concat([gibbs_df, pd.DataFrame(columns=gibbs_df.columns, data=[newline_gibbs])])

                newline = [locus_tag, gene, start_pos, start_c, 0, "in_frame", strand, 0, 1, False, 0, False,
                best_binding_region.__str__(), pval, sequence_start.__str__(), upstream_seq, MFE]
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


                for i in range(len(upstream_seq) - 2):
                    codon = upstream_seq[i:i + 3]
                    found_motif = re.search(re_motiv, codon)
                    if found_motif is not None:
                        fm = found_motif[0]
                        motif_df_alt.loc[i, fm] = motif_df_alt.loc[i, fm] + 1


motif_df_alt.to_csv("../data/motif_df_alt.tsv", sep="\t")


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
        # drop lower one:
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
#sum(reduced_df_ss["peak_found_Hör"])


#meydan_genes = pd.read_csv("../data/wigglefiles/meydan_alt_tis.csv")














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


# create random sequences:
alphabet = ["A", "U", "G", "C"]


with open("../analysis/revision_06_2022/random_seqs_MFE.csv", "w") as file:
    file.write("MFE,peak_type\n")
    with open("../data/shuffled_seqs_esl.txt") as f_shuff:
        for line in f_shuff:
            print(line)
            cmd = " ".join(["echo", "'" + line + "'", "|", "RNAfold", "|", "grep", "-Ev", "'A|U|G|C'", "|",
                            "sed", "-E", "'s/.* \\(([^\\)]+).*$/\\1/'"])
            ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            MFE = float(ps.communicate()[0].decode()[:-2].strip())
            file.write(str(MFE) + ",random sequence" + "\n")


    for oligo in SeqIO.parse("../data/reference_sequences/NC_000913.3.fasta", "fasta"):
        for i in range(10000):
            s = oligo.seq.__str__()
            length = len(s)
            start_rn = random.randint(0, length-46)
            ranseq = s[start_rn:start_rn + 45]

            bits_ranseq = str.encode(ranseq.__str__())
            shuffled = shuffle(bits_ranseq, 2).decode("utf-8")
            print(ranseq, shuffled, "shuffled")
            cmd = " ".join(["echo", "'" + shuffled + "'", "|", "RNAfold", "|", "grep", "-Ev", "'A|U|G|C'", "|",
                            "sed", "-E", "'s/.* \\(([^\\)]+).*$/\\1/'"])
            ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            MFE = float(ps.communicate()[0].decode()[:-2].strip())
            file.write(str(MFE) + ",random sequence" + "\n")



    #for i in range(10000):
    #    sequence = ''.join(random.choice(alphabet) for i in range(45))
    #    cmd = " ".join(["echo", "'" + sequence + "'", "|", "RNAfold", "|", "grep", "-Ev", "'A|U|G|C'", "|",
    #                    "sed", "-E", "'s/.* \\(([^\\)]+).*$/\\1/'"])
    #    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #    MFE = float(ps.communicate()[0].decode()[:-2].strip())
    #    file.write(str(MFE) + ",random sequence" + "\n")


