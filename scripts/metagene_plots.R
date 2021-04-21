library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)



wiggle_R1 <- read.delim(
  "../jj/Documents/riboseq_experiments/oligopool_02_2021/data/wigglefiles/01_no_PNA.bam_3primeend_mod.wig",
  col.names = c("gene", "start", "end", "depth"))
wiggle_R1$rep <- "R1"

wiggle_R2 <- wiggle <- read.delim(
  "../jj/Documents/riboseq_experiments/oligopool_02_2021/data/wigglefiles/07_no_PNA.bam_3primeend_mod.wig",
  col.names = c("gene", "start", "end", "depth"))
wiggle_R2$rep <- "R2"

wiggle_R3 <- wiggle <- read.delim(
  "../jj/Documents/riboseq_experiments/oligopool_02_2021/data/wigglefiles/13_no_PNA.bam_3primeend_mod.wig",
  col.names = c("gene", "start", "end", "depth"))
wiggle_R3$rep <- "R3"

wiggles <- rbind(wiggle_R1, wiggle_R2, wiggle_R3)
  
wiggles$distance_to_start <- wiggles$start - 58

# plot for each replicate:
wiggles$depth_norm <- wiggles$depth/1000000
wiggles %>% group_by(rep, distance_to_start) %>%
  summarise(reads_norm = sum(depth_norm)) %>%
  ggplot(aes(x=distance_to_start, y=reads_norm, group=rep, color=rep)) + 
  scale_x_continuous(limits = c(-10, 50), breaks = seq(-10, 50, by=5)) +
  geom_line(size=1.2) + theme() + theme_pubr()

# plot for whole data:
wiggles$depth_norm <- wiggles$depth_norm/3
wiggles %>% 
  group_by(distance_to_start) %>% 
  summarise(reads_norm = sum(depth_norm)) %>%
  ggplot(aes(x=distance_to_start, y=reads_norm)) + 
  scale_x_continuous(limits = c(-10, 50), breaks = seq(-10, 50, by=5)) +
  geom_line(size=1.2) + 
  geom_line() + theme_pubr() 


  wiggle_paper_F <- read.delim(
  "../jj/Documents/riboseq_experiments/oligopool_02_2021/data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_F_new.wig",
  col.names = c("pos", "depth"), skip = 2)
wiggle_paper_R <- read.delim(
  "../jj/Documents/riboseq_experiments/oligopool_02_2021/data/wigglefiles/wiggles-papercomparison_2/GSM3455900_RET_BWK_U00096_3_R_new.wig",
  col.names = c("pos", "depth"), skip = 2)
head(wiggle_paper_R)

metagene_F <- read.delim(
  "../jj/Documents/riboseq_experiments/oligopool_02_2021/data/wigglefiles/wiggles-papercomparison_2/metagene_F.tab",
  skip = 3, header = FALSE)
head(metagene_F)
dim(metagene_F)
metagene_F[is.na(metagene_F)] = 0 

metagene_R <- read.delim(
  "../jj/Documents/riboseq_experiments/oligopool_02_2021/data/wigglefiles/wiggles-papercomparison_2/metagene_R.tab",
  skip = 3, header = FALSE)
head(metagene_R)
dim(metagene_R)
metagene_R[is.na(metagene_R)] = 0 
metagene_R = - metagene_R


mf_df <- data.frame(distance_to_start=c(seq(-11,-1,1), seq(1,51,1)), metagene_F=colSums(metagene_F),
                    metagene_R = colSums(metagene_R))
mf_df$distance_to_start = mf_df$distance_to_start-1
mf_df$norm_reads = (mf_df$metagene_F + mf_df$metagene_R) / sum(mf_df$metagene_F + mf_df$metagene_R)

ggplot(mf_df, aes(x=distance_to_start, y=norm_reads)) + 
  scale_x_continuous(limits = c(-10, 50), breaks = seq(-10, 50, by=5)) +
  geom_line(size=1.2) + 
  geom_line() + theme_pubr()

#merge 2 dataframes:
exp_tib <- tibble(mf_df[,c(1,4)])
exp_tib$exp <- "Meydan et al."

exp_ours <- wiggles %>% 
  group_by(distance_to_start) %>% 
  summarise(norm_reads = sum(depth_norm))
exp_ours$exp <- "Oligo pool"

total_df <- bind_rows(exp_tib, exp_ours)


ours <- wiggles %>% 
  group_by(distance_to_start) %>% 
  summarise(reads_norm = sum(depth_norm)) %>%
  ggplot(aes(x=distance_to_start, y=reads_norm)) + 
  scale_x_continuous(limits = c(-10, 50), breaks = seq(-10, 50, by=5)) +
  scale_y_continuous(limits = c(0, 0.3))+
  geom_line(size=1.2) + 
  geom_line() + theme_pubr() + labs(x="distance from start codon (nt)", y = "Normalized footprints")

meydan <- exp_tib %>% ggplot(aes(x=distance_to_start, y=norm_reads)) + 
  scale_x_continuous(limits = c(-10, 50), breaks = seq(-10, 50, by=5)) +
  scale_y_continuous(limits = c(0, 0.3))+
  geom_line(size=1.2) + 
  geom_line() + theme_pubr() + labs(x="distance from start codon (nt)", y = "Normalized footprints")


both <- total_df %>% ggplot(aes(x=distance_to_start, y=norm_reads, group=exp, color=exp)) + 
  geom_line(size=1.2) +
  scale_x_continuous(limits = c(-10, 50), breaks = seq(-10, 50, by=5)) +
  scale_y_continuous(limits = c(0, 0.3))+
  geom_line() + theme_pubr() + labs(x="distance from start codon (nt)", y = "Normalized footprints")


ours_rep <- wiggles %>% group_by(rep, distance_to_start) %>%
  summarise(reads_norm = sum(depth_norm)*3) %>%
  ggplot(aes(x=distance_to_start, y=reads_norm, group=rep, color=rep)) + 
  scale_x_continuous(limits = c(-10, 50), breaks = seq(-10, 50, by=5)) +
  scale_y_continuous(limits = c(0, 0.3))+
  geom_line(size=1.2) +  theme_pubr() + labs(x="distance from start codon (nt)", y = "Normalized footprints")



grid.arrange(both +scale_color_grey()+ theme_pubr(),nrow=1)
grid.arrange(meydan +scale_color_grey()+ theme_pubr(),nrow=1)
grid.arrange(ours +scale_color_grey()+ theme_pubr(),nrow=1)
grid.arrange(ours_rep +scale_color_viridis_d()+ theme_pubr(),nrow=1)

svg("../jj/Documents/riboseq_experiments/oligopool_02_2021/analysis/metagene_analysis/comparison.svg")
grid.arrange(both +scale_color_grey()+ theme_pubr(),nrow=1)
dev.off()

svg("../jj/Documents/riboseq_experiments/oligopool_02_2021/analysis/metagene_analysis/meydan.svg")
grid.arrange(meydan +scale_color_grey()+ theme_pubr(),nrow=1)
dev.off()

svg("../jj/Documents/riboseq_experiments/oligopool_02_2021/analysis/metagene_analysis/oligopool_avg.svg")
grid.arrange(ours +scale_color_grey()+ theme_pubr(),nrow=1)
dev.off()

svg("../jj/Documents/riboseq_experiments/oligopool_02_2021/analysis/metagene_analysis/oligopool_all_replicates.svg")
grid.arrange(ours_rep +scale_color_viridis_d()+ theme_pubr(),nrow=1)
dev.off()










