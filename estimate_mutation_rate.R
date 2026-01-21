# ---
# title: "TR mutation rate estimates"
# date: "2025-05-12"
# author: Sudmant, Peter
# ---

library(tidyr)
library(dplyr)
library(ggplot2)

# Load TR genotype/summary data
# gt = copy number
t <- read.table("human_chimp_sumstats_subset.txt.gz", header = F, sep="\t", col.names = c("tr_id", "indiv", "motif_len", "bp_len", "motif", "gt", "species"))

# Load Porubsky data (tab "TR - CHM13" from the Supplemental Table 10)
#t_hsa <- read.table("/Users/petersudmant/Documents/science/sudmantlab/projects/panpan/analyses/TR/input/41586_2025_8922_MOESM12_ESM_TR.txt", header = TRUE, sep="\t")
t_hsa <- read.table("41586_2025_8922_MOESM12_ESM_TR.txt", header=T, sep="\t")

# Load denominator information per sample for mutation rate calculations
t_hsa_denom <- read.table("K1463.CHM13v2.denominator_per_sample.tsv", header = T, sep="\t")

# Quick look at the total passing loci per sample
t_hsa_denom %>%
  group_by(sample_id) %>%
  summarize(n_passing_loci = sum(n_passing_loci))

# Calculate mutation rates from Porubsky dataset
# Filter TRs passing all quality filters, then calculate bp, allele, and count changes per sample
porubsky_rates <- t_hsa %>%
  filter(pass_all_filters=="YES") %>%
  dplyr::select(sample_id, likely_denovo_size, min_motif_size_in_locus) %>%
  group_by(min_motif_size_in_locus, sample_id) %>%
  summarize(bp_change = sum(abs(likely_denovo_size)),
            allele_change = sum(abs(likely_denovo_size/min_motif_size_in_locus)),
            count_change = n()) %>%
  # Join denominator info to normalize rates
  right_join(t_hsa_denom, by=c("sample_id","min_motif_size_in_locus")) %>%
  mutate(bp_change = ifelse(is.na(bp_change),0,bp_change),
         allele_change = ifelse(is.na(allele_change),0,allele_change),
         count_change = ifelse(is.na(count_change),0,count_change)) %>%
  # Sum changes across samples per motif size
  group_by(min_motif_size_in_locus) %>%
  summarize(bp_change = sum(bp_change),
            count_change = sum(count_change),
            allele_change = sum(allele_change),
            n_passing_loci = sum(n_passing_loci)) %>%
  # Calculate mutation rates per locus, dividing by 2 for diploid genomes
  mutate(mu_count = (count_change/n_passing_loci)/2,
         mu_allele = (allele_change/n_passing_loci)/2,
         mu_bp = (bp_change/n_passing_loci)/2) %>%
  filter(mu_count != 0) %>%
  dplyr::select(min_motif_size_in_locus, mu_count, mu_allele, mu_bp) %>%
  mutate(mu = mu_count, 
         type = "Porubsky",
         motif_len = min_motif_size_in_locus) %>%
  dplyr::select(motif_len, mu, mu_allele, type)

# Pivot Porubsky rates to long format for plotting
porubsky_rates <- porubsky_rates %>%
  mutate(Porubsky_per_event = mu,
         Porubsky_per_allele = mu_allele) %>%
  dplyr::select(motif_len, Porubsky_per_event, Porubsky_per_allele) %>%
  pivot_longer(cols = c("Porubsky_per_event","Porubsky_per_allele"), 
               names_to = "type", values_to = "mu")

# Subset data for testing
top_k_trs <- unique(t$tr_id)[1:100000]

t_data <- t %>% 
  filter(tr_id %in% top_k_trs) %>% 
  group_by(tr_id, species) %>%
  # Compute number of individuals, variance, and mean TR copy number per TR
  mutate(n_indivs = n(), v = var(gt), mu_gt = mean(gt), mu_bp = mean(bp_len)) %>% 
  ungroup() %>%
  # Summarize allele frequencies per TR
  group_by(tr_id, species, gt, n_indivs, v, mu_gt, mu_bp, motif_len) %>%
  summarize(n = n()) %>%
  mutate(f = n/n_indivs)

# Calculate within species variance metrics (V1, V1_dsw)
# For chimps
t_PTR <- t_data %>% filter(species=="pantro")
t_PTR <- full_join(t_PTR, t_PTR %>% ungroup() %>% dplyr::select(tr_id, f, gt), by="tr_id") %>%
  # V_1 = Nei's within species variance weighted by the allele frequencies
  mutate(V_i = ((gt.x-gt.y)^2)*f.x*f.y,
         V_dsw_i = (abs(gt.x-gt.y))*f.x*f.y) %>%
  group_by(tr_id) %>%
  summarize(V1 = sum(V_i),
            V1_dsw = sum(V_dsw_i)) %>%
  inner_join(t_PTR, by="tr_id")

# For humans
t_HSA <- t_data %>% filter(species=="homo")
t_HSA <- full_join(t_HSA, t_HSA %>% ungroup() %>% dplyr::select(tr_id, f, gt), by="tr_id") %>%
  mutate(V_i = ((gt.x-gt.y)^2)*f.x*f.y,
         V_dsw_i = (abs(gt.x-gt.y))*f.x*f.y) %>%
  group_by(tr_id) %>%
  summarize(V1 = sum(V_i),
            V1_dsw = sum(V_dsw_i)) %>%
  inner_join(t_HSA, by="tr_id")

# Join species data for between species comparisons
t_joined <- full_join(t_PTR, t_HSA, by = "tr_id", suffix = c("_PTR", "_HSA"))

# Compute between-species variance metrics (ASD, du2, dsw)
du_stats <- t_joined %>%
  # ASD_1 = allele-specific squared diff. between species, weighted by allele freq.
  mutate(ASD_i = ((gt_PTR-gt_HSA)^2)*f_PTR*f_HSA,
         ASD_dsw_i = (abs(gt_PTR-gt_HSA))*f_PTR*f_HSA) %>%
  group_by(tr_id, v_PTR, v_HSA, V1_PTR, V1_HSA, V1_dsw_PTR, V1_dsw_HSA, motif_len_HSA, mu_gt_HSA, mu_gt_PTR, mu_bp_HSA, mu_bp_PTR) %>%
  # Summing over all alleles - total divergence between species per TR
  summarize(ASD = sum(ASD_i),
            ASD_dsw = sum(ASD_dsw_i)) %>%
  # du2 = measure of between species divergence in copy numbers, correcting for within species variance (Goldstein et al. 1995)
  mutate(du2 = ASD - (v_PTR + v_HSA),
         du2_biased = ASD - ((V1_HSA + V1_PTR)/2),
         # dsw = distance-based version of du2.
         dsw = ASD_dsw - ((V1_dsw_PTR + V1_dsw_HSA)/2),
         delta_mu = abs(mu_gt_HSA - mu_gt_PTR),
         delta_bp = abs(mu_bp_HSA - mu_bp_PTR))

# Nei's genetic distance between species (Ds)- measures the overall allele freq. divergence between species
# Ds ignores the length of the TR
nei_stats <- full_join(t_PTR %>% ungroup() %>% dplyr::select(tr_id,species,f,gt), 
                      t_HSA %>% ungroup() %>% dplyr::select(tr_id,species,f,gt,motif_len), 
                      by = c("tr_id","gt"), suffix = c("_PTR", "_HSA")) %>%
  mutate(f_PTR = ifelse(is.na(f_PTR),0,f_PTR),
         f_HSA = ifelse(is.na(f_HSA),0,f_HSA),
         # sum of squared allele freq. within chimps
         Jptr2_i = f_PTR*f_PTR,
         # sum of squared allele freq. within human
         Jhsa2_i = f_HSA*f_HSA,
         # sum of allele freq. between species
         Jhsa_ptr_i = f_HSA*f_PTR) %>%
  group_by(tr_id) %>%
  summarize(Jptr2 = sum(Jptr2_i),
            Jhsa2 = sum(Jhsa2_i),
            Jhsa_ptr = sum(Jhsa_ptr_i),
            motif_len = first(na.omit(motif_len)))

# summarise
nei_stats_sum <- nei_stats %>% 
                group_by(motif_len) %>%
                summarize(Jxy = mean(Jhsa_ptr), 
                          Jx = mean(Jptr2),
                          Jy = mean(Jhsa2)) %>%
                mutate(Ds = -log(Jxy/sqrt(Jx*Jy)))

# Mutation rate estimates normalized by branch lengths
branch_len_HSA = 6.2e6/28
branch_len_PTR = 6.2e6/25
branch_len_total = branch_len_PTR + branch_len_HSA

filt_du_stats <- du_stats %>% filter(!is.na(ASD))

# Summarize mutation rates per motif length for plotting
summarized_du_stats <- filt_du_stats %>% 
                        group_by(motif_len_HSA) %>%
                        summarize(mu = mean(du2)/branch_len_total) %>%
                        mutate(motif_len = motif_len_HSA, type="du2") %>%
                        dplyr::select(motif_len, mu, type)

summarized_dsw_stats <- filt_du_stats %>% 
                        group_by(motif_len_HSA) %>%
                        summarize(mu = mean(dsw)/branch_len_total) %>%
                        mutate(motif_len = motif_len_HSA, type="dsw") %>%
                        dplyr::select(motif_len, mu, type)

summarized_delta_mu_stats <- filt_du_stats %>% 
                        group_by(motif_len_HSA) %>%
                        summarize(mu = mean(delta_bp)/branch_len_total) %>%
                        mutate(motif_len = motif_len_HSA, type="delta_bp") %>%
                        dplyr::select(motif_len, mu, type)

summarized_nei_stats <- nei_stats_sum %>%
                        mutate(mu = Ds / branch_len_total,
                               type="Ds") %>%
                        dplyr::select(motif_len, mu, type)

all_rates <- rbind(summarized_du_stats, summarized_dsw_stats, summarized_delta_mu_stats, summarized_nei_stats, porubsky_rates)

all_rates_filtered <- all_rates %>% filter(!type %in% c("dsw", "Ds", "delta_bp"))

# Plot mutation rates per motif length
g <- ggplot(all_rates %>% filter(motif_len<200))
g + geom_point(aes(x=motif_len, y=mu, color=type)) +
    geom_line(aes(x=motif_len, y=mu, color=type, size=grepl("Porubsky", type), linetype=grepl("Porubsky", type))) +
theme_bw() +
    scale_y_log10("Mutation rate / generation") +
    scale_x_continuous("Motif length (bp)") +
    scale_color_viridis_d() +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_size_manual(values = c(0.5,1.5)) +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=18),
axis.title.x=element_text(size= 16),
axis.title.y=element_text(size = 16),
axis.text.x=element_text(colour="black", size = 16),
axis.text.y=element_text(colour="black", size = 16))

# Example for Huntington’s disease TR
HTT = t %>% 
  filter(tr_id=="HTT") %>%
  arrange(species, gt) %>%
  mutate(y=row_number())

g = ggplot(HTT)

g+geom_segment(aes(x=0,xend=gt,y=y,yend=y,color=species),size=2) +
theme_bw() +
theme(axis.line = element_line(linewidth=1, colour = "black"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(), panel.background = element_blank(),
text = element_text(size=21),
axis.title.x=element_text(size= 20),
axis.title.y=element_text(size = 20),
axis.text.x=element_text(colour="black", size = 20),
axis.text.y=element_text(colour="black", size = 20))

d = 0.018588
s = 1-d
p_len = 20
s^(2*p_len)
