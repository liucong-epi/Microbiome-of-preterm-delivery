# ==============================================================================
# Preterm/Term Birth Microbiota Analysis Pipeline
# Author: [Liu Cong]
# Date: [2026-03]
# Description: Analysis of vaginal/placental microbiota composition, diversity, 
#              and differential abundance between preterm and term birth groups
# ==============================================================================

# 1. Environment Setup ==============================

# Clear workspace and free up memory
rm(list = ls()); gc()

# Load required R packages (install missing packages automatically)
pacman::p_load(
  dplyr, purrr, tidyverse, psych, stringr, reshape2, pheatmap, openxlsx,
  vegan, qiime2R, tibble, ggpubr, tableone, phyloseq, readxl, tidyr, forcats,
  MicrobiotaProcess, patchwork, data.table, magrittr, plyr, circlize, ComplexHeatmap,
  VennDiagram, eulerr, Hmisc, corrplot, ggsignif, ggrepel, cowplot, egg
)

# Define custom color palette for visualization
colors <- c(
  "#E41A1C", "#1E90FF", "#FF8C00", "#4DAF4A", "#40E0D0", "#984EA3", "#FFC0CB",
  "#00BFFF", "#FFDEAD", "#90EE90", "#EE82EE", "#00FFFF", "#F0A3FF", "#0075DC",
  "#993F00", "#4C005C", "#2BCE48", "#FFCC99", "#808080", "#94FFB5", "#8F7C00",
  "#9DCC00", "#426600", "#FF0010", "#5EF1F2", "#00998F", "#740AFF", "#990000", "#FFFF00"
)

# 2. Data Preprocessing ==============================

# 2.1 Vaginal Metagenomic Data ------------------------------

# Import raw vaginal metagenomic data
dat_kra <- read.csv("data/preterm_study_vaginal_metagenomic_data.csv", row.names = 1, header = TRUE)

# Filter low-abundance features 
# Keep features present in ≥5% samples (4 samples) with ≥0.01% relative abundance
tmp <- t(dat_kra)
s_filter <- qiime2R::filter_features(tmp, minsamples = 4, minreads = sum(tmp)*0.0001) %>% 
  t() %>% as.data.frame()

# Normalization
# 1. Relative abundance normalization
s_relative <- funrar::make_relative(as.matrix(s_filter)) %>% as.data.frame()
# write.csv(s_relative,"intermediate_file/preterm_study_vaginal_metagenomic_data_relative.csv",row.names = TRUE)

# Import group metadata
group <- read.csv("data/group.csv")  

# 2.2 Placenta 16S rRNA Data ------------------------------
# Import placenta 16S data (genus/phylum/OTU level)
genus <- read.csv("data/placenta/data_placenta_16S_genus.csv", row.names = 1)
phylum <- read.csv("data/placenta/data_placenta_16S_phylum.csv", row.names = 1)
otu <- read.csv("data/placenta/data_placenta_16S_otu.csv", row.names = 1)

# Import placenta metadata and taxonomy
meta_mat <- read.csv("data/placenta/data_placenta_meta_mat.csv", row.names = 1) 
tax_mat <- read.csv("data/placenta/data_placenta_16S_tax_mat.csv", row.names = 1)

# Quality control for placenta data
# Remove low-abundance ASVs (considered as contamination, < 100 reads)
genus2 <- genus
genus2[genus2 < 100] <- 0
genus2 %<>% .[, colSums(.) > 0]  # Remove genera with zero abundance across all samples
genus2 <- genus2 %>% t() %>% 
  filter_features(minsamples = 4) %>%  # Keep genera present in ≥5% samples
  t() %>% as.data.frame()

# Relative abundance for placenta genus data
g_relative <- funrar::make_relative(as.matrix(genus2))
g_relative[g_relative == "NaN"] <- 0

# Raw genus data for quality control comparison
genus_raw <- genus
genus_raw_relative <- funrar::make_relative(as.matrix(genus_raw))
genus_raw_relative[genus_raw_relative == "NaN"] <- 0
genus_raw_relative %<>% .[rowSums(.) > 0, ]  # Remove samples with zero total abundance

# 3. Core figure Analysis (fig 1) ==============================

# 3.1 fig 1A: Top 10 Species Composition (Circular Bar Plot) ==============================

# Extract top 10 species by total abundance
top10_species <- apply(s_relative, 2, sum) %>% 
  sort(decreasing = TRUE) %>% head(n = 10) %>% names()

genus_top <- s_relative %>% subset(select = all_of(top10_species))
Other <- (rep(1, nrow(genus_top)) - apply(genus_top, 1, sum))
g_top10 <- cbind(genus_top, as.data.frame(Other))

# Reorder samples by Lactobacillus abundance
g_top10 <- g_top10 %>% 
  dplyr::arrange(desc(`Lactobacillus.crispatus`), desc(`Lactobacillus.iners`))

# Reshape data for plotting
g_top10 <- g_top10 %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("species") %>% 
  reshape2::melt(id = "species")

# Merge with metadata
g_top10_2 <- merge(g_top10, group, by.x = "variable", by.y = "SampleID")

# Plot for term birth group
p1_s_term <- ggplot(subset(g_top10_2, group == "term"),
                    aes(variable, value, fill = species)) +
  scale_fill_manual(values = colors) +
  geom_col(position = 'fill', width = 0.9) +
  labs(x = '', y = '', fill = "species", title = "Term birth") +
  guides(colour = guide_legend(ncol = 2)) +
  theme_classic(base_size = 18) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = -5),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) + ylim(-0.50, 1)

# Plot for preterm birth group
p1_s_preterm <- ggplot(subset(g_top10_2, group == "preterm"),
                       aes(variable, value, fill = species)) +
  scale_fill_manual(values = colors) +
  geom_col(position = 'fill', width = 0.9) +
  labs(x = '', y = '', fill = "species", title = "Preterm birth") +
  guides(colour = guide_legend(ncol = 2)) +
  theme_classic(base_size = 18) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = -5),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) + ylim(-0.50, 1)

# Convert to polar coordinates (circular plot)
p1_s_term <- p1_s_term + coord_polar() 
p1_s_preterm <- p1_s_preterm + coord_polar() 

# Combine and save plots
cowplot::plot_grid(p1_s_term, p1_s_preterm, ncol = 2)
ggsave("fig 1A.pdf", width = 12, height = 7)

# 3.2 fig 1B: Alpha Diversity Analysis ------------------------------
# Import normalized vaginal data
aa2 <- read.csv("intermediate_file//preterm_study_vaginal_metagenomic_data_relative.csv", row.names = 1)

# Calculate alpha diversity indices
Shannon1 <- diversity(aa2, index = "shannon", MARGIN = 1, base = exp(1))
Simpson1 <- diversity(aa2, index = "simpson", MARGIN = 1, base = exp(1))

# Reshape data for plotting
index1 <- as.data.frame(cbind(Shannon1, Simpson1)) %>% 
  dplyr::rename("Shannon" = "Shannon1", "Simpson" = "Simpson1")
index1 <- merge(index1, group, by.x = 'row.names', by.y = 'SampleID')
index1.1 <- index1 %>% 
  tidyr::pivot_longer(cols = Shannon:Simpson, names_to = "diversity", values_to = "value")

# Boxplot with Wilcoxon test
p_alpha_s <- ggplot(index1.1, aes(x = group, y = value)) +
  stat_boxplot(geom = "errorbar", width = 0.1, size = 0.6) +
  geom_boxplot(aes(fill = group), outlier.colour = "white", size = 0.6) +
  theme_classic() +
  geom_jitter(width = 0.2, size = 0.3) +
  scale_fill_manual(values = c("#FF6A6A", "#00D7F5")) +
  geom_signif(
    comparisons = list(c("preterm", "term")),
    map_signif_level = TRUE,
    test = wilcox.test,
    size = 0.5, color = "black", textsize = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.09))) + 
  facet_grid(~diversity) +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(colour = "black", size = 8),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.title.x = element_text(colour = "black", size = 8),
    axis.title.y = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black")
  )

# Save plot
p_alpha_s
ggsave("fig 1B_alpha diversity.pdf", p_alpha_s, width = 2.5, height = 2.7)

# 3.3 fig 1C: PCoA (Beta Diversity) ------------------------------
# PCoA with Bray-Curtis distance and Adonis test
pca1.1 <- s_relative 
rownames(pca1.1) == group$SampleID

set.seed(123)
p_ado <- adonis2(pca1.1 ~ group$group, permutations = 999)
p_ado2 <- print(p_ado$`Pr(>F)`[1])
p_ado3 <- round(print(p_ado$R2[1]), 3)

# Calculate Bray-Curtis dissimilarity
sub_beta <- vegdist(pca1.1, method = 'bray')

# Perform PCoA
pcoa <- cmdscale(sub_beta, k = 3, eig = TRUE)
points <- as.data.frame(pcoa$points) 
eig <- pcoa$eig

# Merge with metadata
points <- cbind(points, group)
colnames(points)[1:3] <- c("PC1", "PC2", "PC3") 

# Plot PCoA with confidence ellipse
p_pcoa <- ggplot(points, aes(x = PC1, y = PC2, color = group)) +
  geom_point(alpha = 0.7, shape = 20) +
  labs(
    x = paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep = ""),
    y = paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep = ""),
    title = paste("PCoA Adonis p =", p_ado2, "R2 =", p_ado3)
  ) + 
  theme_bw() + 
  theme(
    legend.title = element_text(hjust = 0.2, size = 8),
    legend.key.height = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(colour = "black", size = 8),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.title.x = element_text(colour = "black", size = 8),
    axis.title.y = element_text(size = 8, colour = "black"),
    plot.title = element_text(size = 8, colour = "black")
  )

# Add 65% confidence ellipse
p_pcoa <- p_pcoa + stat_ellipse(
  aes(fill = group),
  alpha = 0.3,
  linetype = 'dashed',
  level = 0.65,
  type = "norm",
  geom = "polygon",
  show.legend = FALSE,
  linewidth = 0.8
)

# Save plot
p_pcoa
p_width <- 3.5
ggsave("fig 1C_PCoA.pdf", p_pcoa, width = p_width, height = p_width * 0.7)

# 3.4 fig 1D: BIRDMAn Differential Abundance ------------------------------
# Note: aa4 needs to be defined (import from BIRDMAn output)
aa4 <- read.csv("data/birdman_merged_result.csv")  # Uncomment and adjust path

# Effect size plot
custom_colors <- c("blue" = "#00BFC4", "red" = "#F8766D", "grey" = "grey")
ggplot(aa4, aes(x = reorder(feature, beta_var_mean), y = beta_var_mean, group = 1)) +
  geom_line(color = "gray") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.1, alpha = 0.6) +
  geom_pointrange(aes(ymin = beta_low, ymax = beta_up, color = lables), fatten = 0.2, size = 0.1) +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  labs(
    title = "Effect Sizes of Features",
    x = "Feature",
    y = expression(ln(frac(Preterm, Term)) + K)
  ) + 
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.title.x = element_text(size = 6, colour = "black"),
    axis.title.y = element_text(size = 6, colour = "black"),
    plot.title = element_text(size = 6, colour = "black", hjust = 0.5)
  ) +
  geom_text_repel(aes(label = lables2), na.rm = TRUE, size = 1.8, box.padding = 0.5)

# Save plot
ggsave("fig 1D_plot_preterm_result_ggplot.pdf", width = 3, height = 2)

# 3.5 fig 1E: Placenta Microbiota Composition (Relative Abundance) ------------------------------
# Prepare summary data
summary_dat <- genus2 %>% 
  t() %>% 
  colSums() %>% 
  as.data.frame() %>% 
  merge(meta_mat, by = "row.names") %>% 
  dplyr::rename("总计" = ".", 'SampleID' = 'Row.names')

summary_dat$Sampletype %<>% factor(
  levels = c("control", "true_sample"),
  labels = c("blank control", "CAM")
)
summary_dat$总计_log10 <- log10(summary_dat$总计 + 1)


# Prepare metadata for plotting
meta_mat2 <- meta_mat %>% 
  rownames_to_column(var = 'SampleID') %>% 
  merge(summary_dat[c(1:2)], by = 'SampleID') %>% 
  arrange(desc(总计))  # "总计" = total count (keep original for compatibility)

# Rename genera (extract genus name from g__ prefix)
g_relative1 <- g_relative %>% as.data.frame()
names(g_relative1) <- str_extract(names(g_relative1), "(?<=g__).+")
names(g_relative1)[1] <- 'Unknow_genus'

# Reshape data
g_relative1 <- g_relative1 %>%
  rownames_to_column(var = "SampleID") %>% 
  pivot_longer(!'SampleID', names_to = "genus", values_to = "count") %>% 
  as.data.frame()

# Merge with metadata and add Ureaplasma counts
data2 <- merge(g_relative1, meta_mat2, by = "SampleID") %>% 
  dplyr::arrange(desc(总计))
ureaplasma_counts <- data2 %>% 
  filter(genus == "Ureaplasma") %>% 
  select(SampleID, Ureaplasma_count = count)
data2 <- left_join(data2, ureaplasma_counts, by = "SampleID")

# Add blank placeholder for visualization clarity
data2_blank <- data2[data2$SampleID == 'FA14', ]
data2_blank[, c(1, 3:12)] <- NA
data2_blank$SampleID <- '1'
data2_blank$count <- 0

# Factorize genus levels
data2$genus %<>% factor(
  levels = c("Lactobacillus", "Prevotella", "Ureaplasma", "Unkonw_genus"),
  labels = c("Lactobacillus", "Prevotella", "Ureaplasma", "Unkonw_genus")
)

# Plot for preterm group
data2_preterm <- subset(data2, group == "preterm")
data2_preterm <- rbind(data2_preterm, data2_blank)
data2_preterm <- data2_preterm %>%
  arrange(desc(Ureaplasma_count), desc(总计 != 0), SampleID)
data2_preterm$SampleID <- factor(data2_preterm$SampleID, levels = unique(data2_preterm$SampleID))

p_preterm <- ggplot(data2_preterm, aes(SampleID, y = -count, fill = genus)) +
  scale_fill_manual(values = colors) +
  geom_col(position = 'stack', width = 0.95) +
  labs(x = '', y = '', fill = "genus", title = "Preterm") +
  theme_classic(base_size = 10) +
  scale_x_discrete(labels = c(1:40, '', '')) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 5),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) 

# Convert to polar coordinates
p_preterm <- p_preterm + coord_polar() 
p_preterm

# Save preterm plot
p_width <- 3
ggsave("fig 1E_species_composition_relative_preterm.pdf", width = p_width, height = p_width * 1)

# Plot for term group
data2_term <- subset(data2, group == "term")
data2_term <- rbind(data2_term, data2_blank)
data2_term <- data2_term %>%
  arrange(desc(Ureaplasma_count), desc(总计 != 0), SampleID)
data2_term$SampleID <- factor(data2_term$SampleID, levels = unique(data2_term$SampleID))

p_term <- ggplot(data2_term, aes(reorder(SampleID, -count), y = -count, fill = genus)) +
  scale_fill_manual(values = colors) +
  geom_col(position = 'stack', width = 0.9) +
  labs(x = '', y = '', fill = "genus", title = "Term birth") +
  theme_classic(base_size = 10) +
  scale_x_discrete(labels = c(1:39, '', '')) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 5),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

# Convert to polar coordinates
p_term <- p_term + coord_polar()
p_term

# Save term plot
p_width <- 3
ggsave("fig 1E_species_composition_relative_term.pdf", width = p_width, height = p_width * 1)

# Extract and save legend
legend <- cowplot::get_legend(
  ggplot(data2_term, aes(reorder(SampleID, -count), y = -count, fill = genus)) +
    scale_fill_manual(values = colors) +
    geom_col(position = 'stack', width = 0.9) +
    theme_classic(base_size = 10) + 
    guides(col = guide_legend(title = "Legend", nrow = 2)) +
    theme(
      legend.title = element_text(size = 10),
      legend.key.height = unit(0.5, "cm"),
      legend.text = element_text(size = 10),
      legend.direction = "horizontal"
    )
)
ggsave("fig 1E_species_composition_relative_legend.pdf", legend, width = 4, height = 3)

# 3.6 fig 1F: Key Taxa Abundance Comparison (Boxplot) ------------------------------
# Prepare genus-level data for comparison
meta_mat1 <- rownames_to_column(meta_mat, var = 'SampleID')
genus_com <- bind_cols(as.data.frame(g_relative), meta_mat1) %>% 
  subset(group != 'blank')
names(genus_com)[1:4] <- str_extract(names(genus_com)[1:4], "(?<=g__).+")
names(genus_com)[1] <- 'Unkonw_genus'

# 1. Prevotella comparison
tmp_plot <- ggplot2::ggplot(genus_com, aes(x = group, y = Prevotella * 100, fill = group)) +
  geom_boxplot(outlier.colour = "white", outlier.size = 1) +
  geom_jitter(aes(fill = group), width = 0.2, shape = 21, size = 2) +
  scale_fill_manual(values = c('#FF6A6A', '#00D7F5')) +
  scale_color_manual(values = c("black", "black")) +
  theme_bw() +
  labs(
    x = "", 
    y = expression("%" ~ italic("Prevotella") ~ "abundance"), 
    title = "", 
    fill = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 10, hjust = 0.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_signif(
    comparisons = list(c("preterm", "term")),
    map_signif_level = TRUE,
    test = wilcox.test,
    y_position = 110,
    size = 0.5, color = "black"
  ) +  
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.09))
  )

# Save Prevotella plot
tmp_plot
p_width <- 2
ggsave("fig 1F_relative_between_group_prevotella.pdf", width = p_width, height = p_width * 1.5)

# 2. Lactobacillus comparison
tmp_plot <- ggplot2::ggplot(genus_com, aes(x = group, y = Lactobacillus * 100, fill = group)) +
  geom_boxplot(outlier.colour = "white", outlier.size = 1) +
  geom_jitter(aes(fill = group), width = 0.2, shape = 21, size = 2) +
  scale_fill_manual(values = c('#FF6A6A', '#00D7F5')) +
  scale_color_manual(values = c("black", "black")) +
  theme_bw() +
  labs(
    x = "", 
    y = expression("%" ~ italic("Lactobacillus") ~ "abundance"), 
    title = "", 
    fill = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 10, hjust = 0.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_signif(
    comparisons = list(c("preterm", "term")),
    map_signif_level = TRUE,
    test = wilcox.test,
    vjust = -0.2,
    y_position = 110,
    size = 0.5, color = "black"
  ) +  
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.09))
  )

# Save Lactobacillus plot
tmp_plot
p_width <- 2
ggsave("fig 1F_relative_between_group_Lactobacillus.pdf", width = p_width, height = p_width * 1.5)

# 3. Ureaplasma comparison
tmp_plot <- ggplot2::ggplot(genus_com, aes(x = group, y = Ureaplasma * 100, fill = group)) +
  geom_boxplot(outlier.colour = "white", outlier.size = 1) +
  geom_jitter(aes(fill = group), width = 0.2, shape = 21, size = 2) +
  scale_fill_manual(values = c('#FF6A6A', '#00D7F5')) +
  scale_color_manual(values = c("black", "black")) +
  theme_bw() +
  labs(
    x = "", 
    y = expression("%" ~ italic("Ureaplasma") ~ "abundance"), 
    title = "", 
    fill = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 10, hjust = 0.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_signif(
    comparisons = list(c("preterm", "term")),
    map_signif_level = TRUE,
    test = wilcox.test,
    y_position = 110,
    size = 0.5, color = "black"
  ) +  
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.09))
  )

# Save Ureaplasma plot
tmp_plot
p_width <- 2
ggsave("fig 1F_relative_between_group_Ureaplasma.pdf", width = p_width, height = p_width * 1.5)

# 3.8 fig 1H: qPCR Validation of Ureaplasma parvum ------------------------------
# Import qPCR data
dat1 <- read.csv("data/placenta/data_CAM_Up_qPCR_final.csv")
dat1.1 <- pivot_longer(dat1, term:preterm, names_to = "group", values_to = "value")

# Paired Wilcoxon test
wilcox_result <- wilcox.test(dat1$term, dat1$preterm, paired = TRUE) 
wilcox_result %>% print()

# Plot boxplot
plot_tmp <- ggplot2::ggplot(dat1.1, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.colour = "black", outliers = TRUE, outlier.size = 1) +
  scale_fill_manual(values = c('#FF6A6A', '#00D7F5')) +
  scale_color_manual(values = c("black", "black")) +
  theme_bw() +
  labs(
    x = "", 
    y = expression("qPCR of " ~ italic("U.parvum") ~ " serovar 3"), 
    title = "", 
    fill = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 10, hjust = 0.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_signif(
    comparisons = list(c("preterm", "term")),
    map_signif_level = TRUE,
    test = wilcox.test,
    size = 0.5, color = "black",
    annotations = round(wilcox_result$p.value, 3) # ggplot无法使用paired-Wilcox，使用上面检验的结果
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.09)))

# Save plot with fixed panel size
plot_tmp
p_width <- 2
ggsave(
  filename = "fig 1H_CAM_qPCR_boxplot.pdf",
  plot = egg::set_panel_size(
    p = plot_tmp,
    width = unit(p_width, "in"),
    height = unit(p_width * 1.5, "in")
  ),
  width = p_width * 1.2,
  height = p_width * 1.5 * 1.2,
  units = c("in")
)

# 3.9 fig 1I: FISH Staining Quantification ------------------------------
# Import FISH data
fish <- read.csv("data/placenta/data_FISH_Density_DIO.csv")

# 1. Density mean plot
plot_f1 <- ggplot(aes(x = group, y = Density_mean, fill = group), data = fish) +
  geom_boxplot() + 
  theme_classic() +
  scale_fill_manual(values = c('#FF6A6A', '#00D7F5')) +
  scale_color_manual(values = c("black", "black")) +
  theme_bw() +
  labs(
    x = "", 
    y = expression("Density_mean of " ~ italic("Ureaplasma")), 
    title = "", 
    fill = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 10, hjust = 0.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_signif(
    comparisons = list(c("preterm", "term")),
    map_signif_level = TRUE,
    test = wilcox.test,
    textsize = 6, 
    size = 0.5, color = "black"
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.09)))

plot_f1

# 2. IOD (Integrated Optical Density) plot
plot_f2 <- ggplot(aes(x = group, y = IOD, fill = group), data = fish) +
  geom_boxplot() + 
  theme_classic() +
  scale_fill_manual(values = c('#FF6A6A', '#00D7F5')) +
  scale_color_manual(values = c("black", "black")) +
  theme_bw() +
  labs(
    x = "", 
    y = expression("IOD of " ~ italic("Ureaplasma")), 
    title = "", 
    fill = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 10, hjust = 0.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_signif(
    comparisons = list(c("preterm", "term")),
    map_signif_level = TRUE,
    test = wilcox.test,
    textsize = 6, 
    size = 0.5, color = "black"
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.09)))

plot_f2

# Combine Density and IOD plots
(plot_f1 + plot_f2) 
p_width <- 5.4
ggsave(
  "fig 1I_FISH_Density_mean_and_IOD.pdf",
  width = p_width,
  height = p_width * 0.65,
  units = c("in")
)

# 4. Supplementary figure Analysis (sfig 1) ==============================

# 4.1 sfig 1A: Co-occurrence Network (Cytoscape) ------------------------------
# Note: Network visualization was performed using Cytoscape software
# Output files are stored in the "cytoscape/" directory

# 4.2 sfig 1B: Salmonella Abundance in Placenta ------------------------------
salm <- read.csv("data/placenta//data_Salmonella_abundance.csv", row.names = 1)

ggplot(aes(x = 1, y = (Salmonella)), data = salm) +
  geom_boxplot(outlier.color = "white") +
  geom_jitter(width = 0.2, shape = 21, size = 2) +
  labs(
    x = '16S amplicon\nsequencing', 
    y = expression(italic("Salmonella typhimurium") ~ ' abundance')
  ) +
  theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Save plot
p_width <- 2
ggsave("sfig 1B-Salmonella_abundance.pdf", width = p_width, height = p_width * 1.5)

# 4.3 sfig 1C: ASV Count Comparison (Blank vs True Sample) ------------------------------

# Plot ASV count comparison
plot_sum_blank <- ggplot2::ggplot(summary_dat, aes(x = Sampletype, y = 总计, fill = Sampletype)) +
  geom_jitter(width = 0.25, size = 0.75) +
  labs(x = '', y = 'CAM ASV count') +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 10, hjust = 0.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save plot
plot_sum_blank
p_width <- 2
ggsave("sfig 1C_ASV count_between blank true sample.pdf", width = p_width, height = p_width * 1.25)


# 4.4 sfig 1D: PCoA by Delivery Mode ------------------------------

# Prepare data
tmp <- bind_cols(meta_mat["delivery_mode"], g_relative)
tmp %<>% filter(., rowSums(.[, c(2:ncol(.))]) > 0)

# Adonis test
set.seed(123)
p_ado <- adonis2(tmp[, c(2:ncol(tmp))] ~ tmp$delivery_mode, permutations = 999)
p_ado2 <- print(p_ado$`Pr(>F)`[1])
p_ado3 <- round(print(p_ado$R2[1]), 3)

# Calculate Bray-Curtis distance
sub_beta <- vegdist(tmp[, c(2:ncol(tmp))], method = 'bray')

# Perform PCoA
pcoa <- cmdscale(sub_beta, k = 3, eig = TRUE)
points <- as.data.frame(pcoa$points) 
eig <- pcoa$eig

# Merge with metadata
points <- cbind(points, tmp["delivery_mode"])
colnames(points)[1:3] <- c("PC1", "PC2", "PC3") 

# Plot PCoA
p_pcoa <- ggplot(points, aes(x = PC1, y = PC2, color = delivery_mode)) +
  geom_point(alpha = 0.6, size = 1) +
  labs(
    x = paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep = ""),
    y = paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep = ""),
    title = paste("PCoA Adonis p =", p_ado2, "R2 =", p_ado3)
  ) +
  theme_bw() + 
  theme(
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    axis.text.x = element_text(colour = "black", size = 8),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.title.x = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 10, colour = "black")
  )

# Add confidence ellipse
p_pcoa <- p_pcoa + stat_ellipse(
  aes(fill = delivery_mode),
  alpha = 0.3,
  linetype = 'dashed',
  level = 0.65,
  type = "norm",
  geom = "polygon",
  show.legend = FALSE,
  size = 0.8
)

# Save plot
p_pcoa
p_width <- 4
ggsave("sfig 1D_PCoA_between_delivery_mode.pdf", width = p_width, height = p_width * 0.7)

# 4.5 sfig 1E: Placenta Microbiota (Absolute Abundance) ------------------------------
# Prepare metadata
meta_mat2 <- meta_mat %>% 
  rownames_to_column(var = 'SampleID') %>% 
  merge(summary_dat[c(1:2)], by = 'SampleID') %>% 
  arrange(desc(总计))

# Reshape absolute abundance data
tmp <- genus2 %>% as.data.frame()
names(tmp) <- str_extract(names(tmp), "(?<=g__).+")
names(tmp)[1] <- 'Unknow_genus'

tmp <- tmp %>%
  rownames_to_column(var = "SampleID") %>% 
  pivot_longer(!'SampleID', names_to = "genus", values_to = "count") %>% 
  as.data.frame()

# Merge with metadata
data2 <- merge(tmp, meta_mat2, by = "SampleID") %>% 
  dplyr::arrange(desc(总计))
ureaplasma_counts <- data2 %>% 
  filter(genus == "Ureaplasma") %>% 
  select(SampleID, Ureaplasma_count = count)
data2 <- left_join(data2, ureaplasma_counts, by = "SampleID")

# Add blank placeholder
data2_blank <- data2[data2$SampleID == 'FA14', ]
data2_blank[, c(1, 3:12)] <- NA
data2_blank$SampleID <- '1'
data2_blank$count <- 0

# Factorize genus levels
data2$genus %<>% factor(
  levels = c("Lactobacillus", "Prevotella", "Ureaplasma", "Unkonw_genus"),
  labels = c("Lactobacillus", "Prevotella", "Ureaplasma", "Unkonw_genus")
)

# Plot for preterm group
data2_preterm <- subset(data2, group == "preterm")
data2_preterm <- rbind(data2_preterm, data2_blank)

p_preterm <- ggplot(data2_preterm, aes(reorder(SampleID, -count), y = -count, fill = genus)) +
  scale_fill_manual(values = colors) +
  geom_col(position = 'stack', width = 0.95) +
  labs(x = '', y = '', fill = "genus", title = "Preterm") +
  theme_classic(base_size = 10) +
  scale_x_discrete(labels = c(1:40, '', '')) +
  ggplot2::annotate(
    "text",
    x = rep(42, 5),
    y = c(0, -10000, -20000, -40000, -70000),
    label = c("0", "10000", "20000", "40000", "70000"),
    color = "grey",
    size = 3,
    angle = 0,
    hjust = "right"
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 5),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) + 
  ylim(-80000, 0)

# Convert to polar coordinates
p_preterm <- p_preterm + coord_polar() 
p_preterm

# Save preterm plot
p_width <- 3
ggsave("sfig 1E_species_composition_abundance_preterm.pdf", width = p_width, height = p_width * 1)

# Plot for term group
data2_term <- subset(data2, group == "term")
data2_term <- rbind(data2_term, data2_blank)

p_term <- ggplot(data2_term, aes(reorder(SampleID, -count), y = -count, fill = genus)) +
  scale_fill_manual(values = colors) +
  geom_col(position = 'stack', width = 0.9) +
  labs(x = '', y = '', fill = "genus", title = "Term birth") +
  theme_classic(base_size = 10) +
  scale_x_discrete(labels = c(1:39, '', '')) +
  ggplot2::annotate(
    "text",
    x = rep(42, 5),
    y = c(0, -10000, -20000, -40000, -70000),
    label = c("0", "10000", "20000", "40000", "70000"),
    color = "grey",
    size = 3,
    angle = 0,
    hjust = "right"
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 5),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) + 
  ylim(-80000, 0)

# Convert to polar coordinates
p_term <- p_term + coord_polar()
p_term

# Save term plot
p_width <- 3
ggsave("sfig 1E_species_composition_abundance_term.pdf", width = p_width, height = p_width * 1)

# Extract and save legend
legend <- cowplot::get_legend(
  ggplot(data2_term, aes(reorder(SampleID, -count), y = -count, fill = genus)) +
    scale_fill_manual(values = colors) +
    geom_col(position = 'stack', width = 0.9) +
    theme_classic(base_size = 10) + 
    guides(col = guide_legend(title = "Legend", nrow = 2)) +
    theme(
      legend.title = element_text(size = 10),
      legend.key.height = unit(0.5, "cm"),
      legend.text = element_text(size = 10),
      legend.direction = "horizontal"
    )
)
ggsave("sfig 1E_species_composition_abundance_legend.pdf", legend, width = 4, height = 3)

# 4.6 fig 1G & sfig 1F: Ureaplasma vs Placental Infection Status ------------------------------
aa <- read.csv("intermediate_file//preterm_study_vaginal_metagenomic_data_relative.csv") %>% 
  as.data.frame()
aa <- merge(aa, group, by.x = "X", by.y = "SampleID")
aa %<>% mutate(
  infected = ifelse(pla_urea == "rich", "infected", 
                    ifelse(all_microbiome == "sterile", "sterile", "other"))
)
aa$infected %<>% factor(levels = c("sterile", "infected", "other"))

# Custom boxplot function for Ureaplasma comparison
box_plot2 <- function(data, x, compare, title, sig_max) {
  ggplot2::ggplot(data, aes(x = {{x}}, y = Ureaplasma.parvum, fill = {{x}})) +
    geom_boxplot(outlier.colour = "white", outlier.size = 1) +
    geom_jitter(aes(fill = {{x}}), width = 0.2, shape = 21, size = 2) +
    scale_fill_manual(values = c('#FF6A6A', '#00D7F5')) +
    scale_color_manual(values = c("black", "black")) +
    theme_bw() +
    labs(
      x = "", 
      y = expression("% vaginal " ~ italic("U.parvum")), 
      title = title, 
      fill = ""
    ) +
    scale_y_log10(
      breaks = c(0.000001, 0.01, 1, 10, 100),
      limits = c(0.000001, 100),
      labels = c("1e-6", "0.01", "1", "10", "100")
    ) + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(colour = "black", size = 10),
      axis.text.y = element_text(size = 10, colour = "black"),
      plot.title = element_text(size = 10, hjust = 0.5, colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    geom_signif(
      comparisons = list(compare),
      map_signif_level = TRUE,
      test = wilcox.test,
      y_position = 1.5,
      vjust = 0.15,
      size = 0.5, color = "black"
    )
}

# 1. All samples (exclude "other")
box_plot2(data = filter(aa, infected != "other"), x = infected, 
          compare = c("sterile", "infected"), title = "")
p_width <- 2
ggsave("fig 1G_脲原体比较-胎盘有无微生物.pdf", width = p_width, height = p_width * 1.5)

# 2. Preterm only
box_plot2(data = subset(aa, group == "preterm" & infected != "other"),  
          x = infected, compare = c("sterile", "infected"), title = "preterm")
ggsave("sfig 1F-脲原体比较-胎盘有无微生物-早产.pdf", width = p_width, height = p_width * 1.5)

# 3. Term only
box_plot2(data = subset(aa, group == "term" & infected != "other"),  
          x = infected, compare = c("sterile", "infected"), title = "term")
ggsave("sfig 1F-脲原体比较-胎盘有无微生物-足月.pdf", width = p_width, height = p_width * 1.5)

# 4.7 sfig 1G: Vaginal-Placental Microbiota Correlation ------------------------------
# Import data
dat1.1 <- read.csv("intermediate_file//preterm_study_vaginal_metagenomic_data_relative.csv", row.names = 1)
dat1.2 <- read.csv('data/阴道微生物-胎盘微生物2.csv', row.names = 1)
dat1.2[dat1.2 < 100] <- 0  # Remove ASVs < 100 (considered as contamination)

# Calculate Spearman correlation
cor.data1 <- cbind(dat1.1, dat1.2)
cor.mat <- rcorr(as.matrix(cor.data1), type = "spearman")

# Replace NA p-values with 1 (no statistical significance)
cor.mat$P[is.na(cor.mat$P)] <- 1

# Plot correlation heatmap
col2 <- colorRampPalette(c(
  "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
  "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
))

pdf('sfig 1G-胎盘以及阴道微生物的相关性.pdf', width = 15, height = 15)
corrplot(
  cor.mat$r,
  p.mat = cor.mat$P,
  sig.level = c(0.05),
  insig = 'label_sig',
  pch.cex = 0.8,
  pch.col = 'grey20',
  type = "lower",
  cl.pos = "b",
  tl.pos = 'l',
  tl.col = 'black',
  diag = TRUE,
  tl.cex = 0.75,
  col = col2(1000)
)
dev.off()

# 4.8 sfig 1H: VFDB Virulence Factors Analysis ------------------------------

# Batch plot urease-related factors
vfdb2 <- read.csv("data/data_vfdb_urea.csv")
vfdb2$Up_infected %<>% factor(
  levels = c("sterile", "Up-infected"), 
  labels = c("sterile", "infected")
)
vfdb2 %<>% filter(ureasebetasubunitUreB < 100)

for (var in c("ureasebetasubunitUreB", "ureasealphasubunitUreA")) {
  ggplot2::ggplot(vfdb2, aes(x = Up_infected, y = get(var), fill = Up_infected)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(width = 0.2, shape = 21, size = 2) +
    scale_fill_manual(values = c('#FF6A6A', '#00D7F5')) +
    scale_color_manual(values = c("black", "black")) +
    theme_bw() +
    labs(
      x = "", 
      y = gsub("X\\.", " ", var), 
      fill = ""
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(colour = "black", size = 10),
      axis.text.y = element_text(size = 10, colour = "black"),
      plot.title = element_text(size = 10, hjust = 0.5, colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.09))) +
    geom_signif(
      comparisons = list(c("sterile", "infected")),
      map_signif_level = TRUE,
      test = wilcox.test,
      size = 0.5, color = "black"
    )
  
  # Save plot
  p_width <- 2
  ggsave(
    paste0("sfig 1H-", gsub("/", "", var), ".pdf"),
    width = p_width, height = p_width * 1.5
  )
}

# 5 analysis 孕期尿素 ------------------------------

dat <- openxlsx::read.xlsx("data/胎盘菌群研究临床信息摘录数据_整理.xlsx", sheet = 2)
dat[dat == "."] <- NA
dat <- mutate(dat, across(trimester1_BUN:trimester3_UA, as.numeric) )
dat <- mutate(dat, 
              trimester1_BUN_Cr = trimester1_BUN / trimester1_CREA, 
              trimester2_BUN_Cr = trimester2_BUN / trimester2_CREA, 
              trimester3_BUN_Cr = trimester3_BUN / trimester3_CREA )


dat %>% names

library(rstatix)
library(tidyr)

# 将数据转为长格式：每行一个变量值
dat_long <- dat %>%
  select(group, trimester1_BUN:trimester3_UA) %>%  # 选择 group 和所有指标
  pivot_longer(cols = -group, names_to = "variable", values_to = "value")

dat_long$value <- as.numeric(dat_long$value)
# 按变量分组进行 Wilcoxon 检验
wilcox_results <- dat_long %>%
  group_by(variable) %>%
  wilcox_test(value ~ group)

# 添加效应量（可选）
wilcox_results <- wilcox_results %>% add_effsize(value ~ group)

# 查看结果
wilcox_results


library(tableone)
# 定义要分析的变量（从第4列到第12列，根据你的列名顺序）
myVars <- names(dat)[4:15]   # 包含 trimester1_BUN 到 trimester3_UA

# 定义分类变量（此处没有分类变量，如果需要可在此添加）
catVars <- NULL

# 创建表格对象
# - strata: 分组变量（group）
# - data: 数据框
# - nonnormal: 指定哪些变量视为非正态，这里设为所有 myVars，将使用 Wilcoxon 检验
tab <- CreateTableOne(vars = myVars, 
                      strata = "group", 
                      data = dat)     # 强制所有变量使用非参数检验

# 打印表格，显示 p 值
full_tab <- print(tab, 
      showAllLevels = TRUE,   # 显示分组的所有水平（如果 group 是多分类）
      pvalue = TRUE,          # 显示 p 值
      exact = NULL,           # 精确检验（无需设置）
      nonnormal = myVars,
      formatOptions = list(big.mark = ","))

write.csv(full_tab, "tmp/table1_with_pvalues.csv", na = "")










