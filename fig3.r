# Pseudotime trajectory analysis -------------------------------------------------------------------

library(Seurat)
library(patchwork)
library(monocle)
library(tidyverse)
library(ggplot2)

library(CellChat)
library(igraph)
library(ComplexHeatmap)

# Load dataset-----------------------------------------------

rm(list = ls());gc()

pla_FM <- readRDS('data/pla_FM_data.rds') 

# Rename cell types to publishable labels
pla_FM$cell_type <-
  factor(
    pla_FM$celltype,
    levels = c(
      "B_cell",
      "Decidual",
      "Stromal",
      "fibroblasts",
      "Endothelial",
      "EVT",
      "Mac_cluster2", "Mac_cluster1",
      "Mac_cluster3",
      "Myofibroblasts",
      "Mac_pi",
      "T_NK",
      "Endometrial"
    ),
    labels = c(
      "B_cell",
      "Decidual",
      "Stromal",
      "Fibroblasts",
      "Endothelial",
      "EVT",
      "Mac1", "Mac2",
      "Mac_trans",
      "Myofibroblasts",
      "Mac_pi",
      "T_NK",
      "Endometrial"
    )
  )


pla_FM$regroup4 <-
  factor(
    pla_FM$regroup,
    levels = c("caserich", "casepoor", "controlpoor"),
    labels = c("preterm-infect", "preterm-steri", "term-steri") )


# Construct Monocle object----------------------------------

# Inspect Monocle object
cd <- readRDS('data/cell trajectories_124.rds')

# Visualization of pseudotime trajectory --------------------------------------------------------------

# fig3 A  --------------------------------------------------------------

# Step1: Color mapping by infection status
cd$regroup <- factor(cd$regroup, levels=c("caserich","casepoor","controlpoor"))

library(cowplot)

pdf("02_output/fig3/fig-3A-celltrajector_urea1.pdf",  width = 7,  height = 3.5)
(plot_cell_trajectory(cd, color_by = "celltype", cell_size = 0.5) + 
    ggtitle('') +
    facet_wrap(~ regroup4, nrow = 1) + 
    theme(legend.position = 'none')) 
dev.off()

p1 <- plot_cell_trajectory(
  cd,
  show_cell_names = FALSE,
  color_by = "celltype",
  cell_size = 0.5
) + theme(
  legend.title = element_text(size = 7),
  legend.key.size = unit(0.5, "cm"),
  legend.text = element_text(size = 7)
)

p2 <- plot_cell_trajectory(
  cd,
  show_cell_names = FALSE,
  color_by = "Pseudotime",
  cell_size = 0.5
) + theme(
  legend.title = element_text(size = 5),
  legend.key.size = unit(0.5, "cm"),
  legend.text = element_text(size = 5)
)

pdf("02_output/fig3/fig-3A-celltrajector_urea2.pdf", width = 7, height = 3.5)

plot_grid(
  p1,
  p2,
  nrow = 1
)

dev.off()

# fig3 B  --------------------------------------------------------------

pdf("02_output/fig3/fig-3B-celltrajector_state1.pdf", width = 3, height = 3)

plot_cell_trajectory(
  cd,
  show_cell_names = F,
  color_by = "State",
  cell_size = 0.2
)
dev.off()


# fig3 B  --------------------------------------------------------------

# Plot stacked percentage bar chart of cell states across groups (standard histogram)
cd_meta <- cd@phenoData@data

data_bar2 <- data.frame(table(cd_meta$State, cd_meta$regroup))
data_bar2 <- data_bar2 %>%
  group_by(Var2) %>%
  mutate(percentage = Freq / sum(Freq) * 100)

data_bar2$Var3 <- ifelse(data_bar2$Var2 =="controlpoor", "term-steri", 
                         ifelse(data_bar2$Var2== "casepoor", "preterm-steri", "preterm-infe"))
data_bar2$Var3 <- factor(data_bar2$Var3, levels = c("term-steri", "preterm-steri", "preterm-infe"))
df_final2 <- dplyr::rename(data_bar2, state = Var1, group = Var3)
df_final2$state <- factor(df_final2$state, levels = c("1","2","3","4","5"))

# Standard stacked bar chart 
plot1 <- ggplot(df_final2, aes(x = group, y = percentage, fill = state)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_npg() +
  theme_classic() +
  labs(x = "", y = "Ratio") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "right"
  )

pdf("02_output/fig3/fig-3B-Sankey_plot.pdf", width = 4, height = 3.5)
plot1
dev.off()

# fig 3C --------------------------------------------------------------

# Pseudotime branch analysis 

BEAM_res_1 <- readRDS('data/BEAM_res_branch point 1.rds')

# Visualize all genes at branch point 1/2, runtime ~15 min
tmp1 <-
  plot_genes_branched_heatmap(
    cd[row.names(subset(BEAM_res_1, qval < 1e-4)), ],
    branch_point = 1,
    num_clusters = 5, # Number of gene clusters to divide genes into
    cores = 16,
    # cluster_rows = F, 
    #hmcols = NULL, # Default color palette
    use_gene_short_name = T,
    show_rownames = F,
    return_heatmap = T # Whether to return heatmap metadata
  )

# Export heatmap figure
pdf("02_output/fig3/fig_3C-analyzing_branched_1_cluster_4.pdf", width = 4,height = 6)
tmp1$ph_res
dev.off()


# sfig 3A --------------------------------------------------------------

# Visualize selected marker genes
pdf("02_output/fig3/sfig-3A-genes3-branch-point1.pdf", width = 6, height = 3.2)
plot_genes_branched_pseudotime(cd[c("LGALS3","LAPTM5","TYROBP", "IFI16"),], # These four genes belong to immune suppression GO terms
                               color_by = "cell_type", branch_point = 1, cell_size=0.5, ncol = 2)
dev.off()

# CellChat intercellular communication analysis --------------------------------------------------------------


rm(list = ls());gc()

# Load preterm-infect dataset file
cellchat_1 <- readRDS('data/cellchat_1_Secreted.rds')
cellchat_4 <- readRDS('data/cellchat_4_Secreted.rds')

# Combine multiple CellChat objects
object.list <-
  list(caserich = cellchat_1, 
       controlpoor = cellchat_4)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

groupSize1 <- as.numeric(table(cellchat_1@idents))
groupSize2 <- as.numeric(table(cellchat_4@idents))

# Functional similarity calculation
# 1 Calculate functional similarity, structural similarity not used here
set.seed("123")
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

# fig3 D-------------------------------------------------

pdf("02_output/fig3/fig-3D-overall_information_flow_per_signaling_pathway.pdf",width = 7,height = 4)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)  
dev.off()

# fig3 E-------------------------------------------------

pdf("02_output/fig3/fig-3E-calculate_and_visualize_pathway_distance.pdf",width = 2,height = 3)
rankSimilarity(cellchat, type = "functional", color.use = "#EB746B" )
dev.off() 



# sfig 3B --------------------------------------------------------------

# 3 Compare major sender and receiver cell populations in 2D space, identify cell groups with significantly altered sending/receiving signals between datasets

num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})

weight.MinMax <- c(min(num.link), max(num.link)) # Unified dot size range across datasets

gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <-
    netAnalysis_signalingRole_scatter(object.list[[i]],
                                      title = names(object.list)[i],
                                      weight.MinMax = weight.MinMax)}

pdf("02_output/fig3/sfig-3B-compare_major_senders_and_receivers_2D.pdf",width = 4,height = 5.5)
gg
dev.off()

# sfig 3C --------------------------------------------------------------

# 2.2 Incoming signaling pathways
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[1], width = 5, height = 8.5, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2], width = 5, height = 8.5, color.heatmap = "GnBu")

pdf("02_output/fig3/sfig-3C-two_group_comparison_incoming_signals.pdf", width = 8.5, height = 6)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
dev.off()

# fig 3H --------------------------------------------------------------

pdf('02_output/fig3/fig-3H-MIF_ligand_receptor_contribution_1.pdf', width = 4,height = 3)
netAnalysis_contribution(cellchat_1, signaling = "MIF")
dev.off()

# fig 3I --------------------------------------------------------------

# 5 Visualize signal strength of senders and receivers in MIF pathway
pdf("02_output/fig3/fig-3I-macrophage_sender_MIF_pathway.pdf",width = 6,height = 3)
netAnalysis_signalingRole_network(cellchat_1, signaling = "MIF", width = 8, height = 2.5, font.size = 10)
dev.off()

# fig 3G --------------------------------------------------------------

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("caserich", "controlpoor")) # Set factor order

pdf("02_output/fig3/fig-3G-part6.1-gene_expression_MIF.pdf", width = 8, height = 6)
plotGeneExpression(cellchat, signaling = "MIF", 
                   split.by = "datasets", colors.ggplot = T)
dev.off()

# fig 3F --------------------------------------------------------------

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = "MIF")  # Unified edge weight scale across datasets

pdf(paste0("02_output/fig3/fig-3F-chord_diagram_MIF.pdf"), width = 7,height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_aggregate( object.list[[1]], signaling = "MIF", 
                     layout = "circle", edge.weight.max = weight.max[1], 
                     edge.width.max = 10, 
                     vertex.weight = groupSize1, # Circle size corresponds to cell population proportion
                     signaling.name = paste("MIF", names(object.list)[1]) )

netVisual_aggregate( object.list[[2]], signaling = "MIF", 
                     layout = "circle", edge.weight.max = weight.max[1], 
                     edge.width.max = 10, 
                     vertex.weight = groupSize2,
                     signaling.name = paste("MIF", names(object.list)[2]) )
dev.off()
