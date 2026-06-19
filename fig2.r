# Load packages and preprocessing ---------------------------------------------------------------

if(T) {
  # nolint
  require(Seurat) # nolint
  require(harmony)
  require(magrittr)
  require(ggplot2)
  require(reshape2)
  require(readxl)
  require(dplyr)
  #require(DoubletFinder)
  #require(openxlsx)
  require(grDevices)
  require(ggsignif)
  require(job)
  require(future)
  require(ggsci)
}

setwd('/public/download/lc2021/r/data_audit')

# Import data and data preprocessing------------------------------------------

rm(list = ls());gc()

pla_FM <- readRDS('data/pla_FM_data.rds') # Sourced from "2_normalization_dimension_reduction_clustering.r"

# Rename cells to publishable labels
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

# fig-2A -------------------------------------------------------------------

# miniaixs() is a custom function generating small arrow-style axes, script stored at /Birthcohort/Birthcohort_station/r/miniaixs.r
source("miniaixs.r")

Idents(pla_FM) <- "cell_type"
pdf("02_output/fig2/fig-2A-single_cell_qc_cell_composition.pdf", width = 4.5, height = 4)
miniaixs(DimPlot(pla_FM, reduction = "umap", label = T ))
dev.off()

# fig-2B -------------------------------------------------------------------

# Cell composition proportion
cellnum <- table(pla_FM$cell_type, pla_FM$regroup4)
cell.prop <- as.data.frame(prop.table(cellnum, margin = 2))
colnames(cell.prop) <- c("cell_type", "Group", "Proportion")
cell.prop <- merge(cell.prop,
                   unique(pla_FM@meta.data[c("regroup4", "urea")]),
                   by.x = "Group",
                   by.y = "regroup4")

cell.prop$Group <- factor(cell.prop$Group,
                            levels = c("preterm-infect", 
                                       "preterm-steri",
                                       'term-steri'))

cols <- pal_futurama()(13)
cell.prop$urea <- factor(cell.prop$urea, levels = c("rich", "poor"))

pdf("02_output/fig2/fig-2B-cell_composition_proportion_across_groups.pdf",
    width = 5,
    height = 5.5)

p_com_3 <- ggplot(cell.prop, aes(Group, Proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = cols) + # Custom fill color palette
  ggtitle("cell proportion") +
  theme_bw() +
  theme(axis.ticks.length = unit(0.5, "cm")) +
  guides(fill = guide_legend(title = NULL)) + RotatedAxis()

p_com_3

dev.off()

# fig_2C & E Differentially Expressed Genes---------------------------------------------------------------

# Build custom plotting function
cells <- c("B_cell", "Decidual", "Stromal", "fibroblasts", "Endothelial", 
           "EVT", "Mac_cluster1", "Mac_cluster2", "Mac_cluster3", "Myofibroblasts", 
           "Mac_pi", "T_NK", "Endometrial")

mycol<-c("#E54731","#4dbbd5","#009D83",
         "#3A5287","#F39B7F","#9CA7C3",
         "#91D1C2","#DB0000","#7E6148")

plot_11 <- function(DEGs_PMN0, DEGs_PMN0_sig) {
  ggplot(DEGs_PMN0, aes(x=difference, y=avg_log2FC)) + 
    geom_point(size=0.5, color="grey60") + 
    ggrepel::geom_text_repel(data=DEGs_PMN0_sig, aes(label=label), # Add gene text labels
                             color=ii,fontface="italic", size = 3)+
    geom_point(data=DEGs_PMN0[which(DEGs_PMN0$p_val_adj < 0.05 & DEGs_PMN0$avg_log2FC > 1),],
               aes(x=difference, y=avg_log2FC),
               size=0.5, color = ii) +  
    geom_point(data=DEGs_PMN0[which(DEGs_PMN0$p_val_adj < 0.05 & DEGs_PMN0$avg_log2FC < -1),],
               aes(x=difference, y=avg_log2FC),
               size=0.5, color = ii) +
    labs( x = "Delta Percent", 
          y = "Log-fold Change", 
          colour = "Number of cylinders",
          title = paste0(i,""))+
    ylim(-3, 5) +
    xlim(-0.6, 0.6) +  
    theme(axis.text.x = element_text(colour = 'black',size = 10),
          axis.text.y = element_text(colour = 'black',size = 10),
          axis.title = element_text(colour = 'black',size = 10),
          axis.line = element_line(color = 'black', size = 0.2),
          plot.title = element_text(face = "bold", size = 10),
          panel.border = element_rect(fill = NA, color = "black", size = 0.2),
          panel.background = element_blank())+
    geom_hline(yintercept = 0.2,lty= 2,lwd = 0.4)  +
    geom_hline(yintercept = -0.2,lty= 2,lwd = 0.4)
  
}


## fig_2C EVT  --------------------------------------------------------------------

var_gene <- c("BIRC3", "CXCL8", "NOTUM", "RND3", "ASCL2", "PTGDS", "TCF7L2")

# Load DEG result table
DEGs_PMN0 <- read.csv(paste0("data/DEGs/step3-DEGs_1_4_EVT.csv"), row.names = 1)

# Calculate expression percentage difference and filter significant genes (for labeling)
DEGs_PMN0$difference <- DEGs_PMN0$pct.1 - DEGs_PMN0$pct.2
DEGs_PMN0_sig <- DEGs_PMN0[which(DEGs_PMN0$p_val_adj < 0.05 & abs(DEGs_PMN0$avg_log2FC) > 1),]

# Create gene label column
DEGs_PMN0_sig <- DEGs_PMN0_sig %>%  mutate(label = ifelse(rownames(DEGs_PMN0_sig) %in% var_gene, rownames(DEGs_PMN0_sig), ""))

# Generate plot
j <- 1
i <- "EVT"
ii <- mycol[j]
plot1 <- list()

plot1[[j]] <- plot_11(DEGs_PMN0, DEGs_PMN0_sig)

## Export PDF figure
pdf(
  "02_output/fig2/fig_2C_volcano-EVT.pdf",
  width = 2,
  height = 5
)
plot1[[1]]
dev.off()

## fig_2E-Mac1  --------------------------------------------------------------------
plot1 <- list()

i <- "Mac1"
var_gene <- c("CXCL8", "SOD2", "NAMPT", "CD52", "CCL20", "IL1B", "NOTUM", 
              "CCL4L2", "RNASE1", "CCL3L1", "CXCL1", "CXCL5", "IL1RN", "PAPPA2", 
              "CCL4", "PTGS2", "TNFAIP6")
j <- 1
ii <- mycol[j]

# Load DEG result table
DEGs_PMN0 <- read.csv(paste0("data/DEGs/step3-DEGs_1_4_", i, ".csv"), row.names = 1)

# Calculate expression percentage difference and filter significant genes (for labeling)
DEGs_PMN0$difference <- DEGs_PMN0$pct.1 - DEGs_PMN0$pct.2
DEGs_PMN0_sig <- DEGs_PMN0[which(DEGs_PMN0$p_val_adj < 0.05 & abs(DEGs_PMN0$avg_log2FC) > 1),]

# Create gene label column
DEGs_PMN0_sig <- DEGs_PMN0_sig %>%  mutate(label = ifelse(rownames(DEGs_PMN0_sig) %in% var_gene, rownames(DEGs_PMN0_sig), ""))
# DEGs_PMN0_sig$label <- rownames(DEGs_PMN0_sig)

# Generate plot
plot1[[j]] <- plot_11(DEGs_PMN0, DEGs_PMN0_sig)

## fig_2E-Mac2 --------------------------------------------------------------------

i <- "Mac2"
var_gene <- c("HLA-DRB6", "CD74", "HLA-DPB1", "HLA-DRA", "APOE", "HLA-DRB1", 
              "HLA-DPA1", "G0S2", "CXCL1", "TNFAIP6", "MT1E", "CHI3L1", "NAMPT", 
              "SLC39A8", "IL1B", "SOD2", "MT1H", "MT2A", "FCGBP", "CCL20", 
              "MMP9", "MT1G", "CXCL5")
j <- 2
ii <- mycol[j]

# Load DEG result table
DEGs_PMN0 <- read.csv(paste0("data/DEGs/step3-DEGs_1_4_", i, ".csv"), row.names = 1)

# Calculate expression percentage difference and filter significant genes (for labeling)
DEGs_PMN0$difference <- DEGs_PMN0$pct.1 - DEGs_PMN0$pct.2
DEGs_PMN0_sig <- DEGs_PMN0[which(DEGs_PMN0$p_val_adj < 0.05 & abs(DEGs_PMN0$avg_log2FC) > 1),]

# Create gene label column
DEGs_PMN0_sig <- mutate(DEGs_PMN0_sig, label = ifelse(rownames(DEGs_PMN0_sig) %in% var_gene, rownames(DEGs_PMN0_sig), ""))
# DEGs_PMN0_sig$label <- rownames(DEGs_PMN0_sig)

# Generate plot
plot1[[j]] <- plot_11(DEGs_PMN0, DEGs_PMN0_sig)

## Export combined PDF figure
pdf(
  "02_output/fig2/fig_2E_volcano-Mac1_Mac2.pdf",
  width = 4,
  height = 5
)
plot1[[1]] | plot1[[2]]
dev.off()

## sfig_2E-Fibroblasts  --------------------------------------------------------------------

plot1 <- list()

i <- "Fibroblasts"

var_gene <- c("TNFSF10", "KRT7", "FSTL3", "SERPINE2", "PAPPA2", "ALCAM", 
              "THBS1", #"BIRC3", 
              "AGR3", "SERPINB2")
j <- 2
ii <- mycol[j]

# Load DEG result table
DEGs_PMN0 <- read.csv(paste0("data/DEGs/step3-DEGs_1_4_", i, ".csv"), row.names = 1)

# Calculate expression percentage difference and filter significant genes (for labeling)
DEGs_PMN0$difference <- DEGs_PMN0$pct.1 - DEGs_PMN0$pct.2
DEGs_PMN0_sig <- DEGs_PMN0[which(DEGs_PMN0$p_val_adj < 0.05 & abs(DEGs_PMN0$avg_log2FC) > 1),]

# Create gene label column
DEGs_PMN0_sig <- DEGs_PMN0_sig %>%  mutate(label = ifelse(rownames(DEGs_PMN0_sig) %in% var_gene, rownames(DEGs_PMN0_sig), ""))
# DEGs_PMN0_sig$label <- rownames(DEGs_PMN0_sig)

# Generate plot
plot1[[j]] <- plot_11(DEGs_PMN0, DEGs_PMN0_sig)


## sfig_2E-Stromal  --------------------------------------------------------------------

i <- "Stromal"
var_gene <- c("VCAN", "PTX3", "COMP", "IGF2", "SAA1", "IL1R2", "COL1A2", 
              "COL1A1")
j <- 3
ii <- mycol[j]

# Load DEG result table
DEGs_PMN0 <- read.csv(paste0("data/DEGs/step3-DEGs_1_4_", i, ".csv"), row.names = 1)

# Calculate expression percentage difference and filter significant genes (for labeling)
DEGs_PMN0$difference <- DEGs_PMN0$pct.1 - DEGs_PMN0$pct.2
DEGs_PMN0_sig <- DEGs_PMN0[which(DEGs_PMN0$p_val_adj < 0.05 & abs(DEGs_PMN0$avg_log2FC) > 1),]

# Create gene label column
DEGs_PMN0_sig <- mutate(DEGs_PMN0_sig, label = ifelse(rownames(DEGs_PMN0_sig) %in% var_gene, rownames(DEGs_PMN0_sig), ""))
# DEGs_PMN0_sig$label <- rownames(DEGs_PMN0_sig)

# Generate plot
plot1[[j]] <- plot_11(DEGs_PMN0, DEGs_PMN0_sig)

plot1[[j]]


## sfig_2E Merge subfigures 

pdf(
  "02_output/fig2/sfig_2E_volcano-fibroblast_stromal.pdf",
  width = 5,
  height = 5
)
plot1[[2]]|plot1[[3]]
dev.off()

# sfig-2A -------------------------------------------------------------------

# (nFeature_n, RNA count, mitochondrial percentage) 
meta <- readRDS("data/meta.rds") 

meta$cell_type <-
  factor(
    meta$celltype,
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

# Custom ggplot theme
my_theme <- theme_bw() +
  theme(
    axis.text = element_text(family = "Times", color = "black", size = 10),
    axis.title.y = element_text(family = "Times", size = 10),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# QC violin plot function
plot_qc <- function(metric, ylab) {
  ggplot(meta, aes(x = regroup, y = .data[[metric]], fill = regroup)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
    my_theme +
    labs(x = "", y = ylab)
}

# Combine QC plots
pdf('02_output/fig2/sfig-2A-qc-mitochondrial_rna_hb2_final.pdf', width = 9, height = 3)
plot_qc("nFeature_RNA", "nFeature_RNA") |
  plot_qc("nCount_RNA",   "nCount_RNA") |
  plot_qc("percent.mt",   "percent.mt")
dev.off()


# sfig-2B -------------------------------------------------------------------

pdf('02_output/fig2/sfig-2B-qc-mitochondrial_rna_hb2_by_celltype_final.pdf', width = 5, height = 3)
ggplot(aes(y = nFeature_RNA, x = cell_type), data = meta) + 
  geom_violin(aes(fill = cell_type), trim = FALSE, scale = "width") + # Generate violin plot; set white background (no outline)
  # "trim": If TRUE(default), trim violin tails to data range; FALSE retains full tails
  geom_boxplot(aes(fill = cell_type), width=0.2, position = position_dodge(0.9))+ # Overlay boxplot
  # scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ # Custom fill colors
  theme_bw()+ # White background
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=10), # Rotate x-axis tick labels 45°, right-align, Times font
        axis.text.y=element_text(family="Times",size=10,face="plain"), # Y-axis tick label style
        axis.title.y=element_text(family="Times",size = 10,face="plain"), # Y-axis title style
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), # Remove default grey border, thicken axis lines
        # legend.text=element_text(face="italic", family="Times", colour="black",  # Legend text style
        #                          size=10),
        # legend.title=element_text(face="italic", family="Times", colour="black", # Legend title style
        #                           size=10),
        legend.position = "none",
        panel.grid.major = element_blank(),   # Hide major grid
        panel.grid.minor = element_blank())+  # Hide minor grid
  ylab("nFeature_RNA")+xlab("")## Set X/Y axis titles
dev.off()


# sfig-2C -------------------------------------------------------------------

# Display canonical marker genes for all cell populations
marker1 <- c("CD79A", "IGHG1", "IGHG3", # B_cell
             "DKK1", "IGFBP1", "PRL", # Decidual
             "MME", "COL1A1", "COL3A1", "COL1A2", # Stromal
             "PERP", "KRT6A", "KRT17", # fibroblasts
             "LYVE1", "MMRN1", "CD34", "PECAM1",  # Endothelial
             "PAPPA2", "HLA-G", "NOTUM", "KRT7", # EVT
             "CD68", "CD163", "MRC1", "C1QB", # MAC cluster1
             "IL1B", "CXCL8", "IL1RN" , # MAC cluster2
             "CD86",# MAC cluster3
             "ACTA2", "TAGLN", "HSPB1", # Myofibroblasts
             "TOP2A", "MKI67", "CENPF", "NUSAP1", # Mac_pi
             "CCL5", "GZMB", "GZMA", # T_NK
             "SLPI", "PAEP", "TSPAN1", "DEFB1" # Endometrial
             
)

pdf("02_output/fig2/sfig-2C-celltype_marker_dotplot.pdf", width = 15, height = 6)
DotPlot(pla_FM, features = marker1, group.by = 'cell_type', cols = c("lightgrey", "red")) + 
  RotatedAxis() # + # Rotate x-axis text
ggplot2:::coord_flip() # Flip plot axes
dev.off()


# sfig-2D -------------------------------------------------------------------

pla_FM@meta.data$group2 <- factor(pla_FM@meta.data$split, levels = c("control", "case"), labels = c("term", "preterm"))

Idents(pla_FM) <- "group2"
pdf("02_output/fig2/sfig-2D-single_cell_qc_cell_composition_preterm.pdf", width = 4.5, height = 4)
miniaixs(DimPlot(pla_FM, reduction = "umap", label = F))
dev.off()

Idents(pla_FM) <- "urea"
pdf("02_output/fig2/sfig-2D-single_cell_qc_cell_composition_infection_status.pdf", width = 4.5, height = 4)
miniaixs(DimPlot(pla_FM, reduction = "umap", label = F))
dev.off()


# sfig-2F -------------------------------------------------------------------

plot <- VlnPlot(
  subset(pla_FM, subset = cell_type == "EVT"),
  features =  c("BIRC3", "NFKBIA", "CXCL8", "RND3", "NOTUM", "ASCL2"),
  ncol = 3,
  group.by = 'regroup4',
  pt.size = 0) +labs(x = "", y = "") + theme(legend.position = "none"  )

wid = 8
pdf("02_output/fig2/sfig-2F-EVT_gene_violin.pdf", width = wid, height = 5/7 * wid)
print(plot)
dev.off()

# sfig-2G -------------------------------------------------------------------

wid = 8
pdf("02_output/fig2/sfig-2G-fibroblasts_gene_violin.pdf", width = wid, height = 5/7 * wid)
VlnPlot(
  subset(pla_FM, subset = cell_type == "Fibroblasts"),
  features =  c("ALCAM", "THBS1", "AGR3", 
                "SERPINB2", "TNFSF10","SERPINE2"),
  ncol = 3,
  group.by = 'regroup4',
  pt.size = 0
) +  labs(x = "") + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()

# sfig-2I -------------------------------------------------------------------

pdf("02_output/fig2/sfig-2I/hypoxia_genes-EVT.pdf", width = 3, height = 3)
VlnPlot(
  subset(pla_FM, cell_type == "EVT"),
  features =  c("HIF1A"), #EPAS1 hypoxia-inducible factor 2
  group.by = 'regroup4',
  pt.size = 0
) + NoLegend() + labs(x="", y="", title="EVT")
dev.off()

pdf("02_output/fig2/sfig-2I/hypoxia_genes-Decidual.pdf", width = 3, height = 3)
VlnPlot(
  subset(pla_FM, cell_type == "Decidual"),
  features =  c("HIF1A"), #EPAS1 hypoxia-inducible factor 2
  group.by = 'regroup4',
  pt.size = 0
) + NoLegend() + labs(x="", y="", title="Decidual")
dev.off()

pdf("02_output/fig2/sfig-2I/hypoxia_genes-Stromal.pdf", width = 3, height = 3)
VlnPlot(
  subset(pla_FM, cell_type == "Stromal"),
  features =  c("HIF1A"), #EPAS1 hypoxia-inducible factor 2
  group.by = 'regroup4',
  pt.size = 0
) + NoLegend() + labs(x="", y="", title="Stromal")
dev.off()

pdf("02_output/fig2/sfig-2I/hypoxia_genes-Fibroblasts.pdf", width = 3, height = 3)
VlnPlot(
  subset(pla_FM, cell_type == "Fibroblasts"),
  features =  c("HIF1A"), #EPAS1 hypoxia-inducible factor 2
  group.by = 'regroup4',
  pt.size = 0
) + NoLegend() + labs(x="", y="", title="Fibroblasts")
dev.off()

pdf("02_output/fig2/sfig-2I/hypoxia_genes-Endometrial.pdf", width = 3, height = 3)
VlnPlot(
  subset(pla_FM, cell_type == "Endometrial"),
  features =  c("HIF1A"), #EPAS1 hypoxia-inducible factor 2
  group.by = 'regroup4',
  pt.size = 0
) + NoLegend() + labs(x="", y="", title="Endometrial")
dev.off()

pdf("02_output/fig2/sfig-2I/hypoxia_genes-Myofibroblasts.pdf", width = 3, height = 3)
VlnPlot(
  subset(pla_FM, cell_type == "Myofibroblasts"),
  features =  c("HIF1A"), #EPAS1 hypoxia-inducible factor 2
  group.by = 'regroup4',
  pt.size = 0
) + NoLegend() + labs(x="", y="", title="Myofibroblasts")
dev.off()

# sfig-2J -------------------------------------------------------------------

pdf("02_output/fig2/sfig-2J-Mac2_gene_violin.pdf", width = 6, height = 2.5)
VlnPlot(
  subset(pla_FM, subset = cell_type == "Mac2"),
  features =  c("CXCL5", "MMP9", "FCGBP", "CCL20"),
  ncol = 4,
  group.by = 'regroup4',
  pt.size = 0
)  + theme(legend.position = "none", axis.title.x = element_blank())
dev.off()


# Local GSEA analysis ---------------------------------------------------------------

#### GSEA analysis must be run locally !!!!

if(T){
  require('rvcheck')
  require('clusterProfiler')
  require('org.Hs.eg.db')
  require('dplyr')
  require('ggplot2')
  require('enrichplot')
  require('forcats')
  
}

setwd('D:/1_graduate_study/1_ART/13_placenta_project_zhang/placenta_research/microbiome/12_manuscript_writing/code_collection/CHM_short/202601_first_submission/fig2')

# Custom GSEA visualization function
gsea_beatiful <- function(egseKEGG, pathway, var_gene, genes) {
  a <- 
    enrichplot::gseaplot2(
      egseKEGG,
      geneSetID = which(egseKEGG$Description == pathway),
      title = pathway,
      pvalue_table = T, subplots = 1
    )
  
  var_gene <- var_gene
  
  pvalue <- egseKEGG[which(egseKEGG$Description == pathway), c("Description", "pvalue", "p.adjust")][2] 
  padjust <- ifelse(egseKEGG[which(egseKEGG$Description == pathway), c("Description", "pvalue", "p.adjust")][3] < 0.05, "< 0.05", "not significant")
  
  # Custom aesthetic adjustment
  a$data %>%
    # mutate(gene = AnnotationDbi::mapIds(org.Hs.eg.db, a$data$gene, "ENTREZID", "SYMBOL")) %>% #filter(position==1)
    mutate(label=ifelse(gene %in% var_gene, gene, "")) %>% 
    {
      ggplot(.)+
        geom_line(data = .,aes(x=x,y=runningScore))+
        geom_point(data = subset(.,position==1),aes(x=x,y=runningScore),color=rgb(253,174,97,max=255))+
        ggrepel::geom_label_repel	(data = subset(.,position==1),aes(x=x,y=runningScore,label=label),color="#F39B7F",max.overlaps = 30)+
        geom_rect(data = subset(.,position==1),
                  aes(xmin=x-0.001,
                      xmax=x+0.001,
                      ymin=min(runningScore)-0.25,
                      ymax=min(runningScore)-0.125,
                      color=x),
                  show.legend = F)+
        geom_rect(data = .,
                  aes(xmin=x-0.001,
                      xmax=x+0.001,
                      ymin=min(runningScore)-0.375,
                      ymax=min(runningScore)-0.25,
                      color=x
                  ),
                  show.legend = F)+
        geom_hline(data = .,aes(yintercept = min(runningScore)-0.12),size=0.5)+
        geom_hline(yintercept = 0,lty=2,color="grey",size=1)+
        geom_vline(xintercept = median(.$x),lty=2,color="grey",size=0.5)+
        geom_segment(data = .,aes(x=max(x)*(1/4),
                                  xend=max(x)*(1/4)-max(.$x)*(1/16),
                                  y=max(runningScore)+0.1,
                                  yend=max(runningScore)+0.1), 
                     size = 1.2,color=rgb(116,173,209,max=255),
                     arrow = arrow(length = unit(0.3, "cm")))+
        geom_segment(data = .,aes(x=max(x)*(3/4),
                                  xend=max(x)*(3/4)+max(.$x)*(1/16),
                                  y=max(runningScore)+0.1,
                                  yend=max(runningScore)+0.1),
                     size = 1.2,color=rgb(253,174,97,max=255),
                     arrow = arrow(length = unit(0.3, "cm")))+
        annotate("text",x=max(.$x)*(3/8),y=max(.$runningScore)+0.08,label="preterm-inf.",hjust = 0.5,vjust = -0,size=3.8,color=rgb(116,173,209,max=255))+
        annotate("text",x=max(.$x)*(5/8),y=max(.$runningScore)+0.08,label="term-st.",hjust = 0.5,vjust = -0,size=3.8,color=rgb(253,174,97,max=255))+
        # annotate("text",x=max(.$x)*(3/4),y=0+0.05,label="runningScore = 0")+
        # annotate("text",x=max(.$x)*(1/4),y=min(.$runningScore)+0.05,label="P.value=0.05\nP.adj=0.05\nES=0.5")+
        annotate("text",x=max(.$x)*(1/4) + 60,y=min(.$runningScore)+0.2, size = 2.8,label = paste0("P.adj=", padjust))+
        scale_color_gradient2(low = rgb(116,173,209,max=255),
                              midpoint = median(.$x),
                              mid = "white",
                              high = rgb(253,174,97,max=255))+
        scale_x_continuous(limits = c(1,length(.$x)),expand = c(0,0))+
        scale_y_continuous(expand = c(0,0),limits = c(min(.$runningScore)-0.375,max(.$runningScore)+0.15),breaks = seq(round(min(.$runningScore),2),round(max(.$runningScore),2),0.2))+
        labs(x="",title = pathway) + ggprism::theme_prism() + theme(plot.margin = unit(rep(0.5,4),"cm"),
                                                                    #axis.text.x = element_blank(),
                                                                    #axis.ticks.x = element_blank()
                                                                    panel.border = element_rect(size = 1,fill=NA),
                                                                    axis.line = element_blank()
        )
    }
  
}

# fig_2D Local execution ---------------------------------------------------------------

egseGO <- readRDS("01_data/egseGO_1_4_EVT.rds")

gsea_beatiful(
  egseKEGG = egseGO,
  pathway = "apoptotic process",
  var_gene = c("")
)

ggsave("02_output/fig_2D-1.pdf", width = 5, height = 5)

gsea_beatiful(
  egseKEGG = egseGO,
  pathway = "innate immune response-activating signaling pathway",
  var_gene = c("")
)

ggsave("02_output/fig_2D-2.pdf", width = 5, height = 5)


# fig_2F---------------------------------------------------------------
egseKEGG <- readRDS("01_data/egseKEGG_1_4_Mac1.rds")

gsea_beatiful(
  egseKEGG = egseKEGG,
  pathway = "Cytokine-cytokine receptor interaction",
  var_gene = c("")
)
ggsave("02_output/fig_2F-1.pdf", width = 5, height = 5)

gsea_beatiful(
  egseKEGG = egseKEGG,
  pathway = "NF-kappa B signaling pathway",
  var_gene = c("")
)
ggsave("02_output/fig_2F-2.pdf", width = 5, height = 5)


# sfig_2H---------------------------------------------------------------

egseGO <- readRDS("01_data/egseGO_1_4_Fibroblasts.rds")

gsea_beatiful(
  egseKEGG = egseGO,
  pathway = "innate immune response",
  var_gene = c("")
)
ggsave("02_output/sfig_2H-1.pdf", width = 5, height = 5)

gsea_beatiful(
  egseKEGG = egseGO,
  pathway = "pattern recognition receptor signaling pathway",
  var_gene = c("")
)
ggsave("02_output/sfig_2H-2.pdf", width = 5, height = 5)

# sfig_2K---------------------------------------------------------------

egseKEGG <- readRDS( "01_data/egseKEGG_1_4_Mac2.rds")

gsea_beatiful(
  egseKEGG = egseKEGG,
  pathway = "Th1 and Th2 cell differentiation",
  var_gene = c("")
) 
ggsave("02_output/sfig_2K-1.pdf", width = 5, height = 5)

gsea_beatiful(
  egseKEGG = egseKEGG,
  pathway = "IL-17 signaling pathway",
  var_gene = c("")
)
ggsave("02_output/sfig_2K-2.pdf", width = 5, height = 5)

gsea_beatiful(
  egseKEGG = egseKEGG,
  pathway = "TNF signaling pathway",
  var_gene = c("")
)
ggsave("02_output/sfig_2K-3.pdf", width = 5, height = 5)

# MIF expression difference across urea levels---------------------------------------------------------------

# Differential expression of MIF gene---------------------------------------------------------------

dt_emp <- read.csv("data/meta_data/data_vfdb_urea_with_sampleID.csv")

pla_FM2 <- subset(pla_FM, subset = cell_type == "Stromal")
pla_FM2@meta.data$family = gsub("FM", "", pla_FM2@meta.data$orig.ident)

pla_FM2@meta.data <- merge(pla_FM2@meta.data, dt_emp, by.x = "family", by.y = "family", all.x = T)
pla_FM2@meta.data$UreA <- as.factor(pla_FM2@meta.data$UreA)

# dir.create("02_output/fig2/tmp")
plot <- VlnPlot(
  pla_FM2,
  features =  c("MIF"),
  group.by = 'UreA',
  pt.size = 0) +labs(x = "", y = "MIF", title = "Stromal") + theme(legend.position = "none")

wid = 8
pdf("02_output/fig2/tmp/sfig-2F-Stromal_UreA.pdf", width = wid, height = 5/7 * wid)
print(plot)
dev.off()

pla_FM2 <- subset(pla_FM, subset = cell_type == "Stromal")
pla_FM2@meta.data$family = gsub("FM", "", pla_FM2@meta.data$orig.ident)

pla_FM2@meta.data <- merge(pla_FM2@meta.data, dt_emp, by.x = "family", by.y = "family", all.x = T)
pla_FM2@meta.data$UreB <- as.factor(pla_FM2@meta.data$UreB)

# dir.create("02_output/fig2/tmp")
plot <- VlnPlot(
  pla_FM2,
  features =  c("MIF"),
  group.by = 'UreB',
  pt.size = 0) +labs(x = "", y = "MIF", title = "Stromal") + theme(legend.position = "none")

wid = 8
pdf("02_output/fig2/tmp/sfig-2F-Stromal_UreB.pdf", width = wid, height = 5/7 * wid)
print(plot)
dev.off()


# Differential expression of MIF signaling pathway---------------------------------------------------------------

mif_df <- subsetCommunication(cellchat_1, signaling = "MIF")

# Retrieve names of the third dimension (usually interaction names or pathway labels)
interaction_names <- dimnames(cellchat@net$caserich$prob)[[3]]

# Find indexes of interactions belonging to MIF pathway (adjust based on actual label format)
# Common format example: "MIF_ligand_receptor" or pathway stored in separate slot
mif_idx <- grep("^MIF", interaction_names, ignore.case = TRUE)

if (length(mif_idx) > 0) {
  # Extract probability array for MIF pathway (cell group × cell group × MIF interactions)
  mif_prob_array <- cellchat@net$caserich$prob[, , mif_idx, drop = FALSE]
  
  # Convert array to long-format dataframe
  library(reshape2)
  mif_df <- melt(mif_prob_array, varnames = c("source", "target", "interaction"), 
                 value.name = "prob")
  mif_df <- mif_df[!is.na(mif_df$prob) & mif_df$prob > 0, ]
  
  # Add pathway identifier column
  mif_df$pathway <- "MIF"
} else {
  warning("No interactions associated with MIF pathway were found")
}
