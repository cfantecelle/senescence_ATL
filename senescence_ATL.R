################################################################################
################################################################################
##                                                                            ##
##  Author of this script: Carlos Henrique Fantecelle                         ##
##                                                                            ##
##  Analysis of the Gene Expression data (RNA-Seq) from several RNASeq        ##
##  studies in leishmania infections.                                         ##
##  This includes multiple leishmaniasis forms. Authors of studies used are   ##
##  listed below.                                                             ##
##                                                                            ##
##  Amorim et al, 2019. (Cutaneous form)                                      ##
##  Christensen et al, 2016 and 2019. (Cutaneous and Diffuse forms)           ##
##  Maretti-Mira et al, 2012. (Cutaneous and Mucosal forms)                   ##
##                                                                            ##
################################################################################
################################################################################

##----- Dependencies -----------------------------------------------------------

# Set system language
Sys.setenv(LANG = "en")

# Data Import, Manipulation & Visualisation
library(tidyverse)
library(RColorBrewer)
library(circlize)
library(xlsx)
library(readxl)
library(writexl)
library(reshape2)
library(ggpubr)
library(ggcorrplot)
library(ggbeeswarm)
library(ggcorrplot)
library(CompBioTools)
library(volcano3D)
library(pROC)
library(plotly)

# Bioconductor Packages
library(AnnotationHub)
library(ensembldb)
library(ComplexHeatmap)
library(sva)
library(DESeq2)
library(tximport)
library(GSVA)

# Other packages
library(omnideconv)

## Setting up colors
HC_colour <- "#f8766d"
LCL_colour <- "#7cae00"
MCL_colour <- "#00bfc4"
DCL_colour <- "#c77cff"

## Life, Universe and Everything
set.seed(42)

#@ Improving plotting in plot panel
trace(grDevices:::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)

trace(grDevices::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)


##----- Importing data ---------------------------------------------------------

## Importing samples data

# Phenotypic data
all_pdata <- as.data.frame(read_excel("data/metadata/all_pdata_updated.xlsx"))
all_pdata$group <- factor(all_pdata$group, 
                          levels = c("HC", "DCL", "LCL", "MCL"))
levels(all_pdata$group)
rownames(all_pdata) <- all_pdata$SRA

# List of salmon directories
sample_names <- list.files("data/salmon_seqbias")

# Removing .salmon from name
sample_names <- str_replace_all(sample_names, ".salmon", "")

# Getting quantification files directories
quant_dirs <- list.files("data/salmon_seqbias", pattern = "quant.sf", 
                         recursive = TRUE, full.names = TRUE)
names(quant_dirs) <- sample_names
head(quant_dirs)

## Fetching ensdb from AnnotationHub to summarize to gene level

# Getting AnnotationHub
ah <- AnnotationHub()

# Querying for EnsDB release 101
query(ah, c("EnsDb", "Homo sapiens", "101"))

# Downloading
ensdb <- ah[["AH83216"]]
keytypes(ensdb)

##----- Assembling count data using tximport -----------------------------------

## Preparing tx2 gene data.frame

# Creating complete tx2gene data.frame using ENSDB
k <- keys(ensdb, keytype = "TXNAME")
tx2gene_ensdb <- ensembldb::select(ensdb, k, 
                                   columns = c("GENEID", "SYMBOL", "GENEBIOTYPE",
                                               "DESCRIPTION"),
                                   keytype = "TXNAME")

table(tx2gene_ensdb$GENEBIOTYPE)

# Checking 'unknown' transcripts
novel_tx <- tx2gene_ensdb[str_detect(tx2gene_ensdb$SYMBOL, ".*\\.\\d?"),]

# Filtering to obtain only known protein coding, IG and TR genes
tx2gene_ensdb_biotypes <- tx2gene_ensdb %>%
  dplyr::filter(str_detect(GENEBIOTYPE, 
                           "protein_coding|IG_[:alpha:]{1}_gene|TR_[:alpha:]{1}_gene")) %>%
  dplyr::filter(!str_detect(SYMBOL, ".*\\.\\d?")) %>%
  dplyr::filter(!str_detect(DESCRIPTION, "readthrough"))


table(tx2gene_ensdb_biotypes$GENEBIOTYPE)

## Counting transcripts with tximport

# Summarizing all_leish dataset using gene symbols
txi_all_leish <- tximport(quant_dirs,
                          type = "salmon",
                          tx2gene = tx2gene_ensdb_biotypes[,c(1,3)], 
                          ignoreTxVersion = TRUE)

dim(txi_all_leish$counts)


# Reordering datasets to match
identical(colnames(txi_all_leish$counts), all_pdata$SRA)
all_pdata <- all_pdata[order(all_pdata$SRA),]
identical(colnames(txi_all_leish$counts), all_pdata$SRA)

##----- Estimating and correcting batch effects with sva package ---------------

# Creating initial dds object
dds <- DESeqDataSetFromTximport(txi_all_leish, 
                                colData = all_pdata,
                                design = ~ group)
dds <- DESeq2::estimateSizeFactors(dds)
dim(dds)

# Filtering low expression genes
norm_counts <- counts(dds, normalized = TRUE)
keep <- rowSums(norm_counts >= 10) >= 5
dds <- dds[keep,]
dim(dds)

# VST to plot PCA before sva
vst_dds_b4 <- vst(dds, blind = FALSE)

## Setting up for sva

# Counts matrix
dat0 <- counts(dds, normalized=TRUE)

# Full model
mod1 <- model.matrix(~ group, data = colData(dds))

# Null model
mod0 <- model.matrix(~ 1, data = colData(dds))

# Running svaseq
svseq <- svaseq(dat0, mod1, mod0, n.sv = 3)
plot(svseq$sv, pch = 19, col = all_pdata$group)

# Plotting relationship between SVs and known confounders/groups
par(mfrow = c(2,3))
for (n in 1:3) {
  boxplot(svseq$sv[,n] ~ all_pdata$study)
  points(svseq$sv[,n] ~ jitter(as.numeric(as.factor(all_pdata$study))),
         col=as.numeric(as.factor(all_pdata$study)), pch=16, cex=1.5)
  
}

for (n in 1:3) {
  boxplot(svseq$sv[,n] ~ all_pdata$group)
  points(svseq$sv[,n] ~ jitter(as.numeric(as.factor(all_pdata$group))),
         col=as.numeric(as.factor(all_pdata$group)), pch=16, cex=1.5)
  
}
par(mfrow = c(1,1))

# Checking relationship between known covariates and SVs
summary(lm(svseq$sv ~ all_pdata$group))
summary(lm(svseq$sv ~ all_pdata$study))

##----- Creating dds object with SVs included in the design --------------------

# Preparing SVs
svseq.df <- svseq$sv
colnames(svseq.df) <- c("SV1", "SV2", "SV3")

# Joining pdata with surrogate variables
all_pdata.sva <- cbind(all_pdata, svseq.df)

# Creating new dds object with SVs in the design
dds.sva <- DESeqDataSetFromMatrix(counts(dds), colData = all_pdata.sva,
                                  design = ~ SV1 + SV2 + SV3 + group)

# Running DEseq2
dds.sva <- DESeq(dds.sva)

# Applying vst for further visualization of results
vst_dds.sva <- vst(dds.sva, blind = FALSE)

resultsNames(dds.sva)

# Saving results
sva_res.DCLvsHC <- results(dds.sva, name="group_DCL_vs_HC", alpha = 0.05)
sva_res.LCLvsHC <- results(dds.sva, name="group_LCL_vs_HC", alpha = 0.05)
sva_res.MCLvsHC <- results(dds.sva, name="group_MCL_vs_HC", alpha = 0.05)
sva_res.DCLvsLCL <- results(dds.sva, contrast=c("group","DCL","LCL"), alpha = 0.05)
sva_res.MCLvsLCL <- results(dds.sva, contrast=c("group","MCL","LCL"), alpha = 0.05)
sva_res.DCLvsMCL <- results(dds.sva, contrast=c("group","DCL","MCL"), alpha = 0.05)

# Creating list of results
sva_results <- list(DCLvsHC = sva_res.DCLvsHC,
                    LCLvsHC = sva_res.LCLvsHC,
                    MCLvsHC = sva_res.MCLvsHC,
                    DCLvsLCL = sva_res.DCLvsLCL,
                    MCLvsLCL = sva_res.MCLvsLCL,
                    DCLvsMCL = sva_res.DCLvsMCL)

##----- Plotting PCAs for "clean" datasets -------------------------------------

# Separating vst uncleaned
vst_dds_batches <- vst_dds.sva

# Cleaning vst
assay(vst_dds.sva) <- limma::removeBatchEffect(assay(vst_dds.sva), 
                                               covariates = svseq.df)

## Cleaned dataset

# Run PCA 
pca_sva <- prcomp(t(assay(vst_dds.sva)))

# Rounding values
pcaVars_sva <- round((pca_sva$sdev^2/sum(pca_sva$sdev^2))*100,2)

# Checking if match with pdata
identical(rownames(pca_sva$x), all_pdata$SRA)

# Creating df for plotting with study/group vars
pca_sva.plot <- cbind(as.data.frame(pca_sva$x[,1:2]), Group = all_pdata$group, 
                      Study = all_pdata$study)

# Checking group var
pca_sva.plot$Group

# Saving plot
pca_plotB <- ggplot(pca_sva.plot, aes(x = PC1, y = PC2, colour = Group, 
                                      fill = Group, shape = Study)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  #coord_fixed() +
  labs(x = paste0("PC1: ", pcaVars_sva[1], "%"),
       y = paste0("PC2: ", pcaVars_sva[2], "%"),
       title = "PCA after regression of SVs") +
  theme_bw() +
  scale_colour_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, "DCL" = DCL_colour, "MCL" = MCL_colour),
                      labels = c("HC", "DCL", "LCL", "MLP")) +
  scale_fill_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, "DCL" = DCL_colour, "MCL" = MCL_colour),
                      labels = c("HC", "DCL", "LCL", "MLP")) +
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), 
        axis.title = element_text(size=12), plot.title = element_text(hjust = 0.5, size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(face = "bold"),
        aspect.ratio = 1)

## original dataset

# Run PCA 
pca_vst <- prcomp(t(assay(vst_dds_batches)))

# Rounding values
pcaVars_vst <- round((pca_vst$sdev^2/sum(pca_vst$sdev^2))*100,2)

# Checking if match with pdata
identical(rownames(pca_vst$x), all_pdata$SRA)

# Creating df for plotting with group/study vars
pca_vst.plot <- cbind(as.data.frame(pca_vst$x[,1:2]), Group = all_pdata$group, 
                      Study = all_pdata$study)

# Saving plot
pca_plotA <- ggplot(pca_vst.plot, aes(x = PC1, y = PC2, colour = Group, 
                                      fill = Group, shape = Study)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  #coord_fixed() +
  labs(x = paste0("PC1: ", pcaVars_vst[1], "%"),
       y = paste0("PC2: ", pcaVars_vst[2], "%"),
       title = "PCA before regression of SVs") +
  theme_bw() +
  scale_colour_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                                 "DCL" = DCL_colour, "MCL" = MCL_colour),
                      labels = c("HC", "DCL", "LCL", "MLP")) +
  scale_fill_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                               "DCL" = DCL_colour, "MCL" = MCL_colour),
                    labels = c("HC", "DCL", "LCL", "MLP")) +
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title = element_text(size=12), 
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

## Plotting

# Both plots together
tiff("results/figures/PCA_both.tiff", res = 300, type = "cairo", unit = "in", 
     width = 10, height = 4)

ggarrange(pca_plotA, 
          ggplot()+theme_minimal(), 
          pca_plotB,
          ggplot()+theme_minimal(),
          labels = c("A", " ", "B", " "), widths = c(1,0.2,1,0.15), nrow = 1,
          common.legend = TRUE, legend = "right")

dev.off()

# Only cleaned data
pca_plot_after <- ggplot(pca_sva.plot, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  scale_colour_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                                 "DCL" = DCL_colour, "MCL" = MCL_colour),
                      labels = c("HC", "DCL", "LCL", "MLP")) +
  #coord_fixed() +
  labs(x = paste0("PC1: ", pcaVars_sva[1], "%"),
       y = paste0("PC2: ", pcaVars_sva[2], "%"),
       title = "Principal Component Analysis") +
  theme_bw() +
  theme(axis.text.x = element_text(size=8), 
        axis.text.y = element_text(size=8), 
        axis.title = element_text(size=10),
        legend.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

##----- Preparing data flor plotting -------------------------------------------

# vst matrix (after cleaning batch effects)
sva_cleaned <- assay(vst_dds.sva)

# Centering counts around mean
sva_cleaned.centred <- sva_cleaned - rowMeans(sva_cleaned)

## ----- 3D Volcano ------------------------------------------------------------

# Preparing LRT test
dds.sva_lrt <- DESeq(dds.sva, test = "LRT", reduced = ~ SV1 + SV2 + SV3)
#res_volc3d <- deseq_polar(dds.sva, dds.sva_lrt, "group")

# Checking matches
identical(rownames(dds.sva_lrt), rownames(dds.sva))
identical(rownames(dds.sva_lrt), rownames(sva_res.DCLvsLCL))

# Preparing extra results comparison
sva_res.LCLvsMCL <- results(dds.sva, contrast=c("group","LCL","MCL"), 
                            alpha = 0.05)
sva_res.LCLvsDCL <- results(dds.sva, contrast=c("group","LCL","DCL"), 
                            alpha = 0.05)

# Getting p-values
pvals_volc3d <- data.frame(genes = rownames(dds.sva_lrt),
                           lrt = results(dds.sva_lrt)$padj,
                           LCL_vs_DCL = sva_res.LCLvsDCL$padj,
                           LCL_vs_MCL = sva_res.LCLvsMCL$padj,
                           DCL_vs_MCL = sva_res.DCLvsMCL$padj)

# Adjusting names
rownames(pvals_volc3d) <- pvals_volc3d$genes
pvals_volc3d$genes <- NULL

# Subsetting groups
patients_sra <- all_pdata[all_pdata$group!="HC", "SRA"]
patients_grp <- str_replace_all(all_pdata[all_pdata$group!="HC", "group"],
                                "MCL", "MLP")
patients_grp <- factor(patients_grp, levels = c("LCL", "DCL", "MLP"))

# Getting counts
counts_volc3d <- assay(vst_dds.sva)[,patients_sra]

# Colours
LvM_colour = "#3EB762"
LvD_colour = "#a2808d"
MvD_colour = "#649EE2"

# Saving polar graph
polar_graph <- polar_coords(outcome = patients_grp,
                            data = t(counts_volc3d),
                            pvals = as.matrix(pvals_volc3d),
                            scheme = c("gray70", DCL_colour, MvD_colour,
                                       MCL_colour, LvM_colour, LCL_colour,
                                       LvD_colour))



# Creating 3D view and adding animation
volc3dplot <- volcano3D(polar_graph)
add_animation(volc3dplot)

# Saving HTML widget
htmlwidgets::saveWidget(as_widget(add_animation(volc3dplot)), 
                        "results/volcano3D.html")

# Creating top-view plot
radial_volcano <- radial_ggplot(polar_graph, marker_outline_width = 0, 
                                marker_alpha = 1, axis_label_size = 0) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  scale_fill_manual(labels = c("NS", "DCL", "DCL + MLP",
                               "MLP", "LCL + MLP", "LCL",
                               "LCL + DCL"),
                    values = c("gray70", DCL_colour, MvD_colour,
                               MCL_colour, LvM_colour, LCL_colour,
                               LvD_colour),
                    name = "Group") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 10, face = "bold", color = "black", hjust = 0.5),
        legend.text = element_text(size = 10, color = "black"),
        legend.position = "right",
        legend.justification = "center")


# Saving
tiff("results/figures/3D_volcano_topview.tiff", w = 6, h = 4, res = 300, 
     type="cairo", units = "in")

print(radial_volcano)

dev.off()

# Saving plot for paper
tiff("results/paper/Figure1.tiff", w = 12, h = 4, res = 300, 
     type="cairo", units = "in")

ggarrange(pca_plot_after, radial_volcano, labels = "AUTO", nrow = 1)

dev.off()

##----- Senescence analysis ----------------------------------------------------

## Importing sets, adjusting name, reordering, checking DEGs and prepare matrix

# Senescence genes
genes_sen_filter <- read_excel("data/genes/genes_senescence_filter.xlsx")
genes_sen_filter <- geneNamesForPlotting(genes_sen_filter)
genes_sen_filter <- genes_sen_filter[order(genes_sen_filter$hgnc_symbol),]
genes_sen_filter <- geneCheck(sva_results, genes_sen_filter)
sen_filter_sva <- sva_cleaned[rownames(sva_cleaned) %in% 
                                genes_sen_filter$hgnc_symbol,]
identical(rownames(sen_filter_sva), genes_sen_filter$hgnc_symbol)
rownames(sen_filter_sva) <- genes_sen_filter$names_plot

## Preparing data for plotting violins

# Preparing dataframe
sen_filter_violin <- as.data.frame(t(sen_filter_sva))
identical(rownames(sen_filter_violin), all_pdata$SRA)
sen_filter_violin$Sample <- rownames(sen_filter_violin)
sen_filter_violin$Group <- all_pdata$group
sen_filter_violin$Study <- all_pdata$study

# Melting dataframe
sen_filter_violin <- melt(sen_filter_violin,
                          id.vars = c("Group", "Study", "Sample"), 
                          variable.name = "Gene", 
                          value.name = "Expression")

# Creating annotation df
identical(as.factor(genes_sen_filter$names_plot), 
          unique(sen_filter_violin$Gene))

annotation_sen_violin <- rbind(
  data.frame(
    Gene = genes_sen_filter$names_plot,
    Group1 = "HC",
    Group2 = "LCL",
    padj = genes_sen_filter$LCLvsHC_padj,
    log2fc = genes_sen_filter$LCLvsHC_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_sen_filter$names_plot,
    Group1 = "LCL",
    Group2 = "MCL",
    padj = genes_sen_filter$MCLvsLCL_padj,
    log2fc = genes_sen_filter$MCLvsLCL_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_sen_filter$names_plot,
    Group1 = "HC",
    Group2 = "MCL",
    padj = genes_sen_filter$MCLvsHC_padj,
    log2fc = genes_sen_filter$MCLvsHC_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_sen_filter$names_plot,
    Group1 = "LCL",
    Group2 = "DCL",
    padj = genes_sen_filter$DCLvsLCL_padj,
    log2fc = genes_sen_filter$DCLvsLCL_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_sen_filter$names_plot,
    Group1 = "DCL",
    Group2 = "MCL",
    padj = genes_sen_filter$DCLvsMCL_padj,
    log2fc = genes_sen_filter$DCLvsMCL_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_sen_filter$names_plot,
    Group1 = "HC",
    Group2 = "DCL",
    padj = genes_sen_filter$DCLvsHC_padj,
    log2fc = genes_sen_filter$DCLvsHC_log2FC,
    yposition = 0
  )
)

annotation_sen_violin <- calc_yPosition2(annotation_sen_violin,
                                         sen_filter_violin,
                                         "Expression",
                                         "Gene",
                                         "Gene",
                                         "Group2",
                                         "Group1",
                                         "padj",
                                         k = 0.03)

annotation_sen_violin <- pvalueCodes(annotation_sen_violin, "padj")

# Saving violin plot, but first adjusting order of factors
sen_filter_violin$Gene <- factor(sen_filter_violin$Gene, 
                                 levels = c("ATF5 (ATF)",
                                            "B3GAT1 (CD57)",
                                            "SESN2 (Sestrin 2)",
                                            "CDKN2A (p16)",
                                            "CDKN1A (p21)",
                                            "MAPK11 (p38-beta)",
                                            "NFKB2 (p52)",
                                            "TP53 (p53)"))

annotation_sen_violin$Gene <- factor(annotation_sen_violin$Gene, 
                                     levels = c("ATF5 (ATF)",
                                                "B3GAT1 (CD57)",
                                                "SESN2 (Sestrin 2)",
                                                "CDKN2A (p16)",
                                                "CDKN1A (p21)",
                                                "MAPK11 (p38-beta)",
                                                "NFKB2 (p52)",
                                                "TP53 (p53)"))

tiff("results/paper/Figure2A.tiff", res = 300, type = "cairo", unit = "in", 
     width = 10, height = 4)

ggplot(data = sen_filter_violin, aes(x = Gene, y = Expression, fill = Group)) +
  geom_point(aes(x = Group, y = Expression*1.16), alpha = 0) +
  geom_violin(aes(x = Group, y = Expression), color = "black") +
  facet_wrap(~ Gene, scale = "free", ncol = 4) +
  geom_bracket(inherit.aes = FALSE, 
               data = annotation_sen_violin[annotation_sen_violin$pcode!="",],
               vjust = 0.5, size = 0.3, tip.length = 0.0125, label.size = 3,
               aes(y.position = yposition, label = pcode, 
                   xmin = Group1, xmax = Group2)) +
  scale_fill_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                               "DCL" = DCL_colour, "MCL" = MCL_colour),
                    labels = c("HC", "DCL", "LCL", "MLP")) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.ticks.x = element_blank())

dev.off()

##----- Inflammation analysis --------------------------------------------------

## Importing sets, adjusting name, reordering, checking DEGs and prepare matrix

# Inflammation genes
genes_inf_filter <- read_excel("data/genes/genes_inflammation_filter.xlsx")
genes_inf_filter <- geneNamesForPlotting(genes_inf_filter)
genes_inf_filter <- genes_inf_filter[order(genes_inf_filter$hgnc_symbol),]
genes_inf_filter <- geneCheck(sva_results, genes_inf_filter)
inf_filter_sva <- sva_cleaned[rownames(sva_cleaned) %in% 
                                genes_inf_filter$hgnc_symbol,]
identical(rownames(inf_filter_sva), genes_inf_filter$hgnc_symbol)
rownames(inf_filter_sva) <- genes_inf_filter$names_plot

## Preparing data for plotting violins

# Preparing dataframe
inf_filter_violin <- as.data.frame(t(inf_filter_sva))
identical(rownames(inf_filter_violin), all_pdata$SRA)
inf_filter_violin$Sample <- rownames(inf_filter_violin)
inf_filter_violin$Group <- all_pdata$group
inf_filter_violin$Study <- all_pdata$study

# Melting dataframe
inf_filter_violin <- melt(inf_filter_violin,
                          id.vars = c("Group", "Study", "Sample"), 
                          variable.name = "Gene", 
                          value.name = "Expression")

# Creating annotation df
identical(as.factor(genes_inf_filter$names_plot), 
          unique(inf_filter_violin$Gene))

annotation_inf_violin <- rbind(
  data.frame(
    Gene = genes_inf_filter$names_plot,
    Group1 = "HC",
    Group2 = "LCL",
    padj = genes_inf_filter$LCLvsHC_padj,
    log2fc = genes_inf_filter$LCLvsHC_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_inf_filter$names_plot,
    Group1 = "LCL",
    Group2 = "MCL",
    padj = genes_inf_filter$MCLvsLCL_padj,
    log2fc = genes_inf_filter$MCLvsLCL_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_inf_filter$names_plot,
    Group1 = "HC",
    Group2 = "MCL",
    padj = genes_inf_filter$MCLvsHC_padj,
    log2fc = genes_inf_filter$MCLvsHC_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_inf_filter$names_plot,
    Group1 = "LCL",
    Group2 = "DCL",
    padj = genes_inf_filter$DCLvsLCL_padj,
    log2fc = genes_inf_filter$DCLvsLCL_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_inf_filter$names_plot,
    Group1 = "DCL",
    Group2 = "MCL",
    padj = genes_inf_filter$DCLvsMCL_padj,
    log2fc = genes_inf_filter$DCLvsMCL_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_inf_filter$names_plot,
    Group1 = "HC",
    Group2 = "DCL",
    padj = genes_inf_filter$DCLvsHC_padj,
    log2fc = genes_inf_filter$DCLvsHC_log2FC,
    yposition = 0
  )
)

annotation_inf_violin <- calc_yPosition2(annotation_inf_violin,
                                         inf_filter_violin,
                                         "Expression",
                                         "Gene",
                                         "Gene",
                                         "Group2",
                                         "Group1",
                                         "padj",
                                         k = 0.08)

annotation_inf_violin <- pvalueCodes(annotation_inf_violin, "padj")

# Saving violin plot

tiff("results/figures/Inflammation.tiff", res = 300, type = "cairo",
     unit = "in", width = 8, height = 6)

inf_violin_plot <- ggplot(data = inf_filter_violin, 
                          aes(x = Gene, y = Expression, fill = Group)) +
  geom_violin(aes(x = Group, y = Expression), color = "black") +
  geom_point(aes(x = Group, y = Expression*1.45), alpha = 0) +
  facet_wrap(~ Gene, scale = "free", ncol = 3) +
  geom_bracket(inherit.aes = FALSE, 
               data = annotation_inf_violin[annotation_inf_violin$pcode!="",],
               vjust = 0.5, size = 0.3, tip.length = 0.0125, label.size = 3,
               aes(y.position = yposition, label = pcode, 
                   xmin = Group1, xmax = Group2)) +
  scale_fill_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                               "DCL" = DCL_colour, "MCL" = MCL_colour),
                    labels = c("HC", "DCL", "LCL", "MLP")) +
  scale_y_continuous(breaks = function(x) seq(ceiling(x[1]), 
                                              floor(x[2]), by = 4)) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.ticks.x = element_blank(),
        plot.margin = ggplot2::margin(0.75, 0.25, 0.25, 0.25, "cm"))
print(inf_violin_plot)

dev.off()

##----- Cytoxicity analysis ----------------------------------------------------

## Importing sets, adjusting name, reordering, checking DEGs and prepare matrix

# Cytotoxicity genes
genes_cyto_filter <- read_excel("data/genes/genes_cytotoxicity_filter.xlsx")
genes_cyto_filter <- geneNamesForPlotting(genes_cyto_filter)
genes_cyto_filter <- genes_cyto_filter[order(genes_cyto_filter$hgnc_symbol),]
genes_cyto_filter <- geneCheck(sva_results, genes_cyto_filter)
cyto_filter_sva <- sva_cleaned[rownames(sva_cleaned) %in% 
                                 genes_cyto_filter$hgnc_symbol,]
identical(rownames(cyto_filter_sva), genes_cyto_filter$hgnc_symbol)
rownames(cyto_filter_sva) <- genes_cyto_filter$names_plot

## Preparing data for plotting violins

# Preparing dataframe
cyto_filter_violin <- as.data.frame(t(cyto_filter_sva))
identical(rownames(cyto_filter_violin), all_pdata$SRA)
cyto_filter_violin$Sample <- rownames(cyto_filter_violin)
cyto_filter_violin$Group <- all_pdata$group
cyto_filter_violin$Study <- all_pdata$study

# Melting dataframe
cyto_filter_violin <- melt(cyto_filter_violin,
                          id.vars = c("Group", "Study", "Sample"), 
                          variable.name = "Gene", 
                          value.name = "Expression")

# Creating annotation df
identical(as.factor(genes_cyto_filter$names_plot), 
          unique(cyto_filter_violin$Gene))

annotation_cyto_violin <- rbind(
  data.frame(
    Gene = genes_cyto_filter$names_plot,
    Group1 = "HC",
    Group2 = "LCL",
    padj = genes_cyto_filter$LCLvsHC_padj,
    log2fc = genes_cyto_filter$LCLvsHC_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_cyto_filter$names_plot,
    Group1 = "LCL",
    Group2 = "MCL",
    padj = genes_cyto_filter$MCLvsLCL_padj,
    log2fc = genes_cyto_filter$MCLvsLCL_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_cyto_filter$names_plot,
    Group1 = "HC",
    Group2 = "MCL",
    padj = genes_cyto_filter$MCLvsHC_padj,
    log2fc = genes_cyto_filter$MCLvsHC_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_cyto_filter$names_plot,
    Group1 = "LCL",
    Group2 = "DCL",
    padj = genes_cyto_filter$DCLvsLCL_padj,
    log2fc = genes_cyto_filter$DCLvsLCL_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_cyto_filter$names_plot,
    Group1 = "DCL",
    Group2 = "MCL",
    padj = genes_cyto_filter$DCLvsMCL_padj,
    log2fc = genes_cyto_filter$DCLvsMCL_log2FC,
    yposition = 0
  ),
  data.frame(
    Gene = genes_cyto_filter$names_plot,
    Group1 = "HC",
    Group2 = "DCL",
    padj = genes_cyto_filter$DCLvsHC_padj,
    log2fc = genes_cyto_filter$DCLvsHC_log2FC,
    yposition = 0
  )
)

annotation_cyto_violin <- calc_yPosition2(annotation_cyto_violin,
                                         cyto_filter_violin,
                                         "Expression",
                                         "Gene",
                                         "Gene",
                                         "Group2",
                                         "Group1",
                                         "padj",
                                         k = 0.08)

annotation_cyto_violin <- pvalueCodes(annotation_cyto_violin, "padj")

# Saving violin plot

tiff("results/figures/Cytotoxicity.tiff", res = 300, type = "cairo", 
     unit = "in", width = 8, height = 6)

cyto_violin_plot <- ggplot(data = cyto_filter_violin, 
                           aes(x = Gene, y = Expression, fill = Group)) +
  geom_violin(aes(x = Group, y = Expression), color = "black") +
  geom_point(aes(x = Group, y = Expression*1.5), alpha = 0) +
  facet_wrap(~ Gene, scale = "free", ncol = 3) +
  geom_bracket(inherit.aes = FALSE, 
               data = annotation_cyto_violin[annotation_cyto_violin$pcode!="",],
               vjust = 0.5, size = 0.3, tip.length = 0.0125, label.size = 3,
               aes(y.position = yposition, label = pcode, 
                   xmin = Group1, xmax = Group2)) +
  scale_fill_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                               "DCL" = DCL_colour, "MCL" = MCL_colour),
                    labels = c("HC", "DCL", "LCL", "MLP")) +
  scale_y_continuous(breaks = function(x) seq(ceiling(x[1]), 
                                              floor(x[2]), by = 4)) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.ticks.x = element_blank(),
        plot.margin = ggplot2::margin(0.75, 0.25, 0.25, 0.25, "cm"))
print(cyto_violin_plot)

dev.off()

## ----- Correlation analysis --------------------------------------------------

## Preparing data

# Creating matrices
inf_bind <- t(rbind(sen_filter_sva, inf_filter_sva)) %>% as.data.frame()
cyto_bind <- t(rbind(sen_filter_sva, cyto_filter_sva)) %>% as.data.frame()

# Correlation matrices
inf_cormat <- round(cor(inf_bind, method = "spearman"), 2)
inf_cormat <- inf_cormat[9:17, 1:8]
cyto_cormat <- round(cor(cyto_bind, method = "spearman"), 2)
cyto_cormat <- cyto_cormat[9:17, 1:8]

# P-values
inf_pmat <- cor_pmat(inf_bind, method = "spearman", exact = TRUE)
inf_pmat <- inf_pmat[9:17, 1:8]
cyto_pmat <- cor_pmat(cyto_bind, method = "spearman", exact = TRUE)
cyto_pmat <- cyto_pmat[9:17, 1:8]

# Melting, renaming, coding p
inf_cormat_melt <- melt(inf_cormat)
colnames(inf_cormat_melt) <- c("Gene", "Sen", "Correlation")
inf_pmat_melt <- melt(inf_pmat)
colnames(inf_pmat_melt) <- c("Gene", "Sen", "pvalue")
inf_pmat_melt <- pvalueCodes(inf_pmat_melt, "pvalue")
inf_cormat_melt$Sen <- factor(inf_cormat_melt$Sen, levels = c("ATF5 (ATF)",
                                                           "B3GAT1 (CD57)",
                                                           "SESN2 (Sestrin 2)",
                                                           "CDKN2A (p16)",
                                                           "CDKN1A (p21)",
                                                           "MAPK11 (p38-beta)",
                                                           "NFKB2 (p52)",
                                                           "TP53 (p53)"))

cyto_cormat_melt <- melt(cyto_cormat)
colnames(cyto_cormat_melt) <- c("Gene", "Sen", "Correlation")
cyto_pmat_melt <- melt(cyto_pmat)
colnames(cyto_pmat_melt) <- c("Gene", "Sen", "pvalue")
cyto_pmat_melt <- pvalueCodes(cyto_pmat_melt, "pvalue")
cyto_cormat_melt$Sen <- factor(cyto_cormat_melt$Sen, 
                               levels = c("ATF5 (ATF)",
                                          "B3GAT1 (CD57)",
                                          "SESN2 (Sestrin 2)",
                                          "CDKN2A (p16)",
                                          "CDKN1A (p21)",
                                          "MAPK11 (p38-beta)",
                                          "NFKB2 (p52)",
                                          "TP53 (p53)"))

# Adding p-value to cor df
identical(inf_cormat_melt$Gene, inf_pmat_melt$Gene)
inf_cormat_melt$pcode <- inf_pmat_melt$pcode
identical(cyto_cormat_melt$Gene, cyto_pmat_melt$Gene)
cyto_cormat_melt$pcode <- cyto_pmat_melt$pcode


## Plotting

# Color map
col_corrmaps <- colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
col_corrmaps(seq(-1, 1, length = 9))

# Inflammation
tiff("results/figures/Inflammation_Senescence_CorMap.tiff", units="in", width=6,
     height=6, res=300, type = "cairo")

ggplot(data = inf_cormat_melt, aes(Sen, Gene, fill = Correlation))+
  geom_tile(color = "white") +
  labs(title = "Inflammation", x = "", y = "") +
  scale_fill_stepsn(colours = col_corrmaps(seq(-1, 1, length = 9)),
                    name="Spearman correlation",
                    limits = c(-1,1),
                    n.breaks = 9,
                    show.limits = TRUE) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, 
                                   size = 10, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Sen, Gene, label = pcode), color = "black", size = 3)+
  scale_y_discrete(position = "right")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.ticks = element_blank(),
        legend.position = "left",
        legend.direction = "vertical", 
        legend.text=element_text(size=8),
        legend.title = element_text(size=10, angle = 90), 
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 9), 
        axis.text.y = element_text(size = 9)) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 10,
                               title.position = "left", title.hjust = 0.5))

dev.off()

# Cytotoxicity
tiff("results/figures/Cytotoxicity_Senescence_CorMap.tiff", units="in", width=6,
     height=6, res=300, type = "cairo")

ggplot(data = cyto_cormat_melt, aes(Sen, Gene, fill = Correlation))+
  geom_tile(color = "white") +
  labs(title = "Cytotoxicity", x = "", y = "") +
  scale_fill_stepsn(colours = col_corrmaps(seq(-1, 1, length = 9)),
                    name="Spearman correlation",
                    limits = c(-1,1),
                    n.breaks = 9,
                    show.limits = TRUE) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, 
                                   size = 10, hjust = 1))+
  coord_fixed() +
  geom_text(aes(Sen, Gene, label = pcode), color = "black", size = 3)+
  scale_y_discrete(position = "right")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "left",
        legend.direction = "vertical", 
        legend.text=element_text(size=8),
        legend.title = element_text(size=10, angle = 90), 
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 9), 
        axis.text.y = element_text(size = 9)) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 10,
                               title.position = "left", title.hjust = 0.5))

dev.off()

##----- Preparing data for omnideconv ------------------------------------------

# Reading data
sc_matrix <- read_csv("data/signature/skin_sc_matrix.csv")
sc_matrix <- as.data.frame(sc_matrix)
rownames(sc_matrix ) <- sc_matrix$index
sc_matrix$index <- NULL
sc_matrix <- as.data.frame(t(sc_matrix))
sc_obs <- read_csv("data/signature/skin_sc_obs.csv")

# Checking if matches
identical(colnames(sc_matrix), sc_obs$index)

# Extracting ids
sc_ids <- sc_obs$reduced_celltypes

# Build signature matrix
sig_matrix <- build_model(sc_matrix, sc_ids,
                          method = "scaden",
                          bulk_gene_expression = txi_all_leish$counts,
                          model_path = "data/signature/scaden_model")

# It does not store the matrix, but a path which might not work properly
sig_matrix <- "data/signature/scaden_model"

# Running deconv
deconv_scaden <- deconvolute(txi_all_leish$counts,
                             sig_matrix,
                             method = "scaden")

rowSums(deconv_scaden)
head(deconv_scaden) * 100

## Prepare for plotting

# Converting to df
deconv_scaden <- as.data.frame(deconv_scaden)

# Adding metadata
identical(rownames(deconv_scaden), all_pdata$SRA)
deconv_scaden$Group <- all_pdata$group

# Melting
deconv_scaden_melt <- melt(deconv_scaden,
                           id.vars = c("Group"),
                           value.name = "Proportion",
                           variable.name = "Cell")

deconv_scaden_melt$Proportion <- deconv_scaden_melt$Proportion*100

# Looping to get stats for each cell type correctly
for (cell in unique(deconv_scaden_melt$Cell)) {
  if (match(cell, unique(deconv_scaden_melt$Cell)) == 1) {
    manual_stats_scaden <- compare_means(Proportion ~ Group, 
                                         data = deconv_scaden_melt[deconv_scaden_melt$Cell==cell,],
                                         method = "wilcox.test", 
                                         paired = FALSE, 
                                         group.by = "Cell", 
                                         p.adjust.method = "BH")
  } else {
    temp <- compare_means(Proportion ~ Group, 
                          data = deconv_scaden_melt[deconv_scaden_melt$Cell==cell,],
                          method = "wilcox.test", 
                          paired = FALSE, 
                          group.by = "Cell", 
                          p.adjust.method = "BH")
    manual_stats_scaden <- rbind(manual_stats_scaden, temp)
  }
  
}


# Extracting annotation for plotting brackets
annotation_scaden <- as.data.frame(manual_stats_scaden[,c(1,3,4,6)])
colnames(annotation_scaden) <- c("Cell", "Group1", "Group2", "padj")

# Separating factors for ordering
scaden_cell_factors <- unique(annotation_scaden$Cell)

# Converting columns to character
annotation_scaden$Cell <- as.character(annotation_scaden$Cell)
deconv_scaden_melt$Cell <- as.character(deconv_scaden_melt$Cell)

# Removing not significant padj to avoid plotting
annotation_scaden$padj <- ifelse(annotation_scaden$padj <= 0.05, 
                                 annotation_scaden$padj, 
                                 NA)

# Calculatting bracket positions
annotation_scaden <- calc_yPosition2(annotation_scaden,
                                     deconv_scaden_melt,
                                     "Proportion",
                                     "Cell",
                                     "Cell",
                                     "Group2",
                                     "Group1",
                                     "padj",
                                     k = 0.1)

# Calculating * for pvalues
annotation_scaden <- pvalueCodes(annotation_scaden, "padj")

# Adjusting levels for plotting order
deconv_scaden_melt$Group <- factor(deconv_scaden_melt$Group, 
                                   levels = c("HC", "DCL", "LCL", "MCL"))
deconv_scaden_melt$Cell <- factor(deconv_scaden_melt$Cell, 
                                  levels = scaden_cell_factors)
annotation_scaden$Cell <- factor(annotation_scaden$Cell, 
                                 levels = scaden_cell_factors)
cell_filter <- c("CD4 T Cell", "CD8 T Cell", "NK Cell", 
                 "Plasma", "Pericyte", "Schwann", "Treg")

# Saving total figure
png("results/figures/SCADEN_Deconvolution.png", w = 1250, h = 1250, 
    type = "cairo-png", res = 100)

ggplot(deconv_scaden_melt, aes(x = Group, y = Proportion, fill = Group)) +
  facet_wrap(~ Cell, scales = "free", ncol = 4) +
  geom_bar(position = "dodge", stat = "summary", fun = mean, color = "black") +
  geom_point(aes(y = Proportion*1.3), alpha = 0) +
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.4, width = 0.4) +
  geom_bracket(inherit.aes = FALSE, data = annotation_scaden,
               vjust = 0.5, size = 0.4, tip.length = 0.025, label.size = 4,
               aes(y.position = yposition*0.85, label = pcode, xmin = Group1, 
                   xmax = Group2)) +
  scale_fill_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                               "DCL" = DCL_colour, "MCL" = MCL_colour),
                    labels = c("HC", "DCL", "LCL", "MLP")) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.ticks.x = element_blank())

dev.off()

# Saving supplementary figure 1
tiff("results/paper/SuppFig1.tiff", w = 10, h = 6, type = "cairo", 
     res = 300, units = "in")

ggplot(deconv_scaden_melt[!deconv_scaden_melt$Cell%in%cell_filter,], 
       aes(x = Group, y = Proportion, fill = Group)) +
  facet_wrap(~ Cell, scales = "free", ncol = 4) +
  geom_bar(position = "dodge", stat = "summary", fun = mean, color = "black") +
  geom_point(aes(y = Proportion*1.3), alpha = 0) +
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.4, width = 0.4) +
  geom_bracket(inherit.aes = FALSE, 
               data = annotation_scaden[!annotation_scaden$Cell%in%cell_filter,],
               vjust = 0.5, size = 0.4, tip.length = 0.025, label.size = 4,
               aes(y.position = yposition*0.85, label = pcode, 
                   xmin = Group1, xmax = Group2)) +
  scale_fill_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                               "DCL" = DCL_colour, "MCL" = MCL_colour),
                    labels = c("HC", "DCL", "LCL", "MLP")) +
  scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.ticks.x = element_blank())

dev.off()

## Specific cells for figure
cells_plot <- c("CD4 T Cell", "CD8 T Cell", "NK Cell")

deconv_scaden_melt_b <- deconv_scaden_melt[deconv_scaden_melt$Cell %in% 
                                             cells_plot,]
annotation_scaden_b <- annotation_scaden[annotation_scaden$Cell %in%
                                           cells_plot,]

tiff("results/figures/SCADEN_Deconvolution_subset.tiff", w = 7, h = 2.5, 
     type = "cairo", res = 300, units = "in")

scaden_subset_plot <- ggplot(deconv_scaden_melt_b, 
                             aes(x = Group, y = Proportion, fill = Group)) +
  facet_wrap(~ Cell, scales = "free", ncol = 3) +
  geom_bar(position = "dodge", stat = "summary", fun = mean, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.4, width = 0.4) +
  geom_bracket(inherit.aes = FALSE, data = annotation_scaden_b,
               vjust = 0.5, size = 0.4, tip.length = 0.0175, label.size = 3,
               aes(y.position = yposition*0.61, label = pcode, 
                   xmin = Group1, xmax = Group2)) +
  scale_fill_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                               "DCL" = DCL_colour, "MCL" = MCL_colour),
                    labels = c("HC", "DCL", "LCL", "MLP")) +
  scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        axis.ticks.x = element_blank(),
        plot.margin = ggplot2::margin(0.75, 0.25, 0.25, 0.25, "cm"))
print(scaden_subset_plot)

dev.off()


## Saving paper figure 3
tiff("results/paper/Figure3.tiff", w = 12, h = 8, type = "cairo", 
     res = 300, units = "in")

ggarrange(ggarrange(inf_violin_plot + theme(legend.position = "none"), 
          cyto_violin_plot + theme(legend.position = "none"), nrow = 1, 
          labels = c("A", "B")), 
          ggarrange(scaden_subset_plot, empty_plot, widths = c(1, 0.75)), 
          common.legend = FALSE, labels = c("", "C"), 
          nrow = 2, heights = c(1,0.5))

dev.off()

## ----- GSVA of inflammation, cytotoxicity and senescence ---------------------

## Preparing data

# Creating list of genesets
genesets_list <- list("Inflammation" = genes_inf_filter$hgnc_symbol,
                      "Cytotoxicity" = genes_cyto_filter$hgnc_symbol,
                      "Senescence" = genes_sen_filter$hgnc_symbol)

# Performing GSVA
gsva_results <- gsva(
  sva_cleaned,
  genesets_list,
  method = "gsva",
  kcdf = "Gaussian",
  min.sz = 5,
  max.sz = 500,
  mx.diff = TRUE,
  verbose = FALSE
)

# Adjusting table
gsva_results <- as.data.frame(t(gsva_results))

# Checking if orders match
identical(rownames(gsva_results), all_pdata$SRA)

# Adding extra columns
gsva_results$Sample <- rownames(gsva_results)
gsva_results$Group <- all_pdata$group
gsva_results$LesionSize <- all_pdata$lesion_size_mm2
gsva_results$Age <- all_pdata$age_years

# Cor testing
cor.test(gsva_results$Inflammation, gsva_results$Senescence, 
         method = "spearman")
cor.test(gsva_results$Cytotoxicity, gsva_results$Senescence, 
         method = "spearman")
cor.test(gsva_results$Cytotoxicity, gsva_results$Inflammation, 
         method = "spearman")

## Saving plots for final figure

# Correlations
corplot_sencyto <- ggplot(data = gsva_results, 
                          aes(y = Cytotoxicity, x = Senescence)) +
  geom_point(size = 2, aes(color = Group)) +
  geom_smooth(method = lm, color = "#FD4646", fill = "#FD4646", alpha = 0.2) +
  stat_cor(size = 3, method = "spearman",
           label.x.npc = "left", label.y.npc = "top",
           digits = 3, cor.coef.name = "rho",
           aes(label = paste(..r.label.., 
                             cut(..p..,
                                 breaks = c(-Inf, 0.0001, 0.001, 
                                            0.01, 0.05, Inf),
                                 labels = c("'****'", "'***'", 
                                            "'**'", "'*'", "'ns'")), 
                             sep = "~`, `~"))) +
  scale_color_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                                "DCL" = DCL_colour, "MCL" = MCL_colour),
                    labels = c("HC", "DCL", "LCL", "MLP")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold", color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.ticks.x = element_blank(),
        plot.margin = ggplot2::margin(0.75, 0.25, 0.25, 0.25, "cm")) +
  labs(x = "Cytotoxicity score", y = "Senescence score") +
  lims(x = c(-1,1), y = c(-1,1)) 

corplot_seninfl <- ggplot(data = gsva_results, 
                          aes(y = Inflammation, x = Senescence)) +
  geom_point(size = 2, aes(color = Group)) +
  geom_smooth(method = lm, color = "#FD4646", fill = "#FD4646", alpha = 0.2) +
  stat_cor(size = 3, method = "spearman",
           label.x.npc = "left", label.y.npc = "top",
           digits = 3, cor.coef.name = "rho",
           aes(label = paste(..r.label.., 
                             cut(..p..,
                                 breaks = c(-Inf, 0.0001, 0.001,
                                            0.01, 0.05, Inf),
                                 labels = c("'****'", "'***'", 
                                            "'**'", "'*'", "'ns'")), 
                             sep = "~`, `~"))) +
  scale_color_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                                "DCL" = DCL_colour, "MCL" = MCL_colour),
                     labels = c("HC", "DCL", "LCL", "MLP")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold", color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.ticks.x = element_blank(),
        plot.margin = ggplot2::margin(0.75, 0.25, 0.25, 0.25, "cm")) +
  labs(x = "Inflammation score", y = "Senescence score") +
  lims(x = c(-1,1), y = c(-1,1))

# Genes
corplot_cytoinfl <- ggplot(data = gsva_results, 
                           aes(x = Cytotoxicity, y = Inflammation)) +
  geom_point(size = 2, aes(color = Group)) +
  geom_smooth(method = lm, color = "#FD4646", fill = "#FD4646", alpha = 0.2) +
  stat_cor(size = 3, method = "spearman",
           label.x.npc = "left", label.y.npc = "top",
           digits = 3, cor.coef.name = "rho",
           aes(label = paste(..r.label.., 
                             cut(..p.., 
                                 breaks = c(-Inf, 0.0001, 0.001, 
                                            0.01, 0.05, Inf),
                                 labels = c("'****'", "'***'",
                                            "'**'", "'*'", "'ns'")), 
                             sep = "~`, `~"))) +
  scale_color_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                                "DCL" = DCL_colour, "MCL" = MCL_colour),
                     labels = c("HC", "DCL", "LCL", "MLP")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold", color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.ticks.x = element_blank(),
        plot.margin = ggplot2::margin(0.75, 0.25, 0.25, 0.25, "cm")) +
  labs(x = "Inflammation score", y = "Cytotoxicity score") +
  lims(x = c(-1,1), y = c(-1,1))

## ----- PLotting score differences --------------------------------------------

## Preparing data

# Melting
gsva_melt <- melt(gsva_results[,c(1:5)],
                  id.vars = c("Sample", "Group"),
                  variable.name = "Signature",
                  value.name = "Score")

# Looping to get stats correctly
for (signature in unique(gsva_melt$Signature)) {
  if (match(signature, unique(gsva_melt$Signature)) == 1) {
    manual_stats_gsva <- compare_means(Score ~ Group, 
                                       data = gsva_melt[gsva_melt$Signature==signature,],
                                       method = "wilcox.test", paired = FALSE, 
                                       group.by = "Signature", 
                                       p.adjust.method = "BH")
  } else {
    temp <- compare_means(Score ~ Group, 
                          data = gsva_melt[gsva_melt$Signature==signature,],
                          method = "wilcox.test", paired = FALSE, 
                          group.by = "Signature", 
                          p.adjust.method = "BH")
    manual_stats_gsva <- rbind(manual_stats_gsva, temp)
  }
  
}

# Extracting annotation for plotting brackets
annotation_gsva <- as.data.frame(manual_stats_gsva[,c(1,3,4,6)])
colnames(annotation_gsva) <- c("Signature", "Group1", "Group2", "padj")

# Separating factors for ordering
gsva_signature_factors <- unique(annotation_gsva$Signature)

# Converting columns to character
annotation_gsva$Signature <- as.character(annotation_gsva$Signature)
gsva_melt$Signature <- as.character(gsva_melt$Signature)

# Removing not significant padj to avoid plotting
annotation_gsva$padj <- ifelse(annotation_gsva$padj <= 0.05, 
                               annotation_gsva$padj, 
                               NA)

# Calculatting bracket positions
annotation_gsva <- calc_yPosition2(annotation_gsva,
                                   gsva_melt,
                                   "Score",
                                   "Signature",
                                   "Signature",
                                   "Group2",
                                   "Group1",
                                   "padj",
                                   k = 0.3)

# Calculating * for pvalues
annotation_gsva <- pvalueCodes(annotation_gsva, "padj")

# Adjusting levels for plotting order
gsva_melt$Group <- factor(gsva_melt$Group, 
                          levels = c("HC", "DCL", "LCL", "MCL"))
gsva_melt$Signature <- factor(gsva_melt$Signature, 
                              levels = gsva_signature_factors)
annotation_gsva$Signature <- factor(annotation_gsva$Signature, 
                                    levels = gsva_signature_factors)

tiff("results/figures/GSVA_scores_difference.tiff", w = 8, h = 2.5, 
     type = "cairo", res = 300, units = "in")

gsva_stats_plot <- ggplot(gsva_melt, aes(x = Group, y = Score, fill = Group)) +
  facet_wrap(~ Signature, scales = "free", ncol = 3) +
  geom_violin() +
  geom_point(aes(y = Score+1.4), alpha = 0) +
  geom_bracket(inherit.aes = FALSE, data = annotation_gsva,
               vjust = 0.5, size = 0.4, tip.length = 0.025, label.size = 4,
               aes(y.position = yposition, label = pcode, 
                   xmin = Group1, xmax = Group2)) +
  scale_fill_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                               "DCL" = DCL_colour, "MCL" = MCL_colour),
                    labels = c("HC", "DCL", "LCL", "MLP")) +
  scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        axis.ticks.x = element_blank(),
        plot.margin = ggplot2::margin(0.75, 0.25, 0.25, 0.25, "cm"))
print(gsva_stats_plot)

dev.off()

## ----- SCADEN correlation with gene set scores ------------------------------

## SCADEN correlations

# Preparing deconvolution data.frame
deconv_scaden_100 <- deconv_scaden[,1:19]*100
identical(rownames(deconv_scaden_100), all_pdata$SRA)
deconv_scaden_100$Study <- all_pdata$study
deconv_scaden_100$Group <- all_pdata$group
identical(rownames(gsva_results), rownames(deconv_scaden_100))
scores_scaden_bind <- as.data.frame(cbind(deconv_scaden_100[,-21], 
                                          gsva_results))

# Melting twice for facet grid
scores_scaden_melt <- melt(scores_scaden_bind[,c(1,2,14,20:27)],
                           id.vars = c(cells_plot, "Study", "Group", 
                                       "LesionSize", "Age", "Sample"),
                           variable.name = "Signature",
                           value.name = "Score")

scores_scaden_melt <- melt(scores_scaden_melt,
                           id.vars = c("Study", "Group", "LesionSize", "Age", 
                                       "Sample", "Signature", "Score"),
                           variable.name = "Cell",
                           value.name = "Proportion")

# Plotting

tiff("results/figures/Correlations_Cells_Scores.tiff", w = 8, h = 6, 
     type = "cairo", res = 300, units = "in")

cell_scores_corplot <- ggplot(scores_scaden_melt, 
                              aes(y = `Proportion`, x = `Score`)) +
  facet_grid(rows = vars(Cell), cols = vars(Signature), scale = "free") +
  #geom_point(aes(y = adjustment), alpha = 0) +
  geom_point(size = 1, aes(colour=Group)) +
  geom_smooth(method = lm, color = "#FD4646", fill = "#FD4646", alpha = 0.2) +
  stat_cor(size = 3, method = "spearman",
           label.x.npc = "left", label.y.npc = "top",
           digits = 2, cor.coef.name = "rho",
           aes(label = paste(after_stat(r.label), 
                             cut(..p.., 
                                 breaks = c(-Inf, 0.0001, 0.001, 
                                            0.01, 0.05, Inf),
                                 labels = c("'****'", "'***'", 
                                            "'**'", "'*'", "'ns'")), 
                             sep = "~`, `~"))) +
  scale_colour_manual(values = c("HC" = HC_colour, "LCL" = LCL_colour, 
                                 "DCL" = DCL_colour, "MCL" = MCL_colour),
                      labels = c("HC", "DCL", "LCL", "MLP")) +
  scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 22, size = 5, 
                                                  color = "black", 
                                                  fill = c(HC_colour, 
                                                           DCL_colour, 
                                                           LCL_colour, 
                                                           MCL_colour)))) +
  theme(axis.title = element_text(size = 12, face = "bold", color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, face = "bold", color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        axis.ticks.x = element_blank(),
        plot.margin = ggplot2::margin(0.75, 0.25, 0.25, 0.25, "cm")) +
  labs(x = "GSVA Score", y = "Cell proportion")
print(cell_scores_corplot)

dev.off()

## Figure for paper
empty_plot <- ggplot() + theme_void()

scores_cor_plots <- ggarrange(corplot_sencyto, corplot_seninfl, 
                              corplot_cytoinfl, nrow = 1, 
                              legend = "none", align = "h")
gsva_stats_plot
cell_scores_corplot

tiff("results/paper/Figure4.tiff", w = 12, h = 5, 
     type = "cairo", res = 300, units = "in")

ggarrange(ggarrange(gsva_stats_plot, empty_plot, scores_cor_plots, nrow = 3, 
                    legend = "none", heights = c(0.8,0.1,1), 
                    labels = c("A", "", "B")),
          empty_plot, cell_scores_corplot, ncol = 3, common.legend = TRUE, 
          legend = "top", labels = c("", "", "C"), widths = c(1,0.05,0.8))

dev.off()

## ----- Scores ROC AUC --------------------------------------------------------

# Setting up df
gsva_mcl <- gsva_results[gsva_results$Group%in%c("LCL", "MCL"),c(1:3, 5)]
gsva_mcl$Group <- droplevels(gsva_mcl$Group)

# ROC
scores_roc <- roc(formula = Group ~ ., data = gsva_mcl)

# Extracting data for plotting
scores_roc_df <- data.frame()

# Looping to create dataframe with each score
for (i in 1:length(scores_roc)) {
  df_temp <- data.frame(Score = names(scores_roc)[i],
                        y = scores_roc[[i]]$sensitivities,
                        x = scores_roc[[i]]$specificities)
  scores_roc_df <- rbind(scores_roc_df, df_temp)
  rm(df_temp)
}

## Legend AUCs

# Creating df
scores_leg_df <- data.frame()

# Looping through scores
for (i in 1:length(scores_roc)) {
  
  auc <- auc(scores_roc[[i]])
  ci <- ci.auc(scores_roc[[i]])
  ci_l <- round(ci[1], 2)
  ci_u <- round(ci[3], 2)
  
  df_temp <- data.frame(Score = names(scores_roc)[i],
                        AUC = paste0(" AUC = ", round(auc, 3), 
                                     "\n(95% CI = ", ci_l, " - ", ci_u, ")"))
  scores_leg_df <- rbind(scores_leg_df, df_temp)
  rm(df_temp)
}

# Saving figure
tiff("results/figures/Figure5A.tiff", w = 7, h = 3, res = 300, 
     type="cairo", units = "in")

roc_scores_plot <- ggplot(scores_roc_df[order(scores_roc_df$x),], 
                          aes(x=1-x, y=y)) +
  facet_wrap(~Score, nrow=1) +
  geom_path(colour = "#FF2345", size = 0.75) +
  geom_polygon(fill = "#FF2345", alpha=0.1, colour=0) +
  labs(x = "FPR (1 - Specificity)", y = "TPR (Sensitivity)") +
  geom_segment(aes(y = 0, yend = 1, x = 0, xend = 1), color = "grey50", 
               linetype="dashed", size=0.25) +
  geom_text(data = scores_leg_df, x = 0.7, y = 0.1, size = 2.5, 
            aes(label=AUC)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, 
                                   hjust = 1.2, vjust = 1.2), 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 12, face = "bold"),
        plot.margin = ggplot2::margin(0.75, 0.25, 0.25, 0.25, "cm"))
print(roc_scores_plot)

dev.off()

## ----- ROC AUC of senescence genes -------------------------------------------

# Sen genes
genes_sen_filter

# LCL MCL pacients
lclmcl_pct <- all_pdata[all_pdata$group %in% c("MCL", "LCL"), "SRA"]

# Sva cleaned data
sva_cleaned

# Subsetting
sva_roc_df <- as.data.frame(
  t(sva_cleaned[genes_sen_filter$hgnc_symbol,lclmcl_pct]))
identical(colnames(sva_roc_df), genes_sen_filter$hgnc_symbol)
colnames(sva_roc_df) <- genes_sen_filter$names_plot

# Adding grouping
identical(rownames(sva_roc_df), all_pdata[all_pdata$SRA %in% lclmcl_pct, "SRA"])
sva_roc_df$Group <- relevel(all_pdata[all_pdata$SRA %in% lclmcl_pct, "group"], 
                            ref = "LCL") %>% droplevels()

## Preparing data for plotting

# Creating dataframes
roc_df = data.frame()
leg_df = data.frame()

# Looping to get values for each gene
for (gene in colnames(sva_roc_df[,1:8])) {
  
  roc_formula <- formula(paste0("Group ~ `", gene, "`"))
  roc_res <- roc(formula = roc_formula, data = sva_roc_df)
  
  auc <- auc(roc_res)
  ci <- ci.auc(roc_res)
  ci_l <- round(ci[1], 2)
  ci_u <- round(ci[3], 2)
  
  temp <- data.frame(x = roc_res$specificities,
                     y = roc_res$sensitivities,
                     gene = gene)
  
  roc_df <- rbind(roc_df, temp)
  temp2 <- data.frame(gene = gene,
                      auc = paste0(" AUC = ", round(auc, 3), 
                                   "\n(95% CI = ", ci_l, " - ", ci_u, ")"))
  leg_df <- rbind(leg_df, temp2)
  
}

# Adjusting order
roc_df$gene <- factor(roc_df$gene, levels = c("ATF5 (ATF)",
                                              "B3GAT1 (CD57)",
                                              "SESN2 (Sestrin 2)",
                                              "CDKN2A (p16)",
                                              "CDKN1A (p21)",
                                              "MAPK11 (p38-beta)",
                                              "NFKB2 (p52)",
                                              "TP53 (p53)"))

leg_df$gene <- factor(leg_df$gene, levels = c("ATF5 (ATF)",
                                              "B3GAT1 (CD57)",
                                              "SESN2 (Sestrin 2)",
                                              "CDKN2A (p16)",
                                              "CDKN1A (p21)",
                                              "MAPK11 (p38-beta)",
                                              "NFKB2 (p52)",
                                              "TP53 (p53)"))

tiff("results/figures/Figure5B.tiff", w = 8, h = 4, res = 300, 
     type="cairo", units = "in")

roc_genes_plot <- ggplot(roc_df[order(roc_df$x),], aes(x=1-x, y=y)) +
  facet_wrap(~gene, nrow=2) +
  geom_path(colour = "#FF2345", size = 0.75) +
  geom_polygon(fill = "#FF2345", alpha=0.1, colour=0) +
  labs(x = "FPR (1 - Specificity)", y = "TPR (Sensitivity)") +
  geom_segment(aes(y = 0, yend = 1, x = 0, xend = 1), color = "grey50", 
               linetype="dashed", size=0.25) +
  geom_text(data = leg_df, x = 0.7, y = 0.1, size = 2.5, aes(label=auc)) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 10, angle = 45, 
                                 hjust = 1.2, vjust = 1.2), 
        axis.title.x=element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        strip.text.x = element_text(size = 12, face = "bold"),
        plot.margin = ggplot2::margin(0.75, 0.25, 0.25, 0.25, "cm"))
print(roc_genes_plot)

dev.off()


## Saving figure for paper

tiff("results/paper/Figure5.tiff", w = 8, h = 8, res = 300, 
     type="cairo", units = "in")

ggarrange(roc_scores_plot, empty_plot, roc_genes_plot, nrow = 3, 
          heights = c(0.7, 0.05, 1),
          labels = c("A", "", "B"))

dev.off()

## -----------------------------------------------------------------------------

## Saving DEG tables

# Subsetting
sva_results_sig <- lapply(sva_results, subset, padj < 0.05) %>% 
  lapply(as.data.frame)

# Rownames
for (i in 1:length(sva_results_sig)) {
  sva_results_sig[[i]][, "gene"] <- rownames(sva_results_sig[[i]])
}

# Writing file
writexl::write_xlsx(sva_results_sig, "results/tables/DEGs.xlsx")