# This script is used to identify genes whose expression are regulated by IFNG in cancer cells.
# We examined 48 human melanoma cell lines with or without IFNG treatment (RNA-seq data)

##••-- Library --••##
library(tidyverse)
library(Seurat)
library(anndata)
library(ggsci)
library(reticulate)
library(SeuratDisk)
library(limma)
library(ggthemes)
library(glue)
library(patchwork)
library(Matrix)
library(data.table)
library(rcna)
library(edgeR)
library(DESeq2)
library(ggpubr)
library(lme4)
library(lmerTest)
library(effectsize)
library(emmeans)
library(ggpubr)
use_python("/home/liubaoli/miniconda3/bin/")

##••-- Library --••##
bda = readr::read_tsv("/home/pauling/projects/01.melanoma/09.data.new/18.bulk.RNA.ifng/GSE154996_FPKM_Matrix.txt.gz")
bda.sample = tibble(sample = colnames(bda)[-c(1:6)]) %>%
  dplyr::mutate(label = ifelse(stringr::str_detect(sample,"KO"), "KO", "notKO")) %>%
  dplyr::mutate(group = ifelse(stringr::str_detect(sample,"Ctrl"), "Ctrl", "Expt")) 

bda <- bda[,-c(2:6)]
bda.pro <- bda[,-1]
bda.pro <- as.matrix(bda.pro)
rownames(bda.pro) <- bda$GeneId


#---- DE

bda.sample.filt <- bda.sample %>%
  dplyr::filter(label == "notKO")

counts1 = bda.pro[,bda.sample.filt$sample]

dge <- DGEList(counts=counts1)
dge <- calcNormFactors(dge)

keep <- filterByExpr(dge, group=bda.sample.filt$group)
dge <- dge[keep,]
normalized_counts <- cpm(dge, log=TRUE, prior.count=5) 

bda.sample.filt <- bda.sample.filt %>%
  dplyr::mutate(Sample = stringr::str_remove(sample, "\\.Ctrl.FPKM")) %>%
  dplyr::mutate(Sample = stringr::str_remove(Sample, "\\.Expt.FPKM")) %>%
  dplyr::mutate(Sample = stringr::str_remove(Sample, "\\.Ctrl.Counts")) %>%
  dplyr::mutate(Sample = stringr::str_remove(Sample, "\\.Expt.Counts"))

bda.sample.filt$Sample <- as.factor(bda.sample.filt$Sample)
bda.sample.filt$group <- as.factor(bda.sample.filt$group)

design <- model.matrix(~ Sample + group, bda.sample.filt)  # includes an intercept
colnames(design)
y <- voom(counts1, design, plot = F)

fit  <- lmFit(y, design)
fit  <- eBayes(fit)

ifng.up.gene <- topTable(fit, coef = "groupExpt", number = Inf) %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(logFC > 0.3 & adj.P.Val < 0.01) 

ifng.down.gene <- topTable(fit, coef = "groupExpt", number = Inf) %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(logFC < -0.3 & adj.P.Val < 0.01) 


## Boxplot
as.data.frame(t(bda.pro[c('MAP3K4',"AGO2", 'ZCCHC9', 'MAP2K6', 'ZMYND11',"ARID2"),])) %>%
  dplyr::mutate(
    group = bda.sample$group,
    sample = bda.sample$sample) %>%
  dplyr::mutate(Sample = stringr::str_remove(sample, "\\.Ctrl.FPKM")) %>%
  dplyr::mutate(Sample = stringr::str_remove(Sample, "\\.Expt.FPKM")) %>%
  dplyr::mutate(Sample = stringr::str_remove(Sample, "\\.Ctrl.Counts")) %>%
  dplyr::mutate(Sample = stringr::str_remove(Sample, "\\.Expt.Counts")) %>%
  tidyr::gather(key = "gene", value = "expr", -Sample, -sample, -group) %>%
  dplyr::filter(gene %in% ifng.down.gene$gene) %>%
  ggplot(aes(gene, log2(expr + 1))) +
  geom_boxplot(aes(color = group),outlier.size = -1) +
  geom_point(aes(color = group), position = position_jitterdodge(jitter.width = 0.1)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(color="black"),
    strip.background = element_blank(),
    #legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 11,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(color="black", size = 14),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1)
    #axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#B4EEB4","#3CB371")) +
  labs(x = "")
