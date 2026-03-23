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
library(doParallel)
library(ComplexHeatmap)
use_python("/home/liubaoli/miniconda3/bin/")


##••-- Load data (ICB outcome-associated cell populations and causal genes) --••##
load("projects/01.melanoma/09.data.new/19.causal.feature.heatmap.da.rda")

##••-- Clustering and heatmap --••##
mat = Heatmap(t(final.causal.matr[,-1]), clustering_distance_rows = "pearson", clustering_distance_columns = "pearson", clustering_method_columns = "complete")

causal.feature = colnames(final.causal.matr[,-1])[row_order(mat)]
sample.order = final.causal.matr$sample[column_order(mat)]

final.causal.matr %>%
  tidyr::gather(key = "feature", value = "score", -sample) %>%
  dplyr::inner_join(final.meta.pro, by = "sample") %>%
  dplyr::mutate(score = ifelse(score > 3, 3, score)) %>%
  dplyr::mutate(score = ifelse(score < -3, -3, score)) %>%
  ggplot(aes(factor(sample, levels = rev(sample.order)), factor(feature, levels = rev(causal.feature)))) +
  geom_tile(aes(fill = score), color = "white") +
  #scale_fill_distiller(palette = "RdBu") +
  scale_fill_gradientn(colours = rev(gplots::redblue(200)), na.value = "grey80") +
  facet_grid(. ~ response, scales = "free", space = "free") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "grey80"),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(color="black"),
    strip.background = element_blank(),
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 11,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(color="black", size = 11),
    axis.text.x = element_text(color="black", size = 0)
  )

