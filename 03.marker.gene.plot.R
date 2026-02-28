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
library(presto)
use_python("/home/liubaoli/miniconda3/bin/")


##••-- Load data and plot --••##
tcell.sce <- readr::read_rds("projects/01.melanoma/02.data/01.adata/07.tcell.rds.gz") # T cells as an example
markers <- FindAllMarkers(tcell.sce, only.pos = T)

use.markers <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice(1:5)

receptor.expr = DotPlot(tcell.sce, features = rev(unique(use.markers$gene)))
receptor.expr$data %>%
  #dplyr::mutate(avg.exp.scaled = ifelse(avg.exp.scaled > 1.5, 1.5, avg.exp.scaled)) %>%
  ggplot(aes(features.plot, factor(id, levels = unique(use.markers$cluster)))) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_distiller(palette = "RdBu") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_line(linetype = "dashed", color = "grey"),
        legend.position = "bottom",
        axis.title = element_text(size = 13,color="black"),
        axis.text = element_text(size = 11,color="black"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black", angle = 45, hjust = 1, face = "italic")
  ) +
  labs(x = "", y = "")
