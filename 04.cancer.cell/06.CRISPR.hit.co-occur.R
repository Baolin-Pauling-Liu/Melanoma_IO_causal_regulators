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

##••-- Load dada --••##
screen.da <- readr::read_rds("projects/01.melanoma/09.data.new/00.screen/20.screen.data.rds.gz")

screen.da <- screen.da.new %>%
  dplyr::mutate(rank = ifelse(rank_m == "neg_pos_cbind", ifelse(regulator == "Negative", `neg|rank`, `pos|rank`), rank))

neg.top.genes <- screen.da %>%
  dplyr::filter(rank < 150 & regulator == "Negative") %>%
  #dplyr::arrange(rank) %>%
  dplyr::count(gene) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::filter(n > 6)

##••-- Generate matrix --••##
screen.matr <- screen.da %>%
  dplyr::filter(regulator == "Negative") %>%
  dplyr::filter(gene %in% neg.top.genes$gene) %>%
  dplyr::select(gene, group, rank, model, size) %>%
  dplyr::mutate(score = ifelse(rank < 150, 1, 0)) %>%
  dplyr::select(-size) %>%
  dplyr::mutate(group = paste0(group, model)) %>%
  tidyr::complete(gene, group, fill = list(score = 0))

screen.matr.pro <- screen.matr[,c(1,2,5)] %>%
  dplyr::distinct(gene, group, .keep_all = T) %>%
  tidyr::spread(key = group, value = score)

screen.matr <- as.matrix(screen.matr.pro[,-1])
rownames(screen.matr) <- screen.matr.pro$gene


##••-- Co-occurrence --••##
library(discover)
events = discover.matrix(screen.matr)
sim.gene.co = discover::pairwise.discover.test(events[neg.top.genes$gene,], alternative = "greater")

for (r in 1:nrow(sim.gene.co$q.values)) {
  for (col in 1:ncol(sim.gene.co$q.values)) {
    if(is.na(sim.gene.co$q.values[r,col])){
      sim.gene.co$q.values[r,col] = sim.gene.co$q.values[col,r]
    }
  }
}


as.data.frame(sim.gene.co$q.values) %>%
  tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(key = "gene2", value = "fdr", -gene)  %>%
  dplyr::filter(fdr < 0.05) %>%
  dplyr::arrange(gene, fdr) %>%
  readr::write_tsv("projects/01.melanoma/09.data.new/07.cancer.cell/37.gene.co.occur.all.tsv")

aa = as.data.frame(sim.gene.co$q.values) %>%
  tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(key = "gene2", value = "fdr", -gene)  %>%
  dplyr::filter(fdr < 0.05) %>%
  dplyr::arrange(gene, fdr) %>%
  dplyr::filter(gene %in% genes_of_interest) #### Input genes of interest, e.g., NR-neg genes


##••-- Plot --••##
as.data.frame(sim.gene.co$q.values) %>%
  tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(key = "gene2", value = "fdr", -gene) %>%
  dplyr::filter(dr <= 0.05 & gene2 %in% aa$gene2 & gene %in% aa$gene) %>%
  dplyr::filter(!is.na(fdr)) %>%
  dplyr::mutate(sig = ifelse(fdr <= 0.05, ifelse(fdr < 0.01, ifelse(fdr < 0.001,"***","**"), "*"),NA)) %>%
  dplyr::mutate(sig = ifelse(fdr > 0.05 & fdr < 0.11, "`", sig)) %>%
  ggplot(aes(factor(gene2, levels = unique(aa$gene2)), factor(gene, levels = rev(unique(aa$gene))))) +
  geom_tile(aes(fill = -log10(fdr)), color = "white") +
  geom_text(aes(label = sig), color = "white") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(linetype = "dashed", color = "grey80"),
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 13,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = "", y = "") +
  scale_fill_gradientn(colours = c(ggsci::pal_material(palette = "teal")),  na.value = "white")
