
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
library(Hmisc)
use_python("/home/liubaoli/miniconda3/bin/")

##••-- Load data --••##
#save(comb.module, all.module, file = "projects/01.melanoma/09.data.new/14.inter.cell/02.all.module.rda")
load("projects/01.melanoma/09.data.new/14.inter.cell/02.all.module.rda")

##••-- Process module data --••##
comb.module.pro <- comb.module %>%
  dplyr::select(-sample) %>%
  as.matrix()

rownames(comb.module.pro) <- comb.module$sample

result <- rcorr(comb.module.pro, type = "pearson")
cor_matrix <- result$r   # Correlation coefficients
pval_matrix <- result$P  

fdr_matrix <- apply(pval_matrix, 2, p.adjust, method = "fdr")  ### Adjusted P value

cor.da <- as.data.frame(cor_matrix) %>%
  tibble::rownames_to_column(var = "m1") %>%
  tidyr::gather(key = "m2", value = "cor", -m1) %>%
  as.tibble()

pval.da <- as.data.frame(fdr_matrix) %>%
  tibble::rownames_to_column(var = "m1") %>%
  tidyr::gather(key = "m2", value = "pvalue", -m1) %>%
  as.tibble()


##••-- Cytokines and chmokines: correlation with the frequency of CXCL13+ CD8 T cells within immune cells --••##
cytokines = tibble(module = colnames(comb.module)[-1]) %>%
  dplyr::filter(!stringr::str_detect(module, "M"))

chemo.da <- all.module %>%
  dplyr::filter(stringr::str_detect(gene, "CXCL") | stringr::str_detect(gene, "CCL") | stringr::str_detect(gene, "CCR") | stringr::str_detect(gene, "CXCR")) %>%
  dplyr::distinct(gene, module)

#--- cytokines and chmokines --- correlation CXCL13 CD8
ccs = c(cytokines$module[1:203], chemo.da$module, "CD8.CXCL13+")

result <- rcorr(comb.module.pro[,unique(ccs)], type = "pearson")
cor_matrix <- result$r   # Correlation coefficients
pval_matrix <- result$P  

fdr_matrix <- apply(pval_matrix, 2, p.adjust, method = "fdr")

cor.da <- as.data.frame(cor_matrix) %>%
  tibble::rownames_to_column(var = "m1") %>%
  tidyr::gather(key = "m2", value = "cor", -m1) %>%
  as.tibble()

pval.da <- as.data.frame(fdr_matrix) %>%
  tibble::rownames_to_column(var = "m1") %>%
  tidyr::gather(key = "m2", value = "pvalue", -m1) %>%
  as.tibble()

cor.da %>%
  dplyr::mutate(fdr = pval.da$pvalue) %>%
  dplyr::filter(m1 == "CD8.CXCL13+") %>%    ### Replce this with cell populations of interest
  dplyr::arrange(cor) %>%
  dplyr::filter(!is.na(fdr)) %>%
  dplyr::mutate(fdr = ifelse(fdr == 0, 1e-15, fdr)) %>%
  dplyr::mutate(name = ifelse(cor > 0.6 | cor < -0.4, m2, NA)) %>%
  ggplot(aes(1:nrow(.), cor)) +
  geom_point(aes(size = -log10(fdr), color = cor)) +
  theme_bw() +
  geom_text_repel(aes(label = name), box.padding = 0.8, max.overlaps = 100) +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "grey80"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(color="black"),
    strip.background = element_blank(),
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 13,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(color="black", size = 11),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_gradientn(colours = c(rev(ggsci::pal_material(palette = "blue")(10)),"white",ggsci::pal_material(palette = "purple")(10)))


