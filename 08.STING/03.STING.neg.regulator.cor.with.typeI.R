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

##••-- Load data: type I interferon response activity in each cell type across samples --••##
file.path <- "/home/pauling/projects/01.melanoma/09.data.new/16.typeI/"
files <- list.files(file.path)

typei.score <- tibble(file = files) %>%
  dplyr::mutate(
    data = purrr::map(
      .x = file,
      .f = function(.x){
        tmp.da <- readr::read_rds(file.path(file.path, .x)) %>%
          dplyr::select(-cellcounts)
        
        tmp.da1 <- tmp.da %>%
          dplyr::select(sample, label) %>%
          dplyr::mutate(label = ifelse(label == "MigDC", "mregDC", label)) %>%
          dplyr::mutate(label = ifelse(label == "ILC", "NK", label))
        
        tmp.da1 <- tmp.da %>%
          dplyr::select(sample, label) %>%
          dplyr::mutate(label = ifelse(label == "MigDC", "mregDC", label)) %>%
          dplyr::mutate(label = ifelse(label == "ILC", "NK", label)) %>%
          dplyr::mutate(label = ifelse(label == "Plasma", "plasma", label)) %>%
          dplyr::mutate(label = ifelse(label == "cDC1", "DC1", label)) %>%
          dplyr::mutate(label = ifelse(label == "cDC2", "DC2", label)) %>%
          dplyr::mutate(label = ifelse(label == "B cell", "Bcell", label)) 
        
        tmp.da2 <- tmp.da %>%
          dplyr::select(-sample, -label)
        
        colnames(tmp.da2) <- paste0(unique(tmp.da1$label), "_", colnames(tmp.da2))
        cbind(tmp.da2, tmp.da1) %>%
          dplyr::select(-label)
      }
    )
  )

##••-- Data process: type I interferon response activity in each cell type across samples --••##
comb.typei <- typei.score$data[[1]]
for (i in 2:nrow(typei.score)) {
  comb.typei <- comb.typei %>%
    dplyr::full_join(typei.score$data[[i]], by = "sample")
}

#--- STING negative regulators
load("projects/01.melanoma/09.data.new/14.inter.cell/02.all.module.rda")
comb.module.ma <- comb.module[, stringr::str_detect(colnames(comb.module), "M") | colnames(comb.module) == "sample"]
comb.module.ma <- comb.module.ma %>%
  dplyr::inner_join(comb.typei, by = "sample") %>%
  dplyr::select(-sample) %>%
  as.matrix()

result <- rcorr(comb.module.ma, type = "pearson")
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

all.da <- cor.da %>%
  dplyr::mutate(fdr = pval.da$pvalue) %>%
  dplyr::filter(stringr::str_detect(m1, "IFNa1") | stringr::str_detect(m1, "IFNb")) %>%
  dplyr::filter(stringr::str_detect(m2, "_M")) %>%
  dplyr::mutate(label = stringr::str_sub(m1, 1, -6)) %>%
  dplyr::filter(label == stringr::str_sub(m2, 1, stringr::str_length(label))) %>%
  dplyr::inner_join(all.module, by = c("m2" = "module"))


sting.neg <- readr::read_rds("/home/pauling/projects/01.melanoma/09.data.new/16.typeI/STING.neg.regulators.rds.gz")
sting.neg <- sting.gene %>%
  dplyr::filter(gene %in% all.da$gene) %>%
  dplyr::distinct(gene)

cor.da %>%
  dplyr::mutate(fdr = pval.da$pvalue) %>%
  dplyr::filter(stringr::str_detect(m1, "IFNa1") | stringr::str_detect(m1, "IFNb")) %>%
  dplyr::filter(stringr::str_detect(m2, "_M")) %>%
  dplyr::mutate(label = stringr::str_sub(m1, 1, -6)) %>%
  dplyr::filter(label == stringr::str_sub(m2, 1, stringr::str_length(label))) %>%
  dplyr::filter(cor < -0.2 & fdr < 0.05) %>%
  dplyr::inner_join(all.module, by = c("m2" = "module")) %>%
  dplyr::filter(gene %in% sting.neg$gene) %>%
  dplyr::mutate(cytokine = ifelse(stringr::str_detect(m1, "IFNa1"), "IFNa1","IFNb")) %>%
  dplyr::distinct(gene, m1, cor, fdr) %>%
  dplyr::mutate(FDR = ifelse(fdr < 1e-7, 1e-7, fdr)) %>%
  ggplot(aes(factor(gene, levels = sting.neg$gene), m1)) +
  geom_point(aes(fill = cor, size = -log10(FDR)), shape = 21, stroke = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = 13,color="black"),
        axis.text = element_text(size = 11,color="black"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = "", y = "") +
  scale_fill_gradientn(colours = ggsci::pal_material(palette = "purple")(10))



##••-- Permutation test --••##
erm.n <- c()
res <- list()
for (i in 1:1000) {
  sample.gene <- all.module %>%
    dplyr::distinct(gene) %>%
    dplyr::sample_n(123) ### There are 123 negative regulators of STING identified by corresponding CRISPR screens
  
  all.da %>%
    dplyr::filter(cor < -0.2 & fdr < 0.05) %>%
    dplyr::filter(gene %in% sample.gene$gene) %>%
    dplyr::distinct(gene, m1) %>%
    dplyr::count(gene) %>%
    dplyr::filter(n > 1) %>%
    nrow(.) -> perm.n[i]
  
  all.da %>%
    dplyr::filter(cor < -0.2 & fdr < 0.05) %>%
    dplyr::filter(gene %in% sample.gene$gene) %>%
    dplyr::distinct(gene, m1, label) %>%
    dplyr::count(label) -> res[[i]]
  
  print(i)
}

sting.da <- all.da %>%
  dplyr::filter(cor < -0.2 & fdr < 0.05) %>%
  dplyr::filter(gene %in% sting.neg$gene) %>%
  dplyr::distinct(gene, m1, label) %>%
  dplyr::count(label)

bind_rows(res) %>%
  dplyr::inner_join(sting.da, by = "label") %>%
  dplyr::arrange(label) %>%
  dplyr::arrange(desc(n.x)) %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(size = n()) %>%
  dplyr::filter(n.x > n.y) %>%
  dplyr::mutate(outlier = n()) %>%
  dplyr::mutate(p = outlier/size) %>%
  dplyr::distinct(label, p) %>%
  dplyr::mutate(cytokine = ifelse(stringr::str_detect(label, "_"), "IFNb","IFNa1")) %>%
  dplyr::mutate(label = stringr::str_remove(label, "_")) %>%
  dplyr::mutate(sig = ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001,"***","**"), "*"),NA)) %>%
  ggplot(aes(cytokine, label)) +
  theme_bw() +
  geom_tile(aes(fill = -log10(p)), color = "white") +
  geom_text(aes(label = sig), color = "white") +
  theme(panel.grid.major = element_blank(),
        strip.text = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.position = "bottom",
        axis.title = element_text(size = 13,color="black"),
        axis.text = element_text(size = 11,color="black"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black")
  ) +
  labs(x = "", y = "") +
  scale_fill_gradientn(colours = c(ggsci::pal_material(palette = "deep-purple")(10)))
