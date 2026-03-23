
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

Tcell.prop.pro <- readr::read_rds("/home/pauling/projects/01.melanoma/09.data.new/00.screen/03.Tcell/Tcell.prop.rds.gz")

comb.typei.pro <- comb.typei %>%
  dplyr::inner_join(Tcell.prop.pro, by = "sample") %>%
  dplyr::select(-sample) %>%
  as.matrix()


##••-- Correlation with frequencies of CXCL13+ CD4 and CD8 T cells within immune cells across samples --••##
library(Hmisc)
result <- rcorr(comb.typei.pro, type = "pearson")
cor_matrix <- result$r   # Correlation coefficients
pval_matrix <- result$P  

cor.da <- as.data.frame(cor_matrix) %>%
  tibble::rownames_to_column(var = "m1") %>%
  tidyr::gather(key = "m2", value = "cor", -m1) %>%
  as.tibble()

pval.da <- as.data.frame(pval_matrix) %>%
  tibble::rownames_to_column(var = "m1") %>%
  tidyr::gather(key = "m2", value = "pvalue", -m1) %>%
  as.tibble()

cor.da %>%
  dplyr::mutate(pvalue = pval.da$pvalue) %>%
  dplyr::filter(m1 %in% colnames(Tcell.prop.pro) & !(m2 %in% colnames(Tcell.prop.pro))) %>%
  dplyr::mutate(sig = ifelse(pvalue < 0.05, ifelse(pvalue < 0.01, ifelse(pvalue < 0.001,"***","**"), "*"),NA)) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::filter(stringr::str_detect(m2, "IFNa1") | stringr::str_detect(m2, "IFNb")) %>
  ggplot(aes(factor(m1, levels = c("CD8.CXCL13+", "CD4.CXCL13")), m2)) +
  geom_tile(aes(fill = cor), color = "white") +
  geom_text(aes(label = sig), color = "white") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_line(linetype = "dashed", colour = "grey80"),
        axis.title   = element_text(size = 13, colour = "black"),
        axis.text    = element_text(size = 13, colour = "black"),
        legend.title = element_text(size = 13),
        legend.text  = element_text(size = 11),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "") +
  scale_fill_gradientn(colours = c(ggsci::pal_material("teal")(10))) +
  scale_color_manual(values = "black", na.value = "white")


##••-- Quantile regression test --••##
quant.res <- list()
for (i in 1:(ncol(comb.typei.pro)-3)) {
  tmp.da <- comb.typei.pro[,c(i, (ncol(comb.typei.pro)-2):ncol(comb.typei.pro))] %>%
    tidyr::drop_na()
  
  colnames(tmp.da) <- c("score","response","Timepoint","sex")
  coefs <- c()
  pvalues <- c()
  m = 0
  for (p in c(0.05,0.08,0.1,0.15,0.85,0.9,0.92,0.95)) {
    m = m + 1
    tmp.sum <- summary(rq(score ~ response + Timepoint + sex, data = tmp.da, tau = p), se="ker")
    coefs[m] <- tmp.sum$coefficients[2,1]
    pvalues[m] <- tmp.sum$coefficients[2,4]
  }
  quant.res[[i]] <- tibble(coef = coefs, pvalue = pvalues) %>%
    dplyr::mutate(feature = colnames(comb.typei.pro)[i]) %>%
    dplyr::mutate(quantile = c(0.05,0.08,0.1,0.15,0.85,0.9,0.92,0.95))
}

bind_rows(quant.res) %>%
  dplyr::filter(!stringr::str_detect(feature, "IFNg")) %>%  #  IFNa1 or IFNb
  dplyr::mutate(sig = ifelse(pvalue < 0.05, ifelse(pvalue < 0.01, ifelse(pvalue < 0.001,"***","**"), "*"),NA)) %>%
  #tidyr::separate(feature, c("celltype","cytokine"), sep = "_") %>%
  ggplot(aes(as.character(quantile), feature)) +
  geom_tile(aes(fill = coef), color = "white") +
  geom_text(aes(label = sig), color = "white") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_line(linetype = "dashed", colour = "grey80"),
        axis.title   = element_text(size = 13, colour = "black"),
        axis.text    = element_text(size = 13, colour = "black"),
        legend.title = element_text(size = 13),
        legend.text  = element_text(size = 11),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "") +
  #scale_fill_gradientn(colours = c(rev(ggsci::pal_material(palette = "blue")(10)[-c(2,4,6,8)]),"white",ggsci::pal_material(palette = "purple")(10))) +
  scale_fill_gradientn(colours = rev(gplots::redblue(100)), na.value = "grey80")  
  scale_color_manual(values = "black", na.value = "white")
