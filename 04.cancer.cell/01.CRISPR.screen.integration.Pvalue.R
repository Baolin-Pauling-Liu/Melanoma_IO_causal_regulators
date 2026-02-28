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
library(RobustRankAggreg)
use_python("/home/liubaoli/miniconda3/bin/")

##••-- Load data --••##
screen.da.new <- readr::read_rds("projects/01.melanoma/09.data.new/00.screen/21.screen.data.new.rds.gz")

screen.categories <- screen.da.new %>%
  dplyr::distinct(annotation)

screen.da.new <- screen.da.new %>%
  dplyr::mutate(screen_id = paste0(group, annotation, model))

category.res <- list()
for (i in 1:8) {
  print(i)
  tmp.screen <- screen.da.new %>%
    dplyr::filter(annotation == screen.categories$annotation[i])
  
  uni.screen <- unique(tmp.screen$screen_id)
  
  rank_depl <- list()
  rank_enrich <- list()
  
  for (m in 1:length(uni.screen)) {
    ind.screen <- tmp.screen %>%
      dplyr::filter(screen_id == uni.screen[m])
    
    if(unique(ind.screen$rank_m) == "comb"){
      rank_depl[[m]] <- ind.screen %>%
        dplyr::arrange(rank) %>%
        dplyr::pull(gene) %>%
        unique()
      
      rank_enrich[[m]] <- ind.screen %>%
        dplyr::arrange(desc(rank)) %>%
        dplyr::pull(gene) %>%
        unique()
      
    }else if(unique(ind.screen$rank_m) == "neg_pos"){
      rank_depl[[m]] <- ind.screen %>%
        dplyr::filter(regulator == "Negative") %>%
        dplyr::arrange(rank) %>%
        dplyr::pull(gene) %>%
        unique()
      
      rank_enrich[[m]] <- ind.screen %>%
        dplyr::filter(regulator == "Positive") %>%
        dplyr::arrange(rank) %>%
        dplyr::pull(gene) %>%
        unique()
    }else if(unique(ind.screen$rank_m) == "neg_pos_cbind"){
      rank_depl[[m]] <- ind.screen %>%
        dplyr::arrange(`neg|rank`) %>%
        dplyr::pull(gene) %>%
        unique()
      
      rank_enrich[[m]] <- ind.screen %>%
        dplyr::arrange(`pos|rank`) %>%
        dplyr::pull(gene) %>%
        unique()
    }
  }
  
  N_total <- length(unique(unlist(rank_depl)))
  res_depl <- aggregateRanks(glist = rank_depl, N = N_total, method = "RRA")
  res_depl <- res_depl %>%
    dplyr::mutate(fdr = p.adjust(res_depl$Score, method = "fdr")) %>%
    dplyr::mutate(group = screen.categories$annotation[i]) %>%
    dplyr::mutate(regulator = "Negative")
  
  N_total <- length(unique(unlist(rank_enrich)))
  res_enrich <- aggregateRanks(glist = rank_enrich, N = N_total, method = "RRA")
  res_enrich <- res_enrich %>%
    dplyr::mutate(fdr = p.adjust(res_enrich$Score, method = "fdr")) %>%
    dplyr::mutate(group = screen.categories$annotation[i]) %>%
    dplyr::mutate(regulator = "Positive")
  
  category.res[[i]] <- res_depl %>%
    dplyr::full_join(res_enrich, by = "Name") %>%
    dplyr::rename(gene = Name)
}

save(category.res, file = "projects/01.melanoma/09.data.new/00.screen/screen.pvalue.new.rda")
