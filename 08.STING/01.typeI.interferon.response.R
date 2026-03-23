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
use_python("/home/liubaoli/miniconda3/bin/")


##••-- Load data --••##
mac.sce <- readr::read_rds("/home/pauling/projects/01.melanoma/02.data/01.adata/06.mac.sce.rds.gz")  # Here we loaded macrophage data, and can replace this with the obejct of each cell population
mac.bulk <- generate.pseudobulk(.obj = mac.pro, .cluster = levels(mac.pro))

##••-- Define function --••##
cal.cyto.sig <- function(cell_type, .sce, nbin = 24){
  cytokine.genes <- dic.cyto %>% 
    dplyr::filter(Celltype == cell_type) %>% 
    dplyr::rename(Gene = gene)

  #Identify genes specifically induced by IFNa1/IFNb or IFNg
  typei.genes <- cytokine.genes %>%
    dplyr::filter(Cytokine_Str %in% c("IFNa1","IFNb") & Avg_log2FC > 0 & FDR < 1e-5) 
  
  overlap.gene <- cytokine.genes %>%
    dplyr::filter(Cytokine_Str %in% c("IFNg") & Avg_log2FC > 0 & FDR < 1e-5) %>%
    dplyr::filter(Gene %in% typei.genes$Gene)
  
  cytokine.genes.pos <- cytokine.genes %>%
    dplyr::filter(!(Gene %in% c(overlap.gene$gene))) %>%
    dplyr::filter(Avg_log2FC > 0 & FDR < 1e-5) %>%
    dplyr::filter(Cytokine_Str %in% c("IFNa1","IFNb","IFNg"))

  #-- Here we calculate module scores for genes induced by IFNa1/IFNb and IFNg, respestively.
  cytokine.sig.pos <- cytokine.genes.pos %>%
    dplyr::select(Cytokine_Str, Gene) %>%
    dplyr::filter(Gene %in% rownames(.sce)) %>%
    dplyr::select(Cytokine_Str, Gene) %>%
    tidyr::nest(Gene) %>%
    dplyr::mutate(
      genes = purrr::map(
        .x = data,
        .f = function(.x){
          .x$Gene
        }
      )
    )
  
  
  cytokine.sig.pos <- cytokine.sig.pos %>%
    dplyr::mutate(
      n = purrr::map_dbl(
        .x = data,
        .f = function(.x){
          nrow(.x)
        }
      )
    ) %>%
    dplyr::filter(n > 3)
  
  .sce.pos <- AddModuleScore(.sce, features = cytokine.sig.pos$genes, nbin = nbin)
  
  cytokine.sig.score.pos <- .sce.pos@meta.data %>%
    dplyr::select_at(vars(contains("Cluster"))) 
  
  colnames(cytokine.sig.score.pos) <- cytokine.sig.pos$Cytokine_Str
  
  cytokine.sig.score.pos <- cytokine.sig.score.pos %>%
    dplyr::mutate(
      sample = colnames(.sce.pos),
      cellcounts = .sce.pos$cellcounts
    )
  
  cytokine.sig.score <- cytokine.sig.score.pos
  cytokine.sig.score$label <- cell_type
  return(cytokine.sig.score)
}


##••-- Cytokine response activity --••##
dic.cyto <- readr::read_rds("/home/pauling/projects/01.melanoma/02.data/10.cytokine/00.dic.rds.gz")
cyto.score <- cal.cyto.sig("Macrophage", .sce = mac.bulk)

cyto.score %>%
  readr::write_rds('/home/pauling/projects/01.melanoma/09.data.new/16.typeI/05.macrophage.typei.rds.gz', compress = "gz")
