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
sobj.filt <- readr::read_rds("/home/pauling/projects/01.melanoma/02.data/01.adata/14.Bcell.soupx.rds.gz") # Here we loaded B cell data, and can replace this with the obejct of each cell population

mgenes <- readr::read_csv("/home/liubaoli/projects/01.melanoma/02.data/08.hotspot/04.B.cell/Bcell.gene.module.csv") # Co-expressed gene modules in B cells identified by Hotspot
mgenes.pro <- mgenes %>%
  dplyr::filter(!stringr::str_detect(Module, "-1")) %>%
  dplyr::filter(Module != -1) %>%
  tidyr::nest(-Module) %>%
  dplyr::mutate(
    data = purrr::map(
      .x = data,
      .f = function(.x){
        .x$featurekey
      }
    )
  )


##••-- Defining functions --••##
generate.pseudobulk <- function(.obj, .cluster){
  tmp.obj <- .obj[,Idents(.obj) %in% .cluster]
  tmp.bulk <- AggregateExpression(tmp.obj, assays = "RNA", return.seurat = T, group.by = c("Sample"))
  cell.count <- table(tmp.obj$Sample) %>%
    as.data.frame.array()
  
  tmp.bulk$cellcounts = cell.count[colnames(tmp.bulk),1]
  tmp.bulk <- NormalizeData(tmp.bulk)
  
  bulk.matr <- as.matrix(tmp.bulk@assays$RNA$counts)
  
  
  get.norm.matr <- function(input.ma){
    dge <- DGEList(counts=input.ma)
    dge <- calcNormFactors(dge)
    
    normalized_counts <- cpm(dge, log=TRUE, prior.count=5) 
    return(normalized_counts)
  }
  
  bulk.matr.norm <- get.norm.matr(bulk.matr)
  tmp.bulk@assays$RNA$data <- bulk.matr.norm
  
  tmp.bulk <- ScaleData(tmp.bulk)
  return(tmp.bulk)
}

cal.module.sig <- function(cell_type, .sce, mgene.da, nbin = 24){
  
  tmp.mgene <- mgene.da
  
  .sce.pos <- AddModuleScore(.sce, features = tmp.mgene$data, nbin = nbin)
  
  module.score <- .sce.pos@meta.data %>%
    dplyr::select_at(vars(contains("Cluster"))) 
  
  colnames(module.score) <- paste0(cell_type,"_M",tmp.mgene$Module)
  
  module.score <- module.score %>%
    dplyr::mutate(
      sample = colnames(.sce.pos),
      cellcounts = .sce.pos$cellcounts
    )
  return(module.score)
}


##••-- Pseudo-bulk data --••##
bulk.da <- generate.pseudobulk(.obj = sobj.filt, .cluster = levels(sobj.filt)[-5])


##••-- Module score --••##
Bcell.module.score <- cal.module.sig("Bcell", .sce = bulk.da, mgenes.pro)


##••-- Cytokine response --••##
cal.cyto.sig <- function(cell_type, .sce, nbin = 24){
  cytokine.genes <- dic.cyto %>% 
    dplyr::filter(Celltype == cell_type) %>% 
    dplyr::rename(Gene = gene)
  
  cytokine.genes.pos <- cytokine.genes %>%
    dplyr::filter(Avg_log2FC > 0 & FDR < 1e-5)
  
  cytokine.genes.neg <- cytokine.genes %>%
    dplyr::filter(Avg_log2FC < 0 & FDR < 1e-5)
  
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
  
  cytokine.sig.neg <- cytokine.genes.neg %>%
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
    dplyr::filter(n > 2)
  
  cytokine.sig.neg <- cytokine.sig.neg %>%
    dplyr::mutate(
      n = purrr::map_dbl(
        .x = data,
        .f = function(.x){
          nrow(.x)
        }
      )
    ) %>%
    dplyr::filter(n > 2)
  
  .sce.pos <- AddModuleScore(.sce, features = cytokine.sig.pos$genes, nbin = nbin)
  .sce.neg <- AddModuleScore(.sce, features = cytokine.sig.neg$genes, nbin = nbin)
  
  cytokine.sig.score.pos <- .sce.pos@meta.data %>%
    dplyr::select_at(vars(contains("Cluster"))) 
  
  colnames(cytokine.sig.score.pos) <- cytokine.sig.pos$Cytokine_Str
  
  cytokine.sig.score.pos <- cytokine.sig.score.pos %>%
    dplyr::mutate(
      sample = colnames(.sce.pos),
      cellcounts = .sce.pos$cellcounts
    )
  
  cytokine.sig.score.neg <- .sce.neg@meta.data %>%
    dplyr::select_at(vars(contains("Cluster"))) 
  
  
  colnames(cytokine.sig.score.neg) <- cytokine.sig.neg$Cytokine_Str
  
  cytokine.sig.score.neg <- cytokine.sig.score.neg %>%
    dplyr::mutate(
      sample = colnames(.sce.neg)
    )
  
  
  overlap.cyto <- intersect(colnames(cytokine.sig.score.pos), colnames(cytokine.sig.score.neg))
  overlap.cyto <- overlap.cyto[1:(length(overlap.cyto)-1)]
  cytokine.sig.score <- cytokine.sig.score.pos
  cytokine.sig.score[,overlap.cyto] <- cytokine.sig.score[,overlap.cyto] - cytokine.sig.score.neg[,overlap.cyto]
  cytokine.sig.score$label <- cell_type
  return(cytokine.sig.score)
}

dic.cyto <- readr::read_rds("/home/pauling/projects/01.melanoma/02.data/10.cytokine/00.dic.rds.gz")
Bcell.cyto.score <- cal.cyto.sig("B cell", .sce = bulk.da)

save(Bcell.cyto.score, Bcell.module.score, file = "/home/pauling/projects/01.melanoma/09.data.new/09.co.variation/Bcell.module.cytokine.rda")
