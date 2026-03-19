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

##••-- Load pseudo-bulk data --••##
bulk.da <- readr::read_rds("/home/pauling/projects/01.melanoma/09.data.new/00.screen/03.Tcell/CXCL13_CD8_pseudobulk.rds.gz")

use.sample <- tibble(sample = bulk.da$Sample,
       cn = bulk.da$cellcounts) %>%
  dplyr::inner_join(final.meta.pro, by = "sample") %>%
  dplyr::filter(cn > 19)

bulk.da <- bulk.da[,use.sample$sample]

#--- Diff expr ---#
limma.fun <- function(input.ma, input.mda, input.group1, input.group2, gp1.label = "01.g1", gp2.label = "02.g2", filter = F){
  meta.da <- tibble(sample = colnames(input.ma)) %>%
    dplyr::left_join(input.mda, by = "sample") 
  
  meta.da <- meta.da %>%
    dplyr::filter(response %in% c(input.group1, input.group2)) %>%
    dplyr::mutate(response = ifelse(response %in% input.group1, gp1.label, response)) %>%
    dplyr::mutate(response = ifelse(response %in% input.group2, gp2.label, response))
  
  counts1 = input.ma[,meta.da$sample]
  
  dge <- DGEList(counts=counts1)
  dge <- calcNormFactors(dge)
  
  if (filter){
    keep <- filterByExpr(dge, group=meta.da$response)
    dge <- dge[keep,]
  }
  
  meta.da$response <- factor(meta.da$response)
  meta.da$melanome.type <- factor(meta.da$melanome.type)
  meta.da$Timepoint <- factor(meta.da$Timepoint)
  meta.da$sex <- factor(meta.da$sex)
  
  mm <- model.matrix(~ 0 + response + Timepoint + sex, meta.da)
  #y <- voom(counts1, mm, plot = F)
  
  #y <- cpm(dge, log=TRUE, prior.count=1)
  v <- voom(dge, mm)
  fit <- lmFit(v, mm, method = "robust", maxit = 100)
  
  contr <- makeContrasts(response01.g1 - response02.g2, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  
  tmp <- eBayes(tmp)
  topTable(tmp, sort.by = "P", n = Inf)
}

get.norm.matr <- function(input.ma){
  dge <- DGEList(counts=input.ma)
  dge <- calcNormFactors(dge)
  
  normalized_counts <- cpm(dge, log=TRUE, prior.count=5) 
  return(normalized_counts)
}

tcell.bulk <- bulk.da@assays$RNA$counts

Tcell.diff = limma.fun(tcell.bulk, final.meta.pro, 
                        input.group1 = "R",
                        input.group2 = "NR", 
                        filter = T)

save(Tcell.diff, file = "/home/pauling/projects/01.melanoma/09.data.new/00.screen/03.Tcell/CXCL13_CD8_DE.rda")
