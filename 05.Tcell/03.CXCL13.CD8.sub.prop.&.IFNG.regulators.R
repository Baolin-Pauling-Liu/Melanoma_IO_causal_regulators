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

##••-- Load data --••##
cxcl13.cd8.subsets <- readr::read_rds("/home/pauling/projects/01.melanoma/02.data/01.adata/07.cxcl13.tcell.rds.gz")


##••-- Marker genes --••##
library(presto)
subset.markers <- FindAllMarkers(cxcl13.cd8.subsets, only.pos = T)

subset.markers <- subset.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::filter(pct.1 > 0.15)


##••-- Boxplot for genes of interest --••##
tibble(
  expr = cxcl13.cd8.subsets@assays$RNA$data["IFNG",],   ### Replace this with genes of interest, e.g., T cell exhaustion markers
  sample = cxcl13.cd8.subsets$Sample,
  cluster = as.character(Idents(cxcl13.cd8.subsets))
) %>%
  dplyr::group_by(sample, cluster) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::summarise(n = mean(n), expr = mean(expr)) %>%
  ggplot(aes(cluster, expr)) +
  geom_boxplot(aes(color = cluster),outlier.size = -1) +
  geom_jitter(aes(color = cluster),width = 0.1, alpha = 0.7) +
  scale_color_manual(values = ggsci::pal_d3(palette = "category20c")(20)[c(7,1,2,3,8,4,6,9:20)]) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(linetype = "dashed"),
    panel.grid.minor = element_line(linetype = "dashed"),
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 13,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  )

##••-- Assotiation between the abundance of each CXCL13+ CD8 T cell subset and response to ICB --••##

#-- Defining functions: linear model: lm(frequency ~ response + covariates)
cal.es.pvale <- function(cell.prop, cell.type){
  
  beta.coefficients <- c()
  beta.coefficients2 <- c()
  beta.coefficients3 <- c()
  beta.coefficients4 <- c()
  beta.coefficients5 <- c()
  
  pvalues <- c()
  pvalues2 <- c()
  pvalues3 <- c()
  pvalues4 <- c()
  pvalues5 <- c()
  
  for (i in 1:length(cell.type)) {
    tmp.da <- cell.prop %>%
      dplyr::rename(celltype = cluster) %>%
      dplyr::filter(response %in% c("R","NR")) %>%
      dplyr::filter(celltype == cell.type[i]) %>%
      dplyr::mutate(age = as.double(age))
    
    tmp.model <- lm(prop ~ response + sex, tmp.da)
    tmp.model2 <- lm(prop ~ response + sex + age, tmp.da)
    tmp.model3 <- lm(prop ~ response + sex + Timepoint, tmp.da)
    tmp.model4 <- lm(prop ~ response + sex + Site, tmp.da)
    
    tmp.da.cutaneous <- tmp.da %>%
      dplyr::filter(met == 0 & melanome.type == "Cutaneous")
    
    tmp.model5 <- lm(prop ~ response + sex, tmp.da.cutaneous)
    
    model.sum <- summary(tmp.model)
    model.sum2 <- summary(tmp.model2)
    model.sum3 <- summary(tmp.model3)
    model.sum4 <- summary(tmp.model4)
    model.sum5 <- summary(tmp.model5)
    
    standardized_model <- standardize_parameters(tmp.model)
    standardized_model2 <- standardize_parameters(tmp.model2)
    standardized_model3 <- standardize_parameters(tmp.model3)
    standardized_model4 <- standardize_parameters(tmp.model4)
    standardized_model5 <- standardize_parameters(tmp.model5)
    
    beta.coefficients[i] <- standardized_model$Std_Coefficient[2]
    beta.coefficients2[i] <- standardized_model2$Std_Coefficient[2]
    beta.coefficients3[i] <- standardized_model3$Std_Coefficient[2]
    beta.coefficients4[i] <- standardized_model4$Std_Coefficient[2]
    beta.coefficients5[i] <- standardized_model5$Std_Coefficient[2]
    
    pvalues[i] <- model.sum$coefficients[2,4]
    pvalues2[i] <- model.sum2$coefficients[2,4]
    pvalues3[i] <- model.sum3$coefficients[2,4]
    pvalues4[i] <- model.sum4$coefficients[2,4]
    pvalues5[i] <- model.sum5$coefficients[2,4]
  }
  
  tibble(
    celltype = cell.type,
    beta.coefficient = beta.coefficients,
    beta.coefficient2 = beta.coefficients2,
    beta.coefficient3 = beta.coefficients3,
    beta.coefficient4 = beta.coefficients4,
    beta.coefficient5 = beta.coefficients5,
    pvalue = pvalues,
    fdr = p.adjust(pvalues, method = "fdr"),
    pvalue2 = pvalues2,
    pvalue3 = pvalues3,
    fdr3 = p.adjust(pvalues3, method = "fdr"),
    pvalue4 = pvalues4,
    pvalue5 = pvalues5
  ) %>%
    dplyr::arrange(pvalue)
}

#-- Load meta data
final.meta <- readr::read_rds("/home/pauling/projects/01.melanoma/09.data.new/02.sample.info/sample.clinial.meta.new.rds.gz")
major.meta <- readr::read_csv("/home/liubaoli/projects/01.melanoma/02.data/01.adata/01.round1.new.clustering.csv")

major.meta.filt <- major.meta %>%
  dplyr::filter(doublet == "FALSE")

final.meta.pro <- final.meta %>%
  dplyr::filter(response %in% c("R","NR")) %>%
  dplyr::mutate(Site = ifelse(site == "LN","LN","other")) %>%
  dplyr::mutate(Timepoint = ifelse(Timepoint == "pre", "pre", "post/on"))

#CD45 cell count per sample
sample.count <- major.meta.filt %>%
  dplyr::filter(!(new_clusters %in% c("Endothelial","Fibroblast","Malignant cell/melanocyte/epithelial cell"))) %>%
  dplyr::count(Sample) %>%
  dplyr::rename(sample = Sample) 


#Frequency calculation
tibble(
  sample = cxcl13.cd8.subsets$Sample,
  cluster = as.character(Idents(cxcl13.cd8.subsets))
) %>%
  dplyr::count(sample, cluster) %>%
  dplyr::inner_join(sample.count, by = "sample") %>%
  dplyr::mutate(prop = n.x/n.y) %>%
  dplyr::ungroup() %>%
  tidyr::complete(sample, cluster, fill = list(prop = 0)) %>%
  dplyr::inner_join(final.meta, by = "sample") -> pda

#Linear model
cxcl13.sub.diff <- cal.es.pvale(pda, unique(pda$cluster))

#plot
cxcl13.sub.diff %>%
  ggplot(aes(fct_reorder(celltype, beta.coefficient3), beta.coefficient3)) +
  geom_col(aes(fill = celltype), color = "black") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(linetype = "dashed"),
    panel.grid.minor = element_line(linetype = "dashed"),
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 13,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  coord_flip() +
  labs(y = "Effect size (standardized beta coefficient)", x = "") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_fill_manual(values = ggsci::pal_d3(palette = "category20c")(20)[c(7,1,2,3,8,6,4,9:20)])



#IFNG-low vs IFNG-hi CXCL13+ CD8 T cells
tibble(
  sample = as.character(cxcl13.cd8.subsets$Sample),
  cluster = as.character(Idents(cxcl13.cd8.subsets))
) %>%
  dplyr::filter(cluster %in% c("IFNG-hi","IFNG-low")) %>%
  dplyr::count(sample, cluster) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::mutate(prop = log((n+1)/(sum(n)-n+1))) %>%
  dplyr::ungroup() %>%
  tidyr::complete(sample, cluster, fill = list(prop = 0)) %>%
  dplyr::filter(total > 29) %>%
  dplyr::inner_join(final.meta.pro, by = "sample") %>%
  dplyr::filter(cluster == "IFNG-low") %>%
  ggplot(aes(response, prop)) +
  geom_boxplot(aes(color = response), outlier.size = -1) +
  geom_point(aes(color = response, shape = Timepoint, group = response), position = position_jitterdodge(jitter.width = 0.3), size = 3, alpha = 0.8) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linetype = "dashed"),
    panel.grid.minor = element_line(linetype = "dashed"),
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 13,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(x = "", y = "Log-ratio of IFNG-low to\nIFNG-high cell abundance") +
  scale_color_manual(values = c("#00868B","#CD2626")) +
  scale_shape_manual(values = c(20,17))



#------ Differential expression and causal genes: IFNG-hi vs IFNG-low cells ------------#

library(presto)
DimPlot(cxcl13.cd8.subsets, label = T)
ifng.subsets <- cxcl13.cd8.subsets[,Idents(cxcl13.cd8.subsets) %in% c("IFNG-hi","IFNG-low")]
ifng.markers <- FindMarkers(ifng.subsets, ident.1 = "IFNG-low")

ifng.markers <- ifng.markers %>%
  tibble::rownames_to_column(var = "gene")

ifng.markers <- ifng.markers %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(p_val_adj < 1e-5)

ifng.screen <- readr::read_tsv("/home/pauling/projects/01.melanoma/02.data/09.pertubation/07.tcell/01.IL2.ifng")

ifng.screen.cd8 <- ifng.screen %>%
  dplyr::filter(Cytokine == "IFNG") %>%
  dplyr::filter(CD4_or_CD8 == "CD8") %>%
  dplyr::filter(FDR < 0.05)

ifng.screen.cd8 <- ifng.screen.cd8 %>%
  dplyr::arrange(desc(abs(zscore))) %>%
  dplyr::group_by(Hit_Type, CRISPRa_or_i) %>%
  dplyr::mutate(rank = 1:n())

#-- Negative regulators of IFNG that are up-regulated in IFNG-low cells
nr.neg <- ifng.screen.cd8 %>%
  dplyr::arrange(desc(abs(zscore))) %>%
  dplyr::group_by(Hit_Type, CRISPRa_or_i) %>%
  dplyr::mutate(rank = 1:n()) %>%
  dplyr::filter(Hit == "TRUE") %>%
  dplyr::inner_join(ifng.markers, by = c("Gene" = "gene")) %>%
  dplyr::filter(p_val_adj < 1e-4 & abs(avg_log2FC) >= 0.15) %>%
  dplyr::ungroup() %>%
  dplyr::filter((avg_log2FC > 0 & Hit_Type == "Negative Hit")) %>%
  dplyr::distinct(Gene, avg_log2FC, p_val_adj) %>%
  dplyr::mutate(p_val_adj = ifelse(p_val_adj < 1e-20, 1e-20, p_val_adj)) %>%
  dplyr::arrange(desc(avg_log2FC))  
