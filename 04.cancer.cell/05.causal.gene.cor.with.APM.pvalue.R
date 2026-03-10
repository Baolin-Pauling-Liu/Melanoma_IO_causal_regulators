
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

##••-- Meta data --••##
final.meta <- readr::read_rds("/home/pauling/projects/01.melanoma/09.data.new/02.sample.info/sample.clinial.meta.new.rds.gz")

final.meta.pro <- final.meta %>%
  dplyr::filter(response %in% c("R","NR")) %>%
  dplyr::mutate(Timepoint = ifelse(Timepoint == "pre", "pre", "post/on"))


##••-- Load data --••##
##••-- Refine MHC I and II modules for each tumor; some genes may be not expressed or do not correlate with other MHC genes due to certain reasons
##••-- Including these genes in the analysis will introduce bias in identifying causal genes that directly or indirectly suppress MHC gene modules
inp.dir = "/home/liubaoli/projects/01.melanoma/02.data/11.gene.program/02.cancer.cell/57.APM/"
samples = list.files(inp.dir)

MHC_I_genes = c("HLA-A","HLA-B","HLA-C","HLA-E","B2M","TAP1","PSMB9")
MHC_II_genes = c('HLA-DRA', 'HLA-DRB1', 'CD74', 'HLA-DRB5', 'HLA-DPA1', 'HLA-DQA1', 'HLA-DPB1', 'HLA-DQB1', 'HLA-DMA')

MHC_gene = tibble(
  gene = c(MHC_I_genes, MHC_II_genes),
  group = c(rep("MHC_I",length(MHC_I_genes)), rep("MHC_II",length(MHC_II_genes)))
)

samples = samples[-123]
module.list = list()
for (i in 1:length(samples)) {
  module.list[[i]] = readr::read_csv(file.path("/home/liubaoli/projects/01.melanoma/02.data/11.gene.program/02.cancer.cell/57.APM/",samples[i],paste0(samples[i],".module.APM.csv"))) %>%
    dplyr::mutate(sample = samples[i]) 
}

module.list.pro <- module.list
for (i in 1:length(module.list)) {
  m.count = module.list[[i]] %>%
    dplyr::filter(Module != -1) %>%
    dplyr::count(Module) %>%
    dplyr::filter(n > 3)
  
  if(nrow(m.count) == 2){
    filt.mhc = module.list[[i]] %>%
      dplyr::inner_join(MHC_gene,  by = c("featurekey" = "gene")) %>%
      dplyr::filter(Module != -1) %>%
      dplyr::count(Module, group) %>%
      dplyr::group_by(Module) %>%
      dplyr::mutate(frac = n/sum(n)) %>%
      dplyr::filter(frac > 0.74) %>%
      dplyr::ungroup()
    
    if(nrow(filt.mhc) > 0){
      module.list[[i]] = module.list[[i]] %>%
        dplyr::filter(Module != -1) %>%
        dplyr::inner_join(MHC_gene,  by = c("featurekey" = "gene")) %>%
        dplyr::inner_join(filt.mhc, by = "Module") %>%
        dplyr::filter(group.x == group.y) %>%
        dplyr::select(featurekey, sample, group.x) %>%
        dplyr::rename(group = group.x)
    }else if(nrow(filt.mhc) == 0){
      filt.mhc = module.list[[i]] %>%
        dplyr::inner_join(MHC_gene,  by = c("featurekey" = "gene")) %>%
        dplyr::filter(Module != -1) %>%
        dplyr::count(Module, group) %>%
        dplyr::group_by(Module) %>%
        dplyr::arrange(desc(n)) %>%
        dplyr::slice(1:2) %>%
        dplyr::filter(n > 3)
      
      module.list[[i]] = module.list[[i]] %>%
        dplyr::filter(Module != -1) %>%
        dplyr::inner_join(MHC_gene,  by = c("featurekey" = "gene")) %>%
        dplyr::inner_join(filt.mhc, by = "Module") %>%
        dplyr::filter(group.x == group.y) %>%
        dplyr::select(featurekey, sample, group.x) %>%
        dplyr::rename(group = group.x)
    }
  }else if(nrow(m.count) == 1){
    
    print(samples[i])
    
    filt.m = module.list[[i]] %>%
      dplyr::inner_join(MHC_gene,  by = c("featurekey" = "gene")) %>%
      dplyr::filter(Module != -1) %>%
      dplyr::count(Module, group) %>%
      dplyr::filter(n > 3)
    
    module.list[[i]] = module.list[[i]] %>%
      dplyr::filter(Module != -1) %>%
      dplyr::filter(Module %in% filt.m$Module) %>%
      dplyr::inner_join(MHC_gene,  by = c("featurekey" = "gene")) %>%
      dplyr::filter(group %in% filt.m$group) %>%
      dplyr::select(featurekey, sample, group)
  }
}

##••-- MHC I and II module summary --••##
bind_rows(module.list) %>%
  dplyr::count(sample, group) %>%
  dplyr::left_join(final.meta.pro, by = "sample") %>%
  dplyr::filter(response %in% c("R","NR")) %>%
  ggplot(aes(sample, group)) +
  geom_tile(aes(fill = "1"), color = "white", lwd = 0.5, fill = c("#00868B")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(linetype = "dashed", color = "grey80"),
    #legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 11,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 90, vjust = 0.5)
  ) +
  labs(
    x = "", y = ""
  ) +
  facet_grid(. ~ factor(response, levels = c("NR","R")), scales = "free", space = "free")


##••-- Causal gene correlation with MHC I and II modules --••##

sample.module = bind_rows(module.list) %>%
  dplyr::select(-Module) %>%
  dplyr::mutate(sm = paste0(sample,"_", featurekey)) %>%
  dplyr::select(-sample)

cor.list = list()
for (i in 1:length(samples)) {
  cor.list[[i]] = readr::read_csv(file.path("/home/liubaoli/projects/01.melanoma/02.data/11.gene.program/02.cancer.cell/57.APM/",samples[i],paste0(samples[i],".causal.cor.csv"))) %>%
    dplyr::mutate(sample = samples[i]) 
}

for (i in 1:length(cor.list)) {
  cor.list[[i]] = cor.list[[i]] %>%
    tidyr::gather(key = "gene", value = "cor", -featurekey, -sample) %>%
    dplyr::rename(gn = featurekey) %>%
    dplyr::mutate(sm = paste0(sample,"_", gene)) %>%
    dplyr::inner_join(sample.module, by = c("sm"))
}

##••-- MHC I and II module gene weights --••##
gene.weights = bind_rows(cor.list) %>%
  dplyr::filter(gn %in% MHC_gene$gene & gene %in% MHC_gene$gene) %>%
  dplyr::filter(gn != gene) %>%
  dplyr::inner_join(MHC_gene, by = c("gn" = "gene")) %>%
  dplyr::filter(group.x == group.y) %>%
  dplyr::select(-group.x) %>%
  dplyr::rename(group = group.y) %>%
  tidyr::nest(-group, -sample)

mhc.weights <- list()
for (i in 1:nrow(gene.weights)) {
  mhc.weights[[i]] <- gene.weights$data[[i]] %>%
    dplyr::group_by(gn) %>%
    dplyr::summarise(weight = mean(cor)) %>%
    dplyr::mutate(sample = gene.weights$sample[i], group = gene.weights$group[i])
}

mhc.weights <- bind_rows(mhc.weights) %>%
  dplyr::mutate(MHC_ID = paste0(sample,"_" ,group,"_",gn))

##••-- Pvalues --••##

acat <- function(p, w = NULL, lower = 1e-300, upper = 1 - 1e-15) {
  p <- as.numeric(p)
  if (anyNA(p)) stop("p has NA")
  p <- pmin(pmax(p, lower), upper)
  
  if (is.null(w)) w <- rep(1/length(p), length(p)) else {
    w <- as.numeric(w); stopifnot(length(w)==length(p), all(w>=0))
    w <- w / sum(w)
  }
  
  cotpi <- function(x) 1 / tanpi(x)  # 若无 tanpi，可用 1/tan(pi*x)
  T <- sum(w * cotpi(p))
  p_acat <- 0.5 - atan(T) / pi
  max(min(p_acat, 1), 0)
}

causal_gene_sta = bind_rows(cor.list) %>%
  dplyr::filter(!(gn %in% apm.module$featurekey)) %>%
  dplyr::mutate(pvalue = pnorm(-cor, lower.tail = FALSE)) %>%
  dplyr::mutate(pvalue2 = pnorm(cor, lower.tail = FALSE)) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::mutate(ID = paste(sample, group, gene, sep = "_"), Group = group, Sample = sample) %>%
  tidyr::nest(-sample, -group)

sample_gene = list()
for (i in 1:nrow(causal_gene_sta)) {
  sample_gene[[i]] = causal_gene_sta$data[[i]] %>%
    dplyr::mutate(MHC_ID = paste0(Sample,"_" ,Group, "_", gene)) %>%
    dplyr::inner_join(mhc.weights[,c(2,4,5)], by = "MHC_ID") %>%
    tidyr::nest(-Sample, -Group, -gn)
  
  sample_gene[[i]] <- sample_gene[[i]] %>%
    dplyr::mutate(
      pvalue = purrr::map_dbl(
        .x = data,
        .f = function(.x){
          acat(.x$pvalue, .x$weight)
        }
      )
    ) %>%
    dplyr::mutate(
      pvalue2 = purrr::map_dbl(
        .x = data,
        .f = function(.x){
          acat(.x$pvalue2, .x$weight)
        }
      )
    ) %>%
    dplyr::mutate(
      cor = purrr::map_dbl(
        .x = data,
        .f = function(.x){
          mean(.x$cor)
        }
      )
    ) %>%
    dplyr::mutate(
      num_pos = purrr::map_dbl(
        .x = data,
        .f = function(.x){
          length(.x$cor[.x$cor > 0])
        }
      )
    ) %>%
    dplyr::mutate(
      num_neg = purrr::map_dbl(
        .x = data,
        .f = function(.x){
          length(.x$cor[.x$cor < 0])
        }
      )
    )

  sample_gene[[i]]$FDR = p.adjust(sample_gene[[i]]$pvalue, method = "fdr")
  sample_gene[[i]]$FDR2 = p.adjust(sample_gene[[i]]$pvalue2, method = "fdr")
}

apm.sample <- bind_rows(module.list) %>%
  dplyr::count(sample, group) %>%
  dplyr::inner_join(final.meta.pro, by = "sample") %>%
  dplyr::filter(response %in% c("R","NR")) %>%
  dplyr::distinct(sample)

final.meta.pro.filt <- final.meta.pro %>%
  dplyr::filter(sample %in% apm.sample$sample)

causal.neg.immune.regulators <- readr::read_rds("projects/01.melanoma/09.data.new/00.screen/negative_regulators.rda")

neg.causal = bind_rows(sample_gene) %>%
  dplyr::rename(sample = Sample, group = Group) %>%
  dplyr::right_join(final.meta.pro.filt[,c("sample","response")], by = "sample") %>%
  #dplyr::right_join(final.meta[,c("sample")], by = "sample") %>%
  dplyr::mutate(sample = as.factor(sample)) %>%
  dplyr::filter(FDR < 0.01 & cor < -2 & num_pos < 2) 

neg.causal.pro <- neg.causal %>%
  dplyr::distinct(sample, gn) %>%
  dplyr::count(gn) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::filter(gn %in% causal.neg.immune.regulators) %>%
  dplyr::filter(n > 9)

# Causal genes exhibiting inverse correlation with APM genes in each sample
neg.causal %>%
  dplyr::filter(gn %in% neg.causal.pro$gn) %>%
  dplyr::filter(sample %in% final.meta$sample) %>%
  dplyr::group_by(gn, response, sample) %>%
  dplyr::mutate(num = n()) %>%
  dplyr::mutate(label = ifelse(num == 2, "Both", NA)) %>%
  dplyr::ungroup() %>%
  dplyr::select(gn, sample, group,response, num, label, cor) %>%
  dplyr::arrange(cor) %>%
  dplyr::distinct(gn, sample, num, response, label, .keep_all = T) %>%
  dplyr::mutate(label = ifelse(is.na(label), group, label)) %>%
  dplyr::filter(gn %in% neg.causal.pro$gn) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(gn) %>%
  dplyr::mutate(total = n()) %>%
  ggplot(aes(fct_reorder(sample, n, .desc = T), factor(gn, levels = rev(unique(gene.rank$gene))))) +
  geom_tile(aes(fill = label), color = "white",lwd = 0.2) +
  #geom_tile(aes(fill = label), color = "white", lwd = 0.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(linetype = "dashed", color = "grey80"),
    #legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 11,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 90, vjust = 0.5)
  ) +
  labs(
    x = "", y = ""
  ) +
  facet_grid(. ~ factor(response, levels = c("NR","R")), scales = "free", space = "free") +
  scale_fill_manual(values = c("#F08080", "#FFA500", "#F5DEB3")) 


#-- CRISPR screen data
pda <- readr::read_rds("projects/01.melanoma/09.data.new/00.screen/screen_sum_gene.rda")
Ï
pda$annotation <- factor(pda$annotation, 
                         levels = c("T cell killing", "ICB regulators", "ICB regulators (vs.output)",
                                    "Tumor cell killing by NK cells",
                                    "Regulators of MHC class I proteins",
                                    "STING regulation",
                                    "IFNG/TNF.induced.cell.death",
                                    "RNA sensing",
                                    "Phagocytosis by macrophages",
                                    "tumor growth regulators"))

pda %>%
  dplyr::filter(annotation != "ICB regulators (vs.output)") %>%
  dplyr::filter(!is.na(regulator)) %>%
  tidyr::complete(gene, annotation, fill = list(score = 0)) %>%
  dplyr::filter(annotation == "T cell killing") %>%
  dplyr::arrange(rank) -> gene.rank

pda %>%
  dplyr::group_by(annotation, regulator) %>%
  dplyr::ungroup() %>%
  ggplot(aes(-log10(fdr), factor(gene, levels = rev(unique(gene.rank$gene))))) +
  geom_col(aes(fill = regulator), position = position_dodge2()) +
  theme_bw() +
  facet_grid(. ~ annotation, scales = "free") +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "grey80"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(color="black"),
    strip.background = element_blank(),
    #legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 11,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(color="black", size = 14),
    axis.text.x = element_text(color="black")
    #axis.text.x = element_text(color="black")
  ) +
  labs(x = "CRISPR-Based Functional Impact Score
", y = "") +
  scale_fill_manual(values = c("#00868B", "#CD2626"))
