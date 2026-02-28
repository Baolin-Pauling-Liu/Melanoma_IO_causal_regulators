
##••-- load data --••##
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

##••-- DE functions --••##
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
  v <- voomWithQualityWeights(dge, mm, plot = TRUE)
  fit <- lmFit(v, mm)
  
  contr <- makeContrasts(response01.g1 - response02.g2, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  
  tmp <- eBayes(tmp)
  #return(v)
  topTable(tmp, sort.by = "P", n = Inf)
}

##••-- load data --••##
final.meta <- readr::read_rds("/home/pauling/projects/01.melanoma/09.data.new/02.sample.info/sample.clinial.meta.new.rds.gz")

malignant.cellid <- readr::read_csv("/home/pauling/projects/01.melanoma/09.data.new/03.gene.module/01.malignant.cellid.csv")
malignant.cellid %>%
  dplyr::count(Channel) %>%
  dplyr::arrange(n) %>%
  dplyr::filter(Channel %in% final.meta$sample) %>%
  dplyr::filter(n > 99) -> filt.sample

mean.expr <- readr::read_rds("/home/pauling/projects/01.melanoma/09.data.new/01.gene.expr.objects/01.malignant.cell.pseudo.bulk.rds.gz")
get.norm.matr <- function(input.ma){
  dge <- DGEList(counts=input.ma)
  dge <- calcNormFactors(dge)
  
  normalized_counts <- cpm(dge, log=TRUE, prior.count=5) 
  return(normalized_counts)
}

mean.expr.norm <- get.norm.matr(t(mean.expr))

mean.expr <- t(mean.expr)

##••-- DE analsyis --••##
final.meta.pro.sub = final.meta %>%
  dplyr::filter(sample %in% filt.sample$Channel)


diff2.limma = limma.fun(mean.expr, final.meta.pro.sub, 
                        input.group1 = "R",
                        input.group2 = "NR", 
                        filter = T)

##••-- DE analsyis --••##
causal.gene <- readr::read_tsv("/home/pauling/projects/01.melanoma/09.data.new/18.causal.gene.new.Dec29.csv/cancer.cell.causal.gene.tsv")

causal.gene.limma <- diff2.limma %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% causal.gene$gene)

causal.gene.limma$fdr <- p.adjust(causal.gene.limma$P.Value, method = "fdr")

hh <- causal.gene.limma %>%
  dplyr::filter(fdr <= 0.1)

##••-- DE gene boxplot --••##
design <- model.matrix(~ response + sex + Timepoint, data=final.meta.pro.sub)
covars <- model.matrix(~ 0 + Timepoint + sex, data = final.meta.pro.sub)

# 2) 
design_keep <- model.matrix(~ response, data = final.meta.pro.sub)

lcpm_adj <- limma::removeBatchEffect(
  x          = mean.expr.norm[,final.meta.pro.sub$sample],
  design     = design_keep,    
  covariates = covars          
)


as.data.frame(lcpm_adj[unique(hh$gene),]) %>%
  tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(key = "sample", value = "expr", -gene) %>%
  dplyr::inner_join(final.meta.pro.sub[,c("sample","response","Timepoint")], by = "sample") %>%
  dplyr::filter(gene != "NFIB") %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(expr = (expr - median(expr))/(mad(expr))) %>%
  ggplot(aes(factor(gene, levels = unique(gene.rank$gene)), expr)) +
  geom_boxplot(aes(color = response), outlier.size = -1) +
  geom_point(aes(color = response, shape = Timepoint, group = response), position = position_jitterdodge(jitter.width = 0.2), size = 1.3, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_line(linetype = "dashed", color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 13,color="black"),
        axis.text = element_text(size = 11,color="black"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  labs(
    x = "",
    y = "Expression"
  ) +
  scale_color_manual(values = c("#00868B", "#CD5555")) 
  
