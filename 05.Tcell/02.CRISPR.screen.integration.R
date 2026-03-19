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

##••-- Load individual CRISPR screens --••##
ifng.da <- readxl::read_xlsx("projects/01.melanoma/09.data.new/00.screen/03.Tcell/IFNG.IL2.xlsx")
ifng.da = ifng.da %>%
  dplyr::filter(Screen_Version == "Primary") %>%
  dplyr::filter(CD4_or_CD8 == "CD8") %>%
  dplyr::mutate(regulator = ifelse(zscore < 0, "Negative", "Positive")) %>%
  dplyr::group_by(CRISPRa_or_i, regulator) %>%
  dplyr::arrange(desc(abs(zscore))) %>%
  dplyr::mutate(rank = 1:n()) %>%
  dplyr::rename(gene = Gene) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(screen = paste0("IFNG screen (",CRISPRa_or_i,")")) %>%
  dplyr::group_by(CRISPRa_or_i) %>%
  dplyr::mutate(size = n()) %>%
  dplyr::ungroup()

#-----------------------
RARS2 <- readxl::read_xlsx("projects/01.melanoma/09.data.new/00.screen/03.Tcell/RASA2 screen.xlsx")
RARS2 <- RARS2 %>%
  dplyr::group_by(Type) %>%
  dplyr::mutate(size = n()) %>%
  dplyr::ungroup() %>%
  #dplyr::mutate(rank = ifelse(pos.rank < neg.rank, pos.rank, neg.rank)) %>%
  dplyr::mutate(regulator = ifelse(pos.rank < neg.rank, "Negative", "Positive")) %>%
  dplyr::mutate(screen = paste0("T cell proliferation (", Type, ")")) %>%
  dplyr::rename(gene = id) 


#-------- T cell exhaustion
Tex <- readr::read_tsv("projects/01.melanoma/09.data.new/00.screen/03.Tcell/Tcell.exhaustion.cancer.cell.rank.zscore.tsv")
Tex <- Tex %>%
  dplyr::mutate(size = n()) %>%
  dplyr::arrange(desc(abs(z))) %>%
  dplyr::mutate(regulator = ifelse(z > 0, "Negative", "Positive")) %>%
  dplyr::group_by(regulator) %>%
  dplyr::mutate(rank = 1:n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(screen = paste0("T cell exhaustion")) %>%
  dplyr::rename(gene = Gene) 


#-------- T cell stimulation
tcell.stimulation <- readr::read_tsv("/home/pauling/projects/01.melanoma/02.data/09.pertubation/07.tcell/07.tcell.stimultation.prolif.replicates")

tcell.stimulation <- tcell.stimulation %>%
  dplyr::mutate(regulator = ifelse(`neg|rank` < `pos|rank`, "Positive", "Negative")) %>%
  dplyr::mutate(screen = "T cell proliferation (stimulation)") %>%
  dplyr::rename(gene = id) %>%
  #dplyr::mutate(rank = ifelse(`pos|rank` < `neg|rank`, `pos|rank`, `neg|rank`)) %>%
  dplyr::mutate(size = n())

#-------- T cell proliferation in vivo
gmap <- readr::read_rds("/home/pauling/projects/01.melanoma/02.data/97.gene/01.gene.map.rds.gz")
tcell.p.vivo <- readr::read_tsv("/home/pauling/projects/01.melanoma/09.data.new/00.screen/03.Tcell/ICR473")
tcell.p.vivo <- tcell.p.vivo %>%
  dplyr::mutate(size = n()) %>%
  dplyr::mutate(rank = ifelse(pos_rank < neg_rank, pos_rank, neg_rank)) %>%
  dplyr::mutate(regulator = ifelse(pos_rank < neg_rank, "Negative", "Positive")) %>%
  dplyr::mutate(screen = "T cell proliferation (in vivo)")

tcell.p.vivo <- tcell.p.vivo %>%
  dplyr::left_join(gmap, by = c("gene" = "gene.x")) %>%
  dplyr::mutate(gene = ifelse(is.na(gene.y), gene, gene.y))

#-------- T cell expansion (acute chronic intense)
all_sheets <- readxl::excel_sheets(path = "/home/pauling/projects/01.melanoma/09.data.new/00.screen/03.Tcell/Tcell.expansion.cancer.cell.xlsx")
read.xl <- function(p, st = "2", label = "intense"){
  tmp.da <- readxl::read_xlsx("/home/pauling/projects/01.melanoma/09.data.new/00.screen/03.Tcell/Tcell.expansion.cancer.cell.xlsx", sheet = st)
  colnames(tmp.da) <- unname(unlist(tmp.da[1,]))
  tmp.da[-1,] %>%
    dplyr::mutate(screen = paste0("T cell proliferation (", label, ")"))
}

da1 = read.xl(p = "/home/pauling/projects/01.melanoma/09.data.new/00.screen/03.Tcell/Tcell.expansion.cancer.cell.xlsx",
        st = "2", label = "intense")

da2 = read.xl(p = "/home/pauling/projects/01.melanoma/09.data.new/00.screen/03.Tcell/Tcell.expansion.cancer.cell.xlsx",
              st = "3", label = "acute")

da3 = read.xl(p = "/home/pauling/projects/01.melanoma/09.data.new/00.screen/03.Tcell/Tcell.expansion.cancer.cell.xlsx",
              st = "4", label = "chronic")

tcell.prolif.3contexs = bind_rows(da1, da2, da3) %>%
  dplyr::mutate(`pos|rank` = as.double(`pos|rank`), `neg|rank` = as.double(`neg|rank`)) %>%
  dplyr::mutate(regulator = ifelse(`pos|rank` < `neg|rank`, "Negative", "Positive")) %>%
  #dplyr::mutate(rank = as.double(ifelse(`pos|rank` < `neg|rank`, `pos|rank`, `neg|rank`))) %>%
  dplyr::group_by(screen) %>%
  dplyr::mutate(size = n()) %>%
  dplyr::ungroup() %>%
  dplyr::rename(gene = id) %>%
  dplyr::left_join(gmap, by = c("gene" = "gene.x")) %>%
  dplyr::mutate(gene = ifelse(is.na(gene.y), gene, gene.y))


#---- hongbo T cell R

hb <- readr::read_tsv("/home/pauling/projects/01.melanoma/02.data/09.pertubation/07.tcell/06.hongbo.long.live.T")
hb <- hb %>%
  dplyr::arrange(desc(abs(`log2 ratio (TIL/input)`))) %>%
  dplyr::mutate(regulator = ifelse(`log2 ratio (TIL/input)` > 0, "Negative", "Positive")) %>%
  dplyr::mutate(rank = 1:n()) %>%
  dplyr::rename(gene = `Gene Symbol`) %>%
  dplyr::left_join(gmap, by = c("gene" = "gene.x")) %>%
  dplyr::mutate(gene = ifelse(is.na(gene.y), gene, gene.y)) %>%
  dplyr::mutate(screen = "T cell expansion (melanoma cell co-culture)") %>%
  dplyr::mutate(size = nrow(.))


#---- RRA and pvalues --  negative regulators

rank_depl <- list()
rank_depl[[1]] <- ifng.da %>%
  dplyr::filter(CRISPRa_or_i == "CRISPRi") %>%
  dplyr::arrange(zscore) %>%
  dplyr::pull(gene)

rank_depl[[2]] <- ifng.da %>%
  dplyr::filter(CRISPRa_or_i == "CRISPRa") %>%
  dplyr::arrange(zscore) %>%
  dplyr::pull(gene)

RARS2.rank <- RARS2 %>%
  dplyr::arrange(Type, pos.rank) %>%
  dplyr::select(Type, gene) %>%
  tidyr::nest(-Type) %>%
  dplyr::mutate(data = purrr::map(.x = data, .f = function(.x){.x$gene}))

rank_depl[3:8] <- RARS2.rank$data

tcell.prolif.3contexs.pro <- tcell.prolif.3contexs %>%
  dplyr::arrange(screen, `pos|rank`) %>%
  dplyr::select(screen, gene) %>%
  tidyr::nest(-screen) %>%
  dplyr::mutate(data = purrr::map(.x = data, .f = function(.x){.x$gene}))

rank_depl[9:11] <- tcell.prolif.3contexs.pro$data
tcell.stimulation <- tcell.stimulation %>%
  dplyr::arrange(`pos|rank`)

rank_depl[[12]] <-  tcell.stimulation$gene

tcell.p.vivo <- tcell.p.vivo %>%
  dplyr::arrange(pos_rank)

rank_depl[[13]] <-  tcell.stimulation$gene

hb <- hb %>%
  dplyr::arrange(desc(`log2 ratio (TIL/input)`))

rank_depl[[14]] <-  hb$gene

N_total <- length(unique(unlist(rank_depl)))
res_depl <- aggregateRanks(glist = rank_depl, N = N_total, method = "RRA")


#---- RRA and pvalues -- Positive regulators

rank_enr <- list()
rank_enr[[1]] <- ifng.da %>%
  dplyr::filter(CRISPRa_or_i == "CRISPRi") %>%
  dplyr::arrange(desc(zscore)) %>%
  dplyr::pull(gene)

rank_enr[[2]] <- ifng.da %>%
  dplyr::filter(CRISPRa_or_i == "CRISPRa") %>%
  dplyr::arrange(desc(zscore)) %>%
  dplyr::pull(gene)

RARS2.rank <- RARS2 %>%
  dplyr::arrange(Type, neg.rank) %>%
  dplyr::select(Type, gene) %>%
  tidyr::nest(-Type) %>%
  dplyr::mutate(data = purrr::map(.x = data, .f = function(.x){.x$gene}))

rank_enr[3:8] <- RARS2.rank$data

tcell.prolif.3contexs.pro <- tcell.prolif.3contexs %>%
  dplyr::arrange(screen, `neg|rank`) %>%
  dplyr::select(screen, gene) %>%
  tidyr::nest(-screen) %>%
  dplyr::mutate(data = purrr::map(.x = data, .f = function(.x){.x$gene}))

rank_enr[9:11] <- tcell.prolif.3contexs.pro$data
tcell.stimulation <- tcell.stimulation %>%
  dplyr::arrange(`neg|rank`)

rank_enr[[12]] <-  tcell.stimulation$gene

tcell.p.vivo <- tcell.p.vivo %>%
  dplyr::arrange(neg_rank)

rank_enr[[13]] <-  tcell.stimulation$gene

hb <- hb %>%
  dplyr::arrange(`log2 ratio (TIL/input)`)

rank_enr[[14]] <-  hb$gene

N_total <- length(unique(unlist(rank_enr)))
rank_enr <- aggregateRanks(glist = rank_enr, N = N_total, method = "RRA")


save(rank_enr, res_depl, file = "projects/01.melanoma/09.data.new/06.Tcell/05.tcell.screen.rda")



#--- Identify R-neg and NR-neg genes
load("/home/pauling/projects/01.melanoma/09.data.new/00.screen/03.Tcell/CXCL13_CD8_DE.rda")  # Load CXCL13 CD8 T cell DE data --> Tcell.diff

res_depl <- res_depl %>%
  dplyr::mutate(fdr = p.adjust(res_depl$Score, method = "fdr")) %>%
  dplyr::mutate(regulator = "Negative")

neg.gene <- Tcell.screen %>%  ### This is an additional filtering in addition to FDR < 0.05. We request ranks of a given negative regulator < 150 in at least 3 CRISPR screens.   
  dplyr::filter(rank < 150) %>%
  dplyr::count(gene) %>%
  dplyr::filter(n > 2)


r.neg = Tcell.diff %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(logFC > 0) %>%
  dplyr::filter(gene %in% res_depl$Name[res_depl$fdr < 0.05] & gene %in% neg.gene$gene)

nr.neg = Tcell.diff %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(logFC < 0) %>%
  dplyr::filter(gene %in% res_depl$Name[res_depl$fdr < 0.05] & gene %in% neg.gene$gene)

res_depl %>%
  dplyr::filter(Name %in% neg.gene$gene) %>%
  dplyr::filter(fdr < 0.05) %>%
  dplyr::rename(gene = Name) %>%
  View(.)

res_depl %>%
  dplyr::filter(Name %in% neg.gene$gene) %>%
  dplyr::filter(fdr < 0.05) %>%
  dplyr::rename(gene = Name) %>%
  dplyr::mutate(rank = 1:nrow(.)) %>%
  ggplot(aes(fct_reorder(gene, -Score), -log10(fdr))) +
  geom_col(aes(fill = -log10(fdr))) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color="black", angle = 90, size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(linetype = "dashed", color = "grey80"),
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 13,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    axis.text.y = element_text(color="black", size = 14),
    axis.text.x = element_text(color="black")
  ) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_fill_gradientn(colours = c(ggsci::pal_material(palette = "purple")(10)))

Tcell.screen %>%
  dplyr::count(screen) %>%
  ggplot(aes(fct_reorder(screen, n), n)) +
  geom_col(fill = c("#CD6889", "#8968CD")[1]) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color="black", angle = 90, size = 11),
    panel.grid.major = element_line(linetype = "dashed", color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 13,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    axis.text.y = element_text(color="black", size = 14),
    axis.text.x = element_text(color="black")
  ) +
  coord_flip() +
  labs(x = "", y = "") +
  scale_fill_gradientn(colours = c(ggsci::pal_material(palette = "purple")(10))) +
  labs(y = "Number of genes")




