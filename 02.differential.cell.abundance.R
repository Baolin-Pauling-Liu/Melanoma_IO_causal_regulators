##••-- Library --••##
library(tidyverse)
library(ggsci)
library(effectsize)
library(emmeans)
library(ggpubr)
library(sccomp)

##••-- Load data --••##
major.meta <- readr::read_csv("/home/liubaoli/projects/01.melanoma/02.data/01.adata/01.round1.new.clustering.csv")
final.meta <- readr::read_rds("/home/pauling/projects/01.melanoma/09.data.new/02.sample.info/sample.clinial.meta.new.rds.gz")

final.meta.pro <- final.meta %>%
  dplyr::filter(response %in% c("R","NR")) %>%
  dplyr::mutate(Site = ifelse(site == "LN","LN","other")) %>%
  dplyr::mutate(Timepoint = ifelse(Timepoint == "pre", "pre", "post/on"))

##••-- sccomp (non-malignant major cell types) --••##
sample.cell.type.count <- major.meta %>% 
  dplyr::filter(!(new_clusters %in% c("Malignant cell/melanocyte/epithelial cell"))) %>%
  dplyr::count(Sample, new_clusters) %>%
  tidyr::complete(Sample, new_clusters, fill = list(n = 0))

sccomp_input = sample.cell.type.count %>%
  dplyr::select(Sample, new_clusters, n) %>%
  dplyr::rename(sample = Sample) %>%
  dplyr::inner_join(final.meta.pro[,c("sample", "response","sex","Timepoint","Site")], by = "sample")

sccomp_input$sample <- as.factor(sccomp_input$sample)
sccomp_input$Timepoint <- as.factor(sccomp_input$Timepoint)
sccomp_input$response <- as.factor(sccomp_input$response)
sccomp_input$new_clusters <- as.factor(sccomp_input$new_clusters)

estimate <- sccomp_estimate(
  sccomp_input,
  ~ response + sex + Timepoint,
  ~1,
  "sample",
  "new_clusters",
  "n",
  cores = 5,
  verbose = F
)

estimate_clean <- estimate |> sccomp_remove_outliers(cores = 5, verbose = F)
results_sub        <- estimate_clean |> sccomp_test()

pred_df <- sccomp_predict(
  fit                   = estimate_clean,
  formula_composition   = ~ response + sex + timepoint,
  #new_data              = new_data,
  summary_instead_of_draws = TRUE 
)

df_plot <- results_sub %>%
  dplyr::filter(parameter == "responseR") %>%         
  dplyr::rename(cluster = label,
         estimate = c_effect,
         lower    = c_lower,
         upper    = c_upper) %>%
  dplyr::arrange(-log10(c_FDR) * estimate) %>%                        
  dplyr::mutate(cluster = factor(cluster, levels = cluster)) %>%
  dplyr::mutate(c_FDR = ifelse(c_FDR == 0, 1e-4, c_FDR))


ggplot(df_plot, aes(x = estimate, y = cluster)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = -log10(c_FDR)),
                 height = 0, lwd = 0.9) +
  geom_point(size = 4, aes(color = -log10(c_FDR))) +
  theme_bw() +
  theme(
    #legend.position = "none",
    panel.grid.major = element_line(linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 15,color="black"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")#, angle = 45, hjust = 1)
  ) +
  labs(
    x = "Composition effect (posterior mean)",
    y = "Cell type",
    title = ""
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_material(palette = "purple") +
  scale_fill_material(palette = "purple")

##••-- sccomp (immune cell subsets) --••##
cell.sub <- readr::read_rds("/home/pauling/projects/01.melanoma/02.data/00.meta/16.cell.subset.rds.gz")

sccomp_input = cell.sub %>%
  dplyr::count(Sample, label) %>%
  dplyr::rename(sample = Sample) %>%
  tidyr::complete(sample, label, fill = list(n = 0)) %>%
  dplyr::inner_join(final.meta.pro[,c("sample", "response","sex","Timepoint","Site")], by = "sample")

sccomp.fun <- function(input){
  estimate_all_sub <- sccomp_estimate(
    input,
    ~ response + sex + Timepoint,
    ~ 1,
    "sample",
    "label",
    "n",
    cores = 5,
    verbose = F
  )
  
  tmp.estimate_all_sub_clean <- estimate_all_sub |> sccomp_remove_outliers(cores = 5, verbose = F)
  tmp.results_all_sub        <- tmp.estimate_all_sub_clean |> sccomp_test()
  
  print("50% done")
  
  tmp.estimate_sub_co <- sccomp_estimate(
    input,
    ~ response + sex + Timepoint + Site,
    ~ 1,
    "sample",
    "label",
    "n",
    cores = 5,
    verbose = F
  )
  
  tmp.estimate_sub_co_clean <- tmp.estimate_sub_co |> sccomp_remove_outliers(cores = 5, verbose = F)
  tmp.resuts_sub_co        <- tmp.estimate_sub_co_clean |> sccomp_test()
  
  list(tmp.results_all_sub, tmp.resuts_sub_co)
}

diff.cell.res = sccomp.fun(sccomp_input)

df_plot <- diff.cell.res[[1]] %>%
  dplyr::filter(parameter == "responseR") %>%         
  dplyr::rename(cluster = label,
                estimate = c_effect,
                lower    = c_lower,
                upper    = c_upper) %>%
  dplyr::arrange(log10(c_FDR) * estimate) %>%                        
  dplyr::mutate(cluster = factor(cluster, levels = cluster)) %>%
  dplyr::mutate(c_FDR = ifelse(c_FDR == 0, 1e-4, c_FDR))

df_plot %>%
  ggplot(aes(x = estimate, y = cluster)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = -log10(c_FDR)),
                 height = 0, lwd = 0.9) +
  geom_point(size = 4, aes(color = -log10(c_FDR))) +
  theme_bw() +
  theme(
    #legend.position = "none",
    panel.grid.major = element_line(linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 15,color="black"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_material(palette = "purple") +
  scale_fill_material(palette = "purple") +
  coord_flip()

