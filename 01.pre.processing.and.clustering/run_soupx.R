
library(Seurat)
library(SoupX)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

sample <- args[1]

inp.path <- "/home/pauling/projects/01.melanoma/02.data/29.soupx/01.input"
out.path <- "/home/pauling/projects/01.melanoma/02.data/29.soupx/02.output"

sc = load10X(file.path(inp.path,sample))

obj <- Read10X_h5(file.path(inp.path,sample,"filtered_feature_bc_matrix.h5"))
obj <- CreateSeuratObject(obj)
obj <- NormalizeData(obj)
obj <- ScaleData(obj, do.center = F, do.scale = F)
obj <- FindVariableFeatures(obj)
obj <- RunPCA(obj, features = VariableFeatures(obj))
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, dims = 1:30, resolution = 0.5)

sc = setClusters(sc, obj@active.ident)
sc = autoEstCont(sc, doPlot = F, forceAccept=TRUE)
out = adjustCounts(sc)

print(paste0("SoupX: ", sample, "done"))

out %>% readr::write_rds(file.path(out.path,paste0(sample,".rds.gz")), compress = "gz")
