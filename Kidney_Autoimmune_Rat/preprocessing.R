library(Seurat)
library(tidyverse)
options(stringsAsFactors = F)

load("F:/biotech/projects/rat_annoed/seurat/outputs/seurat_preped.rda")

cluster <- readxl::read_xlsx("data/T&NK细分亚群.xlsx")
metadt <- sr_preped@meta.data

metadt$ident1 <- ""
cluster <- as.data.frame(cluster)
rownames(cluster) <- cluster$Cell
head(cluster)
length(cluster$Cell)
cellid <- intersect(rownames(metadt), cluster$Cell)
length(cellid)
length(cluster$Cell)
cluster <- cluster[cellid, ]
metadt[rownames(cluster), "ident1"] <- cluster$Cluster
table(metadt$ident1)
View(metadt)

id0 <- as.list(metadt[, c("manual.ident", "CD4", "CD8", "ident1"), drop = F] )

# ident2
metadt$ident2 <- ""
id2 <- vector(mode = "character", length = length(metadt$ident2))

l <- which(id0$CD8 == "CD8T-MAIT" & id0$manual.ident == "T" & id0$ident1 == "CD8T")
id2[l] <- "MAIT"
table(id2)

l <- which(id0$ident1 == "CD8T" & id0$CD8 != "" &id0$CD8 != "CD8T-MAIT")
id2[l] <- "classical-CD8T"
table(id2)

l <- which(id0$ident1 %in% c("ProliferativeT", "NK") & id0$manual.ident == "T")
id2[l] <- id0$ident1[l]
table(id2)

l <- which(id0$ident1 == "CD4T" & id0$CD4 != "")
id2[l] <- id0$ident1[l]
table(id2)

l <- which(id2 == "")
id2[l] <- id0$manual.ident[l]
table(id2)

l <- which(id2 == "T")
id2[l] <- "Unclassical-T"
table(id2)

metadt$ident2 <- id2

# ident1
id1 <- vector(mode = "character", length = length(metadt$ident1))

l <- which(id0$CD8 == "CD8T-MAIT" & id0$manual.ident == "T" & id0$ident1 == "CD8T")
id1[l] <- "MAIT"
table(id1)

l <- which(id0$ident1 == "CD8T" & id0$CD8 != "" &id0$CD8 != "CD8T-MAIT")
id1[l] <- id0$CD8[l]
table(id1)

l <- which(id1 %in% c("CD8T-other", "Proliferative-CD8T"))
id1[l] <- "CD8T-Other"
table(id1)


l <- which(id0$ident1 %in% c("ProliferativeT", "NK") & id0$manual.ident == "T")
id1[l] <- id0$ident1[l]
table(id1)

l <- which(id0$ident1 == "CD4T" & id0$CD4 != "")
id1[l] <- id0$CD4[l]
table(id1)

l <- which(id1 == "CD4-TEM_TH1like")
id1[l] <- "CD4-TH1"
table(id1)

l <- which(id1 == "CD4-TEMRA_TEFF")
id1[l] <- "CD4-TEM"
table(id1)

l <- which(id1 == "Proliferative-CD4T")
id1[l] <- "ProliferativeT"
table(id1)

l <- which(id1 == "")
id1[l] <- id0$manual.ident[l]
table(id1)

l <- which(id1 == "T")
id1[l] <- "Unclassical-T"
table(id1)

metadt$ident1 <- id1

# macrophage
macro <- readxl::read_xlsx("data/monocytic-reclustersID.xlsx") %>% as.data.frame()
head(macro)
rownames(macro) <- macro$Cell

metadt$macro <- ""
head(macro)
length(macro$Cell)
cellid <- intersect(rownames(metadt), macro$Cell)
length(cellid)
length(macro$Cell)
macro <- macro[cellid, ]
metadt[rownames(macro), "macro"] <- macro$Cluster
table(metadt$macro)
table(metadt$macro, metadt$manual.ident)

id3 <- vector(mode = "character", length = length(metadt$macro))
l <- which(metadt$macro == "0")
id3[l] <- "M1"
table(id3)

l <- which(metadt$macro == "1")
id3[l] <- "cDC2"
table(id3)

l <- which(metadt$macro == "2")
id3[l] <- "M2a"
table(id3)

l <- which(metadt$macro == "3")
id3[l] <- "M2b"
table(id3)

l <- which(metadt$macro == "4")
id3[l] <- "Monocyte"
table(id3)

l <- which(metadt$macro == "5")
id3[l] <- "others"
table(id3)

l <- which(metadt$macro == "6")
id3[l] <- "pDC"
table(id3)

l <- which(metadt$macro == "7")
id3[l] <- "cDC1"
table(id3)

l <- which(metadt$macro == "8")
id3[l] <- "Proliferative-M"
table(id3)
metadt$macro <- id3
table(metadt$macro)


l <- which(metadt$macro %in% c("cDC1", "pDC", "cDC2"))
metadt$ident1[l] <- "DCs"
table(metadt$ident1)

l <- which(metadt$macro %in% c("others", "Proliferative-M"))
metadt$ident1[l] <- "M-Other"
table(metadt$ident1)

l <- which(metadt$macro %in% c("M1", "M2a", "M2b", "Monocyte"))
metadt$ident1[l] <- id3[l]
table(metadt$ident1)

table(metadt$ident1)
table(metadt$ident2)

sr_annoed@meta.data <- metadt
saveRDS(sr_annoed, file = "objs/seurat_annoed.rds")

# 准备1
Idents(sr_annoed) <- sr_annoed$ident1
levels(sr_annoed)
sr_mait_vs_others <- subset(sr_annoed, idents = c("Unclassical-T", "M-Other"), invert = T)
table(sr_mait_vs_others$ident1)
colnames(sr_mait_vs_others@meta.data)
sr_mait_vs_others@meta.data <- sr_mait_vs_others@meta.data[, c(2,3,7,8,9)]
names(sr_mait_vs_others@meta.data)[names(sr_mait_vs_others@meta.data) == "ident1"] <- "celltype"
colnames(sr_mait_vs_others@meta.data)
saveRDS(sr_mait_vs_others, file = "objs/seurat_mait_vs_others.rds")

# 准备2
Idents(sr_annoed) <- sr_annoed$ident2
levels(sr_annoed)
sr_mait_early_vs_lately <- subset(sr_annoed, idents = c("Unclassical-T"), invert = T)
table(sr_mait_early_vs_lately$ident2)
levels(sr_mait_early_vs_lately)
colnames(sr_mait_early_vs_lately@meta.data)
sr_mait_early_vs_lately@meta.data <- sr_mait_early_vs_lately@meta.data[, c(2,3,7,8,10)]
names(sr_mait_early_vs_lately@meta.data)[names(sr_mait_early_vs_lately@meta.data) == "ident2"] <- "celltype"
colnames(sr_mait_vs_others@meta.data)
saveRDS(sr_mait_early_vs_lately, file = "objs/seurat_mait_early_vs_lately.rds")

sr_early <- subset(sr_mait_early_vs_lately, stage == "early")
sr_lately <- subset(sr_mait_early_vs_lately, stage == "lately")
levels(sr_early)
levels(sr_lately)
saveRDS(sr_early, file = "objs/seurat_mait_early.rds")
saveRDS(sr_lately, file = "objs/seurat_mait_lately.rds")
