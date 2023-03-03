# title: "cellchat_mait"
# author: "thathat"
# date: '2023-02-16'
# output: html_document

# 0. 准备工作 ####
# devtools::install_github("renozao/NMF@devel", upgrade = F)
# devtools::install_github("jokergoo/circlize", upgrade = F)
# devtools::install_github("jokergoo/ComplexHeatmap", upgrade = F)
# devtools::install_github("sqjin/CellChat", upgrade = F)
# chooseBioCmirror()
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(NMF)

options(stringsAsFactors = FALSE)


## 0.1 创建CellChat对象1 ####
rm(list = ls())
# Create a CellChat object from a Seurat object.
{
  # load("F:/biotech/projects/rat_annoed/seurat/outputs/seurat_preped.rda")
  # 
  # metadt$ident <- ""
  # id <- vector(mode = "list", length = length(metadt$ident))
  # for (i in 1:length(id)) {
  #   ifelse(metadt$CD4[i] != "", id[i] <- metadt$CD4[i],
  #          ifelse(metadt$CD8[i] != "", id[i] <- metadt$CD8[i], id[i] <- metadt$celltype[i]))
  # }
  # id <- unlist(id)
  # metadt$ident <- id
  # # View(metadt)
  # sr_preped@meta.data <- metadt
  # saveRDS(sr_preped, file = "F:/biotech/projects/rat_annoed/seurat/outputs/seurat_annoed.rds")
  # table(sr_preped$ident)
  # table(sr_preped$celltype)
  # # 准备1
  # celltype <- c("CD8T-MAIT", "CD4-TEM_TH1like", "CD4-TEMRA_TEFF", "CD4-TFH", "CD4-TH17", "CD4-TN", "CD4-Treg", "CD8-TEM", "CD8-TEX", "CD8-TN", "CD8T-other", "PT","LOH2")
  # Idents(sr_preped) <- sr_preped$ident
  # sr <- subset(sr_preped, idents = celltype)
  # CD8_others <- c("CD8-TEM", "CD8-TEX", "CD8-TN", "CD8T-other")
  # CD4 <- c("CD4-TEM_TH1like", "CD4-TEMRA_TEFF", "CD4-TFH", "CD4-TH17", "CD4-TN", "CD4-Treg")
  # metadt <- sr@meta.data
  # id <- vector(mode = "list", length = length(metadt$ident))
  # for (i in 1:length(id)) {
  #   ifelse(metadt$ident[i] %in% CD8_others, id[i] <- "CD8-Others",
  #          ifelse(metadt$ident[i] %in% CD4, id[i] <- "CD4", id[i] <- metadt$ident[i]))
  # }
  # id <- unlist(id)
  # sr_mait1 <- sr
  # sr_mait1$ident <- id
  # table(sr_mait1$ident)
  # saveRDS(sr_mait1, file = "objs/seurat_mait_early_vs_lately.rds")
  # # 准备2
  # celltype <- c("CD8T-MAIT", "CD8-TEM", "CD8-TEX", "CD8-TN", "Monocytic","Proliferative-CD4T", "Proliferative-CD8T")
  # sr <- subset(sr_preped, idents = celltype)
  # table(sr$ident)
  # saveRDS(sr, file = "objs/seurat_mait_vs_others.rds")
}

sr_obj <- readRDS("objs/seurat_mait_vs_others.rds")

cc_obj <- createCellChat(object = sr_obj, group.by = "ident", assay = "RNA")

# see the order of idents
levels(cc_obj@idents)
# > levels(cc_obj@idents)
# [1] "Monocytic"          "CD8-TEM"            "Proliferative-CD4T" "CD8T-MAIT"         
# [5] "CD8-TEX"            "Proliferative-CD8T" "CD8-TN" 

table(cc_obj@idents)

# 1. 预处理 ####
## 1.1-1.2 导入配受体数据库 ####

# use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse
# colnames(CellChatDB$interaction)

# showDatabaseCategory(CellChatDB)
# # Show the structure of the database
# dplyr::glimpse(CellChatDB$interaction)

# unique(CellChatDB$interaction$annotation)
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
# cc_obj@DB <- CellChatDB.use # set the used database in the object

cc_obj@DB <- CellChatDB
# subset the expression data of signaling genes for saving computation cost
cc_obj <- subsetData(cc_obj)

future::plan("multisession", workers = 2) # 多线程加速，可能出bug，可不执行
cc_obj <- identifyOverExpressedGenes(cc_obj)
cc_obj <- identifyOverExpressedInteractions(cc_obj)
cc_obj <- projectData(cc_obj, PPI.mouse)


## 1.3 互作推断（计算） ####

# 默认 population.size = F, raw.use = T 。即不考虑细胞亚群细胞数影响，不调用projectData(cc_obj, PPI.mouse)处理后矩阵
# 如果数据集未经处理，可反映现实情况，二参可调用
# 用时较长
cc_obj <- computeCommunProb(cc_obj, population.size = T, raw.use = F)

# 如果亚群细胞数过小可舍弃，该函数可筛选
cc_obj <- filterCommunication(cc_obj, min.cells = 10)

# 在信号通路水平判断细胞互作
cc_obj <- computeCommunProbPathway(cc_obj)
# Ps.受配体对水平、信号通路水平的细胞通讯网络分别存储在slot 'net'和'netP'中
cc_obj@net
cc_obj@netP
cc_obj@netP$pathways

# 整合细胞通讯网络
cc_obj <- aggregateNet(cc_obj)

head(cc_obj@LR$LRsig)

## 1.4 细胞通路可视化前机器学习 ####
cc_obj <- netAnalysis_computeCentrality(cc_obj, slot.name = "netP")

saveRDS(cc_obj, file = "objs/cellchat_mait_vs_others.rds")

# 2 可视化 ####
## 2.1 可视化整合的细胞通信网络 ####
# 例如，使用圆图显示任意两个细胞组之间的相互作用次数或总交互强度（比重）。
cc_obj <- readRDS(file = "objs/cellchat_mait_vs_others.rds")

# 保证不同样本、不同亚群有可比性
# number of cells in each cell group
groupSize <- as.numeric(table(cc_obj@idents))

# cat(paste("max(cc_obj@net$count): ", max(cc_obj@net$count), "\n", sep = ""), file = 'outputs/for.max.txt', append = T)
# cat(paste("max(cc_obj@net$weight): ", max(cc_obj@net$weight), "\n", sep = ""), file = 'outputs/for.max.txt', append = T)
# cat(paste("max(groupSize): ", max(groupSize), "\n", sep = ""), file = 'outputs/for.max.txt', append = T)

max.count <- max(cc_obj@net$count)
max.weight <- max(cc_obj@net$weight)
pdf(file = "outputs/plots/1.1 netVisual_interaction.pdf")
netVisual_circle(cc_obj@net$count, 
                 vertex.weight = groupSize, 
                 vertex.weight.max = max(groupSize), 
                 edge.weight.max = max.count, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(cc_obj@net$weight, 
                 vertex.weight = groupSize, 
                 vertex.weight.max = max(groupSize), 
                 edge.weight.max = max.weight, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")
dev.off()

netVisual_heatmap(cc_obj, slot.name = "netP", measure = "weight")
dev.off()

### 2.1.1 可视化单对单细胞通信 ####
# 由于细胞通信网络复杂，我们可以检查每个细胞组发送的信号
# 在这里，可控制参数edge.weight.max，以便可以比较不同网络之间的边缘权重
mat <- cc_obj@net$count
pdf("outputs/plots/1.2.1 netVisual_interaction_deploy_counts.pdf")
par(mfrow = c(3,4), xpd = T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   edge.weight.max = max(mat), 
                   weight.scale = T, 
                   title.name = rownames(mat)[i])
}

mat <- cc_obj@net$weight
pdf("outputs/plots/1.2.2 netVisual_interaction_deploy_weight.pdf")
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

## 2.2 信号通路水平可视化 ####
# check the order of cell identity to set suitable vertex.receiver
levels(cc_obj@idents)
vertex.receiver = seq(1,7) # a numeric vector

# 显示显著通路
cc_obj@netP$pathways

cc_obj@LR$LRsig$pathway_name

pwn <- unique(cc_obj@netP$pathways)
# cc_obj@LR$LRsig$antagonist
# ant <- unique(cc_obj@LR$LRsig$antagonist)

netVisual_aggregate(cc_obj, signaling = pwn,  
                    vertex.receiver = vertex.receiver, 
                    vertex.size = groupSize)

pathways.show <- c("SPP1")
for (psi in pathways.show) {
  pdf(paste("1.3 ", psi, "_hierarchy.pdf", sep = ""))
  netVisual_aggregate(cc_obj, signaling = psi, 
                      vertex.receiver = vertex.receiver, 
                      layout = "hierarchy")
  netVisual_aggregate(cc_obj, signaling = psi, 
                      vertex.receiver = vertex.receiver, 
                      layout = "circle")
  netVisual_aggregate(cc_obj, signaling = psi, 
                      vertex.receiver = vertex.receiver, 
                      layout = "chord")
  dev.off()
}


# 计算和可视化每个配体-受体对整个信号通路的贡献度。
netAnalysis_contribution(cc_obj, signaling = pwn)

# 可视化指定单一通路
pathways.show.select <- "SPP1"
# 经典的配受体圈图：
netVisual_aggregate(cc_obj, 
                    signaling = pathways.show.select, 
                    layout = "circle", 
                    vertex.weight = groupSize, 
                    pt.title=20, vertex.label.cex = 1.7)


# 识别细胞群的信号转导作用，通过计算每个细胞群的网络中心性指标，CellChat允许随时识别细胞间通信网络中的主要发送者、接收者、调解者和影响者。

cc_obj <- netAnalysis_computeCentrality(cc_obj, slot.name = "netP")

netAnalysis_signalingRole_heatmap(cc_obj, pattern = "all")
netAnalysis_signalingRole_heatmap(cc_obj, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cc_obj, pattern = "incoming")

pdf(file = 'pathway.pdf', height = 6, width = 12)
par(mfrow = c(3,1), xpd=TRUE)
netAnalysis_signalingRole_network(cc_obj, slot.name = "netP")
dev.off()

netAnalysis_signalingRole_scatter(cc_obj, pattern = "all")
netAnalysis_signalingRole_scatter(cc_obj, pattern = "outgoing")
netAnalysis_signalingRole_scatter(cc_obj, pattern = "incoming")



# 识别特定细胞群的全局通信模式和主要信号。除了探索单个通路的详细通讯外，一个重要的问题是多个细胞群和信号通路如何协调运作。CellChat采用模式识别方法来识别全局通信模式以及每个小群的关键信号。

# 识别分泌细胞外向交流模式。随着模式数量的增加，可能会出现冗余的模式，使得解释通信模式变得困难。我们选择了5种模式作为默认模式。一般来说，当模式的数量大于2时就可以认为具有生物学意义。

# 运行selectK推断模式的数量
selectK(cc_obj, pattern = "outgoing")

nPatterns = 3

cc_obj <- identifyCommunicationPatterns(cc_obj, pattern = "outgoing", k = nPatterns)

# Visualize the communication pattern using river plot
netAnalysis_river(cc_obj, pattern = "outgoing")

# Visualize the communication pattern using dot plot
netAnalysis_dot(cc_obj, pattern = "outgoing")



selectK(cc_obj, pattern = "incoming")
nPatterns = 3

cc_obj <- identifyCommunicationPatterns(cc_obj, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cc_obj, pattern = "incoming")


# dot plot
netAnalysis_dot(cc_obj, pattern = "incoming")


# 信号网络的多重和分类学习分析

# 此外，CellChat 能够量化所有重要信号通路之间的相似性，然后根据其CellChat 网络的相似性对其进行分组。分组可以基于功能或结构相似性进行。

# 功能相似性：功能相似度高表示主要发送器和接收器相似，可解释为两个信号通路或两个配体受体对具有相似的作用。功能相似性分析要求两个数据集之间的细胞群组成相同。

# 结构相似性：结构相似性用于比较其信号网络结构，而不考虑发送器和接收器的相似性。

# 根据信号组的功能相似性识别信号组
cc_obj <- readRDS(file = "objs/cellchat_mait_vs_others_analysed.rds")
cc_obj <- computeNetSimilarity(cc_obj, type = "functional")
cc_obj <- netEmbedding(cc_obj, type = "functional")

cc_obj <- netClustering(cc_obj, type = "functional")

# Visualization in 2D-space
netVisual_embedding(cc_obj, type = "functional", label.size = 3.5)

netVisual_embeddingZoomIn(cc_obj, type = "functional", nCol = 2)


# 基于结构相似性识别信号组
cc_obj <- computeNetSimilarity(cc_obj, type = "structural")
cc_obj <- netEmbedding(cc_obj, type = "structural")

cc_obj <- netClustering(cc_obj, type = "structural")

# Visualization in 2D-space
netVisual_embedding(cc_obj, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cc_obj, type = "structural", nCol = 2)


# 保存cellchat对象
saveRDS(cc_obj, file = "objs/cellchat_mait_vs_others_analysed.rds")


