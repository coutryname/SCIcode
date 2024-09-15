#remotes::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratData)
# AvailableData(
#InstallData("pbmc3k")
imm_anno <- readRDS("imm_anno.RDS")
pbmc3k <- imm_anno
options(stringsAsFactors = FALSE)
data("pbmc3k")
#pbmc3k里的seurat_annotations有一些NA注释，过滤掉
data.input = pbmc3k@assays$RNA@data
meta.data =  pbmc3k@meta.data
meta.data = meta.data[!is.na(meta.data$SCT_snn_res.0.4),]
data.input = data.input[,row.names(meta.data)]

#设置因子水平
meta.data$SCT_snn_res.0.4 = factor(meta.data$SCT_snn_res.0.4,
                                      levels = c("CD4+ T_cell", "CD8+ T_cell", "Macrophage", "B_cell", "NK_cell", 
                                                 "NKT_cell", "Regulatory_T_cell", "Neutrophil", "Paneth_cell", "Plasma_cell", "Stem_cell"))

### 1.2 Create a CellChat object
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "cell_type_imm")
### 1.3 可在cellchat对象的meta插槽中添加表型信息
# 添加meta.data信息
cellchat <- addMeta(cellchat, meta = meta.data)

# 设置默认的labels
levels(cellchat@idents) # show factor levels of the cell labels
cellchat <- setIdent(cellchat, ident.use = "cell_type_imm") 
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
### 1.4 加载CellChat受配体数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
### 1.5 对表达数据进行预处理，用于细胞间通讯分析
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 1) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat,population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
### cellchat取子集
barcode.use = sample(row.names(cellchat@meta),100)
cellchat.subset = subsetCellChat(cellchat,cells.use = barcode.use)
#获取所有的配受体对以及其通讯概率
df.net <- subsetCommunication(cellchat)
head(df.net)
#以通路为单位提取通讯信息
df.pathway = subsetCommunication(cellchat,slot.name = "netP")
df.net <- subsetCommunication(cellchat, sources.use = c(1), targets.use = c(2,3))
head(df.net)
df.net <- subsetCommunication(cellchat, signaling = c("MIF", "TNF"))
head(df.net)
#在信号通路水平推断细胞通讯
#CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通讯概率，
#计算信号通路层面的通讯概率。
#注意：推断出的每个配体-受体对的细胞间通讯网络和
#每个信号通路分别存储在' net '和' netP '插槽中。
cellchat <- computeCommunProbPathway(cellchat)
head(cellchat@net)
head(cellchat@netP)
#Step7. 计算加和的cell-cell通讯网络
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,file='cellchat.RDS')
cellchat <- readRDS("cellchat.RDS")
#然后，我们还可以可视化加和的细胞间通讯网络。
#例如，使用circle plot显示任意两个细胞亚群之间的通讯次数或总通讯强度(权重)：
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weights/strength")
dev.off()
#由于细胞间通讯网络的复杂性，我们可以对每个细胞亚群发出的信号进行检测。
#这里我们还控制参数edge.weight.max，以便我们可以比较不同网络之间的边权值：
mat <- cellchat@net$weight
par(mfrow = c(4,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize,
                   weight.scale = T, edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}
dev.off()

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
pathways.show <- c("TNF") 
vertex.receiver = c(1,2,3,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()
# Circle plot show pathway
par(mfrow=c(1,2))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",label.edge= T)
# Circle plot show L-R pairs
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
LR.show
#"MIF_CD74_CXCR4"
# Vis
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = "MIF_CD74_CXCR4", layout = "circle")
par(mfrow = c(1,2), xpd=TRUE)
# Chord diagram 下面这两幅图等价
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord",title.name = "Chord diagram  1")
netVisual_chord_cell(cellchat, signaling = pathways.show,title.name = "Chord diagram  2: show cell type")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
### Bubble plot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:7), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:5), 
                 signaling = c("MIF","TNF"), remove.isolate = FALSE)
