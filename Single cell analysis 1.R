library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(readr)
library(stringr)
library(progeny)
library(scales)
library(ggsci)
#设定阈值
nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 30
pHB_lower <- 0
pHB_upper <- 5

theme_set(theme_cowplot())

#设定配色方案####
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467")


samples <- read_excel("./metadata/patients_metadata.xlsx", range = cell_cols("A:A")) %>% .$sample_id


for (i in seq_along(samples)[15:20]){
  assign(paste0("scs_data", i), Read10X(data.dir = paste0("./cellranger/", samples[i], "/filtered_feature_bc_matrix")))
}

### 创建seurat对象，也是创建6个对象
for (i in seq_along(samples)[15:20]){
  assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = samples[i], min.cells = 3))
}

### merge一键融合，从15-20个样本
seu_obj <- merge(seu_obj15, y = c(seu_obj16, seu_obj17, seu_obj18, seu_obj19, seu_obj20), add.cell.ids = samples[15:20], project = "lung")
seu_obj <- Read10X(data.dir = "GSE202642_RAW")
seu_obj <- CreateSeuratObject(counts = seu_obj, project = "samples", min.cells = 5, min.features = 500)

### 计算线粒体等基因的比例
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^MT-", col.name = "pMT")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^HBA|^HBB", col.name = "pHB")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^RPS|^RPL", col.name = "pRP")

qcparams <- c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")
for (i in seq_along(qcparams)){
  print(VlnPlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident", pt.size = 0))
}
for (i in seq_along(qcparams)){
  print(RidgePlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident"))
}
### 批量看计算结果，很不错
VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "pMT"), pt.size = 0, group.by = "orig.ident", ncol = 1, log = T)
ggsave2("QC.pdf", path = "figure", width = 20, height = 20, units = "cm")
图片
### 清空环境

remove(seu_obj15)
remove(seu_obj16)
remove(seu_obj17)
remove(seu_obj18)
remove(seu_obj19)
remove(seu_obj20)


remove(scs_data15)
remove(scs_data16)
remove(scs_data17)
remove(scs_data18)
remove(scs_data19)
remove(scs_data20)


# 过滤
# 定义两个过滤作图函数
qc_std_plot_helper <- function(x) x + 
  scale_color_viridis() +
  geom_point(size = 0.01, alpha = 0.3)

qc_std_plot <- function(seu_obj) {
  qc_data <- as_tibble(FetchData(seu_obj, c("nCount_RNA", "nFeature_RNA", "pMT", "pHB", "pRP")))
  plot_grid(
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pMT))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pHB))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pRP))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pMT, color = nFeature_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pHB, color = nFeature_RNA))) + 
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pRP, color = nFeature_RNA))) + 
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pMT, color = nCount_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pHB, color = nCount_RNA))) + 
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pRP, color = nCount_RNA))) + 
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nCount_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nFeature_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    
    
    ggplot(gather(qc_data, key, value), aes(key, value)) +
      geom_violin() +
      facet_wrap(~key, scales = "free", ncol = 5),
    
    ncol = 3, align = "hv"
  )
}


## 过滤前

seu_obj_unfiltered <- seu_obj

# 作图，有点慢，但是一般可以出来
qc_std_plot(seu_obj_unfiltered)

ggsave2("before_filtering.pdf", path = "figure", width = 30, height = 30, units = "cm")
## 过滤
seu_obj <- subset(seu_obj_unfiltered, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & pMT < pMT_upper & pHB < pHB_upper)

#有点慢，同样可以出
qc_std_plot(seu_obj)
ggsave2("after_filtering.pdf", path = "figure", width = 30, height = 30, units = "cm")

# 看过滤了多少
seu_obj_unfiltered
seu_obj
# SCTransform提供了去除深度的策略###############################
#利用seurat的SCTransform消除深度差异，非常推荐的标准化的方法，之后保存下次用！
seu_obj <- SCTransform(seu_obj, verbose = T, vars.to.regress = c("nCount_RNA", "pMT"), conserve.memory = T)

saveRDS(seu_obj, file = "scRNA_SCTransform.RDS")

seu_obj=readRDS('scRNA_SCTransform.RDS')
seu_obj=readRDS('scRNA_main_annotated.RDS')

#降维
seu_obj <- RunPCA(seu_obj)
###根据折线何时平缓选择PC数，一般10~20
ElbowPlot(seu_obj, ndims = 50)

seu_obj <- RunUMAP(seu_obj, dims = 1:15, verbose = T)
seu_obj <- RunTSNE(seu_obj, dims = 1:15, verbose = T)
seu_obj <- FindNeighbors(seu_obj, dims = 1:15)
library(ggplot2)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  seu_obj <- FindClusters(seu_obj, resolution = i)
  print(DimPlot(seu_obj, reduction = "tsne") + labs(title = paste0("resolution: ", i)))
}

###简单做图看下
DimPlot(seu_obj, group.by = "orig.ident",reduction = 'tsne')
###简单做图看下，但是选resolution等于几，手工注释并不重要
DimPlot(seu_obj, group.by = "SCT_snn_res.0.4", label = T,reduction = 'umap')
DimPlot(seu_obj, group.by = "SCT_snn_res.0.4", label = T,reduction = 'tsne')
#主要的细胞类型注释——区分上皮、间质、免疫细胞
### 选定marker,这边纳入了上皮marker如EPCAM，间质如VWF，免疫如PTPRC，KIT是间叶特征marker
mainmarkers <- c("PTPRC","EPCAM",'PECAM1', "CD79A"
                 , "KRT19", "ACTA2")
mainmarkers <- c( "CD3D", "CD8A",  "FOXP3", "TRDC", "NKG7"
                 ,  "CD14"
                 , "CD163", "CD1C",  "LAMP3", "TPSAB1", "CSF3R"
                 ,  "S100A8",  "EPCAM",    "VWF",  "COL1A1",  "CLEC10A")
dir.create('SBHX')
dev.off()

### 对每个marker基因作图并观察,这里建议画umap更清楚一些scale_color_viridis()
for (i in seq_along(mainmarkers)) {
  FeaturePlot(seu_obj, features = mainmarkers[i], coord.fixed = T, order = T, cols = viridis(10))
  ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".tiff"), path = "./SBHX", width = 10, height = 10, units = "cm")
}

### 输出下面两个图，做对比，同时结合一下前面每个基因的featureplot，即通过手头的三个资料去分析每个cluster的细胞大类
### 手工注释，resolution并不太重要，我们就选择0.2为例，共19个cluster，但是很多cluster实际上一样
### 但是建议大一点的resolution，否则可能混淆。
DotPlot(seu_obj, features = mainmarkers, group.by = "SCT_snn_res.0.4") + 
  coord_flip() + 
  scale_color_binned()
ggsave2("DotPlot_mainmarkers.pdf", path = "./SBHX", width = 10, height = 5)

DimPlot(seu_obj, group.by = "SCT_snn_res.0.4", label = T, label.size = 5)
ggsave2("DimPlot_all_clusters.pdf", path = "./SBHX", width = 7, height = 7)


### 接下来小编去文件夹对比我的marker去了，你也赶紧去对比吧
### 仔细观察，发现0，10，15，17上皮，1-9，11，12，18免疫细胞，13,14,16间质，你需要同时考虑到Umap的接近程度和marker的表达。
immuu <- immu
Idents(seu_obj) <- seu_obj$SCT_snn_res.0.4
annotation_curated_main <- read_excel("./annotation/annotation_main.xlsx")
new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(seu_obj)
#加入seurat对象中
seu_obj <- RenameIdents(seu_obj, new_ids_main)
seu_obj@meta.data$Main_cell_type <- Idents(seu_obj)

DimPlot(seu_obj, label = T, label.size = 5,group.by="Main_cell_type" ,
        cols= c("darkgoldenrod2","seagreen","steelblue","#B40F20"))
DimPlot(seu_obj, label = T, label.size = 5,group.by="Main_cell_type" ,
        cols= c("seagreen","darkgoldenrod2","steelblue"))
dev.off()
#添加标签图
DimPlot(seu_obj, reduction = "tsne", group.by = "Main_cell_type", label = T, label.box = T)
DimPlot(seu_obj, reduction = "umap", group.by = "Main_cell_type", label = T, label.box = T)

#出正式图

ggsave2("DimPlot_main_cell_type2.pdf", path = "./SBHX", width = 7, height = 6)
saveRDS(seu_obj, file = "scRNA_main_annotated.RDS")

#### 增加进病人信息
seu_obj=readRDS('scRNA_main_annotated.RDS')
metatable <- read_excel("./metadata/patients_metadata.xlsx")
### 从seurat对象里提取orig.ident(病人)
metadata <- FetchData(seu_obj, "orig.ident")
metatable$cell_id <- rownames(metatable)
metatable$sample_id <- metatable$orig.ident
rownames(metatable) <- metatable$...1
metatable <- metatable[,-1]
metadata <-metatable
write.csv(metadata, file = "metadata.csv")
metatable <- read_excel("./metadata2.xlsx")
### dplyr中的left_join合并
metadata <- left_join(x = metadata, y = metatable, by = "sample_id")
### 重新加行名
rownames(metadata) <- metadata$cell_id

### 通过AddMetaData直接补入15列
seu_obj <- AddMetaData(seu_obj, metadata = metadata)

# 细胞周期运算，简单过

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

score_cc <- function(seu_obj) {
  seu_obj <- CellCycleScoring(seu_obj, s.genes, g2m.genes)
  seu_obj@meta.data$CC.Diff <- seu_obj@meta.data$S.Score - seu_obj@meta.data$G2M.Score
  return(seu_obj)
}

seu_obj <- score_cc(seu_obj)

#看上去还可以
FeatureScatter(seu_obj, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) +
  coord_fixed(ratio = 1)
# 提取子集，重归一化

### saveRDS(seu_obj, file = "scRNA_main_anno_add_pt.RDS")

### Idents设置最主要的区分因子
Idents(seu_obj) <- seu_obj@meta.data$Main_cell_type
epi <- subset(seu_obj, idents = "Endothelial cell")
imm <- subset(seu_obj, idents = "Immune cells")
str <- subset(seu_obj, idents = "Fibroblast")

#也可以这样
#Cells.sub <- subset(seu_obj@meta.data, Main_cell_type =='Immune',)
#scRNAsub <- subset(seu_obj, cells=row.names(Cells.sub))


### 因为总量变了，所以需要重新归一化，否则分布将被打破
epi <- ScaleData(epi)
imm <- ScaleData(imm)
str <- ScaleData(str)

saveRDS(epi, file = "epi.RDS")
saveRDS(imm, file = "imm2次新分群.RDS")
saveRDS(str, file = "str.RDS")
#作图：不同注释方式的Umap
### 修改group.by
DimPlot(seu_obj, group.by = "tissue_type", cols = use_colors, pt.size = 0.1)
ggsave2("区分肿瘤正常.pdf", path = "./figure/", width=7 ,height=6)

###区分病人但不区分N/T，32号病人有些特别？我们现按下不表
DimPlot(seu_obj, group.by = "patient_id", cols = use_colors, pt.size = 0.1)
ggsave2("区分病人.png", path = "./result_3v3/",  width=7 ,height=6)

DimPlot(seu_obj, group.by = "Main_cell_type", cols = use_colors, pt.size = 0.1)
ggsave2("区分细胞类群.png", path = "./result_3v3/",  width=7 ,height=6)

theme_set(theme_cowplot())

use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467")


#加载数据

epi <- readRDS("epi.RDS")
imm2 <- readRDS("imm.RDS")
str <- readRDS("str.RDS")
imm <-imm2

# 主要亚群重聚类####

### 上皮重聚类

epi <- RunPCA(epi)
ElbowPlot(epi,  ndims = 50)


epi <- RunUMAP(epi, dims = 1:20)
epi <- FindNeighbors(epi, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  epi <- FindClusters(epi, resolution = i)
  print(DimPlot(epi, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i)))
}

### 我们依然选择比较大的resolution，手动注释没有关系
### 控制在20个cluster左右即可
Idents(epi) <- epi@meta.data$SCT_snn_res.1


### 免疫重聚类

imm <- RunPCA(imm)
ElbowPlot(imm,  ndims = 50)


imm <- RunUMAP(imm, dims = 1:20)
imm <- FindNeighbors(imm, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  imm <- FindClusters(imm, resolution = i)
  print(DimPlot(imm, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}
### 控制在20个cluster左右即可
Idents(imm) <- imm@meta.data$SCT_snn_res.0.4


### 基质重聚类

str <- RunPCA(str)
ElbowPlot(str, ndims = 50)


str <- RunUMAP(str, dims = 1:20)
str <- FindNeighbors(str, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  str <- FindClusters(str, resolution = i)
  print(DimPlot(str, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}
### 细胞不多，可以少点，随意
Idents(str) <- str@meta.data$SCT_snn_res.2

# 疾病特异性细胞亚群
### 正常和肿瘤的上皮对比
### 发现确实多了一群
DimPlot(epi, group.by = "SCT_snn_res.1", label = T, repel = T, split.by = "tissue_type")

ggsave2("正常-肿瘤上皮.pdf", path = "./figure", width = 30, height = 15, units = "cm")

### 比较正常和肿瘤细胞比例柱形
epi_clusters <- FetchData(epi, vars = c("SCT_snn_res.1", "tissue_type"))

count_tumor <- epi_clusters %>% filter(tissue_type == "Tumor") %>% count() %>% as.numeric()
count_normal <- epi_clusters %>% filter(tissue_type == "Normal") %>% count() %>% as.numeric()
### 每个cluster计数总数
epi_counts <- epi_clusters %>% group_by(tissue_type) %>% count(SCT_snn_res.1)
### 除以总count数，是朴素的比例计算
proportion_tumor <- epi_counts %>% filter(tissue_type == "Tumor") %>% mutate(proportion = n/count_tumor)
proportion_normal <- epi_counts %>% filter(tissue_type == "Normal") %>% mutate(proportion = n/count_normal)

### 组合，最后一行很有意思，认为哪群比例高就属于哪群
proportion_epi <- full_join(proportion_normal, proportion_tumor, by = "SCT_snn_res.1") %>% 
  mutate(proportion.x = ifelse(is.na(proportion.x), 0,  proportion.x)) %>%  
  mutate(proportion.y = ifelse(is.na(proportion.y), 0,  proportion.y)) %>%
  mutate(tissue_type.x = "Normal") %>%
  mutate(tissue_type.y = "Tumor") %>%
  mutate(cluster_type = ifelse(proportion.x > proportion.y, "Normal", "Tumor"))

### 每个细胞的详细归属cluster
cluster_type_data <- left_join(x = epi_clusters, y = proportion_epi, by = "SCT_snn_res.1")
rownames(cluster_type_data) <- rownames(epi_clusters)

### 加入metadata中
epi <- AddMetaData(epi, select(cluster_type_data, cluster_type))


### 柱状图

n1 <- select(proportion_epi, c(tissue_type.x, SCT_snn_res.1, proportion.x)) %>%
  mutate(tissue_type = tissue_type.x) %>% 
  mutate(proportion = proportion.x) %>%
  mutate(tissue_type.x = NULL) %>%
  mutate(proportion.x = NULL)
t1 <- select(proportion_epi, c(tissue_type.y, SCT_snn_res.1, proportion.y)) %>%
  mutate(tissue_type = tissue_type.y) %>% 
  mutate(proportion = proportion.y) %>%
  mutate(tissue_type.y = NULL) %>%
  mutate(proportion.y = NULL)

proportion_epi2 <- rbind(n1, t1)

ggplot(proportion_epi2, aes(fill = tissue_type, y = proportion, x = SCT_snn_res.1)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = use_colors)

ggsave2("正常-肿瘤上皮比较.pdf", path = "figure/", width = 40, height = 20, units = "cm")

### 根据下面的图，我们可以明确地识别出肿瘤上皮中哪些类群是新增加的
### 也就是异型增生的上皮细胞——癌细胞
DimPlot(epi, group.by = "cluster_type",split.by = 'tissue_type' ,cols = use_colors, pt.size = 0.1)
ggsave2("新增的肿瘤细胞.pdf", path = "figure/", width = 8, height = 4)
# 详细细胞类型注释
# 先画一些landscape，了解一下再聚类的质量
### 我们需要了解到的是，不同细胞注释当然是要结合起来看的
### 另外，我们对免疫细胞加入一些阴性marker也很重要
### 髓系
myeloid_markers <- c("S100A12", "FCN1", "S100A8", "S100A9", "CD14", "CTSS","VCAN", "LYZ", 
                     "MARCO", "FCGR1A", "C1QA", "APOC1", "LGMN", "CTSB", "FCGR3A", 
                     "MAFB", "MAF", "CX3CR1", "ITGAM", "CSF1R",
                     "FABP4", "MCEMP1", 
                     "IL1B", "CXCL8", 
                     "APOE", "CD163", "C1QB", "C1QC", 
                     "FCER1A", "CD1C", "CLEC9A", 
                     "LILRA4", "CLEC4C", "JCHAIN", "IL3RA", "NRP1", 
                     "CLEC10A", "PTCRA", "CCR7", "LAMP3", 
                     "ITGAX", "CD68", "MKI67", "CDK1", "EPCAM")

### T细胞
tcell_nk_markers <- c("CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", 
                      "CD8A", "CCL5", "NCR1", "NKG7", "GNLY", "NCAM1", 
                      "KLRD1", "KLRB1", "CD69", "KLRG1", 
                      "MKI67", "CDK1", "EPCAM")


### B细胞
bcell_plasma_mast_markers <- c("MS4A1", "CD19", "CD79A", "JCHAIN", "IGHA1", 
                               "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", 
                               "IGKC", "IGLC2", "IGLC3", "CPA3", "KIT", "MS4A2", 
                               "GATA2",  "MKI67", "CDK1", "EPCAM")
### Immu细胞
immucell_plasma_mast_markers <- c("CD3D","CD4","CD8A","IL7R",
"ICOS", "GZMB","PRF1","NKG7",
"CD79A", "SLAMF7", "FCRL5",
'CD14','CD56', 'CD163', 'C1QA',"FOXP3","IKZF2","CTLA4","CLEC10A")
dev.off()
### Dotplot是最关键的手动注释的工具coord_flip() + 
#scale_color_viridis()
DotPlot(imm, features = myeloid_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("髓系标记热图.pdf", path = "./figure/", width = 20, height = 30, units = "cm")

DotPlot(imm, features = immucell_plasma_mast_markers , group.by = "SCT_snn_res.0.4") + coord_flip()+scale_color_viridis()
ggsave2("免疫细胞标记热图.pdf", path = "./SBHX", width = 20, height = 10, units = "cm")

DotPlot(imm, features = tcell_nk_markers, group.by = "SCT_snn_res.0.4") + coord_flip()
ggsave2("DotPlot_T_NK_markers.pdf", path = "./SBHX/", width = 20, height = 20, units = "cm")

DotPlot(imm, features = bcell_plasma_mast_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DotPlot_B_Plasma_markers.pdf", path = "./figure/", width = 20, height = 20, units = "cm")

DimPlot(imm,, reduction = "umap", group.by = "SCT_snn_res.0.4", label.box = T,label = T, split.by = "tissue_type")
ggsave2("正常肿瘤病人免疫细胞分群.pdf", path = "./SBHX", width =10, height =5)

DimPlot(imm, group.by = "patient_id", cols = use_colors)
ggsave2("病人免疫细胞分群.png", path = "./result_3v3/", width = 10, height = 5)
#添加标签图
DimPlot(seu_obj, reduction = "tsne", group.by = "Main_cell_type", label = T, label.box = T)
DimPlot(imm, reduction = "umap", 
        group.by = "SCT_snn_res.0.4", label = T, 
        label.box = T,split.by = "tissue_type")

# 重新降维聚类

### 上皮重聚类

epi <- RunPCA(epi)
ElbowPlot(epi,  ndims = 50)


epi <- RunUMAP(epi, dims = 1:20)
epi <- FindNeighbors(epi, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  epi <- FindClusters(epi, resolution = i)
  print(DimPlot(epi, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i)))
}

### 我们依然选择比较大的resolution，手动注释没有关系
### 控制在20个cluster左右即可
Idents(epi) <- epi@meta.data$SCT_snn_res.1


### 免疫重聚类
imm <-pbmc
imm <- RunPCA(imm)
ElbowPlot(imm,  ndims = 50)


imm <- RunUMAP(imm, dims = 1:20)
imm <- FindNeighbors(imm, dims = 1:20)
for (i in c(0.01,0.2, 0.3, 0.4, 0.5, 1)) {
  imm <- FindClusters(imm, resolution = i)
  print(DimPlot(imm, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}
### 控制在20个cluster左右即可
Idents(imm) <- imm@meta.data$SCT_snn_res.0.4


### 基质重聚类

str <- RunPCA(str)
ElbowPlot(str, ndims = 50)


str <- RunUMAP(str, dims = 1:20)
str <- FindNeighbors(str, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  str <- FindClusters(str, resolution = i)
  print(DimPlot(str, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}
#细胞不多，可以少点
Idents(imm) <- imm@meta.data$cell_type_immT

# 细胞类型注释的marker寻找
### 我们在之前的推文已经说过，一定要记得去找原文！
### 根据文章，使用Habermann et al.参考基因，
#https://www.biorxiv.org/content/10.1101/753806v1

habermann_Endl <- c("PLVAP", "VWF", "ACKR1", "CDH5")

habermann_imm <- c("CD2", "TRAC", "MS4A1", "IGHG1", "CD14",
                   "CD79A", "CD163", "CD86", "MS4A7", "CSF1R", "NKG7")

habermann_ito <- c("COLA2", "PDPN", "DCN", "LUM", "KRT19", "alpha-smooth muscle actin", "CD36", "F2")
habermann_imm <- c("CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", "KLRB1", 
                   "NKG7", "LYZ", "CD68", "ITGAX", "MARCO", "FCGR1A", "FCGR3A", "C1QA", "APOC1", "S100A12",
                   "FCN1", "S100A9", "CD14", "FCER1A", "CD1C", "CLEC9A", "LILRA4", "CLEC4C", 
                   "JCHAIN", "IGHG1", "MS4A1", "CD19", "CD79A", "CPA3", "KIT", "MKI67",
                   "CDK1", "CTLA4")
habermann_imm <- c("CD4", "FOXP3", "IL2RA", "CCR7", "CXCR5",
                   "CD8A", "CD8B",  "PRF1", "IFNG", 
                    "GZMB",  "TNF", "CD44",
                   "PDCD1", "CTLA4", "LAG3", "TIGIT","EOMES", "TOX")
habermann_immCD4 <- c("A2M", "ABCA2", "ABI3", "ADGRG1", "CD4",
                   "ADGRG5", "ADRB2", "AGPAT4", "AHNAK", "AKR1C3", "APMAP")
habermann_immCD8 <- c("ABCG1", "ACP5", "ACSL4", "ACTN4", "CD8A",
                      "AFAP1L2", "AHI1", "AIP", "AKAP5", "ANKRD10", "ANXA5")
habermann_immTrg <- c("ABI2", "ACOT9", "ACP5", "ACSL4", "ACTA2",
                      "ACTR3C", "ADAM10", "ADAT2", "ADCY3", "ADI1")
imm <- imm_lympho
### 免疫marker基因定细胞类型（Habermann et al，配合参考文献中Tata等人）
### "SCT_snn_res.0.5"到底选多少呢？
### 这需要根据后面已知亚群是否被混淆
### 例如Treg Foxp3在0.5时才被清晰地显现出来，0.4却不行(包含在大群的细胞里面了)
### 因此我需要放弃0.4，选择0.5
### 你可以自己featureplot去体会一下0.4和0.5带来的Treg的微妙差异
DotPlot(imm, features = habermann_immCD4, group.by = "SCT_snn_res.0.4") + coord_flip() + 
  scale_color_viridis()
ggsave2(filename = 'Trg免疫marker基因的dotplot为了注释.pdf',path = "./figure/", width =10, height =15)

DotPlot(imm, features = habermann_imm, group.by = "SCT_snn_res.0.4") + coord_flip() + 
  scale_color_viridis()
ggsave2(filename = '免疫marker基因的dotplot为了注释2.pdf',path = "./figure/", width =10, height =15)
DotPlot(epi, features = habermann_Endl, group.by = "SCT_snn_res.0.4") + coord_flip()1
ggsave2(filename = '内皮marker基因的dotplot为了注释.pdf',path = "./figure/", width =10, height =15)
DotPlot(str, features = habermann_ito, group.by = "SCT_snn_res.0.4") + coord_flip()
ggsave2(filename = '其他marker基因的dotplot为了注释.pdf',path = "./figure/", width =10, height =15)


DimPlot(imm,group.by =  "SCT_snn_res.0.4",label = T,split.by = 'tissue_type')
ggsave2(filename = '正常肿瘤病人免疫细胞分群为了注释.pdf',path = "./figure/", width =14, height =7)

### FOXP3帮助你确定一下resolution
FeaturePlot(imm, features =c('CD3E'),split.by = 'tissue_type')
# 详细细胞类型注释（正式）###########################################
### 先上免疫细胞,思路就是marker基因再加肿瘤正常对比
###0，4，12特征明显，普通T细胞
###1，3，10，13，14，20 正常里多，肿瘤里少 结合CD68 所以是正常肺泡巨噬（需要思维）
###2 明显NK
###5，7 普通T，分两群是因为res太高，看图便知，其实是一群
###6，18 CD8T，肿瘤飙升
###8，19 肿瘤升，CD14阳，Monocyte-derived macrophages
###9 明显Myeloid dendritic cells
###11 结合marker明显Monocytes，正常多，因为肿瘤里会刺激分化
###15 明显B
###16 明显Treg
###17 浆细胞
###21 明显Plasmacytoid dendritic cells
### 载入注释
annotation_curated_imm <- read_excel("./annotation/annotation_imm_T总.xlsx")
new_ids_imm <- annotation_curated_imm$cell_type_immT
imm_anno <- imm_anno修改后
Idents(imm_anno) <- imm_anno@meta.data$SCT_snn_res.0.4
names(new_ids_imm) <- levels(imm_anno)
imm_anno <- RenameIdents(imm_anno, new_ids_imm)
imm_anno@meta.data$cell_type_immT <- Idents(imm_anno)
imm_anno <- subset(imm_anno, subset = cell_type_immT != "Epithelial_contamination")
imm_anno <- ScaleData(imm_anno)
DimPlot(imm_anno,group.by =  "cell_type_immT",label = T,split.by = 'tissue_type', cols = use_colors)
ggsave2("DimPlot_imm_Normal_Tumor细胞注释最新- T细胞已分群.pdf", path = "./SBHX", width = 14, height = 7)
saveRDS(imm_anno, file="imm_anno修改后T细胞分群.rds")
theme_set(theme_cowplot())

### 设计配色，细胞类型注释已经好了，所以根据的是上辑第五篇的注释结果Treg

use_colors <- c(
  `Tumor `= "brown2",
  `Normal` = "deepskyblue2",
  G1 = "#46A9C2",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  `CD4 T` = "#FD6467",
  `CD8 T`= "#fff500",
  `TAMs` = "#FA9FB5",
  `DC`= "#DD3497",
  `Treg`= "#7A0177",
  `Exhausted T Cells`= "#c2e699",
  `Effector T Cells`= "#bcbddc",
  `B cells`= "#4a1486",
  `T cells`= "#46ACC2",
  `NK cells`= "#6bAEd6")
#`T_cell`= "#006837",  `CD4+ _cytotoxic T_cell` = "#6bAEd6",
#载入数据
### 数据在阅读原文里有，当然更建议你自行从第一步制作下来
#saveRDS(imm,file='imm.RDS')
imm_anno <- readRDS("imm_anno修改后.rds")

imm_anno@meta.data$cell_type_imm <- ordered(imm_anno@meta.data$cell_type_imm, levels = c("T cells", "Macrophage", "B_cell","NK_cell","NKT_cell","Regulatory_T_cell","Neutrophil","Paneth_cell",
  "Stem_cell","Plasma_cell"))


### 画一些图帮我们回忆一下这个数据集,例如我们看到肿瘤正常确实差的大
DimPlot(imm_anno, group.by = "tissue_type", cols = use_colors)
dir.create('./figure3')
ggsave2("DimPlot_imm_Normal_Tumor.pdf", path = "./figure3", width = 8, height = 7)


DimPlot(imm_anno, group.by = "cell_type_imm", split.by = "tissue_type", cols = use_colors, pt.size = 0.5)
ggsave2("DimPlot_umap.pdf", path = "./figure3", width = 10, height = 5)
dev.off()
#这一步还需要改改#
DotPlot(imm_anno, features = c("CD2", "LYZ", "TRAC", "CD79A", "IL-7R", "KLRB1", "CD14", "S100A9", 
                               "NKG7", "MS4A7", "MKI67", "CDK1", "CD163", "CSF1R", "IGHG1", 
                               "EPCAM", "CD4", "FOXP3", "IL2RA", "CD8A"), group.by = "cell_type_imm") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip() + 
  scale_color_viridis()
ggsave2("注释后dotplot_可以展示.pdf", path = "./figure2", width = 8, height = 7)
### 画一下好看的最终UMAP图
DimPlot(imm_anno, group.by = "cell_type_imm", split.by = "tissue_type", cols = use_colors, pt.size = 0.5)
ggsave2("DimPlot_umap.pdf", path = "./figure3", width = 10, height = 5)
### dotplot，已经有细胞了，所以这步是帮助我们确认一下结果
DotPlot(imm_anno, features = c("CD2", "CCL5", "TRAC", "CD69", "SRSF7"
                               , "DUSP2", "CCL20", "LTB", "CSF1R", "CD79A"
                               , "MS4A1", "IGHG1", "PTPRC", "EPCAM", "KRT19"
                               , "RPS26", "TNFRSF4", "TNFRSF18", "MKI67"
                               , "NKG7", "TUBB", "HMGN2"
                              ), group.by = "cell_type_imm") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip() + 
  scale_color_viridis()
### 发现非常漂亮！
ggsave2("注释后dotplot_可以展示.pdf", path = "./figure3", width = 8, height = 7)
#细胞比例图面面观####
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(stringr)
library(cowplot)
library(scales)
library(readr)
library(progeny)
library(gplots)
library(tibble)
library(grid)
library(rlang)
theme_set(theme_cowplot())

#载入数据
### 数据在阅读原文里有，当然更建议你自行从第一步制作下来
 #saveRDS(imm_anno,file='imm_anno.RDS')
imm_anno2 <- readRDS("imm2次新分群.RDS")
imm_anno <- imm_anno修改后
imm_anno@meta.data$cell_type_imm <- ordered(imm_anno@meta.data$cell_type_immT,
                                            levels = c("T cells","DC", 
                                                       "B celsl","NK cells","Treg","Neutrophil","undefined"
                                                    ))
# 分群分析############
# 免疫细胞分淋巴系和髓系
###T淋巴细胞亚群
imm_lympho <- subset(imm_anno, subset = cell_type_immTyaquan %in% c("Exhausted T Cells","Effector T Cells","CD8 T","CD4 T"))
###淋巴细胞亚群
imm_lympho <- subset(imm_anno, subset = cell_type_immT %in% c("T cells", 
                                                             "B cells","NK cells",
                                                             "Treg","Neutrophil",
                                                             "undefined","DC"))
 ### 记住重新scale！
imm_lympho <- ScaleData(imm_lympho)
          
          
### 髓系亚群
imm_myelo <- subset(imm_anno, subset = cell_type_imm %in% c("Alveolar_Macrophages",
                                                                      "Monocyte-derived macrophages",
                                                                      "Monocytes",
                                                                      "Myeloid dendritic cells",
                                                                      "Plasmacytoid_dendritic_cells"))
### 记住重新scale
imm_myelo <- ScaleData(imm_myelo)
          
### FetchData快速获得附加信息
### 淋巴系
lympho_counts <- FetchData(imm_lympho, vars = c("tissue_type", "cell_type_immT", "sample_id")) %>%  
mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")))
# 细胞比例图
### count根据patient_id，和tissue_type来统计总数，相当于excel分类汇总了
lympho_counts_tbl <- lympho_counts %>%
dplyr::count(cell_type_immT, tissue_type)
write_csv(lympho_counts_tbl, path = "./SBHX/lympho_counts_t细胞.csv")
          
          
 ### 髓系也一样
myelo_counts <- FetchData(imm_myelo, vars = c("tissue_type", "cell_type_imm", "sample_id", "patient_id")) %>%  
mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")))
          
myelo_counts_tbl <- myelo_counts %>%
dplyr::count(cell_type_imm, patient_id, tissue_type)
write_csv(myelo_counts_tbl, path = "./analysis/lmyelo_counts_tbl.csv")
#淋巴系比例tu          
ggplot(data = lympho_counts, aes(x = tissue_type, fill = cell_type_immT)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("barplot_lymphoid细胞亚群比例图.pdf", path = "./SBHX/", width = 8, height = 2)
#髓系比例图
ggplot(data = myelo_counts, aes(x = tissue_type, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("barplot_myeloid.pdf", path = "./analysis/", width = 8, height = 2)

#病人的异质性可见一斑，这是只有细胞比例图才可以揭示的现象！
lympho_counts %>%
  filter(tissue_type == "Tumor") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("lymphoid_per_patient.pdf", path = "./result_3v3/analysis/", width = 8, height = 2.5)

myelo_counts %>%
  filter(tissue_type == "Tumor") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("barplot_myeloid_per_patient.pdf", path = "./result_3v3/analysis/", width = 8, height = 2.5)

lympho_counts %>%
  filter(tissue_type == "Normal") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("normal_lymphoid.pdf", path = "./result_3v3/analysis/", width = 8, height =2.5)

myelo_counts %>%
  filter(tissue_type == "Normal") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("normal_myeloid.pdf", path = "./result_3v3/analysis/", width =8, height =2.5)

### 提出T_cell来源来分析（mdm）
scRNA_tcell <- subset(imm_anno, subset = cell_type_imm %in% c("CD4+T"
                                                            ,"CD8+T"
                                                            ,"Treg"
                                                            ))
scRNA_mdm=ScaleData(scRNA_tcell)

### 重新降维聚类
scRNA_mdm <- RunPCA(scRNA_mdm)
ElbowPlot(scRNA_mdm,  ndims = 50)

library(ggplot2)
scRNA_mdm<- RunUMAP(scRNA_mdm, dims = 1:20)
scRNA_mdm <- FindNeighbors(scRNA_mdm, dims = 1:20)
for (i in c(0.01,0.1,0.2, 0.3, 0.4, 0.5, 1)) {
  scRNA_mdm <- FindClusters(scRNA_mdm, resolution = i)
  print(DimPlot(scRNA_mdm, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}
saveRDS(scRNA_tcell,file='scRNA_tcell.RDS')
### cluster数自己决定，不要太多即可,否则相似的群多
table(scRNA_mdm$SCT_snn_res.0.1)
Idents(scRNA_mdm) <-scRNA_mdm@meta.data$SCT_snn_res.0.1
saveRDS(scRNA_mdm,file='scRNA_tcell.RDS')
DimPlot(scRNA_mdm,reduction = 'umap',group.by = 'SCT_snn_res.0.1',label = T,split.by = 'tissue_type')
ggsave2(filename = 'T_cell亚群-2.pdf',path = './T_cell fenxi',width = 8,height = 4)
### T细胞来源的巨噬细胞间具有较大异质性，所以考虑再分别注释

tcell_markers <- FindAllMarkers(scRNA_mdm, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)

# group_by和top_n组合
top_mdm_markers <- tcell_markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)
write.csv(top_mdm_markers, file = "top_mdm_markers2.csv")
DoHeatmap(scRNA_mdm,features=top_mdm_markers$gene,group.by = 'SCT_snn_res.0.1')+scale_fill_viridis()
ggsave2(filename = '热图细分.pdf',path = './T_cell fenxi',width = 10,height =9)
### 选定marker,这边纳入了CD8Tcell
mainmarkers <- c( "A2M", "ABCA2",  "ABI3", "ARGLU1", "ADGRG5"
                  ,  "ABCG1", "ACP5",  "ACSL4"
                  , "ACTN4", "ADAM19",  "AHI1", "AKAP5", "AKR1A1"
                  ,  "ARHGAP9",  "ABI2",  "ACOT9",  "ACP5",  "ACSL4",  "ADAM12")
mainmarkers <- c( "CD4", "SELL",  "CCR7", "ITK", "CD8A"
                  ,  "CD8B", "CD27",  "CRTAM")
#dir.create('annotation')
dev.off()

### 对每个marker基因作图并观察,这里建议画umap更清楚一些scale_color_viridis()
for (i in seq_along(mainmarkers)) {
  FeaturePlot(scRNA_mdm, features = mainmarkers[i], coord.fixed = T, order = T, cols = cividis(10))
  ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".tiff"), path = "./T_cell fenxi", width = 10, height = 10, units = "cm")
}

### 载入注释
annotation_curated_imm <- read_excel("./T_cell fenxi/annotation_imm.xlsx")
imm_anno <- scRNA_mdm
new_ids_imm <- annotation_curated_imm$cell_type_imm
Idents(imm_anno) <- imm_anno@meta.data$SCT_snn_res.0.01
names(new_ids_imm) <- levels(imm_anno)
imm_anno <- RenameIdents(imm_anno, new_ids_imm)
imm_anno@meta.data$cell_type_imm <- Idents(imm_anno)
imm_anno <- subset(imm_anno, subset = cell_type_imm != "Epithelial_contamination")
imm_anno <- ScaleData(imm_anno)
DimPlot(imm_anno,group.by =  "cell_type_imm",label = T,split.by = 'tissue_type', cols = use_colors)
ggsave2("Tcell亚群注释Normal_Tumor细胞注释图.pdf", path = "./figure2", width = 14, height = 7)
saveRDS(imm_anno, file="imm_anno.rds")
theme_set(theme_cowplot())
#UBE2S，RTKN2，STMN，SLC38A1，S100A10，PTTG1，CITED2， CD69
### dotplot，已经有细胞了，所以这步是帮助我们确认一下结果
DotPlot(imm_anno, features = c("UBE2S", "RTKN2",  "STMN1","SLC38A1"
                               ,"S100A10", "PTTG1", "CD69",
                      "CITED2"), group.by = "cell_type_imm") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip() + 
  scale_color_viridis()
### 发现非常漂亮！
rm(imm_anno)
ggsave2("注释后dotplot_可以展示.pdf", path = "./figure2", width = 8, height = 7)

### T细胞差异基因
devtools::install_github("elliefewings/DoMultiBarHeatmap")
#T细胞亚群
imm_T <- subset(imm_anno, subset = cell_type_imm %in% c("CD4+T_cell",
                                                        "CD8+T_cell"))

imm_T <- ScaleData(imm_T)

Idents(imm_T) <- imm_T@meta.data$cell_type_imm

markers_T <- FindAllMarkers(imm_T, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)

top_markers_T <- markers_T %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)

DoMultiBarHeatmap(imm_T, features = top_markers_T$gene, group.by = "cell_type_imm", additional.group.by = "tissue_type", draw.lines = F) +
  scale_fill_viridis()

write.csv(top_markers_T, file = "top_markers_T.csv")
ggsave2("T_cell_marker_gene.pdf", path = "./result_3v3/", width = 10, height = 9)
