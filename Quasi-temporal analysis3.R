#安装需要的R包
#重新安装获得此功能，如果以前没有安装可以忽略，直接安装就可以
#install.packages("devtools")
#devtools::install_github("junjunlab/ClusterGVis")
#BiocManager::install(“org.Hs.eg.db”)
#install.packages(“ggplot2”)
########----滚打包收，差异、聚类、拟时序和富集分析----####
## 加载R包，读取数据
library(Seurat)
library(ClusterGVis)
library(org.Hs.eg.db)
library(ggplot2)
af=readRDS("scRNA_tcell.rds")
##随机抽取10%个细胞
N=length(colnames(af))/4
N=round(N)
scobj<-af[,sample(x=colnames(af),size = N,replace=F)]
##注释结果可视化
DimPlot(scobj, reduction = "umap", group.by ="cell_type_imm", label = T,raster = T)
new.cluster.ids <- c("CD4+T", "CD8+T", "Treg")
names(new.cluster.ids) <- levels(scobj)
scobj <- RenameIdents(scobj, new.cluster.ids)

#寻找marker基因
pbmc.markers.all <- Seurat::FindAllMarkers(scobj,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)

#挑选表达量前5的marker基因
pbmc.markers <- pbmc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC)

# [1] 以细胞群的均值作为输入,且默认对基因进行去重
st.data1 <- prepareDataFromscRNA(object = scobj,
                                 diffData = pbmc.markers,
                                 showAverage = TRUE)
# 富集分析，取TOP通路
# 这里使用该包自带的GO富集方法
enrich <- enrichCluster(object = st.data1,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 3,
                        seed = 5201314)
head(enrich, 3)
#挑选需要展示的marker基因
markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,
                                             replace = T)]

#绘制cluster基因表达折线图
visCluster(object = st.data1,
           plot.type = "line")

#绘制热图
visCluster(object = st.data1,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:4))

#绘制热图并添加富集注释和分组折线图
pdf('tmp.pdf',height = 10,width = 12,onefile = F)
visCluster(object = st.data1,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:4),
           col = rep(jjAnno::useMyCol("stallion",n = 9),each = 5),
           add.bar = T)
dev.off()

pdf('tmp1.pdf',height = 10,width = 14,onefile = F)
visCluster(object = st.data1,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:4),
          col = rep(jjAnno::useMyCol("stallion",n = 9),each = 5),
           add.bar = T)
dev.off()
