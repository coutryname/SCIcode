#设置工作目录
setwd("GSE185839")
library(dplyr)
library(Seurat)
library(patchwork)
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "GSE185839_RAW")
pbmc <-scRNA_harmony
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "scRNAlist2", min.cells = 3, min.features = 200)
pbmc

#使用pbmc.data函数可以读取表达矩阵，这里我们简单看一下这三个基因在表达矩阵中的样子
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
#可以看到下面的是稀疏矩阵的形式
#然后我们简单探索一下如果不是稀疏矩阵，运算速度会是怎样的？
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
## 709591472 bytes
sparse.size <- object.size(pbmc.data)
sparse.size
## 29905192 bytes
dense.size/sparse.size
## 23.7 bytes
#这里可以很清楚的看到稀疏矩阵确实可以明显的帮助我们减少内存的消耗

#上面构建完Seurat对象后，我们接下来第一步就是对我们的数据进行质控
#首先我们计算一下线粒体Gene的比例（这里可以思考一下，为什么我们要计算线粒体基因的比例）
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#计算完线粒体基因的比例后，我们可视化一下
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##第一张图为不同样本上基因的表达数量，
#其中一个点代表一个细胞；第二张图为不同细胞RNA表达的数量；
#第三个图则代表线粒体的百分比，如果线粒体基因占比过高则说明有问题，
#细胞要么是坏死要么就是特殊细胞（心肌细胞或者神经细胞）根据需求往往是10%一下，
#但是根据特殊情况25%以下也有可能，所以需要我们认清楚我们的样本来源~
#（总而言之，这个可视化的目的就是为了评估我们的样本，也就是这些细胞的质量如何）

#上面这么看可能不够直观，下面再进行相关性的可视化
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#第一张图，是样本中RNA的表达数量和线粒体基因的比例，
#两者之间没什么关系，这个才是正常的样本，因为线粒体基因的数量大约是300多个，
#而且十分稳定，基本上细胞坏死，RNA降解掉，线粒体基因的数量都不会大幅度改变，
#所以如果RNA数量降低，而线粒体基因升高，则提示细胞可能存在坏死；
#第二章图，基因的数量和RNA的数量呈正相关，这个也是很好理解的，
#所以根据这个可视化评估，我们的细胞是正常且合格的
#（其实这个图就是看趋势，上面的那三个小提琴图则是看数量，一般来说线粒体基因比例在10%一下，
#然后基因的数量和RNA的数量只要大部分比较密集即可，因为这个根据样本的变化，
#这两个数值也会有不同的变化，比如说有的要求基因的数量是4000左右，
#RNA的数量为20000-25000左右，和这个不一样则只能说明样本的不同）

##然后这里我们设置如下标准，作为我们的质控指标
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#当然这里我们也可以简化代码，变成以下的样子
pbmc <- NormalizeData(pbmc)
#然后接下来我们就要寻找高可变基因了，关于高可变基因的解读我们后面会进行解答
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#这里默认寻找2000个高可变基因，当然自己也可以人为设置，一般是1000-5000之间
# 这里我们展示一下前10个改变最大的基因
top10 <- head(VariableFeatures(pbmc), 10)
# 对高可变基因进行可视化
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#下面可视化的解读则是，左边和右边的图表达的是一个东西，
#只不过右边的图把表达差异在各个样本中最显著的基因表示了出来，然后高可变基因用红色表示

#前面我们提供如果对数据进行了挑选和改变就要对数据进行标准化来平衡数据之间的改动~
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#这里默认对所有基因进行线性归一化，这一步如果没有服务器往往运行比较慢，
#我们也可以使用其默认形式，只对高可变基因进行数据标准化的操作
pbmc <- ScaleData(pbmc)#默认形式

#接下来我们就要进行PCA降维，具体解释我们会在后续解答~
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#然后针对降维的结果表示为PC的形式，针对展示PC的可视化这里也提供了两种，但是为了教程的连贯性这里就不展现了，感兴趣的同学可以照着运行即可得到可视化结果
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
#当然一个PC内包含的基因，我们可以用热图的形式来展现出来
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
dev.off()
#接下来我们既然已经对庞大的数据进行了降维（也就是聚堆）的形式，那么我们究竟要选择几个PC来代表这么庞大的数据呢，肯定不能都选，否则我的数据量还是这些，就没有我们前面一直渗透的缩小缩小再缩小的含义了
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#上面两个代码是通过不同的方式来帮助我们选择PC的数目，并且分别都对应不同的可视化
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)

#既然我们确定了我们的PC数，也就是我们用这些PC来指代我们的整体
#为什么这里我们要确定降维的维数→可以简单理解为只有确定了这个维数，我们才可以知道后续要把哪些主成分进进行聚类
#下面我们开始进行细胞聚类，这里是把我们看似一个整体的大量细胞分出亚群
pbmc <- FindNeighbors(pbmc, dims = 1:11)#这里的参数10就是我们前面确定的维数
pbmc <- FindClusters(pbmc, resolution = 0.5)#这里的0.5决定了细胞亚群的多少，Seurat官网推荐0.4-1.2都可以进行尝试
#这里的resolution参数的选择有一个小技巧，比如我们一开始设置0.5，然后我们进行后续的操作，
#当我们tSNE可视化细胞亚群结果后，发现有的细胞亚群太大，
#也就是细胞太多，这个时候我们可以回过头来把这个参数调大点，这样来回修改操作,最后得到我们满意的tSNE可视化结果~
#让我们看一下细胞亚群分类的结果
head(Idents(pbmc), 5)
#前面既然已经把亚群分好了，但是光秃秃的数据肯定是没有感召力的，我们就要把我们的亚群进行可视化
#这里介绍两种方法， 一种是UMAP，一种是tSNE，方法没有好坏之分，文章中也可以都放上去，总的来说一句话，挑你看的顺眼而且结果好的就可以（Ps：一般来说UMAP运算速度快，tSNE比较常用）
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

#另一种可视化方法-tSNE
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne", label = TRUE)
save(pbmc,file = "pbmcy-8sample.Rdata")
#细胞注释我粗略的整理出来三种方法，可以供大家选择
#1.首先是依赖标记基因的细胞注释即该标记基因可以代表细胞亚群的功能和属性
#2.通过R包（SingleR包进行自动细胞注释）
#3.人工注释（当然前提你需要知道你研究的组织总各种细胞亚群的marker基因是什么，然后通过数据库查找判断细胞类型，这个方法其实包含在第一种之中）
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)#这个比较是基于cluster2和其它所有亚群进行对比
head(cluster2.markers, n = 5)
#接下来找一下cluster5和cluster0、1、2、3的marker
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
#寻找每一个cluster里与其它所有细胞相比之后的差异marker，并且每个亚群只展示FC最大的两位
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

VlnPlot(pbmc, features = c("BHLHE40", "  CD27  "))#改基因名就行
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))
top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
#手动注释
new.cluster.ids <- c("Endothelial cells (LEC.CFA)", " iPS_cells:PDB_1lox-17Puro-5", "DC:monocyte-derived:immature", "Neutrophil:inflam", "Macrophages (MF.103-11B+)", "B cells", 
                     "Monocyte:CD16+", "Endothelial cells (LEC)", "Epithelial cells (Ep.8wk.MEClo)", "Endothelial cells (LEC)", "Epithelial cells (Ep.5wk.MEC.Sca1+)", "Tgd (Tgd.VG3+24AHI)", "Tgd (Tgd.VG4+24AHI)", "Stromal cells"
                     , "Epithelial cells (Ep.5wk.MEC.Sca1+)")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
saveRDS(pbmc, file = "pbmc3k_final.rds")
