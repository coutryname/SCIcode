library(monocle3)
library(scales)
library(ggsci)
#show_col(pal_aaas("default")(7))  #图一
#show_col(pal_aaas("default", alpha = 0.6)(7))   #图二
#show_col(pal_aaas("default", alpha = 0.2)(4))   #图三
library(Seurat)
library(ggplot2)
#加载数据
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
#InstallData('pbmc3k')
#pbmc3k = LoadData("pbmc3k")
##Seurat版本更迭，需要更新一下才能正常读取
#pbmc3k = UpdateSeuratObject(object = pbmc3k)
#pbmc3k <- NormalizeData(pbmc3k, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc3k <- FindVariableFeatures(pbmc3k, selection.method = "vst", nfeatures = 2000) 

#pbmc3k <- ScaleData(pbmc3k, features = VariableFeatures(pbmc3k))
#pbmc3k <- RunPCA(pbmc3k, features = VariableFeatures(pbmc3k)) 

#ElbowPlot(pbmc3k, ndims=50, reduction="pca") 

#pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10) 
#pbmc3k <- FindClusters(pbmc3k, resolution = 0.5)

#pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)

af=readRDS("scRNA_tcell.rds")
##随机抽取10%个细胞
N=length(colnames(af))/4
N=round(N)
scRNA_tpm<-af[,sample(x=colnames(af),size = N,replace=F)]
af<- scRNA_tpm
##查看一下数据自带的细胞类型注释信息
colnames(pbmc3k@meta.data)
#[1] "orig.ident"         "nCount_RNA"         "nFeature_RNA"       "seurat_annotations"
#[5] "RNA_snn_res.0.5"    "seurat_clusters"  
DimPlot(pbmc3k, reduction = "umap",group.by = 'cell_type_imm',label = T,
        cols = c(pal_jco("default", alpha = 0.6)(9)))
expression_matrix = pbmc3k@assays$RNA@data
cell_metadata = data.frame(pbmc3k@meta.data)
gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
colnames(gene_annotation)=c("gene_short_name")
##构建Monocle3 cds对象
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#预处理，相当于Seurat流程normalize过程
#counts数据选择norm_method = c("log")
#data数据选择norm_method = c("none")
cds <- preprocess_cds(cds, num_dim = 50,norm_method = c("none"))

#去除批次效应,多个样本时可以通过此方式去除批次
cds <- align_cds(cds, alignment_group = "orig.ident")

## 降维，默认是"Umap"方式
cds <- reduce_dimension(cds,cores=5)

## 聚类分群，分辨率调小，是为了让细胞是一群可以更好展示
cds <- cluster_cells(cds,resolution = 0.0000001)

## 拟时序
cds <- learn_graph(cds)

##选择特定细胞作为起点，代码比交互式页面方便许多，且不会出现选中自己不想要的细胞
##这里我们假定以B细胞为起点
myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes}

cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'cell_type_imm',my_select = "NKT_cell"))

##使用Seurat的UMAP信息，这样可以与Seurat对象的细胞分布保持一致
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA_harmony_NKT_cell, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

##不同细胞类型拟时序数值
##拟时序值越高表示细胞分化程度越高，这里仅为演示，并非真实分化情况
plot_cells(cds, color_cells_by = "pseudotime",
           show_trajectory_graph=T) + plot_cells(cds,
                                                 color_cells_by = "cell_type_imm",
                                                 label_cell_groups=FALSE,
                                                 label_leaves=FALSE,
                                                 label_branch_points=FALSE,
                                                 graph_label_size=1)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(9))
#差异基因展示
Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=6)
#按莫兰指数选择TOP基因
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by="cell_type_imm",min_expr=0.5, ncol= 2,cell_size=1.5) + scale_color_manual(values = pal_jco("default", alpha = 0.6)(9)) 
##当然我们也可以选择我们感兴趣的基因进行展示
genes <- c('CD69','CITED2','CKLF','PTTG1','S100A10','SLC38A1','STMN1','RTKN2','UBE2S')
genes_cds <- cds[rowData(cds)$gene_short_name %in% genes, ]
plot_genes_in_pseudotime(genes_cds,color_cells_by="cell_type_imm",min_expr=0.5,ncol= 2,cell_size=1.5)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(9))
#monocle2拟时序_通路得分
