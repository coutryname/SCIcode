
library(monocle)
library(cowplot)
cds=readRDS("cds.RDS")
af=readRDS("scRNA_tcell.rds")
af<- scRNA_tpm
##随机抽取10%个细胞
N=length(colnames(af))/4
N=round(N)
scRNA_tpm<-af[,sample(x=colnames(af),size = N,replace=F)]
scRNA_tpm.markers <- FindAllMarkers(scRNA_tpm, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
logFCfilter=1           
adjPvalFilter=0.05
scRNA_tpm.markers=scRNA_tpm.markers[(abs(as.numeric(as.vector(scRNA_tpm.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(scRNA_tpm.markers$p_val_adj))<adjPvalFilter),]
monocle.matrix=as.matrix(scRNA_tpm@assays$RNA@counts, 'sparseMatrix')
monocle.sample=scRNA_tpm@meta.data[,8,drop=F]
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
marker=scRNA_tpm.markers
#将Seurat结果转换为monocle
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
#给其中一列数据重命名
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#添加细胞聚类数据
pData(cds)$Type=scRNA_tpm$tissue_type
pData(cds)$cellType=scRNA_tpm$cell_type_imm
pData(cds)$Cluster=scRNA_tpm$SCT_snn_res.0.4
pData(cds)$T_cell=scRNA_tpm$cell_type_imm$T_cell
#伪时间分析流程
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#monocle选择高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
#plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)
cds <- setOrderingFilter(cds, marker$gene)
saveRDS(cds, file = "cds.RDS")
#plot_ordering_genes(cds)
pdf(file="cluster.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
dev.off()
pdf(file="af_N-T Type.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Type")
dev.off()
pdf(file="af_cellType.trajectory.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cellType")
dev.off()
pdf(file="af_trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()
#输出拟时序差异基因
groups=subset(pData(cds),select='State')
pbmc=AddMetaData(object=af, metadata=groups, col.name="group")
geneList=list()
for(i in levels(factor(groups$State))){
  pbmc.markers=FindMarkers(pbmc, ident.1 = i, group.by = 'group')
  sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
  sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
  write.table(sig.markers,file=paste0("05.monocleDiff.", i, ".txt"),sep="\t",row.names=F,quote=F)
  geneList[[i]]=row.names(sig.markers)
}
#保存交集基因
unionGenes=Reduce(union,geneList)
write.table(file="monocleDiff.union.txt",unionGenes,sep="\t",quote=F,col.names=F,row.names=F)
#得到这些基因后我们可以随便查看一个基因在时序中的位置

pdf(file = "GNLY_af_cds_geneselect.pdf",width =8,height = 7)
pData(cds)[,'GNLY'] = scRNA_tpm@assays$RNA@scale.data['GNLY',]
plot_cell_trajectory(cds, color_by = 'GNLY') + scale_color_viridis()()
dev.off()
#差异基因展示
Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=6)
#按莫兰指数选择TOP基因
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[scRNA_tpm.markers,],color_cells_by="cell_type",min_expr=0.5, ncol= 2,cell_size=1.5) + scale_color_manual(values = pal_jco("default", alpha = 0.6)(9)) 

genes <- c('CD19','CD4')
genes_cds <- cds[rowData(cds)$gene_short_name %in% genes, ]
plot_genes_in_pseudotime(genes_cds,color_by="cellType",min_expr=0.5,cell_size=1.5)+ scale_color_manual(values = pal_jco("default", alpha = 0.6)(9))
