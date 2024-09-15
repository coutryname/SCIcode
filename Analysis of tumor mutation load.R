# 安装并加载所需的R包
library(TCGAbiolinks)
library(maftools)
query <- GDCquery(
  project = "TCGA-LIHC", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(query)

GDCprepare(query, save = T,save.filename = "TCGA-LIHC_SNP.Rdata") # 这里的Rdata文件是一个数据框，可直接用maftools读取使用

load(file = "TCGA-LIHC_SNP.Rdata")
maf.stad <- data

## 查看数据
class(maf.stad)
dim(maf.stad)
maf.stad[1:5, 1:10]
maf <- read.maf(maf.stad)
# 绘制MAF文件的整体结果图
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
# 绘制oncoplot图(只展示前 15 个基因)
vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
pdf('maf.pdf',height = 5,width = 5,onefile = F)
oncoplot(maf = maf, colors = vc_cols, top = 15)
dev.off()
# 使用 somaticInteractions 函数检测互斥突变或同时突变的基因
pdf('maf2.pdf',height = 14,width = 16,onefile = F)
somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))
dev.off()
# 绘制Transitions Transversions汇总图
pdf('maf4.pdf',height = 5,width = 5,onefile = F)
laml.titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
dev.off()
#2. 计算 TMB
stad.tmb <- tmb(maf, captureSize = 38, logScale = T)
dim(stad.tmb)
# [1] 371   4

head(stad.tmb)
# 根据TMB平均值进行分组
library(dplyr)
stad.tmb <- stad.tmb %>% mutate(group = if_else(total_perMB_log > mean(total_perMB_log), "TMB_high","TMB_low"))
# 只取前12位患者编号
stad.tmb$patient <- substr(stad.tmb$Tumor_Sample_Barcode, 1, 16)
stad.tmb$Tumor_Sample_Barcode <- gsub("-",".",stad.tmb$Tumor_Sample_Barcode)
write.table(stad.tmb, file = "stad.tmb.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#整理完毕#
# 加载自身临床数据
clin_info[1:4,1:4]
#   patient         entity_submitter_id status time
# 1 TCGA-BR-A4J4 TCGA-BR-A4J4-01A-12R-A251-31      0   16
# 2 TCGA-BR-A4IZ TCGA-BR-A4IZ-01A-32R-A251-31      1  273
# 3 TCGA-RD-A7C1 TCGA-RD-A7C1-01A-11R-A32D-31      1  507
# 4 TCGA-BR-6852 TCGA-BR-6852-01A-11R-1884-13      0 1367
TMB <- merge(clin, stad.tmb, by = "patient")

fit <- survfit(Surv(time, status) ~ group, data = TMB)

ggsurvplot(fit, data = TMB,
           pval = T,
           risk.table = T,
           surv.median.line = "hv", #添加中位生存曲线
           palette = c("red", "blue"),  #更改线的颜色
           legend.labs = c("High risk","Low risk"), #标签
           legend.title = "TMB Score",
           title = "Overall survival", #标题
           ylab = "Cumulative survival (percentage)", xlab = " Time (Days)", #更改横纵坐标
           censor.shape = 124, censor.size = 2, conf.int = FALSE, #删失点的形状和大小
           break.x.by = 500 #横坐标间隔
)