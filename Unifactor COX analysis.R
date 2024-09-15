td <- read.table("counts.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv <- group
group <-surv 
rownames(surv) <- surv$X
surv <- surv[,-(1:26)]
comgene <- intersect(rownames(a),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
a <- a[comgene,]
td <- surv
pFilter=0.05
outResult=data.frame()
sigGenes=c("OS","OS.time")
for(i in colnames(td[,3:ncol(td)])){
  tdcox <- coxph(Surv(OS.time, OS) ~ td[,i], data = td)
  tdcoxSummary = summary(tdcox)
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"]
  if(pvalue<pFilter){
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.table(outResult,file="UniCoxSurvival1.txt",sep="\t",row.names=F,quote=F)
UniCoxSurSigGeneExp=td[,sigGenes]
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)
write.table(UniCoxSurSigGeneExp,file="UniCoxSurSigGeneExp1.txt",sep="\t",row.names=F,quote=F)
tducs <- read.table("UniCoxSurvival1.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(tducs)
hr <- sprintf("%.3f",tducs$"HR")
hrLow  <- sprintf("%.3f",tducs$"L95CI")
hrHigh <- sprintf("%.3f",tducs$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tducs$pvalue<0.001, "<0.001", sprintf("%.3f", tducs$pvalue))


pdf(file="UniCoxSurForestPlot1.pdf", width = 7,height =5)
n <- nrow(tducs)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(2.5,2))


xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.2-0.5*0.2,n:1,pValue,adj=1,cex=text.cex);text(1.2-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)
axis(1)


dev.off()
