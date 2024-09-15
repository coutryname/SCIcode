td=read.table("UniCoxSurSigGeneExp.txt",header=T,sep="\t",row.names=1)   #读入td数据       
clinicData=read.table("Clinical.txt",sep="\t",header=T,row.names=1)    #读入临床数据
sameSample=intersect(row.names(clinicData),row.names(td)) #临床数据可能和td数据样品不完全一致，将样品名取个交集，获得两组数据中都有的样品
td=td[sameSample,c("surtime","surstat","gene39","gene87")] #td数据取子集，我们只需要"surtime","surstat","gene39","gene87"就可以了
clinicData=clinicData[sameSample,] #按照相同样品取临床数据子集
tdcl=cbind(td,clinicData) 

#####
td <- surv
td <- tdcl # 偷懒不想改下面代码，遂把tdcl还是命名为td
tdmultiCox=coxph(Surv(OS.time, OS) ~ ., data = td) #这时候结果已经有了
tdmultiCox #直接查看

tdmultiCoxSum=summary(tdmultiCox)
outResult=data.frame()
outResult=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CIH=tdmultiCoxSum$conf.int[,"upper .95"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"])
outResult=cbind(id=row.names(outResult),outResult)
write.table(outResult,file="multiCoxClinical3.txt",sep="\t",row.names=F,quote=F)
######
tdmcs <- read.table("multiCoxclinical3.txt",header=T,sep="\t",row.names=1)
gene <- rownames(tdmcs)
hr <- sprintf("%.3f",tdmcs$"HR")
hrLow  <- sprintf("%.3f",tdmcs$"L95CI")
hrHigh <- sprintf("%.3f",tdmcs$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tdmcs$pvalue<0.001, "<0.001", sprintf("%.3f", tdmcs$pvalue))


pdf(file="multiCoxSurForestPlot3.pdf", width = 6,height = 5)
n <- nrow(tdmcs)
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
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)


dev.off()

