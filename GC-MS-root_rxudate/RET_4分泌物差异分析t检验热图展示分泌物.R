
library("gplots")
library("RColorBrewer")
library("ggplot2")
library("vegan")
#清空内存
rm(list=ls()) 
#########下面来分类统计分泌物的信息
##############全部峰值数据来做统计#########去除磷酸之后进行分析
# 读入mapping文件
design = read.table("RET_map_R.txt", header=T, row.names= 1, sep="\t") 
# 读取OTU表，这里我选择的是整个otu表格，但是一般没有必要全部做差异的啊，相对丰度高的做做就可以了
otu_table = read.delim("RET_a3筛选有化学式的分泌物并做一个分析.txt", row.names= 1, sep="\t",header=T,check.names=F)
head(otu_table)
idx = rownames(design) %in% colnames(otu_table) 
sub_design = design[idx,]
count = otu_table[, rownames(sub_design)]
head(count)
# 转换原始数据为百分比，
norm = t(t(count)/colSums(count,na=T)) * 100 # normalization to total 100
head(norm)
norm=as.data.frame(norm)

a=norm
head(a)
#预生成2个长度与输入文件行数相同的全为0的向量，将用于存储p value和差异倍数（log2FC）
Pvalue<-c(rep(0,nrow(a)))
fdr<-c(rep(0,nrow(a)))
log2_FC<-c(rep(0,nrow(a)))
head(a)
a[is.na(a)] <- 0
###########开始运行脚本
for(i in 1:nrow(a)){
  if(sd(a[i,1:3])==0&&sd(a[i,4:6])==0){
    Pvalue[i] <-"NA"
    log2_FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(a[i,1:3]),as.numeric(a[i,4:6]))
    Pvalue[i]<-y$p.value
    log2_FC[i]<-log2((mean(as.numeric(a[i,1:3]))+0.001)/(mean(as.numeric(a[i,4:6]))+0.001)) 
    fdr[i]=p.adjust(Pvalue[i], "BH") 
  }
}
# 在原文件后面加入log2FC，p value和FDR,共3列；
out<-cbind(a,log2_FC,Pvalue,fdr)

out$tax=otu_table$classification
head(out)
dim(out)

WT <-subset(out,fdr < 0.05 )
dim(WT)
head(WT)

library(dplyr)# desc()表示从大到小
WT1<- arrange(WT, desc(tax))
row.names(WT1)=row.names(WT)
head(WT1)
write.table(out,file="RET_a4化合物t检验结果CSF1和CRF1.xls",quote=FALSE,sep="\t",row.names=T)
################这里直接使用R语言做批量t检验代码如下###################
###########开始开始做另外一组
for(i in 1:nrow(a)){
  if(sd(a[i,7:9])==0&&sd(a[i,10:12])==0){
    Pvalue[i] <-"NA"
    log2_FC[i]<-"NA"
  }else{
    y=t.test(as.numeric(a[i,7:9]),as.numeric(a[i,10:12]))
    Pvalue[i]<-y$p.value
    log2_FC[i]<-log2((mean(as.numeric(a[i,7:9]))+0.001)/(mean(as.numeric(a[i,10:12]))+0.001)) 
    fdr[i]=p.adjust(Pvalue[i], "BH") 
  }
}
# 在原文件后面加入log2FC，p value和FDR,共3列；
out1<-cbind(a,log2_FC,Pvalue,fdr)

out1$tax=otu_table$classification
head(out)
dim(out)

WT <-subset(out1,fdr < 0.05 )
dim(WT)
head(WT)

library(dplyr)# desc()表示从大到小
WT2<- arrange(WT, desc(tax))
row.names(WT2)=row.names(WT)
head(WT2)
write.table(out1,file="RET_a4化合物t检验结果CSF2和CRF2.xls",quote=FALSE,sep="\t",row.names=FALSE)
################这里直接使用R语言做批量t检验代码如下###################


wt=WT1[1:6]
head(wt)
######现在来做差异分析，得到差异分泌物
wt =as.data.frame(wt)
library(pheatmap)

#设置颜色梯度"#1B9E77"
color = colorRampPalette(c( "white", "#FFFFE5","#67001F"))(60)
wt2<--log10(wt+0.000001)
wt2<-sqrt(wt)

wt2<-sqrt(wt)
wt2[wt2>0.2]<-0.2


#我们支持分组,这里我随便造一个分组，作为纵向分组
annotation_row = data.frame(tax=WT1$tax)  
rownames(annotation_row) = rownames(wt2)
#我再造一个横向分组
annotation_col = data.frame(design$SampleType[1:6])  
rownames(annotation_col) = colnames(wt)
#install.packages("pheatmap")
#我们开始做分组


p=pheatmap(wt2,fontsize=6,cellwidth = 8, cellheight =3,cluster_rows = FALSE,
           cluster_cols = F,
           color = color,border_color = "white",
           #cutree_col = 2,
           #gaps_row = c(7,14,21),
           annotation_col = annotation_col,annotation_row = annotation_row,
           labels_row = NULL,labels_col = NULL)

ggsave("a7_热图重新CSF1和CRF1差异.pdf", p, width = 10, height = 10)

##########
#########################
wt=WT2[1:6]
head(wt)
######现在来做差异分析，得到差异分泌物
wt =as.data.frame(wt)
library(pheatmap)

#设置颜色梯度"#1B9E77"
color = colorRampPalette(c( "white", "#FFFFE5","#67001F"))(60)
wt2<--log10(wt+0.000001)
wt2<-sqrt(wt)

wt2<-sqrt(wt)
wt2[wt2>0.2]<-0.2


#我们支持分组,这里我随便造一个分组，作为纵向分组
annotation_row = data.frame(tax=WT2$tax)  
rownames(annotation_row) = rownames(wt2)
#我再造一个横向分组
annotation_col = data.frame(design$SampleType[7:12])  
rownames(annotation_col) = colnames(wt)
#install.packages("pheatmap")
#我们开始做分组


p=pheatmap(wt2,fontsize=6,cellwidth = 8, cellheight =3,cluster_rows = FALSE,
           cluster_cols = F,
           color = color,border_color = "white",
           #cutree_col = 2,
           #gaps_row = c(7,14,21),
           annotation_col = annotation_col,annotation_row = annotation_row,
           labels_row = NULL,labels_col = NULL)




ggsave("a7_热图重新CSF1和CRF1差异.pdf", p, width = 10, height = 10)




