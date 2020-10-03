###尝试在不同分类层次数据使用Desep包来做差异分析
library(DESeq2)
library(limma)
library(pasilla)
library(DESeq)
#清空内存
rm(list=ls()) 
#统计差异otu
# 读入实验设计
design = read.table("map_HG_kangbing_R.txt", header=T, row.names= 1, sep="\t") 
head(design)
# 读取OTU表,全部otu表没有抽平
otu_table = read.delim("gg135_otu_table.txt", row.names= 1, sep="\t",header=T,check.names=F)
head(otu_table)
dim(otu_table)
idx = rownames(design) %in% colnames(otu_table) 
design = design[idx,]
sub_design = design[idx,]
count = otu_table[, rownames(sub_design)]
head(count)
taxonomy = read.delim("rep_seqs_tax.txt", row.names= 1,header=F, sep="\t")
head(taxonomy)
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")
library(dplyr)
taxonomy$id=rownames(taxonomy)
head(taxonomy)
tax = taxonomy[row.names(count),]
head(tax)
dim(tax)
dim(tax)
# 手动筛选显著的组
res = count[rownames(tax), ] # reorder according to tax
head(res)
dim(res)
res$phylum = gsub("p__","",tax$phylum,perl=TRUE) 
res$class = gsub("c__","",tax$class,perl=TRUE) 
res$order = gsub("o__","",tax$order,perl=TRUE) 
res$family = gsub("f__","",tax$family,perl=TRUE) 
res$genus = gsub("g__","",tax$genus,perl=TRUE) 
res$species = gsub("s__","",tax$species,perl=TRUE) 
head(res)
######下面就可以开始合并不同门类的数据了，首先我们按照科水平进行合并#######
# 按照Species分组，计算Sepal.Length的平均值和标准差
res1=res[c(1:12,16)]
head(res1)

res1[res1==" "] <- "unknow"
res1$family
# is.na(res1$family) <- "unknow"
library(dplyr)
iris_groups<- group_by(res1, family)
df_summarise<- summarise(iris_groups, sum(`B80-1.fq`),
                         sum(`B80-2.fq`),
                         sum(`B80-3.fq`),
                         sum(`B80-4.fq`),
                         sum(`B80-5.fq`),
                         sum(`B80-6.fq`),
                         sum(`LF-1.fq`),
                         sum(`LF-2.fq`),
                         sum(`LF-3.fq`),
                         sum(`LF-4.fq`),
                         sum(`LF-5.fq`),
                         sum(`LF-6.fq`)
)
df_summarise<- as.data.frame(df_summarise)
dim(df_summarise)
head(df_summarise)
row.names(df_summarise)=paste("family",1:291,sep="")
colnames(df_summarise)=c("family",colnames(count))


#########现在来使用Desep来做差异分析

count=as.matrix(df_summarise[2:13])
head(count)

dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = sub_design,
                              design = ~ SampleType)

dds2 <- DESeq(dds)  ##第二步,标准化
resultsNames(dds2)
# 将结果用results()函数来获取，赋值给res变量
diff <-  results(dds2, contrast=c("SampleType","CSF", "CRF"),alpha=0.05)
# summary一下，看一下结果的概要信息
summary(diff)
# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
WT <-subset(diff,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
WT
dim(WT)
diff$level = as.factor(ifelse(diff$padj < 0.05 & diff$log2FoldChange > 1, "enriched",ifelse(diff$padj < 0.05 & diff$log2FoldChange < -1, "depleted","nosig")))
dim(diff)
head(diff)

# 转换原始数据为百分比
norm = t(t(count)/colSums(count,na=T))# * 100 # normalization to total 100
head(norm)
norm=as.data.frame(norm)
normCSF=norm[1:6]
head(normCSF)
normCSF$meanCSF=apply(normCSF,1,mean)
###
normCRF=norm[7:12]
head(normCRF)
normCSF$meanCRF=apply(normCRF,1,mean)
diff1=as.data.frame(diff)
head(diff1)
dim(res1)
#str()
index = merge(normCSF ,diff1, by="row.names",all=F)
head(index)
dim(index)
row.names(index)=index$Row.names
index[grep(".fq|Row.names",colnames(index))]<-NULL
head(index)
dim(index)
#####现在我们来匹配其他分类等级的注释信息，思路就是对对family去除重复然后在进行匹配合并#########
res_tax=res[13:18]
head(res_tax)
res_family=distinct(res_tax, family , .keep_all = TRUE)  
head(res_family)
dim(res_family)

family_tax=merge(df_summarise,res_family,by = intersect(names(df_summarise)[1],names(res_family)[4]))
dim(family_tax)
head(family_tax)
row.names(family_tax)=row.names(df_summarise)
dim(family_tax)
taxx=family_tax[c(1,14:18)]
head(taxx)
#####现在我们来匹配其他分类等级的注释信息，思路就是对对family去除重复然后在进行匹配合并#########
zuizhong = merge(index,taxx, by="row.names",all=F)
head(zuizhong)
row.names(zuizhong)=zuizhong$Row.names
zuizhong$Row.names=NULL

write.table(zuizhong,"科水平差异统计表格CSF-CRF_DESeq2.txt",quote = FALSE,row.names = T,
            col.names = T,sep = "\t")

#------出图

design = read.table("map_HG_kangbing_R.txt", header=T, row.names= 1, sep="\t") 
head(design)


L6 = read.table("科水平差异统计表格CSF-CRF_DESeq2.txt", header=T, row.names= 1,  sep="\t") 
head(L6)

L6 = 

# L6$genus = gsub("g__","",L6$X.6,perl=TRUE) 
row.names(L6)=L6$X.6
wt=L6[7:18]
head(wt)
wt=t(wt)
str(wt)
wt=as.data.frame(wt)
wt$id2=design$SampleType
library (reshape2)
fengdu <- melt (wt, id="id2")
head(fengdu)
library("ggplot2") 
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")

#使用分面箱线图出图不太好看
##来做分面箱线图#########
p=ggplot(fengdu, aes(x = id2, y = value)) + 
  geom_boxplot(alpha=1, outlier.size=1, size=1, width=0.5,notchwidth=1,aes(fill=id2)) +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="", y="Relative abundance")+
  theme_bw()+
  scale_fill_manual(values = mi)+
  theme_light()+theme(axis.text.x = element_text(colour = "black",size = 14),
                      axis.text.y = element_text(colour = "black",size = 14))


p=p+facet_wrap(~variable,ncol = 10)
p
##来做分面箱线图#########
head(fengdu)
fengdu$genus = gsub("g__","",fengdu$variable,perl=TRUE) 
fengdu$value<-sqrt(fengdu$value)
fengdu$value<-sqrt(fengdu$value)
head(fengdu)
#,notch = T#无法添加槽口,position = position_dodge(0.7)
p=ggplot(fengdu, aes(x = genus, y =value,fill=id2)) + 
  geom_boxplot(alpha=1, outlier.size=.5, size=0.8, width=0.5,notchwidth=1,aes(fill=id2)) +  
  geom_jitter(position=position_jitter(0.15), size=1, alpha=0.7,aes(fill=id2))+
  labs(x="", y="Relative abundance")+
  scale_fill_manual(values = mi)+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))

p
p=p+scale_x_discrete(limits=L6$genus)+theme_classic()+
  theme(axis.text.x = element_text(colour = "black",size = 14,angle = 60,hjust = 1,vjust = 1),
        axis.text.y = element_text(colour = "black",size = 14)
  )
p

ggsave("属水平CSF富集细菌平方开平方数据.pdf", p, width = 10, height = 6)

################科水平出图##########
setwd("E:/Shared_Folder/HG_kangbing/nobac_noqianheti_chuli")
design = read.table("map_HG_kangbing_R.txt", header=T, row.names= 1, sep="\t") 
head(design)

setwd("E:/Shared_Folder/HG_kangbing/nobac_noqianheti_chuli/gg13.5_result/taxa_summary")
L5 = read.table("出图科水平L5水平差异CSF_CRF比较测试t检验.txt", header=T, row.names= 1,  sep="\t") 
head(L5)
L5$genus = gsub("f__","",L5$X.4,perl=TRUE) 
row.names(L5)=L5$X.4
wt=L5[3:14]
head(wt)
wt=t(wt)
str(wt)
wt=as.data.frame(wt)
wt$id2=design$SampleType
library (reshape2)
fengdu <- melt (wt, id="id2")
head(fengdu)
library("ggplot2") 
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")

#使用分面箱线图出图不太好看
##来做分面箱线图#########
p=ggplot(fengdu, aes(x = id2, y = value)) + 
  geom_boxplot(alpha=1, outlier.size=1, size=1, width=0.5,notchwidth=1,aes(fill=id2)) +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="", y="Relative abundance")+
  theme_bw()+
  scale_fill_manual(values = mi)+
  theme_light()+theme(axis.text.x = element_text(colour = "black",size = 14),
                      axis.text.y = element_text(colour = "black",size = 14))


p=p+facet_wrap(~variable,ncol = 10)
p
##来做分面箱线图#########
head(fengdu)
fengdu$genus = gsub("f__","",fengdu$variable,perl=TRUE) 
#做一个数据转化
#,notch = T#无法添加槽口
fengdu$value<-sqrt(fengdu$value)
fengdu$value<-sqrt(fengdu$value)
head(fengdu)

p=ggplot(fengdu, aes(x = genus, y =value,fill=id2)) + 
  geom_boxplot(alpha=1, outlier.size=.5, size=0.8, width=0.5,notchwidth=1,aes(fill=id2)) +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="", y="Relative abundance")+
  scale_fill_manual(values = mi)+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))

#scale_y_continuous(breaks = c(0,0.05,0.1,0.15,0.2,0.25,0.3))
#ylim(0,0.26)


p
head(L6)
p=p+scale_x_discrete(limits=L5$genus)+theme_classic()+
  theme(axis.text.x = element_text(colour = "black",size = 14,angle = 60,hjust = 1,vjust = 1),
        axis.text.y = element_text(colour = "black",size = 14)
  )
p
ggsave("科水平CSF富集细菌.pdf", p, width = 10, height = 6)

ggsave("科水平CSF富集细菌开方再开方.pdf", p, width = 10, height = 6)







