#############对微生物数据使用随机森林做一个生物标记物选择###########
#更改目录
# 读入实验设计
design = read.table("RET_map_R.txt", header=T, row.names= 1, sep="\t") 
head(design)
# 读取OTU表,全部otu表没有抽平
otu_table = read.delim("RET_a3筛选有化学式的分泌物并做一个分析.txt", row.names= 1, sep="\t",header=T,check.names=F)
head(otu_table)
dim(otu_table)
#数据整理形式一######
idx = rownames(design) %in% colnames(otu_table) 
design = design[idx,]
sub_design = design[idx,]
count = otu_table[, rownames(sub_design)]
head(count)
norm = t(t(count)/colSums(count,na=T))# * 100 # normalization to total 100
head(norm)


#数据整理形式一######
# 加载随机森林包
library(randomForest)
#######使用随机森林做分类########
set.seed(315)
iris.rf = randomForest(t(norm), design$SampleType, importance=TRUE, proximity=TRUE,ntree=1000)
print(iris.rf)
#######使用随机森林做分类########


###提取分类变量使用MeanDecreaseAccuracy作为重要性指标选择变量#######
a=round(importance(iris.rf), 2)
head(a)
index11 = merge(norm ,a, by="row.names",all=F)
head(index11)
dim(index11)
write.table(index11,"随机森林分类重要变量输出_norm.txt",quote = FALSE,row.names = T,
            col.names = T,sep = "\t")


varImpPlot(iris.rf)  

########这里做变量权重图的表格######
head(a)
library(dplyr)
str(a)
a=as.data.frame(a)
a$id=row.names(a)
a2<- arrange(a, desc(MeanDecreaseAccuracy))
head(a2)
row.names(a2)=a2$id
a3=head(a2,n=35L)
head(a3)
########这里做变量权重图的表格######

#######开始出图，做火柴图########
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
library("ggplot2") 
p=ggplot(a3, aes(x = MeanDecreaseAccuracy, y = reorder(id,MeanDecreaseAccuracy)))  +
  geom_segment(aes(yend=id),xend=0,size=3,colour = "#1B9E77" )
#geom_point(size=4,pch=20, colour = "#1B9E77")
p 

p = p+theme_bw()+theme(axis.text.x = element_text(colour = "black",size = 20,face = "bold"),
                       axis.text.y = element_text(colour = "black",size = 10,face = "bold"))
p
ggsave("a7_随机森林35个变量重要性火柴棒图.pdf", p, width = 4, height = 6)

