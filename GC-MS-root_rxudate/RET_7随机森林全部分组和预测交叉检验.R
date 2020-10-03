#############对微生物数据使用随机森林做一个生物标记物选择###########
#更改目录
# 读入实验设计
design = read.table("RET_map_R1.txt", header=T, row.names= 1, sep="\t") 
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
#######开始出图，做火柴图########

#######将挑选的OTU做一张热图展示##########
setwd("E:/Shared_Folder/HG_kangbing/nobac_noqianheti_chuli")
write.table(a3,"随机森林分类挑选41个.txt",quote = FALSE,row.names = T,
            col.names = T,sep = "\t")

# 转换原始数据为百分比
norm = t(t(count)/colSums(count,na=T))# * 100 # normalization to total 100
head(norm)
norm=as.data.frame(norm)
norm1 = norm[row.names(a3),]
wt=norm1

#wt2<--log10(wt+0.000001)
wt2<-sqrt(wt)

wt2<-sqrt(wt2)
wt2[wt2>0.25]<-0.25

#我们支持分组,这里我随便造一个分组，作为纵向分组
annotation_row = data.frame(tax=a3$genus)  
rownames(annotation_row) = rownames(wt2)
#我再造一个横向分组
annotation_col = data.frame(design$SampleType)  
rownames(annotation_col) = colnames(wt)
#install.packages("pheatmap")
#我们开始做分组
library(pheatmap)
#设置颜色梯度"#1B9E77"
color = colorRampPalette(c( "white","#FFFFE5","#67001F" ))(60)
pheatmap(wt2,fontsize=6,cellwidth = 8, cellheight = 8,cluster_rows = FALSE,
         cluster_cols = F,
         color = color ,border_color = "white",
         cutree_col = 2,
         #gaps_row = c(3,72,77,79,93,95),
         annotation_col = annotation_col,
         #annotation_row = annotation_row,
         labels_row = NULL,labels_col = NULL)

#######将挑选的OTU做一张热图展示##########


######################使用钢L3数据做一次随机森林##############
#############对微生物数据使用随机森林做一个生物标记物选择###########
#统计差异otu
setwd("E:/Shared_Folder/HG_kangbing/nobac_noqianheti_chuli/taxa_summary")

# 读入实验设计
design = read.table("map_HG_kangbing_R.txt", header=T, row.names= 1, sep="\t") 
head(design)
# 读取OTU表,全部otu表没有抽平
otu_table = read.delim("otu_table_tax_L4.txt", row.names= 1, sep="\t",header=T,check.names=F)
head(otu_table)
row.names(otu_table)=otu_table$id
otu_table$id=NULL
#
idx = rownames(design) %in% colnames(otu_table) 
design = design[idx,]
sub_design = design[idx,]
count = otu_table[, rownames(sub_design)]
head(count)

library(randomForest)
set.seed(315)
iris.rf = randomForest(t(count), design$SampleType, importance=TRUE, proximity=TRUE,ntree=1000)
print(iris.rf)

b=round(importance(iris.rf), 2)
write.table(b,"随机森林分类l4水平变量文件.txt",quote = FALSE,row.names = T,
            col.names = T,sep = "\t")
#将重要变量简单出图查看#
varImpPlot(iris.rf)  

## 交叉验证选择Features#############目的是选择多少个变量合适
set.seed(315) 
# rfcv是随机森林交叉验证函数：Random Forest Cross Validation
result = rfcv(t(count), design$SampleType, cv.fold=10)
result$error.cv
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
######

#############使用随机森林模型来预测分组###############
count1=t(count)
head(count1)
count1=as.data.frame(count1)
count1$id=factor(rep(c("CSF", "CRF"), c(6,6)))
#dim(count1)
# 随机1-2有放回抽样150次，概率分别为0.67和0.33，用于获取训练集
ind=sample(2,nrow(count1),replace=TRUE, prob=c(0.67,0.33))  
# 2/3作预测建模
iris.rf = randomForest(id ~ .,count1[ind==1,] , ntree=1000, 
                       nPerm=10,  proximity=TRUE, importance=TRUE)  
print(iris.rf)  
# 1/3验证预测模型
iris.pred = predict(iris.rf, count1[ind==2,] )  
# 输出预测与观测对应表
table(observed=count1[ind==2,"id"], predicted=iris.pred) 
#############使用随机森林模型来预测分组###############

##提取变量
imp= as.data.frame(iris.rf$importance)
head(imp)
#变量安装重要性降序排列.order函数可以对一个向量进行排序，这里选择第三列，就是表征变量重要性的那列
imp = imp[order(imp[,3],decreasing = T),]
head(imp)
##保存
write.table(imp,file = "importance_class.txt",quote = F,sep = '\t', row.names = T, col.names = T)
# 简单可视化
varImpPlot(iris.rf, main = "Top 23 - Feature OTU importance",n.var = 25, bg = par("bg"),
           color = par("fg"), gcolor = par("fg"), lcolor = "gray" )
###############
# 读取所有feature贡献度
imp = read.table("importance_class.txt", header=T, row.names= 1, sep="\t") 
# 交叉验证后我选择前14个变量
imp = head(imp, n=14)
head(imp)
# 反向排序X轴,为什么要
imp=imp[order(1:14,decreasing = T),]
head(imp)
# imp物种名分解
# 去除公共部分
imp$temp = gsub("k__Bacteria;p__","",rownames(imp),perl=TRUE) 
head(imp)
# 提取门名称,实际上是使用gsub进行查找替换
imp$phylum = gsub(";[\\w-\\[\\]_]+","",imp$temp,perl=TRUE) # rowname unallowed same name
head(imp)
imp$phylum = gsub("[\\[\\]]+","",imp$phylum,perl=TRUE) 
# 提取纲名称
imp$class = gsub("[\\w\\[\\];_]+;c__","",imp$temp,perl=TRUE)  
imp$class = gsub("[\\[\\]]+","",imp$class,perl=TRUE)
##提取目水平
imp$order = gsub("[\\w\\[\\];_]+;c__","",imp$class,perl=TRUE)  
##
imp$order = gsub("[\\w\\[\\];_]+;o__","",imp$order,perl=TRUE)  

head(imp)
# str(imp)
# 添加纲level，制作因子变量,
imp$order=factor(imp$order,levels = imp$order)

# 图1. 绘制物种类型种重要性柱状图

p=ggplot(data = imp, mapping = aes(x= order,y=MeanDecreaseAccuracy,fill=phylum)) + 
  geom_bar(stat="identity")+coord_flip()
p
ggsave(paste("rf_imp_feature",".pdf", sep=""), p, width = 4, height =4)

############交叉验证选择Features###########################





##########下面这部分是时间序列的，可以不做################
#我们先注意一下随机森林回归和随机森林分类的差别：（1）默认mtry是p/3而不是p1/2，其中p表示预测变量数
#（2）默认节点大小为5而不是1（3）只有一个测量变量的重要性。
#mtry参数表征默认在每个节点抽取的变量数
ozone.rf= randomForest(id ~ ., data=count1, mtry=568,
                       importance=TRUE, na.action=na.omit,ntree=1000)
print(ozone.rf)
head(round(importance(ozone.rf), 2))
varImpPlot(ozone.rf)  
# 先清空NA的样本，验证不允许有NA
airquality = na.omit(airquality)
myairquality= cbind(airquality[1:6], matrix(runif(96 * nrow(airquality)), nrow(airquality), 96))
# 交驻验证添加了随机数的训练集，分组，交叉验证的次数
#十折交叉验证：用来测试精度。是常用的精度测试方法。将数据集分成十分，
#轮流将其中9份做训练1份做测试，10次的结果的均值作为对算法精度的估计，
#一般还需要进行多次10倍交叉验证求均值
result= rfcv(myairquality, airquality$Ozone, cv.fold=3)

# 绘制错误率曲线，观查错误率与使用Markers数量的变化
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
# 使用replicate进行多次交叉验证，可选
result= replicate(5, rfcv(myairquality, airquality$Ozone), simplify=FALSE)
error.cv= sapply(result, "[[", "error.cv")
matplot(result[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type="l",
        lwd=c(2, rep(1, ncol(error.cv))), col=1, lty=1, log="x",
        xlab="Number of variables", ylab="CV Error")
##########下面这部分是时间序列的，可以不做