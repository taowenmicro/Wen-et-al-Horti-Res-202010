#清空内存#######
rm(list=ls()) 

library("ggplot2")
library("vegan")

# 读入实验设计和Alpha多样性值
design = read.table("map_HG_kangbing_R.txt", header=T, row.names= 1, sep="\t") 
head(design)
bray_curtis = read.table("./beta/bray_curtis_CSS_otu_table.txt", sep="\t",row.names= 1,header=T, check.names=F)

#过滤数据并且排序
idx = rownames(design) %in% colnames(bray_curtis)
idx
sub_design = design[idx,]
sub_design
bray_curtis = bray_curtis[rownames(sub_design), rownames(sub_design)] 

pcoa = cmdscale(bray_curtis, k=2, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
#?cmdscale k=3数据在其中所表示的空间的最大维数;eig=T特征值是否应该归还
points = as.data.frame(pcoa$points) # 获得坐标点get coordinate string, format to dataframme
head(points)
colnames(points) = c("x", "y") #命名行名
eig = pcoa$eig
#eig特征值得到
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])
#write.table(points,"pcoa_bray_curtis.txt",quote = FALSE,row.names = F,
#col.names = T,sep = "\t")
points


#显著性检测
#做显著性分析ANOSIM相似性分析是一种非参数检验，用来检验组间（两组或多组）差异是否显著大于组内差异，从而判断分组是否有意义利用Bray-Curtis算法计算两两样品间的距离，然后将所有距离从小到大进行排序，并计算R和P值。
#Adonis又称置换多元方差分析或非参数多元方差分析。它利用半度量（如Bray-Curtis）或度量的距离矩阵（如Euclidean）对总方差进行分解，分析不同分组因素对样品差异的解释度，并使用置换检验对其统计学意义进行显著性分析
ado = adonis(bray_curtis ~ design$SampleType,permutations = 999,method="bray")
a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
R2 <- paste("adonis:R ",a, sep = "")
b = as.data.frame(ado$aov.tab[6])[1,1]
p = paste("p: ",b, sep = "")
title = paste(R2," ",p, sep = "")
################################
pairwise.adonis1 <-function(x,factors,p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] ,permutations = 999);
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}
pairwise.adonis1(bray_curtis,design$SampleType, p.adjust.m= "bonferroni")
#可选anosim检测
#anosim.result<-anosim(bray_curtis, design$SampleType,permutations = 999)
#summary(anosim.result)
head(points)

#定义颜色
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")

# 绘图
p <- ggplot(points, aes(x=x, y=y, color=SampleType)) +
  geom_point(colour="grey50", size =2)+
  geom_point(alpha=.7, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(plot.title = element_text(hjust = 0.5))
p
p = p+stat_ellipse( linetype = 2,level = 0.65)+stat_ellipse( linetype = 1,level = 0.8)+theme_bw()+
  scale_colour_manual(values = mi)+
  theme(
    axis.text = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    #legend.position="none",
    axis.title.y =element_text(size = 24,face = "bold"),
    axis.title.x =element_text(size = 24,face = "bold"),
    title = element_text(size = 24,face = "bold")
  )
p
ggsave("./beta/bray_curtis_beta.pdf", p, width = 10, height = 8)



#下面用其他两个距离矩阵做pcoa分析
weighted_unifrac = read.table("./a2_beta/beta/weighted_unifrac_CSS_otu_table.txt", sep="\t",row.names= 1,header=T, check.names=F)
#下面要对刚才导入的矩阵进行排序，按照design进行排序，因此，第一列必须设置成为列名
#过滤数据并且排序
idx = rownames(design) %in% colnames(weighted_unifrac)#前一个数据集的行名和都一个数据集的列名做匹配，输出逻辑变量
idx
sub_design = design[idx,]
sub_design
weighted_unifrac = weighted_unifrac[rownames(sub_design), rownames(sub_design)] 
pcoa = cmdscale(weighted_unifrac, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # 获得坐标点get coordinate string, format to dataframme
colnames(points) = c("x", "y") #命名行名
eig = pcoa$eig
#eig特征值得到
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])
#write.table(points,"pcoa_bray_curtis.txt",quote = FALSE,row.names = F,
#col.names = T,sep = "\t")
points


#显著性检测
#做显著性分析ANOSIM相似性分析是一种非参数检验，用来检验组间（两组或多组）差异是否显著大于组内差异，从而判断分组是否有意义利用Bray-Curtis算法计算两两样品间的距离，然后将所有距离从小到大进行排序，并计算R和P值。
#Adonis又称置换多元方差分析或非参数多元方差分析。它利用半度量（如Bray-Curtis）或度量的距离矩阵（如Euclidean）对总方差进行分解，分析不同分组因素对样品差异的解释度，并使用置换检验对其统计学意义进行显著性分析
ado = adonis(weighted_unifrac~ design$SampleType,permutations = 999,method="bray")
a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
R2 <- paste("anosis:R ",a, sep = "")
b = as.data.frame(ado$aov.tab[6])[1,1]
p = paste("p: ",b, sep = "")
title = paste(R2," ",p, sep = "")
title
################################
pairwise.adonis1 <-function(x,factors,p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] ,permutations = 999);
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}
pairwise.adonis1(bray_curtis,design$SampleType, p.adjust.m= "bonferroni")
#可选anosim检测
#anosim.result<-anosim(bray_curtis, design$SampleType,permutations = 999)
#summary(anosim.result)


# 绘图
p <- ggplot(points, aes(x=x, y=y, color=SampleType)) +
  geom_point(colour="grey50", size =2)+
  geom_point(alpha=.7, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(plot.title = element_text(hjust = 0.5))

p = p+stat_ellipse( linetype = 2,level = 0.65)+stat_ellipse( linetype = 1,level = 0.9)+theme_bw()+
  scale_colour_manual(values = mi)+
  theme(
    axis.text = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    #legend.position="none",
    axis.title.y =element_text(size = 24,face = "bold"),
    axis.title.x =element_text(size = 24,face = "bold"),
    title = element_text(size = 24,face = "bold")
  )
p
ggsave("E:/Shared_Folder/16s_and_RE_data_analyseHG_resistant/No1_16s/result_and_script/a2_beta/beta/weighted_unifrac_beta.pdf", p, width = 10, height = 8)






#下面用其他两个距离矩阵unweighted_unifrac做pcoa分析
unweighted_unifrac = read.table("./a2_beta/beta/unweighted_unifrac_CSS_otu_table.txt", sep="\t",row.names= 1,header=T, check.names=F)
idx = rownames(design) %in% colnames(unweighted_unifrac)#前一个数据集的行名和都一个数据集的列名做匹配，输出逻辑变量
idx
sub_design = design[idx,]
sub_design
unweighted_unifrac = unweighted_unifrac[rownames(sub_design), rownames(sub_design)] 

pcoa = cmdscale(unweighted_unifrac, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # 获得坐标点get coordinate string, format to dataframme
colnames(points) = c("x", "y") #命名行名
eig = pcoa$eig
#eig特征值得到
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)), ])
#write.table(points,"pcoa_bray_curtis.txt",quote = FALSE,row.names = F,
#col.names = T,sep = "\t")

#显著性检测
ado = adonis(unweighted_unifrac~ design$SampleType,permutations = 999,method="bray")
a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
R2 <- paste("anosis:R ",a, sep = "")
b = as.data.frame(ado$aov.tab[6])[1,1]
p = paste("p: ",b, sep = "")
title = paste(R2," ",p, sep = "")
title
################################
pairwise.adonis1 <-function(x,factors,p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] ,permutations = 999);
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}
pairwise.adonis1(bray_curtis,design$SampleType, p.adjust.m= "bonferroni")
#可选anosim检测
#anosim.result<-anosim(bray_curtis, design$SampleType,permutations = 999)
#summary(anosim.result)


# 绘图oint(colour="grey50", size =2)+
  geom_point(alpha=.7, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=title)+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(plot.title = element_text(hjust = 0.5))

p = p+stat_ellipse( linetype = 2,level = 0.65)+stat_ellipse( linetype = 1,level = 0.9)+theme_bw()+
  scale_colour_manual(values = mi)+
  theme(
    axis.text = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    #legend.position="none",
    axis.title.y =element_text(size = 24,face = "bold"),
    axis.title.x =element_text(size = 24,face = "bold"),
    title = element_text(size = 24,face = "bold")
  )
p
ggsave("E:/Shared_Folder/16s_and_RE_data_analyseHG_resistant/No1_16s/result_and_script/a2_beta/beta/unweighted_unifrac_beta.pdf", p, width = 10, height = 8)
