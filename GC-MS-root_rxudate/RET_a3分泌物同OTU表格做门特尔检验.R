#下面开始做两样本相关性检测用otu表格和分泌物数据做
fenmiwu = read.delim("RET_a3筛选有化学式的分泌物CSF1和CRF1做门特尔检验.txt", header=T, row.names= 1,check.names=F) 
head(fenmiwu)
dim(fenmiwu)
fenmiwu1=fenmiwu[1:12]
head(fenmiwu1)
str(fenmiwu1)
design = t(fenmiwu1)
library("vegan")
env.dist1 <- vegdist(scale(design), "euclid")


#?vegdist
env.dist1
otu_table = read.table("gg135_otu_table.txt", row.names= 1, header=T,check.names=F)
head(otu_table)
# 读入实验设计
design = read.table("map_HG_kangbing_R.txt", header=T, row.names= 1, sep="\t") 
head(design)
idx = rownames(design) %in% colnames(otu_table) 
design = design[idx,]
sub_design = design[idx,]
count = otu_table[, rownames(sub_design)]
head(count)

norm = t(t(count)/colSums(count,na=T))# * 100 # normalization to total 100
head(norm)
otu_table=count

otu_table =t(otu_table)
veg.dist1 <- vegdist(otu_table) # Bray-Curtis


veg.dist1
mantel(veg.dist1, env.dist1, method="spear")#

mantel(veg.dist1, env.dist1,method="pearson")#
#Mantel statistic r: 0.5966 
#Significance: 0.003 
mantel(veg.dist1, env.dist1, method="spearman")#
#Mantel statistic r: 0.5573 
#Significance: 0.005 
######到这里结束#########


