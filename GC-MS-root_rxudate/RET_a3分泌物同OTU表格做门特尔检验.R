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
?mantel
mantel(veg.dist1, env.dist1,method="pearson")#
#Mantel statistic r: 0.5966 
#Significance: 0.003 
mantel(veg.dist1, env.dist1, method="spearman")#
#Mantel statistic r: 0.5573 
#Significance: 0.005 
######到这里结束#########






mantel(veg.dist1, env.dist1, method="kendall")#
##########
env.dist1 <- vegdist(scale(otu_table), "euclidean")
env.dist1 <- vegdist(scale(otu_table), "manhattan")
env.dist1 <- vegdist(scale(otu_table), "gower")
env.dist1 <- vegdist(scale(otu_table), "canberra")
env.dist1 <- vegdist(scale(otu_table), "bray")
env.dist1 <- vegdist(scale(otu_table), "kulczynski")
env.dist1 <- vegdist(scale(otu_table), "morisita")
##########




#到此算制作完成
?vegdist()#微生物生态相异指数
?designdist似乎何意自己编辑相似性检验
#尝试用Estimating Morisita-Horn index方法做相似性检验
library("vegan")
setwd("E:/南农大硕士内容/袁军-试验/袁军作图")
Root_exudates = read.table("fmw.txt", header=T, row.names= 1,check.names=F) 
Root_exudates = t(Root_exudates)
Root_exudates
Root.dist1 <- vegdist(scale(Root_exudates), "euclid")
Root.dist1
geo_chip = read.delim("Geochip_data.txt", row.names= 1, sep="\t",header=T,check.names=F)
geo_chip =t(geo_chip)
head(geo_chip)
geo.dist1 <- vegdist(geo_chip) # Bray-Curtis
geo.dist1
#statistic r: 0.2561 
#Significance: 0.07 
geo.dist1 <- designdist(geo_chip, "(2*J)/(A+B)", terms = c("quadratic"))
geo.dist1
#statistic r: -0.2352 
#Significance: 0.928 
mantel(geo.dist1, Root.dist1, method="spearman")#

otu_table = read.delim("16S-OTU.txt", row.names= 1, sep="\t",header=T,check.names=F)
otu_table
otu_table =t(otu_table)
head(otu_table)
otu.dist1 <- vegdist(otu_table) # Bray-Curtis
otu.dist1
#statistic r: 0.7197 
#Significance: 0.003 
otu.dist1 <- designdist(otu_table, "(2*J)/(A+B)", terms = c("quadratic"))
otu.dist1
#statistic r: -0.8948 
#Significance: 1 
mantel(otu.dist1, Root.dist1, method="spearman")#

###################倪兄检验############
setwd("E:/Shared_Folder/NL")
design = read.table("map_R.txt", header=T, row.names= 1, sep="\t") 
head(design)
otu_table = read.table("otu_table20.txt",row.names= 1,header=T, sep="\t",check.names=F)
#过滤数据并且排序
idx = rownames(design) %in% colnames(otu_table)
idx
design = design[idx,]
design
otu_table = otu_table[rownames(design)] 
head(otu_table)
otu_table = t(otu_table)
otu_bray<-vegdist(otu_table, method='bray')
?vegdist
######环境数据
fangfa=c("euclidean")

envyl = design[5:17]
envyl

env_<-vegdist(envyl, method=fangfa)#euclidean
env_
o=mantel(otu_bray, env_, permutations=999)#做999次重复用来计算P值
r=o$statistic
p=o$signif

###########
envyl_Moisture = design[5]
envyl_Moisture

env_Moisture<-vegdist(envyl_Moisture, method=fangfa)#euclidean
env_Moisture
a=mantel(otu_bray, env_Moisture, permutations=999)#做999次重复用来计算P值
r1=a$statistic
p1=a$signif
#################
envyl_NO3.N = design[6]
envyl_NO3.N

env_NO3.N<-vegdist(envyl_NO3.N, method=fangfa)#euclidean
env_NO3.N
b=mantel(otu_bray, env_NO3.N, permutations=999)#做999次重复用来计算P值
r2=b$statistic
p2=b$signif
##################
envyl_NH4.N = design[7]
envyl_NH4.N

env_NH4.N<-vegdist(envyl_NH4.N, method=fangfa)#euclidean
env_NH4.N
c=mantel(otu_bray, env_NH4.N, permutations=999)#做999次重复用来计算P值
r3=c$statistic
p3=c$signif
#############

d=mantel(otu_bray, vegdist(design[8], method=fangfa), permutations=999)#做999次重复用来计算P值
r4=d$statistic
p4=d$signif
#################
e=mantel(otu_bray, vegdist(design[9], method=fangfa), permutations=999)#做999次重复用来计算P值
r5=e$statistic
p5=e$signif
###################
f=mantel(otu_bray, vegdist(design[10], method=fangfa), permutations=999)#做999次重复用来计算P值
r6=f$statistic
p6=f$signif
#############
g=mantel(otu_bray, vegdist(design[11], method=fangfa), permutations=999)#做999次重复用来计算P值
r7=g$statistic
p7=g$signif
#################
h=mantel(otu_bray, vegdist(design[12], method=fangfa), permutations=999)#做999次重复用来计算P值
r8=h$statistic
p8=h$signif
#################
i=mantel(otu_bray, vegdist(design[13], method=fangfa), permutations=999)#做999次重复用来计算P值
r9=i$statistic
p9=i$signif
#############
j=mantel(otu_bray, vegdist(design[14], method=fangfa), permutations=999)#做999次重复用来计算P值
r10=j$statistic
p10=j$signif
####################
k=mantel(otu_bray, vegdist(design[15], method=fangfa), permutations=999)#做999次重复用来计算P值
r11=k$statistic
p11=k$signif
################
l=mantel(otu_bray, vegdist(design[16], method=fangfa), permutations=999)#做999次重复用来计算P值
r12=l$statistic
p12=l$signif
##############
m=mantel(otu_bray, vegdist(design[17], method=fangfa), permutations=999)#做999次重复用来计算P值
r13=m$statistic
p13=m$signif
############
R=c(r,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13)
p=c(p,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13)
envyl

name=colnames(envyl)
name=c("quanbu",name)
web.data<-data.frame(name=name,R_otu=R,P_otu=p)
write.table(web.data,"Mantel testRNA.txt",quote = FALSE,row.names = F,
            col.names = T,sep = "\t")
#################################################
########################################
#################################
###########我们现在来做Mantel test
setwd("E:/Shared_Folder/NL")
design = read.table("map_R.txt", header=T, row.names= 1, sep="\t") 
head(design)
otu_table = read.table("otu_table_L2.txt",row.names= 1,header=T, sep="\t",check.names=F)
#过滤数据并且排序
idx = rownames(design) %in% colnames(otu_table)
idx
design = design[idx,]
design
otu_table = otu_table[rownames(design)] 
head(otu_table)
otu_table = t(otu_table)
otu_bray<-vegdist(otu_table, method='bray')

######环境数据
fangfa=c("euclidean")

envyl = design[5:17]
envyl

env_<-vegdist(envyl, method=fangfa)#euclidean
env_
o=mantel(otu_bray, env_, permutations=999)#做999次重复用来计算P值
r=o$statistic
p=o$signif

###########
envyl_Moisture = design[5]
envyl_Moisture

env_Moisture<-vegdist(envyl_Moisture, method=fangfa)#euclidean
env_Moisture
a=mantel(otu_bray, env_Moisture, permutations=999)#做999次重复用来计算P值
r1=a$statistic
p1=a$signif
#################
envyl_NO3.N = design[6]
envyl_NO3.N

env_NO3.N<-vegdist(envyl_NO3.N, method=fangfa)#euclidean
env_NO3.N
b=mantel(otu_bray, env_NO3.N, permutations=999)#做999次重复用来计算P值
r2=b$statistic
p2=b$signif
##################
envyl_NH4.N = design[7]
envyl_NH4.N

env_NH4.N<-vegdist(envyl_NH4.N, method=fangfa)#euclidean
env_NH4.N
c=mantel(otu_bray, env_NH4.N, permutations=999)#做999次重复用来计算P值
r3=c$statistic
p3=c$signif
#############

d=mantel(otu_bray, vegdist(design[8], method=fangfa), permutations=999)#做999次重复用来计算P值
r4=d$statistic
p4=d$signif
#################
e=mantel(otu_bray, vegdist(design[9], method=fangfa), permutations=999)#做999次重复用来计算P值
r5=e$statistic
p5=e$signif
###################
f=mantel(otu_bray, vegdist(design[10], method=fangfa), permutations=999)#做999次重复用来计算P值
r6=f$statistic
p6=f$signif
#############
g=mantel(otu_bray, vegdist(design[11], method=fangfa), permutations=999)#做999次重复用来计算P值
r7=g$statistic
p7=g$signif
#################
h=mantel(otu_bray, vegdist(design[12], method=fangfa), permutations=999)#做999次重复用来计算P值
r8=h$statistic
p8=h$signif
#################
i=mantel(otu_bray, vegdist(design[13], method=fangfa), permutations=999)#做999次重复用来计算P值
r9=i$statistic
p9=i$signif
#############
j=mantel(otu_bray, vegdist(design[14], method=fangfa), permutations=999)#做999次重复用来计算P值
r10=j$statistic
p10=j$signif
####################
k=mantel(otu_bray, vegdist(design[15], method=fangfa), permutations=999)#做999次重复用来计算P值
r11=k$statistic
p11=k$signif
################
l=mantel(otu_bray, vegdist(design[16], method=fangfa), permutations=999)#做999次重复用来计算P值
r12=l$statistic
p12=l$signif
##############
m=mantel(otu_bray, vegdist(design[17], method=fangfa), permutations=999)#做999次重复用来计算P值
r13=m$statistic
p13=m$signif
############
R=c(r,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13)
p=c(p,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13)
envyl

name=colnames(envyl)
name=c("quanbu",name)
web.data<-data.frame(name=name,R_otu=R,P_otu=p)
write.table(web.data,"Mantel testRNA_L2.txt",quote = FALSE,row.names = F,
            col.names = T,sep = "\t")
##################################
#############################
######################
###############
##########otu_table_L6.txt
setwd("E:/Shared_Folder/NL")
design = read.table("map_R.txt", header=T, row.names= 1, sep="\t") 
head(design)
otu_table = read.table("otu_table_L6.txt",row.names= 1,header=T, sep="\t",check.names=F)
#过滤数据并且排序
idx = rownames(design) %in% colnames(otu_table)
idx
design = design[idx,]
design
otu_table = otu_table[rownames(design)] 
head(otu_table)
otu_table = t(otu_table)
otu_bray<-vegdist(otu_table, method='bray')

######环境数据
fangfa=c("euclidean")

envyl = design[5:17]
envyl

env_<-vegdist(envyl, method=fangfa)#euclidean
env_
o=mantel(otu_bray, env_, permutations=999)#做999次重复用来计算P值
r=o$statistic
p=o$signif

###########
envyl_Moisture = design[5]
envyl_Moisture

env_Moisture<-vegdist(envyl_Moisture, method=fangfa)#euclidean
env_Moisture
a=mantel(otu_bray, env_Moisture, permutations=999)#做999次重复用来计算P值
r1=a$statistic
p1=a$signif
#################
envyl_NO3.N = design[6]
envyl_NO3.N

env_NO3.N<-vegdist(envyl_NO3.N, method=fangfa)#euclidean
env_NO3.N
b=mantel(otu_bray, env_NO3.N, permutations=999)#做999次重复用来计算P值
r2=b$statistic
p2=b$signif
##################
envyl_NH4.N = design[7]
envyl_NH4.N

env_NH4.N<-vegdist(envyl_NH4.N, method=fangfa)#euclidean
env_NH4.N
c=mantel(otu_bray, env_NH4.N, permutations=999)#做999次重复用来计算P值
r3=c$statistic
p3=c$signif
#############

d=mantel(otu_bray, vegdist(design[8], method=fangfa), permutations=999)#做999次重复用来计算P值
r4=d$statistic
p4=d$signif
#################
e=mantel(otu_bray, vegdist(design[9], method=fangfa), permutations=999)#做999次重复用来计算P值
r5=e$statistic
p5=e$signif
###################
f=mantel(otu_bray, vegdist(design[10], method=fangfa), permutations=999)#做999次重复用来计算P值
r6=f$statistic
p6=f$signif
#############
g=mantel(otu_bray, vegdist(design[11], method=fangfa), permutations=999)#做999次重复用来计算P值
r7=g$statistic
p7=g$signif
#################
h=mantel(otu_bray, vegdist(design[12], method=fangfa), permutations=999)#做999次重复用来计算P值
r8=h$statistic
p8=h$signif
#################
i=mantel(otu_bray, vegdist(design[13], method=fangfa), permutations=999)#做999次重复用来计算P值
r9=i$statistic
p9=i$signif
#############
j=mantel(otu_bray, vegdist(design[14], method=fangfa), permutations=999)#做999次重复用来计算P值
r10=j$statistic
p10=j$signif
####################
k=mantel(otu_bray, vegdist(design[15], method=fangfa), permutations=999)#做999次重复用来计算P值
r11=k$statistic
p11=k$signif
################
l=mantel(otu_bray, vegdist(design[16], method=fangfa), permutations=999)#做999次重复用来计算P值
r12=l$statistic
p12=l$signif
##############
m=mantel(otu_bray, vegdist(design[17], method=fangfa), permutations=999)#做999次重复用来计算P值c("quanbu",name)
web.data<-data.frame(name=name,R_otu=R,P_otu=p)
write.table(web.data,"Mantel testRNA_L6.txt",quote = FALSE,row.names = F,
            col.names = T,sep = "\t")
