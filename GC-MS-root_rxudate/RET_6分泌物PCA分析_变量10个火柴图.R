
# 安装包ggbiplot包
#install.packages("devtools", repo="http://cran.us.r-project.org")
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)

fenmiwu = read.delim("RET_a3筛选有化学式的分泌物并做一个分析.txt", header=T, row.names= 1,check.names=F) 
head(fenmiwu)
dim(fenmiwu)
fenmiwu1=fenmiwu[1:12]
head(fenmiwu1)

#导入需要的作图文件
design = read.table("RET_map_R1.txt", header=T, row.names= 1, sep="\t") 

# 过滤数据并排序
idx = rownames(design) %in% colnames(fenmiwu) 
sub_design = design[idx,]
count = fenmiwu[, rownames(sub_design)]
##########去掉全部为0的代谢物
n=ncol(count)
#增加一行，为整列的均值，计算每一列的均值，2就是表示列
count[n+1]=apply(count[c(1:nrow(count)),],1,sum)
#选择sumsqs大于5的otu
count=count[count[n+1] > 0,1:n] 
head(count)
dim(count)
head(count)
##
count1=log10(count)
head(count1)
count[count==0] <- 0.001
norm = t(t(count)/colSums(count,na=T))# * 100 # normalization to total 100
head(norm)
# 进行PCA分析
otu.pca <- prcomp(t(norm), scale. = TRUE)


#现在我们来提取作图坐标
#首先是样品排序轴坐标
predict(otu.pca) 
#也可以使用下面这条命令提取
yangpin<-otu.pca$x
yangpin=as.data.frame(yangpin)
yangpin$SampleType=sub_design$SampleType
#提取荷载坐标
bianliang<-otu.pca$rotation
bianliang=as.data.frame(bianliang)
head(bianliang)

dim(norm)
index = merge(norm ,bianliang, by="row.names",all=F)
head(index)
row.names(index)=index$Row.names
index$Row.names=NULL
dim(index)
wwe=fenmiwu[13]
head(wwe)
index1 = merge(index ,wwe, by="row.names",all=F)
dim(index1)
head(index1)
write.table(index1,"PCA载荷矩阵_norm变换.txt",quote = FALSE,row.names = T,
            col.names = T,sep = "\t")
##手动选择10个最终要的变量 PCA载荷矩阵挑选37个成分提取差异.txt

top10 = read.delim("PCA载荷矩阵挑选前10个成分提取差异.txt", header=T, row.names= 1,check.names=F) 
head(top10)
top10$ID = row.names(top10)
top10$PCone = top10$PC1^2
########这里做变量权重图的表格######

#######开始出图，做火柴图########
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
library("ggplot2") 
p=ggplot(top10, aes(x = PCone, y = reorder(ID,PCone)))  +
  geom_segment(aes(yend=ID),xend=0,size=3,colour = "#1B9E77" )+
  geom_point(size=4,pch=20, colour = "#1B9E77")
p 

p=  p+theme_bw()+theme(axis.text.x = element_text(colour = "black",size = 20,face = "bold"),
                       axis.text.y = element_text(colour = "black",size = 10,face = "bold"))
p
ggsave("a7_PCA10个变量重要性火柴棒图.pdf", p, width = 6, height = 4)



#提取特征根,这里提供的并不是特征值而是标准差，需要求其平方才是特征值
eig=otu.pca$sdev
eig=eig*eig
#############RET_a4化合物t检验结果CSF1和CRF1.xls
chayi = read.delim("RET_a4化合物t检验结果CSF1和CRF1.txt", header=T, row.names= 1,check.names=F) 
head(chayi)
cha=chayi[7:9]
head(cha)


tiaoxuan = read.delim("PCA载荷矩阵挑选37个成分提取差异.txt", header=T, row.names= 1,check.names=F) 
head(tiaoxuan)

index2 = merge(tiaoxuan,cha, by="row.names",all=F)
dim(index2)
head(index2)

write.table(index2,"PCA载荷矩阵挑选37个成分提取差异结果.txt",quote = FALSE,row.names = T,
            col.names = T,sep = "\t")

####
############


#现在我们来提取作图坐标
#首先是样品排序轴坐标
predict(otu.pca) 
#也可以使用下面这条命令提取
yangpin<-otu.pca$x
yangpin=as.data.frame(yangpin)
yangpin$SampleType=sub_design$SampleType
#提取荷载坐标
bianliang<-otu.pca$rotation
bianliang=as.data.frame(bianliang)
#提取特征根,这里提供的并不是特征值而是标准差，需要求其平方才是特征值
eig=otu.pca$sdev
eig=eig*eig
#在这里我设定了随机种子，方便两种形式图形比较
p=ggplot(data =yangpin,aes(x=yangpin$PC1,y=yangpin$PC2,group=SampleType,color=SampleType))+
  geom_point(colour="grey50", size =2)+
  geom_point(alpha=.7, size=5) +
  stat_ellipse(type = "t", linetype = 2)+
  labs(x=paste("PC 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PC 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  labs(title = "PCA") 

p=p+theme_bw()+
  scale_colour_manual(values = mi)+
  theme(
    axis.text = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    #legend.position="none",
    axis.title.y =element_text(size = 24,face = "bold"),
    axis.title.x =element_text(size = 24,face = "bold")
    
  )
ggsave(paste("PCA_分泌物.pdf", sep=""), p, width = 6, height = 4)

#########
# 绘制PCA图，并按组添加椭圆,这个就是个普通pcr图.scale = 1, var.scale = 1,
           groups = sub_design$SampleType, ellipse = TRUE,var.axes = F)+
  geom_point(colour="grey50", size =2)+
  geom_point(alpha=.7, size=5)


p+theme_bw()+
  scale_colour_manual(values = mi)+
  theme(
    axis.text = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    #legend.position="none",
    axis.title.y =element_text(size = 24,face = "bold"),
    axis.title.x =element_text(size = 24,face = "bold")
    
  )




