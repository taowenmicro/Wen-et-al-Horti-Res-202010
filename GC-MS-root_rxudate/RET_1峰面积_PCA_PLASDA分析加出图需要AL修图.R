library(mixOmics)
library(ggplot2)
library(reshape2)
library(ggdendro)
library(grid)
library(ggbiplot)
library (reshape2)
library(dplyr)
# install.packages("BiocManager")
# #安装的软件包可以更新到当前版本
# BiocManager::install("devtools")
# devtools::install_github("vqv/ggbiplot")
# BiocManager::install("mixOmics")

#读取数据力丰-B80分泌物原始数据.txt
design = read.table("RET_map_R.txt", header=T, row.names= 1, sep="\t") 
head(design)
# 读取OTU表,全部otu表没有抽平
Root_exudates = read.csv("RET_a2缺失值处理后分泌物原始数据.csv",row.names = 1)

##有一个出峰值非常离谱，这里去掉Analyte 12，手动去除
#design3<-Root_exudates[grep("Analyte 12",Root_exudates$id),]
#idx=Root_exudates$id%in%design3$id
#Root_exudates=Root_exudates[!idx,]
dim(Root_exudates)
##变成了707个出峰值
# 过滤数据并排序
idx = rownames(design) %in% colnames(Root_exudates) 
idx
sub_design = design[idx,]
count = Root_exudates[, rownames(sub_design)]
head(count)
dim(count)
count[is.na(count)] <- 0.0
##########去掉全部为0的代谢物
n=ncol(count)
#增加一行，为整列的均值，计算每一列的均值，2就是表示列
count[n+1]=apply(count[c(1:nrow(count)),],1,sum)
#选择sumsqs大于5的otu
count=count[count[n+1] > 0,1:n] 
head(count)
dim(count)

###################分泌物每个样品总峰面积图#############
########计算每一列的和作为总的分泌物#
count[n+1,]=apply(count[c(1:nrow(count)),],2,sum)
count1=count[n+1,]  
count1

count1$id=row.names(count1)
fengdu <- melt (count1, id="id")
head(fengdu)
colnames(fengdu)=c("name","break1","mean")

fengdu$breaks = factor(rep(c("B801", "LF1"), c(3,3)))
fengdu$break1=NULL
fengdu$name=NULL
##########下面做差异分析
model<-aov(mean~breaks, data=fengdu)
summary(model)


########计算每一列的和作为总的分泌物

grou11<- group_by(fengdu,breaks)
summarise11<- dplyr::summarise(grou11,  mean(mean), sd(mean))
ummarise11<- as.data.frame(summarise11)
head(ummarise11)
colnames(ummarise11)=c("break1","mean","sd")
ummarise11$mean=ummarise11$mean/100000
ummarise11$sd=ummarise11$sd/100000
ummarise11$break1 = factor(rep(c("CSF1", "CSF2"), c(1,1)))
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
library("ggplot2") 

###################柱状图分泌物每个样品总峰面积图#############
p=ggplot(ummarise11, aes(x = break1, y = mean,fill=break1)) + 
  geom_bar(stat = "identity", width = 0.4,position = "dodge",colour="black")+ 
  geom_errorbar(aes(ymin=ummarise11$mean-ummarise11$sd,
                    ymax=ummarise11$mean+ummarise11$sd),
                colour="black",width=0.1,size=1)+
  scale_x_discrete(limits=c("CSF1","CRF1", "CSF2", "CRF2"))
p
p=p+theme_light()+theme(axis.text.x = element_text(colour = "black",size = 20,face = "bold"),
                        axis.text.y = element_text(colour = "black",size = 20,face = "bold")
                        
)+scale_colour_manual(values = mi)

ggsave("a3_分泌物全部峰总面积四组出图.pdf", p, width = 10, height = 6)
####
###################分泌物每个样品总峰面积图#############



#########做一个PCA分析###############
# 进行PCA分析
otu.pca <- prcomp(t(count), scale. = TRUE)
p=ggbiplot(otu.pca, obs.scale = 1, var.scale = 1,
           groups = sub_design$SampleType, ellipse = TRUE,var.axes = F,circle = TRUE)+
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
#
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
p+theme_bw()+scale_colour_manual(values = mi)+labs(x=paste("PC 1 (", 61.7, "%)", sep=""),
                                                   y=paste("PC 2 (", 14.5, "%)", sep=""))+
  labs(title = "PCA")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
##########
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
set.seed(10)
p=ggplot(data =yangpin,aes(x=yangpin$PC1,y=yangpin$PC2,group=SampleType,color=SampleType))+geom_point(size=5)+
  stat_ellipse(type = "t", linetype = 2)+
  labs(x=paste("PC 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PC 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  labs(title = "PCA") 

mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
p=p+theme_bw()+scale_colour_manual(values = mi)+labs(x=paste("PC 1 (", 61.7, "%)", sep=""),
                                                     y=paste("PC 2 (", 14.5, "%)", sep=""))+
  labs(title = "PCA")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
p
ggsave("a1_PCA 分泌物四组出图.pdf", p, width = 10, height = 6)

#########################
#转换成矩阵
datatm <-as.matrix(count)
XXt<-t(datatm)
#PLS-DA分析，这里也是选取2个主成分

plsda.datatm <-plsda(XXt, as.factor(sub_design$SampleType), ncomp = 2)
unclass(plsda.datatm)
names(plsda.datatm)

plsda.datatm$mat.c
#PLS-DA without centroids
mi=c("#1B9E77" ,"#D95F02")

plotIndiv(plsda.datatm , comp = c(1,2),
          group = sub_design$SampleType, style = 'ggplot2' )
pdf("a2_PLSDA分泌物四组出图.pdf", width = 10, height = 6)
plotIndiv(plsda.datatm , comp = c(1,2),
          group = sub_design$SampleType, style = 'ggplot2',ellipse = TRUE, 
          size.xlabel = 20, size.ylabel = 20, size.axis = 25, pch = 15, cex = 5)

dev.off()

#----提取数据作图

a = unclass(plsda.datatm)
#--提取坐标值
plotdata = as.data.frame(a$variates$X)
plotdata$SampleType = sub_design$SampleType
#-提取解释度
eig = a$explained_variance$X
eig[1]

# library(ggalt)

p = ggplot(data = plotdata,aes(x=comp1,y=comp2,group=SampleType,color=SampleType))+geom_point(size=5)+
  stat_ellipse(type = "t", linetype = 2)+
  # geom_encircle(s_shape=1, expand=0) +
  labs(x=paste("X-variate 1 (", format(100 * eig[1]), "%)", sep=""),
       y=paste("X-variate 2 (", format(100 * eig[2] ), "%)", sep=""))+
  labs(title = "PLS-DA") 
p
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
p=p+theme_bw()+scale_colour_manual(values = mi)+
  labs(title = "PCA")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
p


