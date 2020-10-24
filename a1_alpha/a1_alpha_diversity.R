#清空内存
rm(list=ls()) 
library("ggplot2") 
library("ggpubr") 
# 读入实验设计和Alpha多样性值
design = read.table("map_HG_kangbing_R.txt", header=T, row.names= 1, sep="\t") 
head(design, n=12L)
#QIIME
alpha = read.table("./alpha.txt", header=T, row.names= 1, sep="\t")
head(alpha)
#合并alpha多样性数值和mapping文件
index = merge(alpha, design,by="row.names",all=F)
head(index)

#定义颜色
mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A")
# 绘图代码、预览、保存PDF
p<-ggplot(index, aes(x=SampleType, y=chao1))+
  geom_boxplot(alpha=1, outlier.size=1, size=0.5, width=0.5,notchwidth=0.5,aes(fill=SampleType)) +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="", y="chao1 alpha diversity")+
  theme_bw()+
  scale_fill_manual(values = mi)+
  #labs(title = "toamto hea and dis")+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))
p
p = p+stat_compare_means()+
  theme(axis.text = element_text(size = 20,face = "bold"),
        legend.text = element_text(size = 15,face = "bold")
  )
p
ggsave("./chao1_alpha.pdf", p, width = 10, height = 6)

#显著性检验
#compare_means(chao1~SampleType, data = index, method = "t.test")
alpha_wilcox.test = compare_means(chao1~SampleType, data = index, method = "wilcox.test")


write.table(alpha_wilcox.test,"chao1_alpha_wilcox.test.txt",row.names = T,
            col.names = T,sep = "\t")



p<-ggplot(index, aes(x=SampleType, y=shannon, color=SampleType))+
  geom_boxplot(alpha=1, outlier.size=2, size=2, width=0.5,notchwidth=1, fill="transparent") +  
  geom_jitter( position=position_jitter(0.17), size=2, alpha=0.7)+
  labs(x="", y="shannon alpha diversity")+
  theme_bw()+scale_colour_manual(values = mi)+
  #labs(title = "toamto hea and dis")+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))

p+stat_compare_means()+theme_bw()+
  theme(
    axis.text = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    legend.position="none",
    axis.title.y =element_text(size = 24,face = "bold")
    
  )
ggsave("./shannon_alpha.pdf", p, width = 10, height = 6)

#显著性检验
alpha_wilcox.test = compare_means(shannon~SampleType, data = index, method = "wilcox.test")

write.table(alpha_wilcox.test,"shannon_alpha_wilcox.test.txt",row.names = T,
            col.names = T,sep = "\t")
#################

p<-ggplot(index, aes(x=SampleType, y=observed_otus, color=SampleType))+
  geom_boxplot(alpha=1, outlier.size=0.5, size=0.8, width=0.5, fill="transparent") +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="Groups", y="observed_otus index")+
  theme_classic()+scale_colour_manual(values = mi)+
  #labs(title = "toamto hea and dis")+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(plot.title = element_text(hjust = 0.5))
p+stat_compare_means()+theme_bw()+
  theme(
    axis.text = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    legend.position="none",
    axis.title.y =element_text(size = 24,face = "bold")
    
  )

ggsave("./observed_otus_alpha.pdf", p, width = 10, height = 6)

#显著性检验
alpha_wilcox.test = compare_means(observed_otus~SampleType, data = index, method = "wilcox.test")

write.table(alpha_wilcox.test,"observed_otus_alpha_wilcox.test.txt",row.names = T,
            col.names = T,sep = "\t")



###################
p<-ggplot(index, aes(x=SampleType, y=PD_whole_tree, color=SampleType))+
  geom_boxplot(alpha=1, outlier.size=0.5, size=0.8, width=0.5, fill="transparent") +  
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  labs(x="Groups", y="PD_whole_tree index")+
  theme_classic()+scale_colour_manual(values = mi)+
  #labs(title = "toamto hea and dis")+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(plot.title = element_text(hjust = 0.5))
#显著性检验
p+stat_compare_means()+theme_bw()+
  theme(
    axis.text = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 15,face = "bold"),
    legend.position="none",
    axis.title.y =element_text(size = 24,face = "bold")
    
  )

ggsave("./PD_whole_tree_alpha.pdf", p, width = 10, height = 6)

#显著性检验
alpha_wilcox.test = compare_means(PD_whole_tree ~ SampleType, data = index, method = "wilcox.test")

write.table(alpha_wilcox.test,"PD_whole_tree_alpha_wilcox.test.txt",row.names = T,
            col.names = T,sep = "\t")

#
###当样品重复较少的时候，alpha多样性分析做做柱状图比较好看############
library("ggplot2") 
library("ggpubr") 
library(dplyr)
#首先按照分组求其均值和标准差
# 按照Species分组，计算Sepal.Length的平均值和标准差
iris_groups<- group_by(index, SampleType)
df_summarise<- summarise(iris_groups, mean(chao1), sd(chao1))
# dplyr的计算是以tbl格式返回的，我们可以转成dataframe, 当然这步没必要，但是为了方便显示：
df_summarise<- as.data.frame(df_summarise)
colnames(df_summarise)=c("fenzu","chao1","sd")
df_summarise$chao1=round(df_summarise$chao1,2)
head(df_summarise)

p=ggplot(df_summarise, aes(x = fenzu, y = chao1)) + 
  geom_bar(stat = "identity", width = 0.4,position = "dodge",colour="black",fill="#984EA3")+ 
  geom_errorbar(aes(ymin=chao1-sd,ymax=chao1+sd),colour="#FF7F00",width=0.2,size=2)+
  geom_text(aes(label = chao1, vjust = -4, hjust = 0.5),size=6) +   ## 显示柱条上的数字
  ylim(0, max(df_summarise$chao1)*1.1)+
  labs(x="Groups", y="chao1 index")+
  theme_classic()+
  labs(title = "toamto hea and dis")+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  theme(plot.title = element_text(hjust = 0.5),axis.title = element_text(size = 15),axis.text.x  = element_text(size = rel(2)),axis.ticks.y = element_blank(),axis.text.y = element_blank())
p

ggsave("./chao1_alpha_barplot.pdf", p, width = 10, height = 10)

#写一个循环同时将全部的alpha多样性指标出图
#但是aes_string只能解决ggplot的出图问题，但是在很多的函数使用字符串未加引号，需要转化
head(index)
i =3

for (i in 2:6) {
  name_i <- colnames(index[i])
  p<-ggplot(index, aes_string(x="SampleType", y=name_i))+
    geom_boxplot(alpha=1, outlier.size=1, size=0.5, width=0.5,notchwidth=0.5,aes(fill=SampleType)) +  
    geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
    labs(x="", y="chao1 alpha diversity")+
    theme_bw()+
    scale_fill_manual(values = mi)+
    #labs(title = "toamto hea and dis")+
    guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))
  
  p = p+stat_compare_means()+
    theme(axis.text = element_text(size = 20,face = "bold"),
          legend.text = element_text(size = 15,face = "bold")
    )

  FileName <- paste(name_i,"_boxport", ".pdf", sep = "")
  ggsave(FileName, p, width = 10, height = 10)
  
  
}

