##################################
####################将得到的代谢组数据表格原始表格导入
########待解决问题：本分泌物使用名字，但是公司给我表格有许多名字是一样的，所以下次解决
#清空内存
rm(list=ls()) 
#########下面来分类统计分泌物的信息
# 读取OTU表，这里我选择的是整个otu表格，但是一般没有必要全部做差异的啊，相对丰度高的做做就可以了
otu_table = read.delim("RET_a1原始数据.txt",sep="\t",header=T,check.names=F)
otu_table1=otu_table
otu_table$id=NULL
head(otu_table)
a=otu_table
a=as.data.frame(a)
a[is.na(a)] <- 0.0
str(a)
head(a)
c=data.frame(a=c(rep(0,nrow(a))),b=c(rep(0,nrow(a))),c=c(rep(0,nrow(a))),d=c(rep(0,nrow(a))),e=c(rep(0,nrow(a))),f=c(rep(0,nrow(a))),
             g=c(rep(0,nrow(a))),h=c(rep(0,nrow(a))),i=c(rep(0,nrow(a))),j=c(rep(0,nrow(a))),l=c(rep(0,nrow(a))),m=c(rep(0,nrow(a))))
colnames(c)=colnames(a)
row.names(c)=row.names(a)
c
str(c)
head(c)
for(i in 1:nrow(a)){
  if(sum(a[i,1:3][1,]==0)==1){
    c[i,1:3]= a[i,1:3]
    c[i,1:3][,c[i,1:3][1,]==0]=sum(a[i,1:3]/2)
  }else if(sum(b=a[i,1:3][1,]==0)==2){
    c[i,1:3]= a[i,1:3]
    c[i,1:3][,!c[i,1:3][1,]==0]=0
  }else{
    c[i,1:3]=a[i,1:3]
  }
}

for(i in 1:nrow(a)){
  if(sum(a[i,4:6][1,]==0)==1){
    c[i,4:6]=a[i,4:6]
    c[i,4:6][,c[i,4:6][1,]==0]=sum(a[i,4:6]/2)
  }else if(sum(b=a[i,4:6][1,]==0)==2){
    c[i,4:6]=a[i,4:6]
    c[i,4:6][,!c[i,4:6][1,]==0]=0
  }else{
    c[i,4:6]=a[i,4:6]
  }
}

head(c)
for(i in 1:nrow(a)){
  if(sum(a[i,7:9][1,]==0)==1){
    c[i,7:9]=a[i,7:9]
    c[i,7:9][,c[i,7:9][1,]==0]=sum(a[i,7:9]/2)
  }else if(sum(b=a[i,7:9][1,]==0)==2){
    c[i,7:9]=a[i,7:9]
    c[i,7:9][,!c[i,7:9][1,]==0]=0
  }else{
    c[i,7:9]=a[i,7:9]
  }
}



for(i in 1:nrow(a)){
  if(sum(a[i,10:12][1,]==0)==1){
    c[i,10:12]=a[i,10:12]
    c[i,10:12][,c[i,10:12][1,]==0]=sum(a[i,10:12]/2)
  }else if(sum(b=a[i,10:12][1,]==0)==2){
    c[i,10:12]=a[i,10:12]
    c[i,10:12][,!c[i,10:12][1,]==0]=0
  }else{
    c[i,10:12]=a[i,10:12]
  }
}
head(c,n = 50L)
c$id=otu_table1$id
dim(c)
##最后又将0替换为NA，方便excel使用
#c[c==0] = NA
write.table(c,"RET_a2缺失值处理后分泌物原始数据.txt",row.names = T,
            col.names = T,sep = "/t")



