
###在成图之前首先我们要做一个注释文件：#####

library("FishTacoPlot")

require(FishTacoPlot); require(ggplot2); require(scales); require(grid)


#"E:/Shared_Folder/HG_kangbing/nobac_noqianheti_chuli/fishtaco_input/第三次fishtacomulti_taxa"
#########d第四次
mi=c("#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33", "#7c260b",
     "#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A", "#66A61E", "#E6AB02",
     "#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072" ,"#80B1D3", "#FDB462", "#B3DE69",
     "#FCCDE5","#D9D9D9", "#BC80BD","grey60","#E7278A", "#64A61E","#6794a7", "#014d64",
     "#01a2d9" ,"#7ad2f6" ,"#00887d" ,"#76c0c1", "#7c260b", "#ee8f71", "#adadad",
     "#008FD5" ,"#FF2700", "#77AB43", "#c72e29", "#016392", "#be9c2e", "#098154","#fb832d","#000000")

show_col(mi)
p <- MultiFunctionTaxaContributionPlots(input_dir="../第四次结果B80富集通路/",
                                        input_prefix="fishtaco_out_no_inf",
                                        input_taxa_taxonomy="./tax_100.txt", sort_by="list", plot_type="bars",
                                        input_function_filter_list=c("ko00472",
                                                                     "ko00944",
                                                                     "ko03450",
                                                                     "ko04614",
                                                                     "ko02040",
                                                                     "ko00440","ko00591","ko02060", "ko00072","ko00120","ko02030","ko00531","ko00351","ko00643","ko00281"),
                                        add_predicted_da_markers=F, 
                                        add_original_da_markers=F,input_permutation="single_taxa",
                                        #min_contribution=0.05,min_cont_as_separate=0.001,
                                        scale_pos_by_original_pos_sum=T)+
  scale_fill_manual(values =mi)
p





p <- p + scale_x_continuous(breaks=c(1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15), labels=c("D-Arginine and D-ornithine metabolism",
                                                                                      "Flavone and flavonol biosynthesis",
                                                                                      "Non-homologous end-joining",
                                                                                      "Renin-angiotensin system",
                                                                                      "Flagellar assembly",
                                                                                      "Phosphonate and phosphinate metabolism",
                                                                                      "Linoleic acid metabolism",
                                                                                      "Phosphotransferase system (PTS)",
                                                                                      "Synthesis and degradation of ketone bodies",
                                                                                      "Primary bile acid biosynthesis",
                                                                                      "Bacterial chemotaxis",
                                                                                      "Glycosaminoglycan degradation",
                                                                                      "DDT degradation",
                                                                                      "Styrene degradation",
                                                                                      "Geraniol degradation")) +
  guides(fill=guide_legend(nrow=7)) + ylab("Wilcoxon test statistic (W)") +
  theme(plot.title=element_blank(), axis.title.x=element_text(size=12,colour="black",face="plain"),
        axis.text.x=element_text(size=10,colour="black",face="plain"), axis.title.y=element_blank(),
        axis.text.y=element_text(size=10,colour="black",face="plain"), axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major.x = element_line(colour="light gray"), panel.grid.major.y = element_line(colour="light gray"),
        panel.grid.minor.x = element_line(colour="light gray"), panel.grid.minor.y = element_line(colour="light gray"),
        panel.background = element_rect(fill="transparent",colour=NA), panel.border = element_rect(fill="transparent",colour="black"),
        legend.background=element_rect(colour="black"), legend.title=element_text(size=10), legend.text=element_text(size=8,face="plain"),
        legend.key.size=unit(0.8,"line"), legend.margin=unit(0.1,"line"), legend.position="bottom")

p
#####################分泌物中也会用差异的通路#########


p <- MultiFunctionTaxaContributionPlots(input_dir="../第四次结果B80富集通路",
                                        input_prefix="fishtaco_out_no_inf",
                                        input_taxa_taxonomy="./tax_100.txt", sort_by="list", plot_type="bars",
                                        input_function_filter_list=c("ko00020",
                                                                     "ko00250",
                                                                     "ko00290",
                                                                     "ko00620",
                                                                     "ko00720"),
                                        add_predicted_da_markers=F, 
                                        add_original_da_markers=F,input_permutation="single_taxa",
                                        #min_contribution=0.05,min_cont_as_separate=0.001,
                                        scale_pos_by_original_pos_sum=T)+
  scale_fill_manual(values =mi)
p




p <- p + scale_x_continuous(breaks=c(1, 2, 3,4,5), labels=c("Citrate cycle (TCA cycle)",
                                                            "Alanine, aspartate and glutamate metabolism",
                                                            "Valine, leucine and isoleucine biosynthesis",
                                                            "Pyruvate metabolism",
                                                            "Carbon fixation pathways in prokaryotes")) +
  guides(fill=guide_legend(nrow=7)) + ylab("Wilcoxon test statistic (W)") +
  theme(plot.title=element_blank(), axis.title.x=element_text(size=12,colour="black",face="plain"),
        axis.text.x=element_text(size=10,colour="black",face="plain"), axis.title.y=element_blank(),
        axis.text.y=element_text(size=10,colour="black",face="plain"), axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major.x = element_line(colour="light gray"), panel.grid.major.y = element_line(colour="light gray"),
        panel.grid.minor.x = element_line(colour="light gray"), panel.grid.minor.y = element_line(colour="light gray"),
        panel.background = element_rect(fill="transparent",colour=NA), panel.border = element_rect(fill="transparent",colour="black"),
        legend.background=element_rect(colour="black"), legend.title=element_text(size=10), legend.text=element_text(size=8,face="plain"),
        legend.key.size=unit(0.8,"line"), legend.margin=unit(0.1,"line"), legend.position="bottom")

p

#---附件#--这部分代码用于整理fishtaco原始文件-供参考
##
# ################################第四次跑fishtaco前处理######
# setwd("E:/Shared_Folder/HG_kangbing/nobac_noqianheti_chuli/fishtaco_input/第四次运行fishtaco")
# otu = read.table("gg135_otu_table.txt", header=T, row.names= 1, sep="\t") 
# head(otu)
# # 转换原始数据为百分比
# norm = t(t(otu)/colSums(otu,na=T)) * 100 # normalization to total 100
# 
# 
# write.table(norm,"gg13.5_tu_table_norm.txt",quote = FALSE,row.names = T,
#             col.names = T,sep = "\t")
# ##
# metagenome = read.table("metagenome_predictions.txt", header=T, row.names= 1, sep="\t") 
# head(metagenome)
# 
# # 转换原始数据为百分比
# norm = t(t(metagenome)/colSums(metagenome,na=T)) * 100 # normalization to total 100
# 
# 
# write.table(norm,"metagenome_predictions_norm.txt",quote = FALSE,row.names = T,
#             col.names = T,sep = "\t")
# ##
# ################################第四次跑fishtaco前处理######
# ###在成图之前首先我们要做一个注释文件：######
# otu = read.table("otu_table1285.txt", header=T, row.names= 1, sep="\t") 
# head(otu)
# str(otu)
# tax = read.table("97_otu_taxonomy.txt", header=F, row.names= 1, sep="\t") 
# head(tax)
# tax$tox=rownames(tax)
# str(tax)
# idx =  rownames(tax)%in%rownames(otu) 
# tax1 = tax[idx,]
# rownames(tax1)=tax1$tox
# tax1$tox<-NULL
# head(tax1)
# write.table(tax1,"tax_100.txt",quote = FALSE,row.names = T,
#             col.names = F,sep = "\t")
# ##注释功能
# ###########
# setwd("E:/Shared_Folder/HG_kangbing/nobac_noqianheti_chuli/fishtaco_input/chongxinyunxign")
# 
# zhushi = read.table("TAXFUNpassway差异T格GC5-GF5.txt", header=T, row.names= 1, sep="\t") 
# head(zhushi)#zhushitonglu.txt
# pathway= read.table("zhushitonglu.txt", header=F, row.names= 1, sep="\t") 
# head(pathway)
# dim(pathway)
# str(pathway)
# index = merge(pathway,zhushi, by="row.names",all.x =T)
# head(index)
# dim(index)
# wwtt=data.frame(pathway=index$V2,zhushi=index$tonglu)
# write.table(wwtt,"差异通路注释文件.txt",quote = FALSE,row.names = T,
#             col.names = T,sep = "\t")
