

library("FishTacoPlot")

require(FishTacoPlot); require(ggplot2); require(scales); require(grid)

mi=c("#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33", "#7c260b",
     "#1B9E77", "#D95F02" ,"#7570B3" ,"#E7298A", "#66A61E", "#E6AB02",
     "#8DD3C7", "#FFFFB3", "#BEBADA","#FB8072" ,"#80B1D3", "#FDB462", "#B3DE69",
     "#FCCDE5","#D9D9D9", "#BC80BD","grey60","#E7278A", "#64A61E","#6794a7", "#014d64",
     "#01a2d9" ,"#7ad2f6" ,"#00887d" ,"#76c0c1", "#7c260b", "#ee8f71", "#adadad",
     "#008FD5" ,"#FF2700", "#77AB43", "#c72e29", "#016392", "#be9c2e", "#098154","#fb832d","#000000")
show_col(mi)




p <- MultiFunctionTaxaContributionPlots(input_dir="./",
                                        input_prefix="fishtaco_out_no_inf",
                                        input_taxa_taxonomy="./tax_100.txt", sort_by="list", plot_type="bars",
                                        input_function_filter_list=c("ko00020",
                                                                     "ko00053",
                                                                     "ko00410",
                                                                     "ko00330",
                                                                     "ko00630",
                                                                     "ko00260",
                                                                     "ko00770",
                                                                     "ko00620",
                                                                     "ko00520",
                                                                     "ko00561",
                                                                     "ko00052",
                                                                     "ko00400",
                                                                     "ko00240",
                                                                     "ko00290",
                                                                     "ko00270",
                                                                     "ko00250",
                                                                     "ko00480"
                                        ),
                                        add_predicted_da_markers=F, 
                                        add_original_da_markers=F,input_permutation="single_taxa",
                                        #min_contribution=0.05,min_cont_as_separate=0.001,
                                        scale_pos_by_original_pos_sum=T)+
  scale_fill_manual(values =mi)
p

p <- p + scale_x_continuous(breaks=c(1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), labels=c("Citrate cycle (TCA cycle)",
                                                                                      "Ascorbate and aldarate metabolism",
                                                                                      "beta-Alanine metabolism",
                                                                                      "Arginine and proline metabolism",
                                                                                      "Glyoxylate and dicarboxylate metabolism",
                                                                                      "Glycine, serine and threonine metabolism",
                                                                                      "Pantothenate and CoA biosynthesis",
                                                                                      "Pyruvate metabolism",
                                                                                      "Amino sugar and nucleotide sugar metabolism",
                                                                                      "Glycerolipid metabolism",
                                                                                      "Galactose metabolism",
                                                                                      "Phenylalanine, tyrosine and tryptophan biosynthesis",
                                                                                      "Pyrimidine metabolism",
                                                                                      "Valine, leucine and isoleucine biosynthesis",
                                                                                      "Cysteine and methionine metabolism",
                                                                                      "Alanine, aspartate and glutamate metabolism",
                                                                                      "Glutathione metabolism"
)) +
  guides(fill=guide_legend(nrow=7)) + ylab("Wilcoxon test statistic (W)") +
  theme(plot.title=element_blank(), axis.title.x=element_text(size=12,colour="black",face="plain"),
        axis.text.x=element_text(size=10,colour="black",face="plain"), axis.title.y=element_blank(),
        axis.text.y=element_text(size=10,colour="black",face="plain"), axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major.x = element_line(colour="light gray"), panel.grid.major.y = element_line(colour="light gray"),
        panel.grid.minor.x = element_line(colour="light gray"), panel.grid.minor.y = element_line(colour="light gray"),
        panel.background = element_rect(fill="transparent",colour=NA), panel.border = element_rect(fill="transparent",colour="black"),
        legend.background=element_rect(colour="black"), legend.title=element_text(size=10), legend.text=element_text(size=8,face="plain"),
        legend.key.size=unit(0.8,"line"),legend.margin=unit(0.1,"line"),  legend.position="bottom")

p
