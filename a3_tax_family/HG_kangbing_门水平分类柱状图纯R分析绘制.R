 
library("DESeq2")
library("phyloseq")
library("ggplot2")
library("dplyr")
China_phyloseq <- import_biom(BIOMfilename = "gg13.5_otu_table_tax.biom")

colnames(tax_table(China_phyloseq))<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

metadata <- import_qiime_sample_data("map_HG_kangbing.txt")
China_phyloseq_complete <- merge_phyloseq(China_phyloseq, metadata)
China_phyloseq_complete 

#change col names in tax_table
#head(tax_table(China_phyloseq_complete))
colnames(tax_table(China_phyloseq_complete)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
head(tax_table(China_phyloseq_complete))



###Taxonomy cumulative bar plot, fig 2A
#Taxonomy barblots all replicates shown, got this code somewhere online..
Taxonomies <- China_phyloseq_complete %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)
head(Taxonomies)

#Set colors for classes
Phylum_colors <- c(
  "#89C5DA", "#74D944", "#DA5724", "#CE50CA", "#D3D93E", "#C0717C", "#CBD588", "#5F7FC7", 
  "#673770", "#3F4921", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
  "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
  "#8A7C64", "#00FF66", "#339933", "#993366"
)
#get rid of rank numbers and other ugly things in taxonomic ranks
#Final_table_replicates_class$Kingdom <- gsub("D_0__","", Final_table_replicates_class$Kingdom)
Taxonomies$Family <- gsub("f__","", Taxonomies$Family)
#Taxonomies$Phylum1 <- gsub("_"," ", Taxonomies$Phylum1)
#make abundances per group sum to 1 then make percentages
Taxonomies$Abundance = Taxonomies$Abundance / 6
Taxonomies$Abundance = Taxonomies$Abundance * 100
head(Taxonomies)
dim(Taxonomies)
#########Taxonomy barplots summarized per Group
p = china_barplots <- ggplot(Taxonomies, aes(x = SampleType, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Phylum_colors) +
  theme(axis.title.x = element_blank()) +
  theme(legend.text=element_text(size=6)) +
  scale_y_continuous(name = "Abundance (%)")
print(china_barplots)

ggsave("科水平水平柱状图.pdf", p, width = 6, height =8 )

