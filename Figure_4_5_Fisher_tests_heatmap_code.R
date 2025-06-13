###Code for Fisher tests and Fisher test Figures 

setwd("C:/Users/jacqu/Documents/Berkeley/Project Info")

library(tidyverse)

#read in data
gro = read.csv("Tab11_Master_Bacterial_Groupings.csv", header = TRUE)
ko = read.table("kofam_annotations_significant.txt", sep = "\t")
scaff = read.table("Final_Genomes_scaff2bin.tsv", sep = "\t")

gro$V2 = gro$Genome
scaffgen = left_join(scaff, gro, by = "V2")

#remove genes with more than one KO classifications 
#order rows by V2 then V6 (p value)
ko_ordered = ko[order(ko$V2, ko$V6), ]
#remove rows with duplicate V2 values (e.g., genes with >1 ko hit)
ko_unique <- ko_ordered[!duplicated(ko_ordered$V2), ]

#create V1 column
ko_unique$V1 = ko_unique$V2
#remove everything after the last underscore to go from gene accession to contig name
ko_unique$V1 <- sub("_[^_]*$", "", ko_unique$V1)
#join data frames
ko_scaffgen = left_join(ko_unique, scaffgen, by = "V1")
#create KO column
ko_unique$KO = ko_unique$V3

#data frame of MAGs based on KOs
mag = ko_scaffgen %>% group_by(Genome, V3) %>% summarize(freq = n()) %>% ungroup()
mag2 = mag %>% pivot_wider(values_from = freq, names_from = V3)
mag2[is.na(mag2)] <- 0

#turn into presence absence matrix
mag3 <- mag2 %>% mutate_if(is.numeric, list(~ ifelse(. > 1, 1, .)))
mag_temp = as.matrix(mag3[,2:ncol(mag3)])

#only include KOs with greater than 10 hits 
mag_temp2 = mag_temp[,(colSums(mag_temp) >10)]
#for smaller categories, only include KOs with greater than 5+ hits - e.g. Specific Host and Source Enrichment
#mag_temp2 = mag_temp[,(colSums(mag_temp) > 4)]

mag_temp3 = data.frame("Genome" = mag3$Genome, mag_temp2)

#add categorical information back in 
mag4 = left_join(mag_temp3, gro, by = "Genome")



#############FISHER TESTS
#Table 1
######################################################
########################## Core vs Everybody else
#create example contingency table using "table" function
#test = mag4 %>% select(K00003, Core_MB)
#test2 = table(test)
#ft = fisher.test(test2)

#pvalue
#ft$p.value
#odds ratio
#ft$estimate


#create for loop that loops through all columns starting with "K"
k_cols <- grep("^K", names(mag4), value = TRUE)

df = data.frame("KO" = character(), "p" = numeric(), "estimate" = numeric())


for (col in k_cols) {
  # Do something with the column
  col2 = col
  test = mag4 %>% select(all_of(col), Core_MB)
  test2 = table(test)
  ft = fisher.test(test2)
  new_row = data.frame(col2, ft$p.value, ft$estimate)
  df = rbind(df, new_row)
}

row.names(df) = NULL
names(df) = c("KO", "p", "estimate")
df2 = left_join(df, ko_unique[c("KO", "V7")], by = "KO") %>% unique()
names(df2) = c("KO", "p", "estimate", "function")
df2$p_adjust = p.adjust(df2$p, method = "fdr")
df3 = df2[order(df2$p_adjust), ]

#write.csv(df3, "Fisher_test1_core_vs_all.csv", row.names = FALSE)


#Table 3
######################################################
########################## Core vs all, taxa restricted to Core_MB, e.g. p__Proteobacteria, p__Spirochaetota, p__Bacteroidota 

#add categorical information back in 
mag_temp = left_join(mag3, gro, by = "Genome")

##Choose only MAGs from CORE MB taxa##### 
mag_temp2 = mag_temp %>% filter(Phylum == "p__Proteobacteria"| Phylum == "p__Spirochaetota" | Phylum == "p__Bacteroidota")
mag_temp3 = as.matrix(mag_temp2[,2:ncol(mag3)])
mag_temp4 = mag_temp3[,(colSums(mag_temp3) >10)]
mag4 = data.frame("Genome" = mag_temp2$Genome, mag_temp4, "Core_MB" = mag_temp2$Core_MB)

#create for loop that loops through all columns starting with "K"
k_cols <- grep("^K", names(mag4), value = TRUE)

df = data.frame("KO" = character(), "p" = numeric(), "estimate" = numeric())


for (col in k_cols) {
  # Do something with the column
  col2 = col
  test = mag4 %>% select(all_of(col), Core_MB)
  test2 = table(test)
  ft = fisher.test(test2)
  new_row = data.frame(col2, ft$p.value, ft$estimate)
  df = rbind(df, new_row)
}

row.names(df) = NULL
names(df) = c("KO", "p", "estimate")
df2 = left_join(df, ko_unique[c("KO", "V7")], by = "KO") %>% unique()
names(df2) = c("KO", "p", "estimate", "function")
df2$p_adjust = p.adjust(df2$p, method = "fdr")
df3 = df2[order(df2$p_adjust), ]

#write.csv(df3, "Fisher_test3_core_vs_all_core_taxa.csv", row.names = FALSE)


####In the end, only Table 1 and 3 were used for analysis
#code for other Fisher test comparisons is at the bottom of this code 


############################Fisher test 17
##fisher test of fisher test results - Fisher test of Fisher Test1 functional categories
##changed K00643 5-aminolevulinate synthase [EC:2.3.1.37] from "Amino acid metabolism" to "Cofactor metabolism"
setwd("C:/Users/jacqu/Documents/Berkeley/Project Info/Fisher_tests/")
#read in Fisher test results 
fish = read.csv("Fisher_test1_core_vs_all_for_volcano.csv", header = TRUE)

#remove all "underrepresented" KOs
#add significant and not significant labels 
fish$sig = ifelse(fish$p_adjust < 0.05 & fish$estimate > 1, "Significant", "Not")
fish$Category = gsub("_", " ", fish$Category)
fish$colourif = ifelse(fish$Pvalue > 1.4 & fish$Expression > 1, fish$Category, "Other")

#keep both over and underrepresented KOs
#fish$sig = ifelse(fish$p_adjust < 0.05 , "Significant", "Not")

unicat = fish$Category %>% unique()

#remove other category
#unicat = unicat[-2]

#create empty data frame
df = data.frame("KO" = character(), "p" = numeric(), "estimate" = numeric())

#Fisher test uses alphabetical order for categories, need to use category e.g. "Znot" that will consistently appear at beginning - high estimate means #enriched
for (x in unicat) {
  col2 = x 
  fish$Category2 = ifelse(fish$Category == col2, fish$Category, "00Not")
  test = fish %>% select(Category2, sig)
  test2 = table(test)
  ft = fisher.test(test2)
  new_row = data.frame(col2, ft$p.value, ft$estimate)
  df = rbind(df, new_row)
}

row.names(df) = NULL
names(df) = c("Category", "p", "estimate")
df2 = df
df2$p_adjust = p.adjust(df2$p, method = "fdr")
df3 = df2[order(df2$p_adjust), ]


#write.csv(df3, "Fisher_test17_fisher_test_results1_categories.csv", row.names = FALSE)

##fisher test of fisher test results - Fisher test of Fisher Test1 functional categories
############################Fisher test 18
##changed K00643 5-aminolevulinate synthase [EC:2.3.1.37] from "Amino acid metabolism" to "Cofactor metabolism"
#to remove 0 and Inf estimates, took max and min values for estimate and used those to fill in table 
fish = read.csv("Fisher_test3_core_vs_all_core_taxa_for_volcano.csv", header = TRUE)

colours = c("#FFBC42", "#D81159", "#F0F757", "#218380", "#73D2DE", "#EDADC7", "#C5DECD", "#1C110A","grey95",  "#12ABF8", "#564787", "#F64740", "#C19875", "#FCB07E", "#98473E", "#E4D6A7", "#FFA400") 

#remove all "underrepresented" KOs
#fish = fish %>% filter(estimate > 1)
fish$sig = ifelse(fish$p_adjust < 0.05 & fish$estimate > 1, "Significant", "Not")
#keep both over and underrepresented KOs
#fish$sig = ifelse(fish$p_adjust < 0.05 , "Significant", "Not")
fish$Category = gsub("_", " ", fish$Category)
fish$colourif = ifelse(fish$Pvalue > 1.4 & fish$Expression > 1, fish$Category, "Other")

#volcano plot
xx = ggplot(fish, aes(x= Expression, y = Pvalue)) + geom_point(aes(colour = colourif, size = Pvalue), alpha = 0.8) + theme(panel.background =element_blank(), panel.border = element_rect(colour = "black", fill = NA), legend.key = element_blank(), panel.grid.major = element_line(colour = 'grey95')) + scale_colour_manual(values = colours) + scale_y_continuous(expand = c(0,0), limits = c(0,6)) + scale_x_continuous(limits = c(0,2), expand = c(0,0)) + labs(x = "log(Fold change)", y = "-log(adjusted P value)", colour = "Functional Category") + geom_vline(xintercept = 1, colour = "red", alpha = 0.5) + geom_hline(yintercept = 1.4, colour = "red", alpha = 0.5) + scale_size_continuous( range = c(1,3)) + guides(size = "none")
xx

unicat = fish$Category %>% unique()
df = data.frame("KO" = character(), "p" = numeric(), "estimate" = numeric())

#Fisher test uses alphabetical order for categories, need to use category e.g. "00not" that will consistently appear at beginning - high estimate means #enriched
for (x in unicat) {
  col2 = x 
  fish$Category2 = ifelse(fish$Category == col2, fish$Category, "00Not")
  test = fish %>% select(Category2, sig)
  test2 = table(test)
  ft = fisher.test(test2)
  new_row = data.frame(col2, ft$p.value, ft$estimate)
  df = rbind(df, new_row)
}

row.names(df) = NULL
names(df) = c("Category", "p", "estimate")
df2 = df
df2$p_adjust = p.adjust(df2$p, method = "fdr")
df3 = df2[order(df2$p_adjust), ]


#write.csv(df3, "Fisher_test18_fisher_test_results3_categories.csv", row.names = FALSE)


####################################
####Figure of Fisher category results 


#all significant - ended up using this one!!!!
en = read.csv("Fisher_test17_fisher_test_results1_categories_all_for_figure.csv", header = TRUE)

en2 = en[order(en$Enrichment), ]
en2$Category <- reorder(en2$Category, en2$Enrichment)

xx = ggplot(en2, aes(x = Enrichment, y = Category)) + geom_bar(stat = "identity", position = "dodge", aes(fill = Taxa, alpha = Significant, colour = Significant), width = 0.75, linewidth = 0.25) + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), panel.grid.major = element_line(colour = "grey98"), axis.title = element_text(colour = "grey25")) + geom_vline(xintercept = 0) + scale_alpha_manual(values = c(0.2,1)) + scale_colour_manual(values = c("grey90", "black")) + scale_fill_manual(values = c("#D81159", "#73D2DE")) + guides(alpha = "none", colour = "none") + labs(x = "Enrichment in Core Microbiome") + scale_x_continuous(limits = c(-0.75,0.75))
xx
#ggsave("Fisher_fisher_results_all_signficant_bar.png", height = 6, width = 5.5)


#enriched only 
en = read.csv("Fisher_test17_fisher_test_results_all_enriched_for_figure.csv", header = TRUE)

en2 = en[order(en$Enrichment), ]
en2$Category <- reorder(en2$Category, en2$Enrichment)

xx = ggplot(en2, aes(x = Enrichment, y = Category)) + geom_bar(stat = "identity", position = "dodge", aes(fill = Taxa, alpha = Significant, colour = Significant), width = 0.75, linewidth = 0.25) + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), panel.grid.major = element_line(colour = "grey98"), axis.title = element_text(colour = "grey25")) + geom_vline(xintercept = 0) + scale_alpha_manual(values = c(0.2,1)) + scale_colour_manual(values = c("grey90", "black")) + scale_fill_manual(values = c("#D81159", "#73D2DE")) + guides(alpha = "none", colour = "none") + labs(x = "Enrichment in Core Microbiome")
xx

####
######################################################################

#code to make heatmap figure 

setwd("C:/Users/jacqu/Documents/Berkeley/Project Info")

library(tidyverse)

#read in data
#MAG data
gro = read.csv("MAGS_plasmids.csv", header = TRUE)
#KO data
ko = read.table("kofam_annotations_significant.txt", sep = "\t")
#table linking contigs to MAGs
scaff = read.table("Final_Genomes_scaff2bin.tsv", sep = "\t")

#join tables
gro$V2 = gro$Genome
scaffgen = left_join(scaff, gro, by = "V2")

#remove genes with more than one KO classifications 
#order rows by V2 then V6 (p value)
ko_ordered = ko[order(ko$V2, ko$V6), ]
#remove rows with duplicate V2 values (e.g., genes with >1 ko hit)
ko_unique <- ko_ordered[!duplicated(ko_ordered$V2), ]

#create V1 column
ko_unique$V1 = ko_unique$V2
#remove everything after the last underscore to go from gene accession to contig name
ko_unique$V1 <- sub("_[^_]*$", "", ko_unique$V1)
#join tables by V1
ko_scaffgen = left_join(ko_unique, scaffgen, by = "V1")
#create KO column
ko_unique$KO = ko_unique$V3

#data frame of MAGs based on KOs
mag = ko_scaffgen %>% group_by(Genome, V3) %>% summarize(freq = n()) %>% ungroup()
mag2 = mag %>% pivot_wider(values_from = freq, names_from = V3)
mag2[is.na(mag2)] <- 0

#turn into presence absence matrix
mag3 <- mag2 %>% mutate_if(is.numeric, list(~ ifelse(. > 1, 1, .)))
mag_temp3 = data.frame("Genome" = mag3$Genome, mag3)

#add categorical information back in 
mag4 = left_join(mag_temp3, gro, by = "Genome")


######################################################################################
###individual heatmap code
setwd("C:/Users/jacqu/Documents/Berkeley/Project Info")

library(tidyverse)
library(vegan)

#read in heatmap data
hm = read.csv("heatmap_function_figure4_individual2.csv", header = TRUE)
#subset to presence/absence data
hm2 = data.frame("Genome"= hm$Genome, hm[,6:ncol(hm)])
#rename hm2
ordered = hm2
ordered_long = ordered %>% pivot_longer(-Genome, names_to = "Function", values_to = "Presence")

#add MAG info back in
ordered_long2 = left_join(ordered_long, hm[,1:5], by = "Genome")

#read in function info 
key = read.csv("Function_metabolic_key_plasmid.csv", header = TRUE)

ordered_long3 = left_join(ordered_long2, key, by = "Function")
ordered_long4 = ordered_long3 %>% filter(Keep == "TRUE")

#change order of categories
ordered_long4$Category <- factor(ordered_long4$Category, levels = c("Anoxygenic Photosynthesis", "Carbon Fixation", "C1 Metabolism", "Carbohydrate", "Respiration", "Porphyrin Metabolism", "Nitrogen", "Sulfur", "Phosphate" ))

#change order of categories to match "Order" column
new_levels <- unique(ordered_long4$Category2[order(ordered_long4$Order)])

# Convert x to a factor with the new order of levels
ordered_long4$Category2 <- factor(ordered_long4$Category2, levels = new_levels)


#no order just MAGs by taxa 
xx = ggplot(ordered_long4, aes(x = Category2, y = Genome)) + geom_tile(aes(fill = Presence),colour = "white") + facet_grid(taxa~Category, space = "free", scales = "free" ) + theme(legend.position = "none", panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),panel.spacing = unit(0.1, "lines"), strip.text.x = element_text(angle = 90, size = 8), strip.background = element_rect(fill = "grey90", colour = "black"), strip.text.y = element_text(angle = 0, size = 8), axis.text.x = element_text(angle = 90, size =6, vjust = 0.4, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_fill_gradient(high = "grey20", low = "grey97") + labs(x = "", y = "")
xx

#ggsave("Heatmap_Figure4_individual_plasmid.png", height = 6.5, width = 7.5)

#changed heatmap colour, and order of taxa in inkscape, although this could be easily coded here

