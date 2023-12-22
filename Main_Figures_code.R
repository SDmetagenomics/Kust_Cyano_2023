
### Load Libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zCompositions)
library(lme4)
library(lmerTest)
library(MuMIn)
library(emmeans)
library(multcomp)




####
### Functions for Transforms and Analysis 
####


### Function 1: notin
'%notin%' <- Negate('%in%')

### Function 2: geometric mean
gm.mean <- function(x){
  exp(mean(log(x)))
}

### Function 3: centered log-ratio transform for matrix input
clr.tfm <- function(mat, features = c("rows", "columns")){
  
  if(features == "rows"){
    geo_means <- apply(mat, 2, function(x) gm.mean(x))
    count_ratios <- sweep(mat, 2, geo_means, '/')
    clr_counts <- log(count_ratios)
    
    return(clr_counts)
  }
  
  if(features == "columns"){
    mat <- t(mat)
    geo_means <- apply(mat, 2, function(x) gm.mean(x))
    count_ratios <- sweep(mat, 2, geo_means, '/')
    clr_counts <- log(count_ratios)
    clr_counts <- t(clr_counts)
    
    return(clr_counts)
  }
  
}

### Function 4: additive log-ratio transform for matrix input
alr.tfm <- function(mat, features = c("rows", "columns"), stds = NULL){
  
  if(features == "rows"){
    geo_means <- apply(mat[stds,], 2, function(x) gm.mean(x))
    count_ratios <- sweep(mat, 2, geo_means, '/')
    alr_counts <- log(count_ratios)
    
    return(alr_counts)
  }
  
  if(features == "columns"){
    mat <- t(mat)
    geo_means <- apply(mat[stds,], 2, function(x) gm.mean(x))
    count_ratios <- sweep(mat, 2, geo_means, '/')
    alr_counts <- log(count_ratios)
    alr_counts <- t(alr_counts)
    
    return(alr_counts)
  }
  
}

### Function 5: Summarize Samples - Expects a matrix object where features are ASVs; also can accept a vector of ASVs that are standards 
sample.summary <- function(mat, features = c("rows", "columns"), stds = NULL){
  
  # transpose if features in columns
  if(features == "columns"){
    mat <- t(mat)
  }
  
  # calculate summary statistics 
  tmp_summary <- data.frame(Sample = colnames(mat),
                            Cts = apply(mat, 2, sum),
                            Mean_Cts = apply(mat, 2, mean),
                            SD_Cts = apply(mat, 2, sd),
                            Has_Cts = apply(mat, 2, function(x) sum(x > 0)),
                            Min_Cts = apply(mat, 2, min),
                            Zero_Cts = apply(mat, 2, function(x) sum(x == 0)),
                            Max_Cts = apply(mat, 2, max),
                            Max_Frc = (apply(mat, 2, max) / apply(mat, 2, sum)) * 100,
                            Max_ASV = apply(mat, 2, function(x) names(which(x == max(x)))[1]))
  
  # if std asvs provided calculate their total and fractional counts
  if(is.null(stds) == FALSE){
    
    # identify row numbers that are std ASVs  
    tmp_std_rows <- which(rownames(mat) %in% stds)
    
    # sum std ASV counts
    tmp_stds_ct <- data.frame(Sample = colnames(mat),
                              STD_Cts = apply(mat, 2, function(x) sum(x[tmp_std_rows])))
    
    # merge into tmp_summary
    tmp_summary <- merge(tmp_summary, tmp_stds_ct, by = "Sample")
    
    # add fraction of std asv counts per sample
    tmp_summary$STD_Frc <- (tmp_summary$STD_Cts / tmp_summary$Cts) * 100
    
    # add fraction of rest of asvs - (max + std)
    tmp_summary$Rest_Frc <- ((tmp_summary$Cts - (tmp_summary$STD_Cts + tmp_summary$Max_Cts)) / tmp_summary$Cts) * 100
  }
  
  # return data 
  return(tmp_summary)
  
}  





###
#################################################### Import  Data 1. ####################################################
####

####16S rRNA data

#####Data used for parts of Figure 1

####
### Load All The Data
####

### ASV Table
asvtab <- fread("~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/16S_Analysis_Runs/23_08_08_All_Samples/zotutab.txt")

### Metadata
metadata <- fread("~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/16S_Analysis_Runs/23_08_08_All_Samples/metadata_sd.txt")
good_samples <- fread("~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/16S_Analysis_Runs/23_08_08_All_Samples/good_metagenome_metadata.txt")

## check all good samples are in metadata
test_metadata <- subset(metadata, Tube_No %in% good_samples$tube_no)
rm(test_metadata)

### Taxonomy Data
taxtab <- fread("~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/16S_Analysis_Runs/23_08_08_All_Samples/taxonomy.tsv")

### Predicted Standards
pred_stand <- fread("~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/16S_Analysis_Runs/23_08_08_All_Samples/predicted_standards.txt")
pred_stand <- subset(pred_stand, PRED_STD == T)


### Plot Output Location
plot_out <- "~/Berkeley_Postdoc/Unicom_analyses/Figures/figures_all/"



####
### Generate Full Rank Split Taxonomy
####


### Rename columns 
colnames(taxtab) <- c("ASV_ID", "Full_Tax", "Confidence")


### Reformat Taxonomy Into Ranks by Column

## creating table of taxonomy and setting any that are unclassified as "NA"
asv_split_ranks <- data.frame(ASV_ID = NA, Rank1 = NA, Phylum = NA, Class = NA, Order = NA, Family = NA, Genus = NA, Species = NA)

## run for loop to split taxonomy into columns 
for (i in 1:nrow(taxtab)) {
  asv_tmp <- taxtab$ASV_ID[i]
  names(asv_tmp) <- "ASV_ID"
  tax_tmp <- phyloseq::parse_taxonomy_qiime(taxtab$Full_Tax[i])
  tmp_dat <- c(asv_tmp, tax_tmp)
  
  asv_split_ranks <- dplyr::bind_rows(asv_split_ranks, tmp_dat)
}

## rename columns and merge to tax table
asv_split_ranks <- asv_split_ranks[-1,]
colnames(asv_split_ranks)[2] <- "Domain"
asv_split_ranks[is.na(asv_split_ranks)] <- "Unknown"
taxtab <- merge(taxtab, asv_split_ranks, by = "ASV_ID")
rm(asv_split_ranks)
rm(tmp_dat)


####
### Filter Data + Create ASV Count Matrix + Generate Data Summary
####


### Specify Samples We Want By Filtering Metadata - Only taking samples at bi-weekly intervals
metadata_filt <- subset(metadata, Total_Days %in% c(0,14,28,42,56,70,84,91,105,119))


### Filter ASV table to contain only samples in metadata subset

## Identify columns(samples) numbers to filter 
tmp_select <- which(colnames(asvtab) %in% c("ASV_ID", metadata_filt$Sample))

## filter samples from asvtab
asvtab_filt <- asvtab[,..tmp_select]
rm(tmp_select)


### Convert ASV Table into ASV Count Matrix

## Save row names
tmp_names <- asvtab_filt$ASV_ID

## Make into matrix
asvmat_filt <- as.matrix(asvtab_filt[,-1])
rownames(asvmat_filt) <- tmp_names  


### Create Sample Summaries 

## Create summaries
asvmat_filt_summary <- sample.summary(asvmat_filt, features = "rows", stds = pred_stand$ID)

## Merge with metadata
asvmat_filt_summary <- merge(asvmat_filt_summary, metadata_filt, by = "Sample")


### Filter Metadata Table (Again) Based off Sample Summaries

## Remove low count samples 
metadata_filt <- subset(metadata_filt, Sample != "3055ERFCBW1B")

## Remove samples where max ASV is not main Cyano

# identify samples
tmp_select <- subset(asvmat_filt_summary, Max_ASV %in% c("ASV1", "ASV2", "ASV3", "ASV4"))$Sample

# filter metadata 
metadata_filt <- subset(metadata_filt, Sample %in% tmp_select)

## Remove samples where max ASV does not match expected Cyano

# identify samples
tmp_select <- c("L0902ERSEBW1CPF",
                "L0902DBMBW1BPF",
                "L0902ERFCBW1CPF",
                "L0902SCNW2CPF")

# filter metadata 
metadata_filt <- subset(metadata_filt, Sample %notin% tmp_select)


### Re-Filter ASV Matrix Off Filtered Metadata with bad samples removed 

## Identify samples to keep
tmp_select <- which(colnames(asvmat_filt) %in% metadata_filt$Sample)

## Filter additional samples out of asvmat_filt again
asvmat_filt <- asvmat_filt[,tmp_select]

## Re-generate asvmat sample summary

# identify sample names to keep
tmp_select <- colnames(asvmat_filt)

# subset asvmat_filt_summary to only above samples
asvmat_filt_summary <- subset(asvmat_filt_summary, Sample %in% tmp_select)


### Create Cyanobacteria filtered ASV matrix

## Identify ASVs that correspond to cyanobacteria
tmp_select <- subset(taxtab, Phylum == "Cyanobacteria")$ASV_ID

## Identify rows in asvmat_filt that are not cyanobacteria 
tmp_select <- which(rownames(asvmat_filt) %notin% tmp_select)

## Remove Cyano ASVs from asvmat_filt
asvmat_filt_nocyano <- asvmat_filt[tmp_select,]

## Create sample summaries and assess low count samples
asvmat_filt_nocyano_summary <- sample.summary(asvmat_filt_nocyano, features = "rows", stds = pred_stand$ID)

# ***NOTE: The lowest number of counts is 642 which may be not totally unreasonable for rarefaction of this particular data

## Assess correlation between sample counts with and without cyanobacteria

# plot
plot(asvmat_filt_summary$Cts, asvmat_filt_nocyano_summary$Cts)

# cor
cor.test(asvmat_filt_summary$Cts, asvmat_filt_nocyano_summary$Cts, method = "pearson")
cor.test(asvmat_filt_summary$Cts, asvmat_filt_nocyano_summary$Cts, method = "spearman")


### Re-Create ASV data.frames from Matricies With and Without cyanobacerial ASVs

## w/ Cyano

# Get rownames
tmp_select_row <- which(asvtab_filt$ASV_ID %in% rownames(asvmat_filt))

# Get colnames
tmp_select_col <- which(colnames(asvtab_filt) %in% colnames(asvmat_filt))

# subset asvtab_filt
asvtab_filt <- asvtab_filt[tmp_select_row, c(1, tmp_select_col), with = F]

## w/o Cyano

# Get rownames
tmp_select_row <- which(asvtab_filt$ASV_ID %in% rownames(asvmat_filt_nocyano))

# Get colnames
tmp_select_col <- which(colnames(asvtab_filt) %in% colnames(asvmat_filt_nocyano))

# subset asvtab_filt
asvtab_filt_nocyano <- asvtab_filt[tmp_select_row, c(1, tmp_select_col), with = F]



## w/Cyano

# identify matrix indexes
tmp_drop <- which(rownames(asvmat_filt_rare) %in% pred_stand$ID)

# filter matrix 
asvmat_filt_rare <- asvmat_filt_rare[-tmp_drop,]

## w/o Cyano 

# identify matrix indexes
tmp_drop <- which(rownames(asvmat_filt_rare_nocyano) %in% pred_stand$ID)

# filter matrix 
asvmat_filt_rare_nocyano <- asvmat_filt_rare_nocyano[-tmp_drop,]

########################################################  Figure 1 B ############################################################################

###
### Area Plots for Directional Change and Community Summaries 
####


### Remove STD ASVs from ASV Tables (data.frames)

## w/ Cyano
asvtab_filt_area <- subset(asvtab_filt, ASV_ID %notin% pred_stand$ID)

## w/o Cyano 
asvtab_filt_nocyano_area <- subset(asvtab_filt_nocyano, ASV_ID %notin% pred_stand$ID)


### Melt ASV Data

## w/ Cyano 
asvtab_filt_area_long <- melt(asvtab_filt_area, variable.name = "Sample", value.name = "Counts")

## w/o Cyano 
asvtab_filt_nocyano_area_long <- melt(asvtab_filt_nocyano_area, variable.name = "Sample", value.name = "Counts")


### Merge in Phylum and Order Level Assignments

## Get ASV Taxonomy
tax_small <- taxtab[,c("ASV_ID", "Phylum", "Class", "Order")]

## w/ Cyano

# merge in taxonomy
asvtab_filt_area_long <- merge(asvtab_filt_area_long, tax_small, by = "ASV_ID", all.x = T)

# create "other" group for low abundance taxonomic groups
asvtab_filt_area_long$Phy_Small <- forcats::fct_lump_n(asvtab_filt_area_long$Phylum, 7, w = asvtab_filt_area_long$Counts, other_level = "Other")
asvtab_filt_area_long$Ord_Small <- forcats::fct_lump_n(asvtab_filt_area_long$Order, 12, w = asvtab_filt_area_long$Counts, other_level = "Other")

# replace cyano orders with "Cyanobacteria" and Preserve Armatimonadetes
asvtab_filt_area_long$Ord_Small <- sub("Phormidesmiales", "Cyanobacteria", asvtab_filt_area_long$Ord_Small)
asvtab_filt_area_long$Ord_Small <- sub("Synechococcales", "Cyanobacteria", asvtab_filt_area_long$Ord_Small)
asvtab_filt_area_long$Ord_Small <- sub("Cyanobacteriales", "Cyanobacteria", asvtab_filt_area_long$Ord_Small)
asvtab_filt_area_long$Ord_Small <- ifelse(asvtab_filt_area_long$Order == "Fimbriimonadales", "Fimbriimonadales", asvtab_filt_area_long$Ord_Small)


## w/o Cyano 

# merge in taxonomy
asvtab_filt_nocyano_area_long <- merge(asvtab_filt_nocyano_area_long, tax_small, by = "ASV_ID", all.x = T)

# create "other" group for low abundance taxonomic groups
asvtab_filt_nocyano_area_long$Phy_Small <- forcats::fct_lump_n(asvtab_filt_nocyano_area_long$Phylum, 6, w = asvtab_filt_nocyano_area_long$Counts, other_level = "Other")
asvtab_filt_nocyano_area_long$Ord_Small <- forcats::fct_lump_n(asvtab_filt_nocyano_area_long$Order, 11, w = asvtab_filt_nocyano_area_long$Counts, other_level = "Other")

# replace cyano orders with "Cyanobacteria" and Preserve Armatimonadetes
asvtab_filt_nocyano_area_long$Ord_Small <- sub("Phormidesmiales", "Cyanobacteria", asvtab_filt_nocyano_area_long$Ord_Small)
asvtab_filt_nocyano_area_long$Ord_Small <- sub("Synechococcales", "Cyanobacteria", asvtab_filt_nocyano_area_long$Ord_Small)
asvtab_filt_nocyano_area_long$Ord_Small <- sub("Cyanobacteriales", "Cyanobacteria", asvtab_filt_nocyano_area_long$Ord_Small)
asvtab_filt_nocyano_area_long$Ord_Small <- ifelse(asvtab_filt_nocyano_area_long$Order == "Fimbriimonadales", "Fimbriimonadales", asvtab_filt_nocyano_area_long$Ord_Small)


### Merge in Metadata

## w/ Cyano
asvtab_filt_area_long <- merge(asvtab_filt_area_long, metadata_filt, by = "Sample", all.x = T)

## w/o Cyano
asvtab_filt_nocyano_area_long <- merge(asvtab_filt_nocyano_area_long, metadata_filt, by = "Sample", all.x = T)


######

########################################################  Figure 1 B ############################################################################

######


### Construct Area Plots at various taxonomic levels faceted by different factors

# all together

# summarize data
area_ord_site <- data.frame(asvtab_filt_area_long %>%
                              group_by( Ord_Small, Total_Days) %>%
                              summarize(Total = sum(Counts)))
#area plot 
ggplot(area_ord_site, aes(x = Total_Days, y = Total, group = reorder(Ord_Small, -Total))) +
  geom_area(aes(fill = Ord_Small), stat = "identity", position = "fill", linewidth = 0.2, colour = "black") +
  #geom_vline(aes(xintercept = 84)) +
  #facet_wrap(.~General_Site,nrow = 3) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "top",
        strip.background = element_rect(colour = "black") ,
        legend.key.size = unit(0.3, "cm") ,
        axis.text.x = element_text(colour = "black")) +
  #labs(x = "Days from Start", y = "Relative Abundance (%)", fill = "Order") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,14,28,42,56,70,84,91,105,119)) +
  scale_fill_brewer(palette = "Set3")

ggsave(paste0(plot_out,"Plot8_Area_Order_All.pdf"), width = 7, height = 10, units = "in")   




### Calculate a ton of diversity statistics for the samples + make dataframe

## w/Cyano
asvmat_filt_rare_alpha <- microbiome::alpha(x = asvmat_filt_rare, index = "all", zeroes = TRUE) # inclusion of zeros seems to have no impact

## w/o Cyano 
asvmat_filt_rare_nocyano_alpha <- microbiome::alpha(x = asvmat_filt_rare_nocyano, index = "all", zeroes = TRUE) 


### Assess correlative relatioships of diversity metrics

## w/Cyano
corrplot::corrplot(cor(asvmat_filt_rare_alpha, method = "spearman"))

## w/o Cyano
corrplot::corrplot(cor(asvmat_filt_rare_nocyano_alpha, method = "spearman"))


### Merge in metadata

## w/Cyano
asvmat_filt_rare_alpha_plot <- data.frame(Sample = rownames(asvmat_filt_rare_alpha), asvmat_filt_rare_alpha)
asvmat_filt_rare_alpha_plot <- merge(asvmat_filt_rare_alpha_plot, metadata_filt, by = "Sample")

## w/o Cyano
asvmat_filt_rare_nocyano_alpha_plot <- data.frame(Sample = rownames(asvmat_filt_rare_nocyano_alpha), asvmat_filt_rare_nocyano_alpha)
asvmat_filt_rare_nocyano_alpha_plot <- merge(asvmat_filt_rare_nocyano_alpha_plot, metadata_filt, by = "Sample")


### Perform Linear Mixed Effects Modeling for Repeated Measures Data on Lagged Beta Diversity + Time Series 

## Add factorized data

# w/Cyano
asvmat_filt_rare_alpha_plot$Total_Days_Fac <- as.factor(asvmat_filt_rare_alpha_plot$Total_Days)

# w/o Cyano
asvmat_filt_rare_nocyano_alpha_plot$Total_Days_Fac <- as.factor(asvmat_filt_rare_nocyano_alpha_plot$Total_Days)





######


########################################################  Figure 1 C ############################################################################



####
### Beta-Diversity Analysis
####


### Pre-Filter and perform Zero Imputation

## w/Cyano  

# only keep ASVs with >= 2 postive values across all samples
good_asv <- apply(asvmat_filt, 1 , function(x) sum(x > 0)) >= 2
sum(good_asv)
asvmat_filt_goodasv <- asvmat_filt[good_asv,]
rm(good_asv)

write.table(asvmat_filt_goodasv, file = "~/Berkeley_Postdoc/UniCom_analyses/Tables/asvmat_filt_goodasv.txt", sep = "\t", quote = TRUE, row.names = TRUE, col.names = NA)

# impute zeros using zCompositions (needs samples in rows and ASVs in columns)
asvmat_filt_no0 <- cmultRepl(t(asvmat_filt_goodasv), output = "p-counts")
asvmat_filt_no0 <- t(asvmat_filt_no0)

## w/o Cyano

# only keep ASVs with >= 2 postive values across all samples
good_asv <- apply(asvmat_filt_nocyano, 1 , function(x) sum(x > 0)) >= 2
sum(good_asv)
asvmat_filt_nocyano_goodasv <- asvmat_filt_nocyano[good_asv,]
rm(good_asv)

# impute zeros using zCompositions (needs samples in rows and ASVs in columns)
asvmat_filt_nocyano_no0 <- cmultRepl(t(asvmat_filt_nocyano_goodasv), output = "p-counts")
asvmat_filt_nocyano_no0 <- t(asvmat_filt_nocyano_no0)


### Perform ALR Transformation

## w/Cyano
asvmat_filt_no0_alr <- alr.tfm(asvmat_filt_no0, features = "rows", stds = pred_stand$ID)

## w/o Cyano
asvmat_filt_nocyano_no0_alr <- alr.tfm(asvmat_filt_nocyano_no0, features = "rows", stds = pred_stand$ID)


### Calculate Achinson Distance

## w/Cyano
asvmat_filt_no0_alr_dist <- vegan::vegdist(t(asvmat_filt_no0_alr), method = "euclidian")

## w/o Cyano
asvmat_filt_nocyano_no0_alr_dist <- vegan::vegdist(t(asvmat_filt_nocyano_no0_alr), method = "euclidian")


### Generate Table of Lagged Beta-Diversity Comparisons - w/Cyano

## Create parsable matrix out of distance object
asvmat_filt_no0_alr_dist_mat <- as.matrix(asvmat_filt_no0_alr_dist)

## Get tube numbers
tube_numbers <- unique(metadata_filt$Tube_No)

## Time points for evaluation
#tp_eval <- c(0,14,28,42,56,70,84)
tp_eval <- c(0,14,28,42,56,70,84,91,105,119)

## Create df to hold output
beta.long.out <- data.frame()

## For loop code to pull distances by tube number 
for (i in 1:length(tube_numbers)){
  
  # Get metadata for tube i
  tmp_tube_dat <- subset(metadata_filt, Tube_No == tube_numbers[i])
  
  # if not full dataset then skip to next tube (Dont use this does not make much difference)
  #if(sum(tp_eval %in% tmp_tube_dat$Total_Days) != 7){
  #  next
  #}
  
  # run for loop to pull beta-div data from tube i
  for (j in 1:(length(tp_eval) - 1)){
    
    # grab sample names 
    tmp_sample1 <- subset(tmp_tube_dat, Total_Days == tp_eval[j])$Sample
    tmp_sample2 <- subset(tmp_tube_dat, Total_Days == tp_eval[j+1])$Sample
    
    # check if either sample is missing and if so skip to next iteration otherwise get the data
    if(length(tmp_sample1) == 1 & length(tmp_sample2) == 1){
      
      # get number of comparison
      tmp_compare <- j
      
      # get number of tube
      tmp_tube_no <- unique(tmp_tube_dat$Tube_No)
      
      # get ach dist
      tmp_dist <- asvmat_filt_no0_alr_dist_mat[tmp_sample1,tmp_sample2]
      
      tmp_out <- data.frame(Sample1 = tmp_sample1,
                            Sample2 = tmp_sample2,
                            Tube_No = tmp_tube_no,
                            Comparison = tmp_compare,
                            Ach_Dist = tmp_dist)
      
      # bind data to output df 
      beta.long.out <- rbind(beta.long.out, tmp_out)
      
    } else {
      next
    }
  }
}

## Merge in metadata to beta.long.out
beta.long.out <- merge(beta.long.out, good_samples[,-c(7,8)], by = "Tube_No", all.x = T)
beta.long.out$Tube_No_Char <- as.character(beta.long.out$Tube_No)
beta.long.out$Comparison_Fac <- as.factor(beta.long.out$Comparison)



## w/Cyano

# run LME with factorized time points
beta_mod <- lmer(Ach_Dist ~ General_Site + Strain + Passage_Rate + Comparison_Fac + (1|Tube_No), data = beta.long.out)
summary(beta_mod)
plot(beta_mod)

# assess term significance and model fit 
anova(beta_mod)
write.table(data.frame(anova(beta_mod)), file = "lme_stats/beta_mod_anova.txt", sep = "\t", quote = F)
r.squaredGLMM(beta_mod) 

# get estimated marginal means of each comparison
beta_mod_em <- emmeans(beta_mod, ~ Comparison_Fac)
beta_mod_em

# get statistics of pairwise comparisons
beta_mod_em_cont <- contrast(beta_mod_em, method = "pairwise")
beta_mod_em_cont

# generate compact letter display for statistically different values between comparisons
beta_mod_em_cld <- cld(beta_mod_em, Letters = letters, , reversed = T)
print(beta_mod_em_cld)
write.table(data.frame(beta_mod_em_cld), file = "~/Berkeley_Postdoc/Unicom_analyses/Figures/figures_all/beta_mod_em_cld.txt", sep = "\t", quote = F, row.names = F)



### Plot Data From Lagged Beta Diversity Modeling

## w/ CyAno

# extract emmeans data for plotting
beta_emmeans <- as.data.frame(beta_mod_em)
colnames(beta_emmeans)[1] <- "Comparison"


# plot data by Comparison
ggplot() +
  geom_jitter(data = beta.long.out, aes(x = as.factor(Comparison), y = Ach_Dist), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.5, size = 3) +
  geom_errorbar(data = beta_emmeans, aes(x = Comparison, ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_point(data = beta_emmeans, aes(x = Comparison, y = emmean), shape = 22, size = 3, fill = "#756bb1") +
  geom_line(data = beta_emmeans, aes(x = as.numeric(Comparison), y = emmean), linetype = 2, color = "#756bb1") +
  #geom_smooth() +
  scale_fill_brewer(palette = "Set1") +
  xlab("Comparison") +
  ylab("Delta Achinson Distance") +
  theme_bw()

ggsave(paste0(plot_out,"Plot5_deltaBeta_wPF.pdf"), width = 5, height = 3, units = "in")



### Perform PERMANOVA analysis on beta-diversity as a function of metadata variables 

## w/ Cyano
adonis_out <- vegan::adonis2(asvmat_filt_no0_alr_dist ~ General_Site + Strain + Passage_Rate + as.factor(Total_Days), data = metadata_filt, by = "margin", permutations = 999, parallel = 8)
adonis_out
write.table(data.frame(adonis_out), file = "lme_stats/adonis_out.txt", sep = "\t", quote = F)

## w/o Cyano
adonis_nocyano_out <- vegan::adonis2(asvmat_filt_nocyano_no0_alr_dist ~ General_Site + Strain + Passage_Rate + as.factor(Total_Days), data = metadata_filt, by = "margin", permutations = 999, parallel = 8)
adonis_nocyano_out
write.table(data.frame(adonis_nocyano_out), file = "lme_stats/adonis_nocyano_out.txt", sep = "\t", quote = F)





########################################################  Figure 1 D ############################################################################

## Richness - w/o Cyano



# run LME with factorized time points
alpha_obs_nocyano_mod <- lmer(observed ~ General_Site + Strain + Passage_Rate + Total_Days_Fac + (1|Tube_No), data = asvmat_filt_rare_nocyano_alpha_plot)
summary(alpha_obs_nocyano_mod)
plot(alpha_obs_nocyano_mod)

# assess term significance and model fit 
anova(alpha_obs_nocyano_mod)
write.table(data.frame(anova(alpha_obs_nocyano_mod)), file = "lme_stats/alpha_obs_nocyano_mod_anova.txt", sep = "\t", quote = F)
r.squaredGLMM(alpha_obs_nocyano_mod)

# get estimated marginal means of time
alpha_obs_nocyano_em <- emmeans(alpha_obs_nocyano_mod, ~ Total_Days_Fac)
alpha_obs_nocyano_em_df <- as.data.frame(alpha_obs_nocyano_em)

# get statistics of pairwise comparisons
alpha_obs_nocyano_em_cont <- contrast(alpha_obs_nocyano_em, method = "pairwise")
alpha_obs_nocyano_em_cont

# generate compact letter display for statistically different values within each General_Site
alpha_obs_nocyano_em_cld <- cld(alpha_obs_nocyano_em, Letters = letters, , reversed = F)
print(alpha_obs_nocyano_em_cld)
write.table(data.frame(alpha_obs_nocyano_em_cld), file = "~/Berkeley_Postdoc/Unicom_analyses/Figures/figures_all/alpha_obs_nocyano_em_cld.txt", sep = "\t", quote = F, row.names = F)

#
### Plot Data

## Richness - w/Cyano

ggplot() +
  geom_jitter(data = asvmat_filt_rare_alpha_plot, aes(x = Total_Days_Fac, y = observed), fill = "grey90",height = 0, width = 0.1, shape = 21, alpha = 0.5, size = 3) +
  geom_errorbar(data = alpha_obs_em_df, aes(x = Total_Days_Fac, ymin = lower.CL, ymax = upper.CL), width = 0.2,) +
  geom_point(data = alpha_obs_em_df, aes(x = Total_Days_Fac, y = emmean), shape = 22, size = 3, fill = "steelblue4") +
  geom_line(data = alpha_obs_em_df, aes(x = as.numeric(Total_Days_Fac), y = emmean), linetype = 2, color = "steelblue4") +
  #geom_smooth() +
  xlab("Time From Inoculation (d)") +
  ylab("Richness") +
  theme_bw() 

ggsave(paste0(plot_out,"Plot1_Richness_wPF.pdf"), width = 5, height = 3, units = "in")






##
#################################################### Import  Data 2. ####################################################
####

####singleM analysis data####



#### order enrichment analysis ####


####
### Import All Data
####


### OTU Count Table
otutab_95 <- fread("~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/singleM_Analysis/Analysis_Groups/All_Samples_noPilot/all_samples_nopilot_nr95_OTU.txt")

## Choose otutab
otutab <- otutab_95

## Fix Col Names
correct_names <- read.table("~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/singleM_Analysis/Analysis_Groups/All_Samples_noPilot/correct_col_names.txt", header = T)
colnames(otutab) <- correct_names$otutab_95


### Sample Metadata Table
metadata <- fread("~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/singleM_Analysis/Analysis_Groups/All_Samples_noPilot/All_Full_Sample_Metadata_wCov.txt")

### Create Graphic Output Directory Variable
output <- "~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/singleM_Analysis/Analysis_Groups/All_Samples_noPilot/Plots/"


###
### Filter/Fix Input Data to Specific Marker of Interest (Rpl6)
####


### Subset OTU Table to robust marker set and add Taxonomy

## Specify robust markers + subset OTU table
robust_markers <- c("S3.39.ribosomal_protein_L6_rplF")

otutab_robust <- subset(otutab, marker %in% robust_markers)

## Add unique OTU Names
OTU_ID <- paste0("OTU_ID_", 1:nrow(otutab_robust)) 
otutab_robust$OTU_ID <- OTU_ID


## Create split rank taxonomy

# data frame to hold split rank taxa
tax <- data.frame()

# split taxonomy and bind to tax using for loop (takes a while)
for (i in 1:nrow(otutab_robust)){
  
  tmp_OTU <- otutab_robust$OTU_ID[i]
  names(tmp_OTU) <- "OTU_ID"
  
  tmp_tax <- phyloseq::parse_taxonomy_qiime(otutab_robust$taxonomy[i])
  
  tmp_tax <- c(tmp_OTU, tmp_tax)
  
  tax <- bind_rows(tax, tmp_tax)
}

# make all NA into "Unknown"
tax[is.na(tax)] <- "Unknown"

# remove non-informative taxonomy columns and rename
tax <- tax[,-c(2,10)]
colnames(tax) <- c("OTU_ID","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# merge into otu_table
otutab_robust <- merge(otutab_robust, tax, by = "OTU_ID")


### Finalize Data Frame for L6 OTUs
otu_df_L6 <- subset(otutab_robust, marker == "S3.39.ribosomal_protein_L6_rplF")


####
### Sample Composition Assessment and Relative Abundance Plotting - Non-Compositional 
####


### Make Count Data into Matrices For  Marker

## Matrix for rpL6
otumat_L6 <- subset(otutab_robust, marker == "S3.39.ribosomal_protein_L6_rplF")
otu_names <- otumat_L6$OTU_ID
otumat_L6 <- as.matrix(otumat_L6[,4:137]) #adjust based on numeric columns in otumat
rownames(otumat_L6) <- otu_names


### Calculate Summary Statistics

## Summaries for L6
otumat_L6_summaries <- sample.summary(otumat_L6, features = "rows")

## Output Summary Table for L6
write.table(otumat_L6_summaries, "~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/singleM_Analysis/Analysis_Groups/Source_Samples_noPilot/L6_OTU_Sample_Summaries.txt", quote = F, sep = "\t", row.names = F)


### Assessment of sequencing and marker recovery depth

## Make df with important stats
seq_depth_var <- data.frame(min_seq_depth = min(metadata$Total),
                            max_seq_depth = max(metadata$Total),
                            sd_seq_depth = sd(metadata$Total),
                            spread_seq_depth = max(metadata$Total) / min(metadata$Total),
                            min_L6_depth = min(otumat_L6_summaries$Cts),
                            max_L6_depth = max(otumat_L6_summaries$Cts),
                            sd_L6_depth = sd(otumat_L6_summaries$Cts),
                            spread_L6_depth = max(otumat_L6_summaries$Cts) / min(otumat_L6_summaries$Cts),
                            min_L6_obs = min(otumat_L6_summaries$Has_Cts),
                            max_L6_obs = max(otumat_L6_summaries$Has_Cts),
                            sd_L6_obs = sd(otumat_L6_summaries$Has_Cts),
                            spread_L6_obs = max(otumat_L6_summaries$Has_Cts) / min(otumat_L6_summaries$Has_Cts))

## Transpose df for readabilty and write to disk
seq_depth_var <- t(seq_depth_var)
write.table(seq_depth_var, "~/Dropbox/Banfield_Lab_Files/Projects/Cyano_Moore/Studies/2021_09_16_Com_Construction_Full/singleM_Analysis/Analysis_Groups/Source_Samples_noPilot/L6_OTU_Seq_Depth_Var.txt", quote = F, sep = "\t", row.names = T)


### Filter out bad samples + source (based on metagenomic)

## Sample names to drop 
sample_drop <- c("unicom2021_Eel_S_A7120_W_B",
                 "unicom2021_StrawCreek_N_A7120_W_C",
                 "unicom2021_Dbay_Inner_3055_BW_C",
                 "unicom2021_Dbay_Inner_3055_BW_B",
                 "unicom2021_StrawCreek_S_A7120_W_C",
                 "unicom2021_Eel_Fox_7942_BW_C",
                 "unicom2021_Dbay_Inner_7942_W_A",
                 "unicom2021_Dbay_Inner_3055_BW_A",
                 "unicom2021_Dbay_Mouth_A7120_BW_C",
                 "unicom2021_Eel_Fox_7942_W_C",
                 "unicom2021_Eel_S_7942_BW_B",
                 "unicom2021_Eel_Fox_7942_BW_B",
                 "unicom2021_StrawCreek_S_Inocula")

## Source sample names
source_drop <- subset(metadata, Sample_Type == "Source")$Sample

## Sample names to keep 
sample_keep <- which(colnames(otu_df_L6) %notin% unique(c(sample_drop, source_drop)))

## Create Filtered Count Dataframe
otu_df_L6_filtered_permute <- otu_df_L6[,..sample_keep]

########################################################  Figure 1 B ############################################################################




### Remove Cyanobacterial OTUs 
otu_df_L6_filtered_permute_noCyano <- subset(otu_df_L6_filtered_permute, Phylum != "Cyanobacteria")


### Generate Data That Will Be Sampled + Used as Background Reference (Should still have ALL sample OTUs but Cyanobacterial ones)
otu_permute_data_noCyano <- otu_df_L6_filtered_permute_noCyano[,c("OTU_ID","Phylum","Class","Order","Family","Genus")]


### Quantify Number of Unique OTUs per Sample

## Melt abundance data for easier summary
otu_df_L6_filtered_permute_noCyano_long <- pivot_longer(otu_df_L6_filtered_permute_noCyano,
                                                        cols = starts_with(c("unicom")),
                                                        names_to = "Sample",
                                                        values_to = "Counts")
## Count unique number of OTUs per sample
otu_per_sample_noCyano <- data.table(otu_df_L6_filtered_permute_noCyano_long %>% 
                                       group_by(Sample) %>%
                                       summarise(num_spec = sum(Counts > 0),
                                                 num_counts = sum(Counts)))


### Run For Loop to Generate Permuted Samples 

## Set seed
set.seed(123)

## Set number of permutations
num_perm <- 10000

## Create list to hold final permutation data
permute_otu_master <- vector("list", num_perm)

## Permutation for loop
for(i in 1:num_perm){
  
  ## sample random OTUs
  tmp_otus <- sample(x = otu_permute_data_noCyano$OTU_ID, size = sum(otu_per_sample_noCyano$num_spec), replace = F) 
  
  ## Split into list of vectors
  otu_lists <- split(tmp_otus, rep(1:nrow(otu_per_sample_noCyano), otu_per_sample_noCyano$num_spec))
  
  ## bind to permute_sample_genomes
  tmp_df <- data.frame(Sample = rep(otu_per_sample_noCyano$Sample, otu_per_sample_noCyano$num_spec),
                       OTU_ID = unlist(otu_lists, use.names = F),
                       iteration = i)
  
  permute_otu_master[[i]] <- tmp_df
  
  # print iteration
  print(paste0("iteration: ",i))
}

# Combine all the data.frames stored in the list into one data.frame
permute_otu_master <- dplyr::bind_rows(permute_otu_master)


### Merge in Taxonomic Data to permute_otu_master
permute_otu_master <- merge(permute_otu_master, otu_permute_data_noCyano, by = "OTU_ID")


### Identify Actual Frequencies of Taxa in Culture Samples

## Order Level

# gather presence absence data per sample 
actual_order_frequencies <- data.table(otu_df_L6_filtered_permute_noCyano_long %>%
                                         group_by(Sample, Order) %>%
                                         summarise(total_cts = sum(Counts),
                                                   present = ifelse(total_cts >= 1, 1, 0)))

# summarize presence absence data
actual_order_frequencies <- data.table(actual_order_frequencies %>%
                                         group_by(Order) %>%
                                         summarise(total_obs = sum(present),
                                                   obs_freq = total_obs / 108))


### Summarize Frequencies of Taxa in Permuted Samples

## Order Level

# gather presence absence data per iteration per sample
permute_order_frequencies <- data.table(permute_otu_master %>%
                                          group_by(iteration, Sample, Order) %>%
                                          summarise(total_obs = n(),
                                                    present = ifelse(total_obs >= 1, 1, 0)))
# summarize presence absence data
permute_order_frequencies <- data.table(permute_order_frequencies %>%
                                          group_by(iteration, Order) %>%
                                          summarise(total_obs = sum(present),
                                                    obs_freq = total_obs / 108))


### Summarize Enrichment of Taxa in Cultures vs Permuted

## Order Level

# build table to hold results
order_enrich_stats <- data.frame()

for (i in 1:nrow(actual_order_frequencies)) {
  
  # get tmp order name
  tmp_order_name <- actual_order_frequencies$Order[i]
  
  # get permuted data on order i
  tmp_permuted <- subset(permute_order_frequencies, Order == tmp_order_name)
  
  # add zeros for orders missing in permutation
  if(nrow(tmp_permuted) < num_perm) {
    
    # get which iterations are missing
    tmp_all_iterations <- c(1:num_perm)
    tmp_missing_iterations <- tmp_all_iterations[tmp_all_iterations %notin% tmp_permuted$iteration]
    
    # create fake data to bind
    tmp_fake_iterations <- data.frame(iteration = tmp_missing_iterations,
                                      Order = tmp_order_name,
                                      total_obs = 0,
                                      obs_freq = 0)
    # bind to tmp_permuted
    tmp_permuted <- rbind(tmp_permuted, tmp_fake_iterations)
    
  }
  
  # calculate summary stats
  mean_perm_obs <- mean(tmp_permuted$total_obs)
  sd_perm_obs <- sd(tmp_permuted$total_obs)
  p1_more_extreme <- sum(tmp_permuted$total_obs >= actual_order_frequencies$total_obs[i]) / num_perm
  p2_less_extreme <- sum(tmp_permuted$total_obs <= actual_order_frequencies$total_obs[i]) / num_perm
  tmp_z <- (actual_order_frequencies$total_obs[i] - mean_perm_obs) / sd_perm_obs
  p_val <- min(p1_more_extreme, p2_less_extreme) * 2
  p_val <- ifelse(p_val > 1, 1, p_val)
  
  # make output data
  tmp_output <- data.frame(Order = tmp_order_name,
                           actual_obs = actual_order_frequencies$total_obs[i],
                           mean_perm_obs = mean_perm_obs,
                           sd_perm_obs = sd_perm_obs,
                           Z_score = tmp_z,
                           P_val = p_val)
  
  # bind output 
  order_enrich_stats <- rbind(order_enrich_stats, tmp_output)
}

# Calculate false discovery rate
order_enrich_stats$fdr <- p.adjust(order_enrich_stats$P_val, method = "fdr")

# Identify significant features
order_enrich_stats$sig <- ifelse(order_enrich_stats$fdr <= 0.05, TRUE, FALSE)

# Create column with enriched group name
order_enrich_stats$group <- ifelse(order_enrich_stats$Z_score > 0, "Culture", "Source")


### Plot Data

## Generate Phylum to Order Taxonomy

# identify unique phylum-order combos
phy_order <- unique(tax[,c("Order", "Phylum")])

# remove unknown
phy_order <- subset(phy_order, Order != "Unknown")

## Merge in phyla info
order_enrich_stats_plot <- merge(order_enrich_stats, phy_order, by = "Order", all.x = T)

## Filter down orders to just those with significant enrichment
order_enrich_stats_plot_sig <- subset(order_enrich_stats_plot, fdr < 0.05)

## Count phylum level lineages enriched in culture and source
tmp_check_counts <- plyr::count(order_enrich_stats_plot_sig, vars = c("group","Phylum"))

## Set up colors for plotting
phy_colors <- c("Acidobacteriota"  = "#5EFFA1",
                "Actinobacteriota"  = "#2444FF",
                "Armatimonadota"  = "#7082B2",
                "Bacteroidota"  = "#FDE4A5",
                #"p__Bdellovibrionota" = "#3e6641",
                #"p__Cyanobacteria" = "#4FFF00",
                #"p__Deinococcota" = "#F6CAF2",
                #p__Desulfobacterota_F " = "#AAA430",
                #p__Eremiobacterota" = "#FC4E07",
                "Chlamydiota" = "#F00000",
                #"p__Chloroflexota"  = "#E8861A",
                #"p__CSP1-3" = "#964B00",
                "Gemmatimonadota" = "#EA15D6",
                "Firmicutes"  = "#687a6a",
                "Omnitrophota" = "#79C5D1", 
                #"p__Myxococcota" = "#00AFBB",
                "Patescibacteria"  = "#F6CAF2",
                "Proteobacteria"  = "#B7B180",
                "Planctomycetota"  = "#4DB6FF",
                "Spirochaetota"  = "#F36579",
                "Verrucomicrobiota"  = "#B4FF00",
                "Other" = "grey80")

## generate vector of phyla (not OTHER)
phy_names <- c("Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota", "Chlamydiota",
               "Gemmatimonadota", "Firmicutes", "Omnitrophota", "Patescibacteria", "Proteobacteria",
               "Planctomycetota", "Spirochaetota", "Verrucomicrobiota")

## add column of simplified phylum names to order_enrich_stats_plot_sig
order_enrich_stats_plot_sig$Phylum_simple <- forcats::fct_lump_n(order_enrich_stats_plot_sig$Phylum, n = 12, other_level = "Other")

install.packages("scales")
library(scales)

## Plot Data
ggplot(order_enrich_stats_plot_sig, aes(x = reorder(Order, Z_score), y = Z_score , fill = Phylum_simple)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  #geom_hline(aes(yintercept = 2.995732), linetype = 2) +
  #xlim(c(-25, 25)) +
  #ylim(c(0, 7.5)) +
  #scale_x_continuous(trans = pseudolog10_trans, limits = c(-60,60)) +
  scale_y_continuous(trans = pseudolog10_trans) +
  scale_fill_manual(values = phy_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90)) +
  ylab("Z Score") +
  xlab(NULL)

ggsave(paste0(output,"Plot18_Order_Enrich_Sig.pdf"), width = 12, height = 5)


## Export Permutation Statistical Significance Data
write.table(order_enrich_stats_plot, "Order_Enrich_Stats.txt", quote = F, row.names = F, sep ="\t")

# ## Plot
# ggplot(order_enrich_stats_plot, aes(x = Z_score, y = -log(fdr) , color = sig)) +
#   geom_point(size = 4, alpha = 0.6) +
#   #geom_hline(aes(yintercept = 2.995732), linetype = 2) +
#   #xlim(c(-25, 25)) +
#   #ylim(c(0, 7.5)) +
#   scale_x_continuous(trans = pseudolog10_trans, limits = c(-60,60)) +
#   scale_y_continuous(trans = pseudolog10_trans) +
#   theme_bw()
#   
# 
# ## Plot2
# ggplot(order_enrich_stats_plot, aes(x = reorder(Order, Z_score), y = Z_score , fill = sig)) +
#   geom_bar(stat = "identity") +
#   #geom_hline(aes(yintercept = 2.995732), linetype = 2) +
#   #xlim(c(-25, 25)) +
#   #ylim(c(0, 7.5)) +
#   #scale_x_continuous(trans = pseudolog10_trans, limits = c(-60,60)) +
#   scale_y_continuous(trans = pseudolog10_trans) +
#   theme_bw()





