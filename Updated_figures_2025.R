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
library(vegan)


####
####################################### Dataset 1  - 16S rRNA Amplicon - Load All The Data ######################################
####

#######################################################  Figure 1  and Figure S3 ############################################################# 

### ASV Table
asvtab <- fread("~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/zotutab.txt")

### Metadata
metadata <- fread("~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/metadata_sd.txt")

good_samples <- fread("~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/good_metagenome_metadata.txt")

## check all good samples are in metadata
test_metadata <- subset(metadata, Tube_No %in% good_samples$tube_no)
rm(test_metadata)

### Taxonomy Data
taxtab <- fread("~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/taxonomy.tsv")

### Predicted Standards
pred_stand <- fread("~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/predicted_standards.txt")
pred_stand <- subset(pred_stand, PRED_STD == T)

# ####metadata all samples from table S3
# metadata_TableS3 <- fread("~/Berkeley_Postdoc/UniCom_analyses/Figures/codes/metadata_all_samples_tableS3.txt") 

### Plot Output Location
plot_out <- "~/Berkeley_Postdoc/Unicom_analyses/Figures/figures_all/"
table_out <- "~/Berkeley_Postdoc/Unicom_analyses/Figures/table_out/"



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

asv_tax <- merge(asvtab, taxtab, by ="ASV_ID")


####
### Filter Data + Create ASV Count Matrix + Generate Data Summary
####



### Specify Samples We Want By Filtering Metadata - Only taking samples at bi-weekly intervals
metadata_filt <- subset(metadata, Total_Days %in% c(14,28,42,56,70,84,91,105,119))

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

## Remove low count samples (less than 1000 counts)
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

#  
### Re-Filter ASV Matrix Off Filtered Metadata with bad samples removed 

## Identify samples to keep
tmp_select <- which(colnames(asvmat_filt) %in% metadata_filt$Sample)

## Filter additional samples out of asvmat_filt again
asvmat_filt <- asvmat_filt[,tmp_select]

asvmat_filt.df <- data.frame(ASV_ID = rownames(asvmat_filt), asvmat_filt, row.names = NULL)
asvmat_filt.df <- left_join(asvmat_filt.df, taxtab, by ="ASV_ID")
## Export summary data
write.table(
  data.frame(asvmat_filt.df), 
  file = file.path(table_out, "asvmat_filt.df.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)

## Re-generate asvmat sample summary

# identify sample names to keep
tmp_select <- colnames(asvmat_filt)

# subset asvmat_filt_summary to only above samples
asvmat_filt_summary <- subset(asvmat_filt_summary, Sample %in% tmp_select)

## Export summary data
write.table(
  data.frame(asvmat_filt_summary), 
  file = file.path(table_out, "asvmat_filt_summary.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)


### Create Cyanobacteria filtered ASV matrix

## Identify ASVs that correspond to cyanobacteria
tmp_select <- subset(taxtab, Phylum == "Cyanobacteria")$ASV_ID

## Identify rows in asvmat_filt that are not cyanobacteria 
tmp_select <- which(rownames(asvmat_filt) %notin% tmp_select)

## Remove Cyano ASVs from asvmat_filt
asvmat_filt_nocyano <- asvmat_filt[tmp_select,]

## Create sample summaries and assess low count samples
asvmat_filt_nocyano_summary <- sample.summary(asvmat_filt_nocyano, features = "rows", stds = pred_stand$ID)


## Export summary data
write.table(
  data.frame(asvmat_filt_nocyano_summary), 
  file = file.path(table_out, "asvmat_filt_nocyano_summary.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)

###mean counts
mean(asvmat_filt_nocyano_summary$Cts) #14565.29
mean(asvmat_filt_summary$Cts) #46098

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


####
### Alpha Diversity Analysis
####


### Rarefy data to normalize between samples (needs samples in rows and ASVs in columns)
# Rarefy to the samples with lowest counts number
## Set random seed
set.seed(123)

## Rarefy - wCyano
asvmat_filt_rare <- GUniFrac::Rarefy(t(asvmat_filt), depth = 4746) # this function applies base::sample and samples without replacement
asvmat_filt_rare <- t(asvmat_filt_rare[[1]])

## Rarefy - w/o Cyano
asvmat_filt_rare_nocyano <- GUniFrac::Rarefy(t(asvmat_filt_nocyano), depth = 642)
asvmat_filt_rare_nocyano <- t(asvmat_filt_rare_nocyano[[1]])


### Remove STD ASVs

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


### Calculate a ton of diversity statistics for the samples + make dataframe

## w/Cyano
asvmat_filt_rare_alpha <- microbiome::alpha(x = asvmat_filt_rare, index = "all", zeroes = TRUE) # inclusion of zeros seems to have no impact
write.table(
  data.frame(asvmat_filt_rare_alpha), 
  file = file.path(table_out, "alpha_div_stats_all.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)

## w/o Cyano 
asvmat_filt_rare_nocyano_alpha <- microbiome::alpha(x = asvmat_filt_rare_nocyano, index = "all", zeroes = TRUE) 
write.table(
  data.frame(asvmat_filt_rare_nocyano_alpha), 
  file = file.path(table_out, "alpha_div_stats_nocyano.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)

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


### Perform Linear Mixed Effects Modeling for Time Series + Repeated Measures Data on Lagged Beta Diversity


## Add factorized data

# w/Cyano
asvmat_filt_rare_alpha_plot$Total_Days_Fac <- as.factor(asvmat_filt_rare_alpha_plot$Total_Days)

# w/o Cyano
asvmat_filt_rare_nocyano_alpha_plot$Total_Days_Fac <- as.factor(asvmat_filt_rare_nocyano_alpha_plot$Total_Days)

## Richness - wCyano

# run LME with factorized time points
alpha_obs_mod <- lmer(observed ~ General_Site + Strain + Passage_Rate + Total_Days_Fac + (1|Tube_No), data = asvmat_filt_rare_alpha_plot)
summary(alpha_obs_mod)
plot(alpha_obs_mod)

# Save summary output to a text file
capture.output(summary(alpha_obs_mod), 
               file = file.path(table_out, "alpha_obs_mod_all_lme.txt"))

# assess term significance and model fit 
anova(alpha_obs_mod)
write.table(
  data.frame(anova(alpha_obs_mod)), 
  file = file.path(table_out, "alpha_obs_mod_all_anova.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = TRUE
)

r.squaredGLMM(alpha_obs_mod)

# get estimated marginal means of time
alpha_obs_em <- emmeans(alpha_obs_mod, ~ Total_Days_Fac)
alpha_obs_em_df <- as.data.frame(alpha_obs_em)

# get statistics of pairwise comparisons
alpha_obs_em_cont <- contrast(alpha_obs_em, method = "pairwise")
alpha_obs_em_cont


# generate compact letter display for statistically different values within each General_Site, reload packages if not working
alpha_obs_em_cld <- cld(alpha_obs_em, Letters = letters, , reversed = F)
print(alpha_obs_em_cld)
write.table(
  data.frame(alpha_obs_em_cld), 
  file = file.path(table_out, "alpha_obs_em_cld_all.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)
## Richness - w/o Cyano

# run LME with factorized time points
alpha_obs_nocyano_mod <- lmer(observed ~ General_Site + Strain + Passage_Rate + Total_Days_Fac + (1|Tube_No), data = asvmat_filt_rare_nocyano_alpha_plot)
summary(alpha_obs_nocyano_mod)
plot(alpha_obs_nocyano_mod)
# Save summary output to a text file
capture.output(summary(alpha_obs_nocyano_mod), 
               file = file.path(table_out, "alpha_obs_mod_nocyano_lme.txt"))

# assess term significance and model fit 
anova(alpha_obs_nocyano_mod)
write.table(
  data.frame(anova(alpha_obs_nocyano_mod)), 
  file = file.path(table_out, "alpha_obs_nocyano_mod_anova.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)
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
write.table(
  data.frame(alpha_obs_nocyano_em_cld), 
  file = file.path(table_out, "alpha_obs_nocyano_em_cld.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)
## Shannon Diversity - w/Cyano

# run LME with factorized time points
alpha_shan_mod <- lmer(diversity_shannon ~ General_Site + Strain + Passage_Rate + Total_Days_Fac + (1|Tube_No), data = asvmat_filt_rare_alpha_plot)
summary(alpha_shan_mod)
plot(alpha_shan_mod)

# assess term significance and model fit 
anova(alpha_shan_mod)
write.table(
  data.frame(anova(alpha_shan_mod)), 
  file = file.path(table_out, "alpha_shan_mod_anova_all.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)
r.squaredGLMM(alpha_shan_mod)

# get estimated marginal means of time
alpha_shan_em <- emmeans(alpha_shan_mod, ~ Total_Days_Fac)
alpha_shan_em_df <- as.data.frame(alpha_shan_em)

# get statistics of pairwise comparisons
alpha_shan_em_cont <- contrast(alpha_shan_em, method = "pairwise")
alpha_shan_em_cont

# generate compact letter display for statistically different values within each General_Site
alpha_shan_em_cld <- cld(alpha_shan_em, Letters = letters, , reversed = F)
print(alpha_shan_em_cld)
write.table(
  data.frame(alpha_shan_em_cld), 
  file = file.path(table_out, "alpha_shan_em_cld_all.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)
## Shannon Diversity - w/o Cyano

# run LME with factorized time points
alpha_shan_nocyano_mod <- lmer(diversity_shannon ~ General_Site + Strain + Passage_Rate + Total_Days_Fac + (1|Tube_No), data = asvmat_filt_rare_nocyano_alpha_plot)
summary(alpha_shan_nocyano_mod)
plot(alpha_shan_nocyano_mod)

# assess term significance and model fit 
anova(alpha_shan_nocyano_mod)
write.table(
  data.frame(anova(alpha_shan_nocyano_mod)), 
  file = file.path(table_out, "alpha_shan_nocyano_mod_anova.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)
r.squaredGLMM(alpha_shan_nocyano_mod)

# get estimated marginal means of time
alpha_shan_nocyano_em <- emmeans(alpha_shan_nocyano_mod, ~ Total_Days_Fac)
alpha_shan_nocyano_em_df <- as.data.frame(alpha_shan_nocyano_em)

# get statistics of pairwise comparisons
alpha_shan_nocyano_em_cont <- contrast(alpha_shan_nocyano_em, method = "pairwise")
alpha_shan_nocyano_em_cont

# generate compact letter display for statistically different values within each General_Site
alpha_shan_nocyano_em_cld <- cld(alpha_shan_nocyano_em, Letters = letters, , reversed = F)
print(alpha_shan_nocyano_em_cld)
write.table(
  data.frame(alpha_shan_nocyano_em_cld), 
  file = file.path(table_out, "alpha_shan_nocyano_em_cld.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)


### Plot Data

## Richness - w/Cyano

###Figure 1 D
########################################################  Figure 1 D ############################################################################


ggplot() +
  geom_jitter(data = asvmat_filt_rare_alpha_plot, aes(x = Total_Days_Fac, y = observed), fill = "grey90",height = 0, width = 0.1, shape = 21, alpha = 0.5, size = 3) +
  geom_errorbar(data = alpha_obs_em_df, aes(x = Total_Days_Fac, ymin = lower.CL, ymax = upper.CL), width = 0.2,) +
  geom_point(data = alpha_obs_em_df, aes(x = Total_Days_Fac, y = emmean), shape = 22, size = 3, fill = "steelblue4") +
  geom_line(data = alpha_obs_em_df, aes(x = as.numeric(Total_Days_Fac), y = emmean), linetype = 2, color = "steelblue4") +
  #geom_smooth() +
  xlab("Time From Inoculation (d)") +
  ylab("Richness") +
  theme_bw() 

ggsave(paste0(plot_out,"Figure_1_D_Richness_wPF_noT0.pdf"), width = 5, height = 3, units = "in")

## Richness - w/o Cyano
###Figure S3 G
########################################################  Figure S3 G ############################################################################

ggplot() +
  geom_jitter(data = asvmat_filt_rare_nocyano_alpha_plot, aes(x = Total_Days_Fac, y = observed), fill = "grey90",height = 0, width = 0.1, shape = 21, alpha = 0.5, size = 3) +
  geom_errorbar(data = alpha_obs_nocyano_em_df, aes(x = Total_Days_Fac, ymin = lower.CL, ymax = upper.CL), width = 0.2,) +
  geom_point(data = alpha_obs_nocyano_em_df, aes(x = Total_Days_Fac, y = emmean), shape = 22, size = 3, fill = "steelblue4") +
  geom_line(data = alpha_obs_nocyano_em_df, aes(x = as.numeric(Total_Days_Fac), y = emmean), linetype = 2, color = "steelblue4") +
  #geom_smooth() +
  xlab("Time From Inoculation (d)") +
  ylab("Richness") +
  theme_bw() 

ggsave(paste0(plot_out,"Figure_S3_G_Richness_nocyano_wPF_noT0.pdf"), width = 5, height = 3, units = "in")

## Shannon Diversity - w/Cyano
###Figure S3 E

########################################################  Figure S3 E ############################################################################

ggplot() +
  geom_jitter(data = asvmat_filt_rare_alpha_plot, aes(x = Total_Days_Fac, y = diversity_shannon), fill = "grey90",height = 0, width = 0.05, shape = 21, alpha = 0.5, size = 3) +
  geom_errorbar(data = alpha_shan_em_df, aes(x = Total_Days_Fac, ymin = lower.CL, ymax = upper.CL), width = 0.2,) +
  geom_point(data = alpha_shan_em_df, aes(x = Total_Days_Fac, y = emmean), shape = 22, size = 3, fill = "firebrick4") +
  geom_line(data = alpha_shan_em_df, aes(x = as.numeric(Total_Days_Fac), y = emmean), linetype = 2, color = "firebrick4") +
  #geom_smooth() +
  xlab("Time From Inoculation (d)") +
  ylab("Shannon Diversity") +
  theme_bw() 

ggsave(paste0(plot_out,"Figure_S3_E_Shannon_wPF_noT0.pdf"), width = 5, height = 3, units = "in")  

## Shannon Diversity - w/o Cyano
###Figure S3 F

########################################################  Figure S3 F ############################################################################

ggplot() +
  geom_jitter(data = asvmat_filt_rare_nocyano_alpha_plot, aes(x = Total_Days_Fac, y = diversity_shannon), fill = "grey90",height = 0, width = 0.05, shape = 21, alpha = 0.5, size = 3) +
  geom_errorbar(data = alpha_shan_nocyano_em_df, aes(x = Total_Days_Fac, ymin = lower.CL, ymax = upper.CL), width = 0.2,) +
  geom_point(data = alpha_shan_nocyano_em_df, aes(x = Total_Days_Fac, y = emmean), shape = 22, size = 3, fill = "firebrick4") +
  geom_line(data = alpha_shan_nocyano_em_df, aes(x = as.numeric(Total_Days_Fac), y = emmean), linetype = 2, color = "firebrick4") +
  #geom_smooth() +
  xlab("Time From Inoculation (d)") +
  ylab("Shannon Diversity") +
  theme_bw() 

ggsave(paste0(plot_out,"Figure_S3_F_Shannon_nocyano_wPF_noT0.pdf"), width = 5, height = 3, units = "in")  



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
## Create parsable matrix out of distance object w/o Cyano
asvmat_filt_nocyano_no0_alr_dist_mat <- as.matrix(asvmat_filt_nocyano_no0_alr_dist)


## Get tube numbers
tube_numbers <- unique(metadata_filt$Tube_No)
unique(metadata_filt$Tube_No)

## Time points for evaluation
tp_eval <- c(14,28,42,56,70,84,91,105,119)

## Create df to hold output
beta.long.out <- data.frame()

# Adjusted for-loop for beta diversity comparisons
beta.long.out <- data.frame()

for (i in 1:length(tube_numbers)) {
  
  # Filter tube-specific metadata
  tmp_tube_dat <- subset(metadata_filt, Tube_No == tube_numbers[i])
  
  # Skip tubes with insufficient data
  if (nrow(tmp_tube_dat) < 2) {
    next
  }
  
  # Iterate over time points for evaluation
  for (j in 1:(length(tp_eval) - 1)) {
    
    # Select samples for current and next time points
    tmp_sample1 <- subset(tmp_tube_dat, Total_Days == tp_eval[j])$Sample
    tmp_sample2 <- subset(tmp_tube_dat, Total_Days == tp_eval[j + 1])$Sample
    
    # Skip if either sample is missing or excluded
    if (length(tmp_sample1) == 0 || length(tmp_sample2) == 0) {
      next
    }
    
    # Check if samples exist in the distance matrix
    if (tmp_sample1 %in% rownames(asvmat_filt_no0_alr_dist_mat) &&
        tmp_sample2 %in% colnames(asvmat_filt_no0_alr_dist_mat)) {
      
      # Extract relevant information
      tmp_compare <- j
      tmp_tube_no <- unique(tmp_tube_dat$Tube_No)
      tmp_dist <- asvmat_filt_no0_alr_dist_mat[tmp_sample1, tmp_sample2]
      
      # Create output row
      tmp_out <- data.frame(
        Sample1 = tmp_sample1,
        Sample2 = tmp_sample2,
        Tube_No = tmp_tube_no,
        Comparison = tmp_compare,
        Ach_Dist = tmp_dist
      )
      
      # Append to output
      beta.long.out <- rbind(beta.long.out, tmp_out)
    }
  }
}
# 

# 
## Merge in metadata to beta.long.out
beta.long.out <- merge(beta.long.out, good_samples[,-c(7,8)], by = "Tube_No", all.x = T)
beta.long.out$Tube_No_Char <- as.character(beta.long.out$Tube_No)
beta.long.out$Comparison_Fac <- as.factor(beta.long.out$Comparison)


### Generate Table of Lagged Beta-Diversity Comparisons - w/o Cayno

# Adjusted for-loop for beta diversity comparisons
beta.long.nocyano.out <- data.frame()

for (i in 1:length(tube_numbers)) {
  
  # Filter tube-specific metadata
  tmp_tube_dat <- subset(metadata_filt, Tube_No == tube_numbers[i])
  
  # Skip tubes with insufficient data
  if (nrow(tmp_tube_dat) < 2) {
    next
  }
  
  # Iterate over time points for evaluation
  for (j in 1:(length(tp_eval) - 1)) {
    
    # Select samples for current and next time points
    tmp_sample1 <- subset(tmp_tube_dat, Total_Days == tp_eval[j])$Sample
    tmp_sample2 <- subset(tmp_tube_dat, Total_Days == tp_eval[j + 1])$Sample
    
    # Skip if either sample is missing or excluded
    if (length(tmp_sample1) == 0 || length(tmp_sample2) == 0) {
      next
    }
    
    # Check if samples exist in the distance matrix
    if (tmp_sample1 %in% rownames(asvmat_filt_nocyano_no0_alr_dist_mat) &&
        tmp_sample2 %in% colnames(asvmat_filt_nocyano_no0_alr_dist_mat)) {
      
      # Extract relevant information
      tmp_compare <- j
      tmp_tube_no <- unique(tmp_tube_dat$Tube_No)
      tmp_dist <- asvmat_filt_nocyano_no0_alr_dist_mat[tmp_sample1, tmp_sample2]
      
      # Create output row
      tmp_out <- data.frame(
        Sample1 = tmp_sample1,
        Sample2 = tmp_sample2,
        Tube_No = tmp_tube_no,
        Comparison = tmp_compare,
        Ach_Dist = tmp_dist
      )
      
      # Append to output
      beta.long.nocyano.out <- rbind(beta.long.nocyano.out, tmp_out)
    }
  }
}

## Merge in metadata to beta.long.nocyano.out
beta.long.nocyano.out <- merge(beta.long.nocyano.out, good_samples[,-c(7,8)], by = "Tube_No", all.x = T)
beta.long.nocyano.out$Tube_No_Char <- as.character(beta.long.nocyano.out$Tube_No)
beta.long.nocyano.out$Comparison_Fac <- as.factor(beta.long.nocyano.out$Comparison)


### Perform Linear Mixed Effects Modeling for Time Series + Repeated Measures Data on Lagged Beta Diversity

## w/Cyano
library(emmeans)
# run LME with factorized time points
beta_mod <- lmer(Ach_Dist ~ General_Site + Strain + Passage_Rate + Comparison_Fac + (1|Tube_No), data = beta.long.out)
summary(beta_mod)
plot(beta_mod)

# assess term significance and model fit 
anova(beta_mod)
write.table(
  data.frame(anova(beta_mod)), 
  file = file.path(table_out, "beta_mod_anova.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)
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
write.table(
  data.frame(beta_mod_em_cld), 
  file = file.path(table_out, "beta_mod_em_cld.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)
## w/o Cyano 

# run LME with factorized time points
beta_nocyano_mod <- lmer(Ach_Dist ~ General_Site + Strain + Passage_Rate + Comparison_Fac + (1|Tube_No), data = beta.long.nocyano.out)
summary(beta_nocyano_mod)
plot(beta_nocyano_mod)

# assess term significance and model fit 
anova(beta_nocyano_mod)
write.table(
  data.frame(anova(beta_nocyano_mod)), 
  file = file.path(table_out, "beta_nocyano_mod_anova.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)
r.squaredGLMM(beta_nocyano_mod) 

# get estimated marginal means of each comparison
beta_nocyano_mod_em <- emmeans(beta_nocyano_mod, ~ Comparison_Fac)
beta_nocyano_mod_em

# get statistics of pairwise comparisons
beta_nocyano_mod_em_cont <- contrast(beta_nocyano_mod_em, method = "pairwise")
beta_nocyano_mod_em_cont

# generate compact letter display for statistically different values between comparisons
beta_nocyano_mod_em_cld <- cld(beta_nocyano_mod_em, Letters = letters, , reversed = T)
print(beta_nocyano_mod_em_cld)
write.table(
  data.frame(beta_nocyano_mod_em_cld), 
  file = file.path(table_out, "beta_nocyano_mod_em_cld.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)

### Plot Data From Lagged Beta Diversity Modeling

## w/ Cayno

# extract emmeans data for plotting
beta_emmeans <- as.data.frame(beta_mod_em)
colnames(beta_emmeans)[1] <- "Comparison"

###Figure 1 C

########################################################  Figure 1 C ############################################################################

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

ggsave(paste0(plot_out,"Figure_1C_deltaBeta_wPF_noT0.pdf"), width = 5, height = 3, units = "in")

## w/o Cayno
###Figure S3 H
########################################################  Figure S3 H ############################################################################


# extract emmeans data for plotting
beta_emmeans_nocyano <- as.data.frame(beta_nocyano_mod_em)
colnames(beta_emmeans_nocyano)[1] <- "Comparison"


# plot data by Comparison
ggplot() +
  geom_jitter(data = beta.long.nocyano.out, aes(x = as.factor(Comparison), y = Ach_Dist), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.5, size = 3) +
  geom_errorbar(data = beta_emmeans_nocyano, aes(x = Comparison, ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_point(data = beta_emmeans_nocyano, aes(x = Comparison, y = emmean), shape = 22, size = 3, fill = "#756bb1") +
  geom_line(data = beta_emmeans_nocyano, aes(x = as.numeric(Comparison), y = emmean), linetype = 2, color = "#756bb1") +
  #geom_smooth() +
  scale_fill_brewer(palette = "Set1") +
  xlab("Comparison") +
  ylab("Delta Achinson Distance") +
  theme_bw()

ggsave(paste0(plot_out,"Figure_S3_H_deltaBeta_nocyano_wPF_noT0.pdf"), width = 5, height = 3, units = "in")   


### Perform PERMANOVA analysis on beta-diversity as a function of metadata variables 
#subset metadata
metadata_filt_matched <- metadata_filt[metadata_filt$Sample %in% rownames(asvmat_filt_no0_alr_dist_mat), ]

## w/ Cyano
adonis_out <- vegan::adonis2(asvmat_filt_no0_alr_dist ~ General_Site + Strain + Passage_Rate + as.factor(Total_Days), data = metadata_filt_matched, by = "margin", permutations = 999, parallel = 8)
adonis_out
write.table(
  data.frame(adonis_out), 
  file = file.path(table_out, "adonis_out.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = TRUE
)
## w/o Cyano
metadata_filt_matched_nocyano <- metadata_filt[metadata_filt$Sample %in% rownames(asvmat_filt_nocyano_no0_alr_dist), ]

adonis_nocyano_out <- vegan::adonis2(asvmat_filt_nocyano_no0_alr_dist ~ General_Site + Strain + Passage_Rate + as.factor(Total_Days), data = metadata_filt_matched_nocyano, by = "margin", permutations = 999, parallel = 8)
adonis_nocyano_out
write.table(
  data.frame(adonis_nocyano_out), 
  file = file.path(table_out, "adonis_nocyano_out.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = TRUE
)



####
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


### Construct Area Plots at various taxonomic levels faceted by different factors

## w/ Cyano

# Order by General Site

# summarize data
area_ord_site <- data.frame(asvtab_filt_area_long %>%
                              group_by(General_Site, Ord_Small, Total_Days) %>%
                              summarize(Total = sum(Counts)))
#area plot 

unique(area_ord_site$Ord_Small)
ord_col <-c ("Burkholderiales"= "#8dd3c7",
             "Chitinophagales" = "#ffffb3",
             "Cyanobacteria" = "#35b779" ,
             "Cytophagales" = "#bdb8d7",
             "Fimbriimonadales" = "#80b1d3",
             "Flavobacteriales" = "#fdb462",
             "Rhizobiales" = "#fccde5",
             "Rhodobacterales" = "#d9d9d9",
             "Sphingobacteriales" = "#5300aa",
             "Sphingomonadales" = "#ba80b5",
             "Xanthomonadales" = "#ffed6f",
             "Other" = "#979797")

###Figure S3 A

########################################################  Figure S3 A ############################################################################


ggplot(area_ord_site, aes(x = Total_Days, y = Total, group = reorder(Ord_Small, -Total))) +
  geom_area(aes(fill = Ord_Small), stat = "identity", position = "fill", linewidth = 0.2, colour = "black") +
  #geom_vline(aes(xintercept = 84)) +
  facet_wrap(.~General_Site,nrow = 3) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none",#replace toe "right" or "top" if you want to have legend as part of the column
        strip.background = element_rect(fill = "white", colour = "black"), # White box with black border
        strip.text = element_text(colour = "black") , # Black text inside the box
        #legend.key.size = unit(0.3, "cm") ,
        axis.text.x = element_text(colour = "black")) +
  #labs(x = "Days from Start", y = "Relative Abundance (%)", fill = "Order") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,14,28,42,56,70,84,91,105,119)) +
  scale_fill_manual(values = ord_col)

ggsave(paste0(plot_out,"Figure_S3_A_Area_Order_Site.pdf"), width = 5, height = 7, units = "in")
#ggsave(paste0(plot_out,"Plot8_Area_Order_Site_legend.pdf"), width = 4, height = 7, units = "in")



# Calculate the average abundance across all sites
area_averaged <- area_ord_site %>%
  group_by(Total_Days, Ord_Small) %>%
  summarize(Average_Total = mean(Total), .groups = "drop")

###Figure 1 B

########################################################  Figure 1 B ############################################################################

# Plot the averaged data
ggplot(area_averaged, aes(x = Total_Days, y = Average_Total, group = reorder(Ord_Small, -Average_Total))) +
  geom_area(aes(fill = Ord_Small), stat = "identity", position = "fill", linewidth = 0.2, colour = "black") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    axis.text.x = element_text(colour = "black")
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 14, 28, 42, 56, 70, 84, 91, 105, 119)) +
  scale_fill_manual(values = ord_col) +
  labs(x = "Days from Start", y = "Average Relative Abundance", fill = "Order")


ggsave(paste0(plot_out,"Figure_1_B_Area_Order_Site_Average.pdf"), width = 10, height = 6, units = "in")



## w/o Cyano
# Order by General Site

# summarize data
area_ord_site_nocyano <- data.frame(asvtab_filt_nocyano_area_long %>%
                                      group_by(General_Site, Ord_Small, Total_Days) %>%
                                      summarize(Total = sum(Counts)))
unique(area_ord_site_nocyano$Ord_Small)
ord_col_nocyano <-c ("Acetobacterales"= "#c95d38",
                     "Caulobacterales" = "#7ed9fc",
                     "Burkholderiales"= "#8dd3c7",
                     "Chitinophagales" = "#ffffb3",
                     "Cyanobacteria" = "#35b779" ,
                     "Cytophagales" = "#bdb8d7",
                     "Fimbriimonadales" = "#80b1d3",
                     "Flavobacteriales" = "#fdb462",
                     "Rhizobiales" = "#fccde5",
                     "Rhodobacterales" = "#d9d9d9",
                     "Sphingobacteriales" = "#5300aa",
                     "Sphingomonadales" = "#ba80b5",
                     "Xanthomonadales" = "#ffed6f",
                     "Other" = "#979797")

###Figure S3 B

########################################################  Figure S3 B ############################################################################

#area plot of non cyanobacterial part
ggplot(area_ord_site_nocyano, aes(x = Total_Days, y = Total, group = reorder(Ord_Small, -Total))) +
  geom_area(aes(fill = Ord_Small), stat = "identity", position = "fill", linewidth = 0.2, colour = "black") +
  #geom_vline(aes(xintercept = 84)) +
  facet_wrap(.~General_Site,nrow = 3) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        #strip.background = element_rect(colour = "black") ,
        strip.background = element_rect(fill = "white", colour = "black"), # White box with black border
        legend.key.size = unit(0.3, "cm") ,
        axis.text.x = element_text(colour = "black")) +
  #labs(x = "Days from Start", y = "Relative Abundance (%)", fill = "Order") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0), breaks = c(0,14,28,42,56,70,84,91,105,119)) +
  scale_fill_manual(values = ord_col_nocyano)

ggsave(paste0(plot_out,"Plot9_Area_Order_Site_nocyano.pdf"), width = 5, height = 7, units = "in")

ggsave(paste0(plot_out,"Figure_S3_B_Area_Order_Site_nocyano.pdf"), width = 5, height = 7, units = "in")



####
### Run Differential Abundance Analysis on T0 -- > T84 (BF) && T84 --> T119 (AF)
####


### Load Library
library(Maaslin2)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Maaslin2")


### Parse Data Inputs for Maaslin2

## Subset Metadata of Interest 

# subset metadata
#maaslin2_metadata <- as.data.frame(subset(metadata_filt, Total_Days %in% c(0,14,28,42,56,70,84,91,105,119)))
maaslin2_metadata <- as.data.frame(subset(metadata_filt, Total_Days %in% c(14,28,42,56,70,84,91,105,119)))


# create rows named by sample
rownames(maaslin2_metadata) <- maaslin2_metadata$Sample

# keep only columns we want
maaslin2_metadata <- maaslin2_metadata[,c("Strain","Tube_No","General_Site", "Passage_Rate", "Total_Days")]

# fix stuff in data
maaslin2_metadata$General_Site <- as.factor(sub(" ", "_", maaslin2_metadata$General_Site))
maaslin2_metadata$Strain <- as.factor(maaslin2_metadata$Strain)
maaslin2_metadata$Passage_Rate <- as.factor(maaslin2_metadata$Passage_Rate)
# maaslin2_metadata$Total_Days_Fac <- factor(maaslin2_metadata$Total_Days, levels = c("84","14","28","42","56","70","0","91","105","119"))
maaslin2_metadata$Total_Days_Fac <- factor(maaslin2_metadata$Total_Days, levels = c("84","14","28","42","56","70","91","105","119"))
maaslin2_metadata$Tube_No <- as.factor(maaslin2_metadata$Tube_No)

## Subset ASV Count Matrices of Interest

# get columns
tmp_select <- which(colnames(asvmat_filt) %in% rownames(maaslin2_metadata))

# get data
maaslin2_data <- t(asvmat_filt[,tmp_select]) # wants samples in rows


### Remove low count ASVs

## only keep ASVs with >= 10 postive values across all samples
good_asv <- apply(maaslin2_data, 2 , function(x) sum(x > 0)) >= 10

## count how many ASVs will be retained 
sum(good_asv) # 1682 out of 2126 ASVs
              ###with no T0 1557

## filter maaslin2_data
maaslin2_data <- maaslin2_data[,good_asv]
rm(good_asv)


### Make ASV Matrix into Data_Frame
maaslin2_data <- data.frame(maaslin2_data)


### Run Maaslin2 and Summarize Output 

## Run Maaslin using default norm and model with Total Days == 84 as a reference point  
maaslin2_fit_TD84ref <- Maaslin2(input_data = maaslin2_data,
                                 input_metadata = maaslin2_metadata,
                                 output = "~/Desktop/maaslin2_out/",
                                 min_abundance = 0, 
                                 min_prevalence = 0,
                                 min_variance = 0,
                                 analysis_method = "LM",
                                 max_significance = 0.05,
                                 fixed_effects = c("General_Site", "Strain", "Passage_Rate", "Total_Days_Fac"),
                                 reference = c("General_Site,Discovery_Bay"),
                                 random_effects = c("Tube_No"),
                                 correction ="BH",
                                 cores = 12,
                                 plot_heatmap = F,
                                 plot_scatter = F)

## Extract filtered results of interest


###day 14 -  day 84 comparison  
maaslin2_fit_TD84ref_14d <- subset(maaslin2_fit_TD84ref$results, metadata == "Total_Days_Fac" & value == "14" & qval <= 0.05)


# 84d --> 119d Comparison
maaslin2_fit_TD84ref_119d <- subset(maaslin2_fit_TD84ref$results, metadata == "Total_Days_Fac" & value == "119" & qval <= 0.05)

## Add direction of coefficient relative to 84d time point

# 14d --> 84d Comparison
maaslin2_fit_TD84ref_14d$Coeff_dir <- ifelse(maaslin2_fit_TD84ref_14d$coef < 0, "Increase", "Decrease")

# 84d --> 119d Comparison
maaslin2_fit_TD84ref_119d$Coeff_dir <- ifelse(maaslin2_fit_TD84ref_119d$coef > 0, "Increase", "Decrease")

## Merge in taxonomic information

# 14d --> 84d Comparison
maaslin2_fit_TD84ref_14d <- merge(maaslin2_fit_TD84ref_14d, taxtab, by.x = "feature", by.y = "ASV_ID", all.x = T)

# 84d --> 119d Comparison
maaslin2_fit_TD84ref_119d <- merge(maaslin2_fit_TD84ref_119d, taxtab, by.x = "feature", by.y = "ASV_ID", all.x = T)

## Add in column for streamlined Order taxonomy (Top10)

# 14d --> 84d Comparison
maaslin2_fit_TD84ref_14d$Ord_Small <- forcats::fct_lump_n(maaslin2_fit_TD84ref_14d$Order, 10, other_level = "Other", ties.method = "first")
write.table(maaslin2_fit_TD84ref_14d, "~/Desktop/maaslin2_out/84vs14.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 84d --> 119d Comparison
maaslin2_fit_TD84ref_119d$Ord_Small <- forcats::fct_lump_n(maaslin2_fit_TD84ref_119d$Order, 10, other_level = "Other", ties.method = "first")
write.table(maaslin2_fit_TD84ref_119d, "~/Desktop/maaslin2_out/119vs84.txt", sep = "\t", row.names = FALSE, quote = FALSE)

unique(maaslin2_fit_TD84ref_119d$Ord_Small)
maaslin2_fit_TD84ref_14d
unique(maaslin2_fit_TD84ref_14d$Ord_Small)

### Plot Results as Volcano Scatter

## 14d --> 84d Comparison

# Set Colors
ord_col_volc <-c ("Acetobacterales"= "#8ccfc3",
                  "Alteromonadales" ="#210e45",
                  "Caulobacterales" = "#7ed9fc",
                  "Burkholderiales"= "#8dd3c7",
                  "Chitinophagales" = "#ffffb3",
                  "Cyanobacteriales" = "#35b779" ,
                  "Cytophagales" = "#bdb8d7",
                  "Fimbriimonadales" = "#80b1d3",
                  "Flavobacteriales" = "#fdb462",
                  "Rhizobiales" = "#fccde5",
                  "Rhodobacterales" = "#d9d9d9",
                  "Sphingobacteriales" = "#5300aa",
                  "Sphingomonadales" = "#ba80b5",
                  "Xanthomonadales" = "#ffed6f",
                  "Other" = "#979797")

# Plot 
##Figure_S3_C

########################################################  Figure S3 C ############################################################################

ggplot(maaslin2_fit_TD84ref_14d, aes(x = (-1 * coef), y = -log10(qval), fill = Ord_Small)) +
  geom_point(shape = 21, size = 3) +
  xlim(c(-6,6)) +
  scale_y_log10() +
  scale_fill_manual(values = ord_col_volc) +
  xlab("Coefficent (84 d / 14 d)") +
  
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none",##if you want legen add where
        strip.background = element_rect(colour = "black") ,
        legend.key.size = unit(0.3, "cm") ,
        axis.text.x = element_text(colour = "black")) 


ggsave(paste0(plot_out,"Figure_S3_C_maaslin_84dv14d.pdf"), width = 5, height = 4, units = "in")   
#ggsave(paste0(plot_out,"Plot10_Area_Order_Site_nocyano.pdf"), width = 7, height = 10, units = "in")

## 84d --> 119d Comparison

# Set Colors
##Figure_S3_D

########################################################  Figure S3 D ############################################################################

# Plot
ggplot(maaslin2_fit_TD84ref_119d, aes(x = coef, y = -log10(qval), fill = Ord_Small)) +
  geom_point(shape = 21, size = 3) +
  xlim(c(-4,4)) +
  scale_y_log10() +
  scale_fill_manual(values = ord_col_volc) +
  xlab("Coefficent (119 d / 84 d)") +
  
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        strip.background = element_rect(colour = "black") ,
        legend.key.size = unit(0.3, "cm") ,
        axis.text.x = element_text(colour = "black")) 
ggsave(paste0(plot_out,"Figure_S3_D_maaslin_119dv84d.pdf"),width = 5, height = 4, units = "in")   




####
####################################### Dataset 2  - singleM - Load All The Data ##################################################
####

#######################################################  Figure 2B  and Figure S1, S4, S5A ############################################################# 

### OTU Count Table
otutab_95 <- fread("~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/all_samples_nopilot_nr95_OTU.txt")

## Choose otutab
otutab <- otutab_95

## Fix Col Names
correct_names <- read.table("~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/correct_col_names.txt", header = T)
colnames(otutab) <- correct_names$otutab_95


### Sample Metadata Table
metadata <- fread("~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/All_Full_Sample_Metadata_wCov.txt")


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

### Function 4: Summarize Samples - Expects a matrix object where features are ASVs; also can accept a vector of ASVs that are standards 
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




####
### Pre-Analysis of Observed Markers and Counts  
####


### Summarize Markers by Sample

## Remove un-needed columns and melt
marker_summary <- otutab[,-c("sequence", "taxonomy")]
marker_summary <- melt(marker_summary)
colnames(marker_summary) <- c("marker", "sample", "counts") 

## Summarize markers by sample
marker_summary <- data.table(marker_summary %>% 
                               group_by(marker, sample) %>%
                               summarise(obs = sum(counts > 0),
                                         total = sum(counts)))
marker_summary <- merge(marker_summary, metadata, by.x = "sample", by.y = "Sample")

## Plot Results - Observations
ggplot(marker_summary, aes(x = marker, y = obs)) +
  geom_boxplot(outlier.shape = NA, fill = "steelblue3") +
  geom_hline(aes(yintercept = median(marker_summary$obs)), color = "firebrick", linetype = 2) +
  geom_jitter(shape = 21, fill = "grey50", alpha = 0.5, width = 0.2, height = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10() +
  ggtitle("Unique Observations of Marker per Sample")

ggsave(paste0(plot_out,"Plot1_Markers_per_Sample.pdf"), width = 20, height = 6)

## Plot Results - Total Counts
ggplot(marker_summary, aes(x = marker, y = total)) +
  geom_boxplot(outlier.shape = NA, fill = "steelblue3") +
  geom_hline(aes(yintercept = median(marker_summary$total)), color = "firebrick", linetype = 2) +
  geom_jitter(shape = 21, fill = "grey50", alpha = 0.5, width = 0.2, height = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10() +
  ggtitle("Total Counts of Marker per Sample")

ggsave(paste0(plot_out,"Plot2_Marker_Counts_per_Sample.pdf"), width = 20, height = 6)



## Assessment:
# Most markers have the same observed and count distribution
# There are some differences between source and culture samples likley due to singificantly increased diversity in the source
# This distribution is shared by the markers we are interested in: S9, L6, L3, and S3
# Some markers are poorly identified in this dataset 
# There seems to be a positive trend in both novel observations and increased counts with sequencing depth...we should investigate collectors curves
# with additional sequencing we could prob get more OTUs...investigate collectors curves




####
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


### Make Count Data into Matrices For Each Marker

## Matrix for rpL6
otumat_L6 <- subset(otutab_robust, marker == "S3.39.ribosomal_protein_L6_rplF")
otu_names <- otumat_L6$OTU_ID
otumat_L6 <- as.matrix(otumat_L6[,4:137]) #adjust based on numeric columns in otumat
rownames(otumat_L6) <- otu_names


### Calculate Summary Statistics

## Summaries for L6
otumat_L6_summaries <- sample.summary(otumat_L6, features = "rows")

## Output Summary Table for L6
write.table(
  data.frame(otumat_L6_summaries), 
  file = file.path(table_out, "L6_OTU_Sample_Summaries.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = TRUE
)


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
write.table(
  data.frame(seq_depth_var), 
  file = file.path(table_out, "/L6_OTU_Seq_Depth_Var.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = TRUE
)

### Plot Phylum Relative Abundance of Source Samples 

## Melt data 
otutab_robust_long <- melt(otutab_robust)

## Sum counts by marker, sample, and phylum
otutab_robust_phy_counts <- data.table(otutab_robust_long %>% 
                                         group_by(marker, variable, Phylum) %>%
                                         summarise(counts = sum(value)))

## Assess Phyla that make up < 1% of bulk composition (across all markers)
bulk_phyla <- aggregate(counts ~ Phylum, otutab_robust_phy_counts, sum)
bulk_phyla$frac <- bulk_phyla$counts / sum(bulk_phyla$counts)
bulk_phyla <- bulk_phyla[which(bulk_phyla$frac >= 0.01),]$Phylum

## Create new phylum variable for plotting
otutab_robust_phy_counts$Phylum_small <- ifelse(otutab_robust_phy_counts$Phylum %in% bulk_phyla, otutab_robust_phy_counts$Phylum, "Other")

## merge in metadata with normal sample names 
otutab_robust_phy_counts <- merge(otutab_robust_phy_counts, metadata, by.x = "variable", by.y = "Sample")
#otutab_robust_phy_counts$Phylum_small <- sub("Patescibacteria", "CPR",  otutab_robust_phy_counts$Phylum_small)

otutab_robust_phy_counts <- otutab_robust_phy_counts %>%
  mutate(Phylum_small = recode(Phylum_small, 
                               "Patescibacteria" = "CPR",
                               "Proteobacteria" = "Pseudomonadota"))

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

## Filter out samples that are low quality 
otutab_robust_phy_counts <- subset(otutab_robust_phy_counts, variable %notin% sample_drop)

otutab_robust_phy_counts <- subset(otutab_robust_phy_counts, Sample_Type =="Source" )

## Create phylum colors varaible
phy_cols <- c("Actinobacteriota" = "#2444FF",
              "Bacteroidota" = "#FDE4A5",
              "Bdellovibrionota" = "#3e6641",
              "Cyanobacteria" = "#35B779FF",
              "Omnitrophota" = "#43756B",
              "Patescibacteria" = "#EA15D6",
              "Planctomycetota" = "#4DB6FF",
              "Pseudomonadota" = "#B7B180",
              "Verrucomicrobiota" = "#B4FF00",
              "CPR" = "#79C5D1",
              "Other" = "grey")



# Convert 'Phylum_small' to a factor with levels in the desired order
desired_order = c("Planctomycetota", "Verrucomicrobiota", "Actinobacteriota", "Cyanobacteria",
                  "Pseudomonadota", "Bacteroidota", "CPR")  # specify all levels in order
otutab_robust_phy_counts$Phylum_small <- factor(otutab_robust_phy_counts$Phylum_small, levels = desired_order)

##################################################.  Figure S2 A ####################################################

# rpL6
ggplot(subset(otutab_robust_phy_counts, marker == "S3.39.ribosomal_protein_L6_rplF"), aes(x = variable, y = counts, fill = Phylum_small)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phy_cols) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 9, angle = 60, hjust = 1)) +
  
  ylab("Relative_Abundance") +
  xlab(NULL)

ggsave(paste0(plot_out,"Suplementray_Figure_2A_Phylum_Relative_Abundance_Source_singlem.pdf"), width = 10, height = 10)


####
### Alpha Diversity Analysis
####

# 
# ### Filter out bad samples (based on metagenomic)
# 
# ## Sample names to drop 
# sample_drop <- c("unicom2021_Eel_S_A7120_W_B",
#                  "unicom2021_StrawCreek_N_A7120_W_C",
#                  "unicom2021_Dbay_Inner_3055_BW_C",
#                  "unicom2021_Dbay_Inner_3055_BW_B",
#                  "unicom2021_StrawCreek_S_A7120_W_C",
#                  "unicom2021_Eel_Fox_7942_BW_C",
#                  "unicom2021_Dbay_Inner_7942_W_A",
#                  "unicom2021_Dbay_Inner_3055_BW_A",
#                  "unicom2021_Dbay_Mouth_A7120_BW_C",
#                  "unicom2021_Eel_Fox_7942_W_C",
#                  "unicom2021_Eel_S_7942_BW_B",
#                  "unicom2021_Eel_Fox_7942_BW_B",
#                  "unicom2021_StrawCreek_S_Inocula")

## Sample names to keep 
sample_keep <- which(colnames(otu_df_L6) %notin% sample_drop)

## Create filtered otumat_L6
otu_df_L6_filtered <- otu_df_L6[,..sample_keep]

## Remove OTUs that are all Zeros
otu_df_L6_filtered <- otu_df_L6_filtered[which(rowSums(otu_df_L6_filtered[,4:124]) != 0),]


### Create Data filtered for Cyanobacteria
otu_df_L6_filtered_noCyano <- subset(otu_df_L6_filtered, Phylum != "Cyanobacteria")


### Create Matricies out of filtered data

## L6_filtered
otu_names_tmp <- otu_df_L6_filtered$OTU_ID
otumat_L6_filtered <- as.matrix(otu_df_L6_filtered[,4:124]) #adjust based on numeric columns in otumat
rownames(otumat_L6_filtered) <- otu_names_tmp
rm(otu_names_tmp)

## L6_filtered_noCyano
otu_names_tmp <- otu_df_L6_filtered_noCyano$OTU_ID
otumat_L6_filtered_noCyano <- as.matrix(otu_df_L6_filtered_noCyano[,4:124]) #adjust based on numeric columns in otumat
rownames(otumat_L6_filtered_noCyano) <- otu_names_tmp
rm(otu_names_tmp)


### Assess Summaries of Filtered Samples 

## L6_Filtered
otumat_L6_filtered_summary <- sample.summary(otumat_L6_filtered, features = "rows")

## L6_Filtered_noCyano
otumat_L6_filtered_noCyano_summary <- sample.summary(otumat_L6_filtered_noCyano, features = "rows") ### Prob dont use this as this depletes culture samples too much



### Rarefy data for alpha diversity metric calculation

## Set random seed
set.seed(123)

## Rarefy to lowest sample depth
otumat_L6_filtered_rare <- GUniFrac::Rarefy(t(otumat_L6_filtered), depth = min(otumat_L6_filtered_summary$Cts)) # this function applies the function "base::sample" and samples without replacement
otumat_L6_filtered_rare <- t(otumat_L6_filtered_rare[[1]])


### Calculate Alpha Diversity Stastics

## many alpha-stats from microbiome package
otumat_L6_filtered_rare_alpha <- microbiome::alpha(x = otumat_L6_filtered_rare, index = "all", zeroes = TRUE) # inclusion of zeros seems to have no impact 


### Assess correlative relatioships of diversity metrics
corrplot::corrplot(cor(otumat_L6_filtered_rare_alpha, method = "spearman"))


### Merge in metadata
otumat_L6_filtered_rare_alpha_plot <- data.frame(Sample = rownames(otumat_L6_filtered_rare_alpha), otumat_L6_filtered_rare_alpha)
otumat_L6_filtered_rare_alpha_plot <- merge(otumat_L6_filtered_rare_alpha_plot, metadata, by = "Sample")


### Check numerical distribtuions of data
plot(density(otumat_L6_filtered_rare_alpha_plot$observed)) # unimodal with extreme values
#plot(density(otumat_L6_filtered_rare_alpha_plot$veg_rare)) # unimodal with extreme values
plot(density(otumat_L6_filtered_rare_alpha_plot$chao1)) # unimodal with extreme values
plot(density(otumat_L6_filtered_rare_alpha_plot$diversity_gini_simpson)) # bimodal
plot(density(otumat_L6_filtered_rare_alpha_plot$diversity_inverse_simpson)) # unimodal with extreme values
plot(density(otumat_L6_filtered_rare_alpha_plot$diversity_shannon)) # unimodal with extreme values


### Analysis of Source Samples Individually (Relative to Each Other)

## Subset Data to Source Only
otumat_L6_filtered_rare_alpha_source_plot <- subset(otumat_L6_filtered_rare_alpha_plot, Sample_Type == "Source")

## Summarize group means

# Richness/Observed Diversity
source_rich_summary <- rcompanion::groupwiseMean(observed ~ Location, data = otumat_L6_filtered_rare_alpha_source_plot, bca = T)

# Shannon Diversity
source_shan_summary <- rcompanion::groupwiseMean(diversity_shannon ~ Location, data = otumat_L6_filtered_rare_alpha_source_plot, bca = T)

## Plot Diversity Metrics


############################################## Figure_S4_C ###################################################################

# Richness/Observed Diversity
ggplot() +
  geom_jitter(data = otumat_L6_filtered_rare_alpha_source_plot, aes(x = Location, y = observed), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = source_rich_summary, aes(x = Location, ymin = Bca.lower, ymax = Bca.upper), width = 0.1) +
  geom_point(data = source_rich_summary, aes(x = Location, y = Mean, fill = Location), shape = 22, size = 5) +
  xlab(NULL) +
  ylab("Richness") +
  ylim(c(0,600)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")

ggsave(paste0(plot_out,"Suplementary_Figure_4C_Alpha_Source_Richness.pdf"), width = 2.5, height = 4)

############################################## Figure_S4_D ###################################################################


# Shannon Diversity
ggplot() +
  geom_jitter(data = otumat_L6_filtered_rare_alpha_source_plot, aes(x = Location, y = diversity_shannon), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = source_shan_summary, aes(x = Location, ymin = Bca.lower, ymax = Bca.upper), width = 0.1) +
  geom_point(data = source_shan_summary, aes(x = Location, y = Mean, fill = Location), shape = 22, size = 5) +
  xlab(NULL) +
  ylab("Shannon") +
  ylim(c(3,6)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") 

ggsave(paste0(plot_out,"Suplementary_Figure_4D_Alpha_Source_Shannon.pdf"), width = 2.5, height = 4)

## Perform Statistical Tests

# Richness/Observed Diversity
aov_source_obs <- aov(observed ~ Location, data = otumat_L6_filtered_rare_alpha_source_plot)
summary(aov_source_obs)
tukey_source_obs <- TukeyHSD(aov_source_obs)
tukey_source_obs

# Shannon Diversity
aov_source_shan <- aov(diversity_shannon ~ Location, data = otumat_L6_filtered_rare_alpha_source_plot)
summary(aov_source_shan)
tukey_source_shan <- TukeyHSD(aov_source_shan)
tukey_source_shan


### Analysis of Source vs Culture Samples


## Summarize group means

# Richness/Observed Diversity
sourcecult_rich_summary <- rcompanion::groupwiseMean(observed ~ Sample_Type, data = otumat_L6_filtered_rare_alpha_plot, bca = T)

# Shannon Diversity
sourcecult_shan_summary <- rcompanion::groupwiseMean(diversity_shannon ~ Sample_Type, data = otumat_L6_filtered_rare_alpha_plot, bca = T)  


## Plot Diversity Metrics

############################################## Figure_S4_A ###################################################################


# Observed Diversity
ggplot() +
  geom_jitter(data = otumat_L6_filtered_rare_alpha_plot, aes(x = Sample_Type, y = observed), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = sourcecult_rich_summary, aes(x = Sample_Type, ymin = Bca.lower, ymax = Bca.upper), width = 0.1) +
  geom_point(data = sourcecult_rich_summary, aes(x = Sample_Type, y = Mean, fill = Sample_Type), shape = 22, size = 5) +
  xlab(NULL) +
  ylab("Richness") +
  ylim(c(0,600)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")

ggsave(paste0(plot_out,"Suplementary_Figure_4A_Alpha_SourceCult_Richness.pdf"), width = 2.5, height = 4)

############################################## Figure_S4_B ###################################################################


# Shannon Diversity - more influenced by rare species
ggplot() +
  geom_jitter(data = otumat_L6_filtered_rare_alpha_plot, aes(x = Sample_Type, y = diversity_shannon), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = sourcecult_shan_summary, aes(x = Sample_Type, ymin = Bca.lower, ymax = Bca.upper), width = 0.1) +
  geom_point(data = sourcecult_shan_summary, aes(x = Sample_Type, y = Mean, fill = Sample_Type), shape = 22, size = 5) +
  xlab(NULL) +
  ylab("Shannon") +
  ylim(c(0,6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")

ggsave(paste0(plot_out,"Suplementary_Figure_4B_Alpha_SourceCult_Shannon.pdf"), width = 2.5, height = 4)

## Statistical Testing of Diversity Differences Between Source and Culture

# Observed Diversity
t.test(observed ~ Sample_Type, data = otumat_L6_filtered_rare_alpha_plot)

# Shannon Diversity
t.test(diversity_shannon ~ Sample_Type, data = otumat_L6_filtered_rare_alpha_plot)


### Analysis of Culture Samples Individually (Relative to Each Other)

## Subset data to cultures only
otumat_L6_filtered_rare_alpha_culture_plot <- subset(otumat_L6_filtered_rare_alpha_plot, Sample_Type == "Culture")

## Richness: 

# Statistically evaluate impacting features: Set up model to test global significance 
mod_cult_rich <- lm(observed ~ as.factor(Location) + as.factor(Strain) + as.factor(Passage_Rate), data = otumat_L6_filtered_rare_alpha_culture_plot)
summary(mod_cult_rich) # summarize model
car::Anova(mod_cult_rich, type="III") # use type III (no interactions) accounts for each variance individually while controlling for other variables  

library(relaimpo)
# Calculate factor importance 
calc.relimp(mod_cult_rich, diff = T)

# Perform multiple comparisons on significant main effects with emmeans
em_out_rich <- emmeans::emmeans(mod_cult_rich, pairwise ~ Location)
em_out_rich_em <- as.data.frame(em_out_rich$emmeans)
em_out_rich


############################################## Figure_S4_E ###################################################################


# Plot Richness / Observed Diversity
ggplot() +
  geom_jitter(data = otumat_L6_filtered_rare_alpha_culture_plot, aes(x = Location, y = observed), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = em_out_rich_em, aes(x = Location, ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_point(data = em_out_rich_em, aes(x = Location, y = emmean, fill = Location), shape = 22, size = 5) +
  xlab(NULL) +
  ylab("Richness") +
  #ylim(c(0,60)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")


ggsave(paste0(plot_out,"Suplementary_Figure_4E_Alpha_Cult_Richness.pdf"), width = 2.5, height = 4)

## Shannon Diversity: Statistically evaluate impacting features

# Statistically evaluate impacting features: Set up model to test global significance 
mod_cult_shan <- lm(diversity_shannon ~ as.factor(Location) + as.factor(Strain) + as.factor(Passage_Rate), data = otumat_L6_filtered_rare_alpha_culture_plot)
summary(mod_cult_shan) # summarize model
car::Anova(mod_cult_shan, type="III") # use type III (no interactions) accounts for each variance individually while controlling for other variables

# Calculate factor importance 
calc.relimp(mod_cult_shan, diff = T)

# Perform multiple comparisons on significant main effects with emmeans
em_out_shan <- emmeans::emmeans(mod_cult_shan, pairwise ~ Location)
em_out_shan_em <- as.data.frame(em_out_shan$emmeans)
em_out_shan


## Perform Statistical Tests

# Richness/Observed Diversity
aov_source_obs <- aov(observed ~ Location, data = otumat_L6_filtered_rare_alpha_source_plot)
summary(aov_source_obs)
tukey_source_obs <- TukeyHSD(aov_source_obs)
tukey_source_obs

# Shannon Diversity
aov_source_shan <- aov(diversity_shannon ~ Location, data = otumat_L6_filtered_rare_alpha_source_plot)
summary(aov_source_shan)
tukey_source_shan <- TukeyHSD(aov_source_shan)
tukey_source_shan


############################################## Figure_S4_F ###################################################################



# Shannon Diversity
ggplot() +
  geom_jitter(data = otumat_L6_filtered_rare_alpha_culture_plot, aes(x = Location, y = diversity_shannon), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = em_out_shan_em, aes(x = Location, ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_point(data = em_out_shan_em, aes(x = Location, y = emmean, fill = Location), shape = 22, size = 5) +
  xlab(NULL) +
  ylab("Shannon") +
  #ylim(c(0,4)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")

ggsave(paste0(plot_out,"Suplementary_Figure_4F_Alpha_Cult_Shannon.pdf"), width = 2.5, height = 4)




####
### Enrichment of Taxonomic Groups within Culture Samples
####      


### Filter out bad samples + source (based on metagenomic)


## Source sample names
source_drop <- subset(metadata, Sample_Type == "Source")$Sample

## Sample names to keep
sample_keep <- which(colnames(otu_df_L6) %notin% unique(c(sample_drop, source_drop)))

## Create Filtered Count Dataframe
otu_df_L6_filtered_permute <- otu_df_L6[,..sample_keep]


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
library(data.table)
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

## add column of simplified phylum names to order_enrich_stats_plot_sig
order_enrich_stats_plot_sig$Phylum_simple <- forcats::fct_lump_n(order_enrich_stats_plot_sig$Phylum, n = 15, other_level = "Other")

unique(order_enrich_stats_plot_sig$Phylum_simple)

order_enrich_stats_plot_sig <- order_enrich_stats_plot_sig %>%
  mutate(Phylum_simple = as.character(Phylum_simple))  # Convert to character

order_enrich_stats_plot_sig <- order_enrich_stats_plot_sig %>%
  mutate(Phylum_simple = recode(Phylum_simple, 
                               "Patescibacteria" = "CPR",
                               "Proteobacteria" = "Pseudomonadota"))


## Set up colors for plotting
phy_colors <- c("Acidobacteriota"  = "#5EFFA1",
                "Actinobacteriota"  = "#2444FF",
                "Armatimonadota"  = "#7082B2",
                "Bacteroidota"  = "#FDE4A5",
                "Chlamydiota" = "#F00000",
                "Gemmatimonadota" = "#EA15D6",
                "Firmicutes"  = "#687a6a",
                "Omnitrophota" = "#F6CAF2", 
                "CPR"  = "#79C5D1",
                "Pseudomonadota"  = "#B7B180",
                "Planctomycetota"  = "#4DB6FF",
                "Spirochaetota"  = "#F36579",
                "Verrucomicrobiota"  = "#B4FF00",
                "Other" = "grey80")


## generate vector of phyla (not OTHER)
phy_names <- c("Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota", "Chlamydiota",
               "Gemmatimonadota", "Firmicutes", "Omnitrophota", "CPR", "Pseudomonadota",
               "Planctomycetota", "Spirochaetota", "Verrucomicrobiota")



#install.packages("scales")
library(scales)

############################################## Figure_2_B ###################################################################


## Plot Data
ggplot(order_enrich_stats_plot_sig, aes(x = reorder(Order, Z_score), y = Z_score , fill = Phylum_simple)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +  # ✅ Correct function  scale_fill_manual(values = phy_colors) +
  scale_fill_manual(values = phy_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90)) +
  ylab("Z Score") +
  xlab(NULL)

ggsave(paste0(output,"Figure_2_B_Order_Enrich_Sig.pdf"), width = 12, height = 5)

## Export Permutation Statistical Significance Data
write.table(
  data.frame(order_enrich_stats_plot), 
  file = file.path(table_out, "Order_Enrich_Stats.txt"), 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE
)
# 


### For all non-cyanobacterial genomes (PCA Fig S4G)


## Perform Zero Imputation on Data

## identify genomes with >=2 positive values and only keep these, i adjusted to keep more samples, particualryl source
good_genomes <- apply(otumat_L6_filtered_noCyano, 1, function(x) sum(x > 0)) >= 2
length(good_genomes[which(good_genomes == F)]) # how many genomes to be removed
nrow(otumat_L6_filtered_noCyano) - length(good_genomes[which(good_genomes == F)]) # how many genomes will remain

## remove genomes with < 2 positive values in any sample
otumat_L6_filtered_noCyano_zcompzero <- otumat_L6_filtered_noCyano[good_genomes,]

# impute remaining zeros using zCompositions package
otumat_L6_filtered_noCyano_zcompzero <- zCompositions::cmultRepl(t(otumat_L6_filtered_noCyano_zcompzero), output = "p-counts", z.warning = 1) # Zcompositions needs samples in rows and features in columns
otumat_L6_filtered_noCyano_zcompzero <- t(otumat_L6_filtered_noCyano_zcompzero)


## Perform CLR Transformation (Also lite normalization through CLR)
otumat_L6_filtered_noCyano_zcompzero_clr <- clr.tfm(otumat_L6_filtered_noCyano_zcompzero, features = "rows")


## Estimate impacts of factors on beta-diversity

# Calculate distance matrix
otumat_L6_filtered_noCyano_zcompzero_clr_ach <- vegan::vegdist(t(otumat_L6_filtered_noCyano_zcompzero_clr), method = "euc") 

# get rid of shit samples in metadata
metadata_filt <- subset(metadata, Sample %in% colnames(otumat_L6_filtered_noCyano_zcompzero_clr))

# Run Adonis
vegan::adonis2(otumat_L6_filtered_noCyano_zcompzero_clr_ach ~ Sample_Type + Location, by = "margin", data = metadata_filt, permutations = 999)


## Calculate PCA

# run PCA 
otumat_L6_filtered_noCyano_zcompzero_pca <- prcomp(otumat_L6_filtered_noCyano_zcompzero_clr, center = T, scale = T) # center and scale so value magnitude is negated 
plot(otumat_L6_filtered_noCyano_zcompzero_pca) # plot scree 
summary(otumat_L6_filtered_noCyano_zcompzero_pca)

## merge in metadata 
otumat_L6_filtered_noCyano_zcompzero_pca_plot <- data.frame(Sample = rownames(otumat_L6_filtered_noCyano_zcompzero_pca$rotation), otumat_L6_filtered_noCyano_zcompzero_pca$rotation)
otumat_L6_filtered_noCyano_zcompzero_pca_plot <- merge(otumat_L6_filtered_noCyano_zcompzero_pca_plot, metadata_filt, by = "Sample")


########################################################  Figure S4 G ############################################################################


## plot PC1 and PC2 by location
ggplot(otumat_L6_filtered_noCyano_zcompzero_pca_plot, aes(x = PC1, y = PC2, fill = Location, shape = Sample_Type)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(21,22)) +
  xlab("PC1 (90.3 %)") +
  ylab("PC2 (0.9 %)") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom")



ggsave(paste0(plot_out,"Suplementary_Figure_4G_PCA_Source_v_Culture_NoCyano.pdf"), width = 5, height = 4)    



####
####################################### Dataset 2  - genome resolved - Load All The Data ##################################################
####

#######################################################  Figure 2B  and Figure S1, S2, S4, S5A ############################################################# 



### Input Data
abund_long <- read.table(file = "~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/coverM_Full_Stats_Filt.txt", header = T, sep = "\t") 
counts_wide <- read.table(file = "~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/coverM_Counts_Wide_Filt.txt", header = T, sep = "\t")
abund_wide <- read.table(file = "~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/coverM_Coverage_Wide_Filt.txt", header = T, sep = "\t") 
genome_info <- read.table(file = "~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/Final_Genomes_Summary.tsv", header = T, sep = "\t") 
metadata <- read.table(file = "~/Berkeley_Postdoc/UniCom_analyses/Figures/Data/All_Full_Sample_Metadata_wCov.txt", header = T, sep = "\t") 


### Check sample names are same in all metadata
abund_samp <- unique(abund_long$Sample)
metadata_samp <- unique(metadata$Sample)
abund_samp %in% metadata_samp


### Check all genome names are same in abundance and genome_info
abund_genome <- unique(abund_long$Genome)
genome_info_genome <- unique(genome_info$Genome)
sum(abund_genome %in% genome_info_genome) # should be 537; extra genome in abund_genome is "Unmapped"


### Defining Functions For Analysis

## (1) geometric mean
gm.mean <- function(x){
  exp(mean(log(x)))
}

## (3) sample summary statistics - expects samples in columns and featurs in rows
sample.summary_gen <- function(mat_raw, metaD = NULL, std.asv.names = NULL){
  
  # calculate summary statistics
  Sample_Summaries <- data.frame(Sample = colnames(mat_raw),
                                 Total_Cts = colSums(mat_raw),
                                 Mean_Cts = round(apply(mat_raw, 2, function(x) mean(x)), digits = 0),
                                 Mean_Cts_no0 = round(apply(mat_raw, 2, function(x) mean(x[x>0])), digits = 0),
                                 Min_Cts = apply(mat_raw, 2, function(x) min(x)),
                                 Max_Cts = apply(mat_raw, 2, function(x) max(x)),
                                 Gm_Mean = apply(mat_raw, 2, function(x) gm.mean(x + 1)), # adding 1 so that we can deal with zeros 
                                 Gm_Mean_no0 = apply(mat_raw, 2, function(x) gm.mean(x[x>0])),
                                 Has_Cts = colSums(mat_raw > 0),
                                 SD = apply(mat_raw, 2, FUN = sd),
                                 CV = apply(mat_raw, 2, function(x) sd(x)/mean(x)),
                                 Max_Feature = colnames(t(mat_raw))[max.col(t(mat_raw))],
                                 Max_Feature_Frac = apply(mat_raw, 2, function(x) max(x)) / colSums(mat_raw))
  
  # add total counts w/o standards if std asvs provided
  if(is.null(std.asv.names) == F){
    mat_no_stds <- mat_raw[rownames(mat_raw) %notin% std.asv.names,]
    Sample_Summaries$Cts_NoStd <- colSums(mat_no_stds)
    Sample_Summaries$Frac_NoStd <- Sample_Summaries$Cts_NoStd / Sample_Summaries$Total_Cts
  }
  
  # remove row names
  rownames(Sample_Summaries) <- NULL
  
  # merge in metadata based on "Sample"
  if(is.null(metaD) == F){
    Sample_Summaries <- merge(Sample_Summaries, metaD, by = "Sample")
  }
  
  # return output
  return(Sample_Summaries)
  
}

## NOT IN Operator for Arg Parseing
'%notin%' <- Negate('%in%')

## Create vector of samples to remove due to not passing filtering criteria
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
# 1. Fraction of mapped reads to a culture was unusually low - 1 sample 
# 2. Cyanobacterium expected as dominant is not correct or makes up < 90 % of Cyano fraction - 10 Samples


####
### Create Relative Abundance Plots of Phyla for Sample Groups 
####


# ### Prepare Data for Relative Abundance Plotting
#   
    
## Merge All Info For Plotting (Repeated Code)
abund_long_plot <- merge(abund_long, genome_info, by = "Genome") # leaving out unmapped on purpose (for now?)
abund_long_plot <- merge(abund_long_plot, metadata, by = "Sample")

## Filter Out Samples We Do Not Want 
abund_long_plot_filt <- subset(abund_long_plot, Sample %notin% sample_drop)    

## Subset to groupings of Cultures and Source samples
abund_long_plot_filt_culture <- subset(abund_long_plot_filt, Sample_Type == "Culture")
abund_long_plot_filt_source <- subset(abund_long_plot_filt, Sample_Type == "Source")


## Assess frequency of Phylum Level Taxonomy

# for culutres
culture_phy_cov <- data.table(abund_long_plot_filt_culture %>%
                                group_by(Phylum) %>%
                                summarise(obs = n(),
                                          cov = sum(Coverage_Filt)))
culture_phy_cov$cov_frac <- culture_phy_cov$cov / sum(culture_phy_cov$cov)

# for source material
source_phy_cov <- data.table(abund_long_plot_filt_source %>%
                               group_by(Phylum) %>%
                               summarise(obs = n(),
                                         cov = sum(Coverage_Filt)))
source_phy_cov$cov_frac <- source_phy_cov$cov / sum(source_phy_cov$cov)


### Plot Relative Abundance In Source Samples


## Generate vector of colors for phyla based on detected phyla above
tax_cols_culture <- c(p__Acidobacteriota = "#5EFFA1",
                      p__Actinobacteriota = "#2444FF",
                      p__Armatimonadota = "#7082B2",
                      p__Bacteroidota = "#FDE4A5",
                      #p__Bdellovibrionota = "#3e6641", #(not in cultures)
                      p__Chlamydiota = "#F00000",
                      #p__Chloroflexota = "#E8861A", #(not in cultures)
                      #"p__CSP1-3" = "#964B00", #(not in cultures)
                      p__Cyanobacteria = "#35B779FF",
                      #Desulfobacterota = "#F3A712", #(not in cultures)
                      p__Firmicutes = "#687a6a",
                      p__Gemmatimonadota = "#EA15D6",
                      #p__Myxococcota = "#00AFBB", #(not in cultures)
                      #Ochrophyta = "#90AAE1", #(not in cultures)
                      #Omnitrophota = "#90AAE1", #(not in cultures)
                      #p__Patescibacteria = "#79C5D1", #(not in cultures)
                      p__Planctomycetota = "#4DB6FF",
                      p__Proteobacteria = "#B7B180",
                      #Other = "#000000", 
                      p__Spirochaetota = "#F36579",#EAD637
                      p__Verrucomicrobiota = "#B4FF00")




########################################################  Figure S2 B ############################################################################



## Plot relative abundance stacked bar for source samples
ggplot(abund_long_plot_filt_source, aes(x = Sample, y = Coverage_Filt, fill = reorder(Phylum, -Coverage_Filt))) +
  geom_bar(stat = "identity", position = "fill") + 
  #geom_hline(aes(yintercept = 0.50)) +
  scale_fill_manual(values = tax_cols_culture) +
  labs(fill = "Phylum") +
  xlab(NULL) +
  ylab("Fraction of Filtered Coverage") +
  theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))

ggsave(paste0(plot_out,"Suplementary_Figure_2B_Phylum_Relative_Abundance_Source.pdf"), width = 10, height = 10)


## Get Mean Relative Abundances at Phylum Level + Plot different ways 

# make abundace fraction matrix 
abund_wide_mat_genomes <- abund_wide$Genome
abund_wide_mat <- as.matrix(abund_wide[,-1])
abund_wide_mat_norm <- vegan::decostand(abund_wide_mat, method = "total", MARGIN = 2)

# make abundance fraction df
abund_wide_df_norm <- data.frame(Genome = abund_wide_mat_genomes, abund_wide_mat_norm)

# melt abundance fraction df 
abund_long_df_norm <- pivot_longer(abund_wide_df_norm,
                                   cols = starts_with("unicom"),
                                   names_to = "Sample",
                                   values_to = "Fraction")

# add in sample metadata
abund_long_df_norm <- merge(abund_long_df_norm, metadata, by = "Sample")

# filter out samples we dont want

# Filter Out bad samples
abund_long_df_norm <- subset(abund_long_df_norm, Sample %notin% sample_drop)    

# Remove all source samples
abund_long_df_norm <- subset(abund_long_df_norm, Sample_Type == "Culture")

# Add in genome metadata
abund_long_df_norm <- merge(abund_long_df_norm, genome_info, by = "Genome")


########################################################  Figure 3 A ############################################################################


# Plot stacked bars - FIGURE 3A - Mean Community Compositions w/ Cyanobacteria
ggplot(abund_long_df_norm, aes(x = Strain, y = Fraction, fill = reorder(Phylum, -Fraction))) +
  geom_bar(stat = "identity", position = "fill") + 
  #geom_hline(aes(yintercept = 0.50)) +
  scale_fill_manual(values = tax_cols_culture) +
  labs(fill = "Phylum") +
  xlab(NULL) +
  ylab("Fraction of Coverage") +
  facet_wrap( .~ Location, scales = "free_x") +
  theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))

ggsave(paste0(plot_out,"Figure_3A.pdf"), width = 6, height = 4, units = "in", dpi = 300)



########################################################  Figure 3 B ############################################################################

# Remove Cyanobacteria from plots
abund_long_df_norm_no_cyano <- subset(abund_long_df_norm, Phylum != "p__Cyanobacteria")

# Plot stacked bars - FIGURE 3B - Mean Community Compositions w/o Cyanobacteria
ggplot(abund_long_df_norm_no_cyano, aes(x = Strain, y = Fraction, fill = reorder(Phylum, -Fraction))) +
  geom_bar(stat = "identity", position = "fill") + 
  #geom_hline(aes(yintercept = 0.50)) +
  scale_fill_manual(values = tax_cols_culture) +
  labs(fill = "Phylum") +
  xlab(NULL) +
  ylab("Fraction of Coverage") +
  facet_wrap( .~ Location, scales = "free_x") +
  theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))

ggsave(paste0(plot_out,"Figure_3B.pdf"), width = 6, height = 4, units = "in", dpi = 300)





########################################################  Figure S1 ############################################################################
####
### Visualize the cyanobacterial relative abundance 
####


#####Filter for Cyanobacteria and Culture Samples #####
abund_long_plot_cyano <- abund_long_plot_filt %>%
  filter(Phylum == "p__Cyanobacteria", Sample_Type == "Culture")

##### Aggregate Counts at the Genus Level #####
top_cyano_genera <- abund_long_plot_cyano %>%
  group_by(Genus) %>%
  summarise(Total_Abundance = sum(Counts_Filt, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance)) %>%  # Sort from highest to lowest
  top_n(4, Total_Abundance)  # Select the top 4 genera

##### Modify the Dataset: Keep Top 4, Assign the Rest as "Others" #####
abund_long_plot_cyano <- abund_long_plot_cyano %>%
  mutate(Genus = ifelse(Genus %in% top_cyano_genera$Genus, Genus, "Others"))  # Group non-top genera as "Others"

genus_rename <- c(
  "g__Synechococcus" = "3055/7942",
  "g__Trichormus" = "A7120",
  "g__Synechocystis" = "6803",
  "g__Nodosilinea" = "L0902"
  
)

# Apply name changes
abund_long_plot_cyano <- abund_long_plot_cyano %>%
  mutate(Genus = recode(Genus, !!!genus_rename))

#####  Recalculate Abundance by Sample, Strain, and Genus #####
abund_long_plot_cyano_grouped <- abund_long_plot_cyano %>%
  group_by(Sample, Strain, Genus, Map_Frac_Filt) %>%  # Include Strain to avoid faceting errors
  summarise(Coverage_Filt = sum(Coverage_Filt, na.rm = TRUE)) %>%
  ungroup()

##### Create a Stacked Bar Plot #####
ggplot(abund_long_plot_cyano_grouped, aes(x = reorder(Sample, Map_Frac_Filt), y = Coverage_Filt, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bar plot (relative abundance)
  geom_hline(aes(yintercept = 0.10), color = "purple", linetype = 2, size = 1) +
  theme_minimal() +
  theme(#axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    axis.text.x = element_blank(),  # ❌ Remove x-axis text labels
    axis.ticks.x = element_blank(),  # ❌ Remove x-axis ticks# Rotate x-axis labels
    legend.position = "top") +  # Place legend to the top
  facet_wrap(.~ Strain, nrow = 5, scales = "free_x") +  # Separate by Strain
  xlab(NULL) + ylab("Cyanobacterial Relative Abundance") +  # Label axes
  guides(fill = guide_legend(title = "Genus"))

ggsave(paste0(plot_out,"Suplementary_Figure_1_Cyano_Genera_Composition.pdf"), width = 5, height = 8)



####
############ Calculation of Alpha and Beta Diversity and Statistical Comparisons
####


## Create vectors of samples to remove cyanobacterial genomes
genome_drop <- c("unicom2021_Dbay_Mouth_L0902_BW_C_idba_metabat.8",
                 "unicom2021_Eel_Fox_3055_W_B_spades_maxbin.002",
                 "unicom2021_StrawCreek_S_6803_BW_A_idba_maxbin.001",
                 "unicom2021_StrawCreek_S_A7120_W_A_spades_metabat.3"
)

## Create vector of source samples that should be removed sometimes
#source_samples <- subset(metadata, Sample_Type == "Source")$Sample

## Create matrix of count data for cultures only + bad samples removed

# remove unmapped 
count_wide_mat_all <- subset(counts_wide, Genome != "Unmapped")  

# make into a matrix (counts)
genome_names <- count_wide_mat_all$Genome
count_wide_mat_all <- as.matrix(count_wide_mat_all[,-1])
rownames(count_wide_mat_all) <- genome_names

# remove bad samples
count_wide_mat_all <- count_wide_mat_all[,which(colnames(count_wide_mat_all) %notin% sample_drop)]

# remove cyano genomes
count_wide_mat_nocyano <- count_wide_mat_all[!rownames(count_wide_mat_all) %in% genome_drop, ]

## Create metadata set  for bad removed
all_metadata <- subset(metadata, Sample %notin% sample_drop)
###

######
####
### Alpha Diversity Analysis
####

### Assess Summaries of Filtered Samples 

## mat_all_Filtered
otumat_mat_all_filtered_summary <- sample.summary(count_wide_mat_all, features = "rows")

### Rarefy data for alpha diversity metric calculation

## Set random seed
set.seed(123)

## Rarefy to lowest sample depth
otumat_mat_all_filtered_rare <- GUniFrac::Rarefy(t(count_wide_mat_all), depth = min(otumat_mat_all_filtered_summary$Cts)) # this function applies the function "base::sample" and samples without replacement
otumat_mat_all_filtered_rare <- t(otumat_mat_all_filtered_rare[[1]])


### Calculate Alpha Diversity Stastics

## many alpha-stats from microbiome package
otumat_mat_all_filtered_rare_alpha <- microbiome::alpha(x = otumat_mat_all_filtered_rare, index = "all", zeroes = TRUE) # inclusion of zeros seems to have no impact 


### Assess correlative relatioships of diversity metrics
corrplot::corrplot(cor(otumat_mat_all_filtered_rare_alpha, method = "spearman"))


### Merge in metadata
otumat_mat_all_filtered_rare_alpha_plot <- data.frame(Sample = rownames(otumat_mat_all_filtered_rare_alpha), otumat_mat_all_filtered_rare_alpha)
otumat_mat_all_filtered_rare_alpha_plot <- merge(otumat_mat_all_filtered_rare_alpha_plot, metadata, by = "Sample")


### Check numerical distribtuions of data
plot(density(otumat_mat_all_filtered_rare_alpha_plot$observed)) # unimodal with extreme values
#plot(density(otumat_mat_all_filtered_rare_alpha_plot$veg_rare)) # unimodal with extreme values
plot(density(otumat_mat_all_filtered_rare_alpha_plot$chao1)) # unimodal with extreme values
plot(density(otumat_mat_all_filtered_rare_alpha_plot$diversity_gini_simpson)) # bimodal
plot(density(otumat_mat_all_filtered_rare_alpha_plot$diversity_inverse_simpson)) # unimodal with extreme values
plot(density(otumat_mat_all_filtered_rare_alpha_plot$diversity_shannon)) # unimodal with extreme values


### Analysis of Source Samples Individually (Relative to Each Other)

## Subset Data to Source Only
otumat_mat_all_filtered_rare_alpha_source_plot <- subset(otumat_mat_all_filtered_rare_alpha_plot, Sample_Type == "Source")

## Summarize group means

# Richness/Observed Diversity
source_rich_summary <- rcompanion::groupwiseMean(observed ~ Location, data = otumat_mat_all_filtered_rare_alpha_source_plot, bca = T)

# Shannon Diversity
source_shan_summary <- rcompanion::groupwiseMean(diversity_shannon ~ Location, data = otumat_mat_all_filtered_rare_alpha_source_plot, bca = T)

## Perform Statistical Tests

# Richness/Observed Diversity

library(multcompView)

# Run ANOVA
aov_source_obs <- aov(observed ~ Location, data = otumat_mat_all_filtered_rare_alpha_source_plot)
summary(aov_source_obs)

# Tukey's HSD test
tukey_source_obs <- TukeyHSD(aov_source_obs)
print(tukey_source_obs)

# Extract group significance letters
tukey_letters_obs <- multcompLetters4(aov_source_obs, tukey_source_obs)
tukey_letters_obs

# Extract letters from Tukey results
letters_df_obs <- data.frame(
  Sample_Type = names(tukey_letters_obs$Location$Letters),  # Extract group names
  Letter = tukey_letters_obs$Location$Letters  # Extract corresponding significance letters
)

# Print for verification
print(letters_df_obs)



# Shannon Diversity
# Run ANOVA
aov_source_shan <- aov(diversity_shannon ~ Location, data = otumat_mat_all_filtered_rare_alpha_source_plot)
summary(aov_source_shan)
# Tukey's HSD test
tukey_source_shan <- TukeyHSD(aov_source_shan)
tukey_source_shan
# Extract group significance letters
tukey_letters_shan <- multcompLetters4(aov_source_shan, tukey_source_shan)
tukey_letters_shan

# Extract letters from Tukey results
letters_df_shan <- data.frame(
  Sample_Type = names(tukey_letters_shan$Location$Letters),  # Extract group names
  Letter = tukey_letters_shan$Location$Letters  # Extract corresponding significance letters
)
# Print for verification
print(letters_df_shan)



## Plot Diversity Metrics


############################################## Figure_S5_C ###################################################################

# Richness/Observed Diversity
ggplot() +
  geom_jitter(data = otumat_mat_all_filtered_rare_alpha_source_plot, aes(x = Location, y = observed), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = source_rich_summary, aes(x = Location, ymin = Bca.lower, ymax = Bca.upper), width = 0.1) +
  geom_point(data = source_rich_summary, aes(x = Location, y = Mean, fill = Location), shape = 22, size = 5) +
  
  # Add statistical significance letters
  geom_text(data = letters_df_obs, 
            aes(x = Sample_Type, y = max(otumat_mat_all_filtered_rare_alpha_plot$observed) + 0.2, 
                label = Letter), size = 6) + 
  xlab(NULL) +
  ylab("Richness") +
  #ylim(c(0,250)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1), # Increase x-axis text size
        legend.position = "none",
        axis.text.y = element_text(size = 14),  # Increase y-axis text size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.title.y = element_text(size = 16))  # Increase y-axis title size



ggsave(paste0(plot_out,"Suplementary_Figure_5C_Alpha_Source_Richness_MAG.pdf"), width = 3, height = 5)

############################################## Figure_S5_D ###################################################################

# Shannon Diversity
ggplot() +
  geom_jitter(data = otumat_mat_all_filtered_rare_alpha_source_plot, aes(x = Location, y = diversity_shannon), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = source_shan_summary, aes(x = Location, ymin = Bca.lower, ymax = Bca.upper), width = 0.1) +
  geom_point(data = source_shan_summary, aes(x = Location, y = Mean, fill = Location), shape = 22, size = 5) +
  # Add statistical significance letters
  geom_text(data = letters_df_shan, 
            aes(x = Sample_Type, y = max(otumat_mat_all_filtered_rare_alpha_plot$diversity_shannon) + 0.2, 
                label = Letter), size = 6) + 
  xlab(NULL) +
  ylab("Shannon") +
  #ylim(c(3,6)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1), # Increase x-axis text size
        legend.position = "none",
        axis.text.y = element_text(size = 14),  # Increase y-axis text size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.title.y = element_text(size = 16))  # Increase y-axis title size



ggsave(paste0(plot_out,"Suplementary_Figure_5D_Alpha_Source_Shannon_MAG.pdf"), width = 3, height = 5)



### Analysis of Source vs Culture Samples


## Summarize group means

# Richness/Observed Diversity
sourcecult_rich_summary <- rcompanion::groupwiseMean(observed ~ Sample_Type, data = otumat_mat_all_filtered_rare_alpha_plot, bca = T)

# Shannon Diversity
sourcecult_shan_summary <- rcompanion::groupwiseMean(diversity_shannon ~ Sample_Type, data = otumat_mat_all_filtered_rare_alpha_plot, bca = T)  

## Statistical Testing of Diversity Differences Between Source and Culture

# Run Welch's t-test fro observed and shannon diversity
t_test_richness <- t.test(observed ~ Sample_Type, data = otumat_mat_all_filtered_rare_alpha_plot)
t_test_richness
t_test_shannon <- t.test(diversity_shannon ~ Sample_Type, data = otumat_mat_all_filtered_rare_alpha_plot)
t_test_shannon
# Automatically assign letters based on p-value
assign_letters <- function(p_value) {
  if (p_value < 0.05) {
    return(c("a", "b"))  # Different letters = significant difference
  } else {
    return(c("a", "a"))  # Same letters = no significant difference
  }
}

# Generate letters for both richness and Shannon diversity
richness_letters <- assign_letters(t_test_richness$p.value)
shannon_letters <- assign_letters(t_test_shannon$p.value)

# Create dataframes for ggplot
letters_df_richness <- data.frame(
  Sample_Type = unique(otumat_mat_all_filtered_rare_alpha_plot$Sample_Type),
  Letter = richness_letters
)

letters_df_shannon <- data.frame(
  Sample_Type = unique(otumat_mat_all_filtered_rare_alpha_plot$Sample_Type),
  Letter = shannon_letters
)

# Print to verify
print(letters_df_richness)
print(letters_df_shannon)



## Plot Diversity Metrics

############################################## Figure_S5_A ###################################################################


# Observed Diversity
ggplot() +
  geom_jitter(data = otumat_mat_all_filtered_rare_alpha_plot, aes(x = Sample_Type, y = observed), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = sourcecult_rich_summary, aes(x = Sample_Type, ymin = Bca.lower, ymax = Bca.upper), width = 0.1) +
  geom_point(data = sourcecult_rich_summary, aes(x = Sample_Type, y = Mean, fill = Sample_Type), shape = 22, size = 5) +
  # Add automatically generated significance letters
  geom_text(data = letters_df_richness, 
            aes(x = Sample_Type, y = max(otumat_mat_all_filtered_rare_alpha_plot$observed) + 10, 
                label = Letter), size = 6) + 
  xlab(NULL) +
  ylab("Richness") +
  # ylim(c(0,259)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1), # Increase x-axis text size
        legend.position = "none",
        axis.text.y = element_text(size = 14),  # Increase y-axis text size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.title.y = element_text(size = 16))  # Increase y-axis title size



ggsave(paste0(plot_out,"Suplementary_Figure_5A_Alpha_SourceCult_Richness_MAG.pdf"), width = 3, height = 5)

############################################## Figure_S5_B ###################################################################


# Shannon Diversity - more influenced by rare species
ggplot() +
  geom_jitter(data = otumat_mat_all_filtered_rare_alpha_plot, aes(x = Sample_Type, y = diversity_shannon), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = sourcecult_shan_summary, aes(x = Sample_Type, ymin = Bca.lower, ymax = Bca.upper), width = 0.1) +
  geom_point(data = sourcecult_shan_summary, aes(x = Sample_Type, y = Mean, fill = Sample_Type), shape = 22, size = 5) +
  # Add automatically generated significance letters
  geom_text(data = letters_df_shannon, 
            aes(x = Sample_Type, y = max(otumat_mat_all_filtered_rare_alpha_plot$diversity_shannon) + 0.2, 
                label = Letter), size = 6) + 
  xlab(NULL) +
  ylab("Shannon") +
  #ylim(c(0,6)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1), # Increase x-axis text size
        legend.position = "none",
        axis.text.y = element_text(size = 14),  # Increase y-axis text size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.title.y = element_text(size = 16))  # Increase y-axis title size



ggsave(paste0(plot_out,"Suplementary_Figure_5B_Alpha_SourceCult_Shannon_MAG.pdf"), width = 3, height = 5)


### Analysis of Culture Samples Individually (Relative to Each Other)

## Subset data to cultures only
otumat_mat_all_filtered_rare_alpha_culture_plot <- subset(otumat_mat_all_filtered_rare_alpha_plot, Sample_Type == "Culture")

## Richness: 

# Statistically evaluate impacting features: Set up model to test global significance 
mod_cult_rich <- lm(observed ~ as.factor(Location) + as.factor(Strain) + as.factor(Passage_Rate), data = otumat_mat_all_filtered_rare_alpha_culture_plot)
summary(mod_cult_rich) # summarize model
car::Anova(mod_cult_rich, type="III") # use type III (no interactions) accounts for each variance individually while controlling for other variables  

library(relaimpo)
# Calculate factor importance 
calc.relimp(mod_cult_rich, diff = T)

# Perform multiple comparisons on significant main effects with emmeans
em_out_rich <- emmeans::emmeans(mod_cult_rich, pairwise ~ Location)
em_out_rich_em <- as.data.frame(em_out_rich$emmeans)
em_out_rich


# Run ANOVA
aov_cult_obs <- aov(observed ~ Location, data = otumat_mat_all_filtered_rare_alpha_culture_plot)
summary(aov_cult_obs)

# Tukey's HSD test
tukey_cult_obs <- TukeyHSD(aov_cult_obs)
print(tukey_cult_obs)

# Extract group significance letters
tukey_letters_obs_cult <- multcompLetters4(aov_cult_obs, tukey_cult_obs)
tukey_letters_obs_cult

# Extract letters from Tukey results
letters_df_obs_cult <- data.frame(
  Sample_Type = names(tukey_letters_obs_cult$Location$Letters),  # Extract group names
  Letter = tukey_letters_obs_cult$Location$Letters  # Extract corresponding significance letters
)

# Print for verification
print(letters_df_obs_cult)



############################################## Figure_S5_E ###################################################################


# Plot Richness / Observed Diversity
ggplot() +
  geom_jitter(data = otumat_mat_all_filtered_rare_alpha_culture_plot, aes(x = Location, y = observed), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = em_out_rich_em, aes(x = Location, ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_point(data = em_out_rich_em, aes(x = Location, y = emmean, fill = Location), shape = 22, size = 5) +
  # Add automatically generated significance letters
  geom_text(data = letters_df_obs_cult, 
            aes(x = Sample_Type, y = max(otumat_mat_all_filtered_rare_alpha_culture_plot$observed), 
                label = Letter), size = 6) + 
  xlab(NULL) +
  ylab("Richness") +
  #ylim(c(0,60)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1), # Increase x-axis text size
        legend.position = "none",
        axis.text.y = element_text(size = 14),  # Increase y-axis text size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.title.y = element_text(size = 16))  # Increase y-axis title size




ggsave(paste0(plot_out,"Suplementary_Figure_5E_Alpha_Cult_Richness_MAG.pdf"), width = 3, height = 5)

## Shannon Diversity: Statistically evaluate impacting features

# Statistically evaluate impacting features: Set up model to test global significance 
mod_cult_shan <- lm(diversity_shannon ~ as.factor(Location) + as.factor(Strain) + as.factor(Passage_Rate), data = otumat_mat_all_filtered_rare_alpha_culture_plot)
summary(mod_cult_shan) # summarize model
car::Anova(mod_cult_shan, type="III") # use type III (no interactions) accounts for each variance individually while controlling for other variables

# Calculate factor importance 
calc.relimp(mod_cult_shan, diff = T)

# Perform multiple comparisons on significant main effects with emmeans
em_out_shan <- emmeans::emmeans(mod_cult_shan, pairwise ~ Location)
em_out_shan_em <- as.data.frame(em_out_shan$emmeans)
em_out_shan

# Run ANOVA
aov_cult_shan <- aov(diversity_shannon ~ Location, data = otumat_mat_all_filtered_rare_alpha_culture_plot)
summary(aov_cult_shan)

# Tukey's HSD test
tukey_cult_shan <- TukeyHSD(aov_cult_shan)
print(tukey_cult_shan)

# Extract group significance letters
tukey_letters_shan_cult <- multcompLetters4(aov_cult_shan, tukey_cult_shan)
tukey_letters_shan_cult

# Extract letters from Tukey results
letters_df_shan_cult <- data.frame(
  Sample_Type = names(tukey_letters_shan_cult$Location$Letters),  # Extract group names
  Letter = tukey_letters_shan_cult$Location$Letters  # Extract corresponding significance letters
)

# Print for verification
print(letters_df_shan_cult)



############################################## Figure_S5_F ###################################################################



# Shannon Diversity
ggplot() +
  geom_jitter(data = otumat_mat_all_filtered_rare_alpha_culture_plot, aes(x = Location, y = diversity_shannon), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = em_out_shan_em, aes(x = Location, ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_point(data = em_out_shan_em, aes(x = Location, y = emmean, fill = Location), shape = 22, size = 5) +
  # Add automatically generated significance letters
  geom_text(data = letters_df_obs_cult, 
            aes(x = Sample_Type, y = max(otumat_mat_all_filtered_rare_alpha_culture_plot$diversity_shannon), 
                label = Letter), size = 6) + 
  xlab(NULL) +
  ylab("Shannon") +
  #ylim(c(0,4)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1), # Increase x-axis text size
              legend.position = "none",
              axis.text.y = element_text(size = 14),  # Increase y-axis text size
              axis.title.x = element_text(size = 16),  # Increase x-axis title size
              axis.title.y = element_text(size = 16))  # Increase y-axis title size


ggsave(paste0(plot_out,"Suplementary_Figure_5F_Alpha_Cult_Shannon_MAG.pdf"), width = 3, height = 5)



######
##### beta diversity 
####

### Perform Zero Imputation on Data

## identify genomes with >=2 positive values and only keep these, i adjusted to keep more samples, particualryl source
good_genomes <- apply(count_wide_mat_nocyano, 1, function(x) sum(x > 0)) >= 2
length(good_genomes[which(good_genomes == F)]) # how many genomes to be removed
nrow(count_wide_mat_nocyano) - length(good_genomes[which(good_genomes == F)]) # how many genomes will remain

## remove genomes with < 2 positive values in any sample
count_wide_mat_all_zcompzero_count <- count_wide_mat_nocyano[good_genomes,]
# count_wide_mat_all_zcompzero.t <- t(count_wide_mat_all_zcompzero)

# impute remaining zeros using zCompositions package
count_wide_mat_all_zcompzero <- zCompositions::cmultRepl(t(count_wide_mat_all_zcompzero_count), output = "p-counts", z.warning = 1) # Zcompositions needs samples in rows and features in columns
count_wide_mat_all_zcompzero <- t(count_wide_mat_all_zcompzero)


### Perform CLR Transformation (Also lite normalization through CLR)
count_wide_mat_all_zcompzero_clr <- clr.tfm(count_wide_mat_all_zcompzero, features = "rows")


### Estimate impacts of factors on beta-diversity

## Calculate distance matrix
count_wide_mat_all_zcompzero_clr_ach <- vegdist(t(count_wide_mat_all_zcompzero_clr), method = "euc") 

## Create new metadata without shit samples
metadata_meta_filt <- subset(metadata, Sample %notin% sample_drop)

## Run Adonis
adonis2(count_wide_mat_all_zcompzero_clr_ach ~ Sample_Type + Location, data = metadata_meta_filt)


### Visualize beta-diversity using PCA

## run PCA 
count_wide_mat_all_zcompzero_pca <- prcomp(count_wide_mat_all_zcompzero_clr, center = T, scale = T) # center and scale so value magnitude is negated 
plot(count_wide_mat_all_zcompzero_pca) # plot scree 
summary(count_wide_mat_all_zcompzero_pca)

## merge in metadata 
count_wide_mat_all_zcompzero_pca_plot <- data.frame(Sample = rownames(count_wide_mat_all_zcompzero_pca$rotation), count_wide_mat_all_zcompzero_pca$rotation)
count_wide_mat_all_zcompzero_pca_plot <- merge(count_wide_mat_all_zcompzero_pca_plot, all_metadata, by = "Sample")



########################################################  Figure S5 B ############################################################################

## plot PC1 and PC2 by location
ggplot(count_wide_mat_all_zcompzero_pca_plot, aes(x = PC1, y = PC2, fill = Location, shape = Sample_Type)) +
  geom_point(size = 4, alpha =0.8) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(21,22)) +
  scale_x_reverse() +  # Reverse the x-axis
  xlab("PC1 (51.8 %)") +
  ylab("PC2 (6.4 %)") +
  theme_bw()+
  theme(legend.position = "bottom")

# save plot 
ggsave(paste0(plot_out,"Suplementary_Figure_5B_MAG.pdf"), width = 5, height = 4, units = "in", dpi = 300)

#
####
#### create a dataset for cultures only including cultures for UMAP  visualization####


## Create vector of source samples that should be removed sometimes
source_samples <- subset(metadata, Sample_Type == "Source")$Sample

## Create matrix of count data for cultures only + bad samples removed

# remove unmapped 
count_wide_mat_cult <- subset(counts_wide, Genome != "Unmapped")  

# make into a matrix (counts)
genome_names <- count_wide_mat_cult$Genome
count_wide_mat_cult <- as.matrix(count_wide_mat_cult[,-1])
rownames(count_wide_mat_cult) <- genome_names

# remove bad and source samples
count_wide_mat_cult <- count_wide_mat_cult[,which(colnames(count_wide_mat_cult) %notin% sample_drop)]
count_wide_mat_cult <- count_wide_mat_cult[,which(colnames(count_wide_mat_cult) %notin% source_samples)]

## Create metadata set for cultures only + bad removed
cult_metadata <- subset(metadata, Sample %notin% sample_drop)
cult_metadata <- subset(cult_metadata, Sample %notin% source_samples)


### Perform Zero Imputation on Data

## identify genomes with >=2 positive values and only keep these
good_genomes <- apply(count_wide_mat_cult, 1, function(x) sum(x > 0)) >= 2
length(good_genomes[which(good_genomes == F)]) # how many genomes to be removed
nrow(count_wide_mat_cult) - length(good_genomes[which(good_genomes == F)]) # how many genomes will remain

## remove genomes with < 2 positive values in any sample
count_wide_mat_cult_zcompzero <- count_wide_mat_cult[good_genomes,]

# impute remaining zeros using zCompositions package
count_wide_mat_cult_zcompzero <- zCompositions::cmultRepl(t(count_wide_mat_cult_zcompzero),  z.warning =1, output = "p-counts") # Zcompositions needs samples in rows and features in columns
count_wide_mat_cult_zcompzero <- t(count_wide_mat_cult_zcompzero)
nrow(count_wide_mat_cult_zcompzero)  # how many genomes will remain


### Perform CLR Transformation (Also lite normalization through CLR)
count_wide_mat_cult_zcompzero_clr <- clr.tfm(count_wide_mat_cult_zcompzero, features = "rows")



### Visualize beta-diversity using UMAP

## Load UMAP library
library(umap)

## set seed
set.seed(123)

## run umap

# count_wide_mat_cult_zcompzero_clr<-count_wide_mat_all_zcompzero_clr
count_wide_mat_cult_zcompzero_umap <- umap(t(count_wide_mat_cult_zcompzero_clr))


## create umap component df 
count_wide_mat_cult_zcompzero_umap_df <- count_wide_mat_cult_zcompzero_umap$layout
count_wide_mat_cult_zcompzero_umap_df <- data.frame(Sample = rownames(count_wide_mat_cult_zcompzero_umap_df), count_wide_mat_cult_zcompzero_umap_df)

## merge metadata for plotting
count_wide_mat_cult_zcompzero_umap_plot <- merge(count_wide_mat_cult_zcompzero_umap_df, metadata, by = "Sample")

########################################################  Figure S5 C ############################################################################

## plot umap for Location
ggplot(count_wide_mat_cult_zcompzero_umap_plot, aes(x = X1, y = X2, fill = Location)) +
  geom_point(shape = 21, size = 4, alpha = 0.8) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  #geom_text(aes(label = BioRep), check_overlap = F) +
  #stat_ellipse(aes(group = tax_fam, color = tax_fam), linetype = 2) +
  theme_bw() +
  theme(legend.position = "bottom")+
  xlab("UMAP1") +
  ylab("UMAP2")

## save plot
ggsave(paste0(plot_out,"Suplementary_Figure_55_C_UMAP_by_location.pdf"), width = 5, height = 4)

########################################################  Figure S5 D ############################################################################


## plot umap for Strain
ggplot(count_wide_mat_cult_zcompzero_umap_plot, aes(x = X1, y = X2, fill = Strain)) +
  geom_point(shape = 21, size = 4, alpha = 0.8) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  #geom_text(aes(label = BioRep), check_overlap = F) +
  #stat_ellipse(aes(group = tax_fam, color = tax_fam), linetype = 2) +
  theme_bw() +
  theme(legend.position = "bottom")+
  xlab("UMAP1") +
  ylab("UMAP2")

## save plot 
ggsave(paste0(plot_out,"Suplementary_Figure_55_C_UMAP_by_strain.pdf"), width = 5, height = 4)


########################################################  Figure S6 ############################################################################

#choose specific genome set of core microbiome

genome_order_vector <- c(
  
  # Rhizobiales
  "unicom2021_StrawCreek_S_L0902_W_A_idba_concoct_4" = "Rhizobiales",
  "unicom2021_StrawCreek_S_6803_W_A_idba_concoct_22" = "Rhizobiales",
  "unicom2021_StrawCreek_S_3055_W_B_idba_metabat.8" = "Rhizobiales",
  "unicom2021_Eel_Fox_6803_BW_C_idba_concoct_8" = "Rhizobiales",
  "unicom2021_StrawCreek_S_L0902_W_B_idba_maxbin.005" = "Rhizobiales",
  "unicom2021_Eel_Fox_6803_W_C_idba_vamb_133" = "Rhizobiales",
  "unicom2021_Eel_S_6803_BW_A_idba_concoct_8" = "Rhizobiales",
  
  # Rhodobacterales
  "unicom2021_Dbay_Inner_L0902_BW_A_idba_concoct_19_sub" = "Rhodobacterales",
  "unicom2021_Dbay_Inner_6803_BW_B_idba_concoct_2" = "Rhodobacterales",
  "unicom2021_Dbay_Inner_A7120_BW_A_idba_concoct_16" = "Rhodobacterales",
  "unicom2021_Dbay_Inner_A7120_BW_C_idba_concoct_21" = "Rhodobacterales",
  "unicom2021_Eel_S_L0902_BW_C_idba_vamb_855" = "Rhodobacterales",
  
  # Sphingomonadales
  "unicom2021_StrawCreek_N_L0902_W_C_idba_maxbin.006" = "Sphingomonadales",
  "unicom2021_Eel_Fox_7942_BW_B_idba_concoct_21" = "Sphingomonadales",
  "unicom2021_Eel_Fox_A7120_W_A_idba_metabat.4" = "Sphingomonadales",
  "unicom2021_Eel_S_6803_BW_A_idba_maxbin.004" = "Sphingomonadales",
  "unicom2021_Dbay_Inner_A7120_BW_C_idba_maxbin.043" = "Sphingomonadales",
  
  # Acetobacterales
  "unicom2021_Dbay_Mouth_6803_BW_A_idba_concoct_19_sub" = "Acetobacterales",
  "unicom2021_StrawCreek_S_A7120_W_C_idba_metabat.10" = "Acetobacterales",
  
  # Other Orders
  "unicom2021_Dbay_Inner_A7120_BW_C_idba_metabat.14" = "GCA-2729495",
  "unicom2021_StrawCreek_S_A7120_W_C_idba_concoct_21_sub" = "Pseudomonadales",
  "unicom2021_Eel_Fox_7942_W_A_spades_concoct_4" = "Leptospirales",
  "unicom2021_Eel_Fox_3055_W_C_idba_maxbin.003" = "Chitinophagales",
  "unicom2021_Eel_Fox_3055_W_B_spades_concoct_1" = "Cytophagales",
  "unicom2021_Eel_Fox_6803_W_C_idba_maxbin.011" = "NS11-12g"
)

# get only the culture samples
culture_samples <- metadata$Sample[metadata$Sample_Type == "Culture"]

# convert count data to matrix and set row names
genome_names <- counts_wide$Genome
count_wide_mat <- as.matrix(counts_wide[, -1])  # Remove Genome column
rownames(count_wide_mat) <- genome_names

# keep only culture samples
count_wide_mat_culture <- count_wide_mat[, colnames(count_wide_mat) %in% culture_samples]

# remove bad samples & genomes
count_wide_mat_culture <- count_wide_mat_culture[, !colnames(count_wide_mat_culture) %in% sample_drop]

count_wide_mat_culture_nocyano <- count_wide_mat_culture[rownames(count_wide_mat_culture) != "Unmapped", ]

# convert counts to presence/absence (1 = present, 0 = absent)
count_wide_pa <- count_wide_mat_culture_nocyano
count_wide_pa[count_wide_pa > 0] <- 1

# subset only genomes of interest
selected_genomes <- names(genome_order_vector)
count_wide_pa_selected <- count_wide_pa[rownames(count_wide_pa) %in% selected_genomes, ]

# order genomes according to predefined vector
count_wide_pa_selected <- count_wide_pa_selected[selected_genomes, ]

# create annotation dataframe for genomes (rows)
annotation_df <- data.frame(Order = genome_order_vector)
rownames(annotation_df) <- selected_genomes

# define colors for orders
order_colors <- c(
  "Rhizobiales" = "#fccde5",
  "Rhodobacterales" = "#d9d9d9",
  "Sphingomonadales" = "#ba80b5",
  "Acetobacterales" = "#8ccfc3",
  "GCA-2729495" = "#f37e71",
  "Pseudomonadales" = "#b15928",
  "Leptospirales" = "#80b0d2",
  "Chitinophagales" = "#ffffb3",
  "Cytophagales" = "#bdb8d7",
  "NS11-12g" = "#f9b163"
)

library(dplyr)
# create sample annotation dataframe (columns)
sample_annotation_df <- metadata %>%
  dplyr::filter(Sample %in% colnames(count_wide_pa_selected)) %>%
  dplyr::select(Sample, Location) %>%
  tibble::column_to_rownames("Sample")

# Step 11: Define colors for locations
location_colors <- c(
  "Disc_Bay" = "#1f78b4",
  "Eel_River" = "#6a3d9a",
  "Straw_Creek" = "#ff7f00"
)

# create list of annotation colors
annotation_colors <- list(
  Order = order_colors, 
  Location = location_colors
)
library(pheatmap)

# create heatmap with **both row (order) and column (location) annotations**
pheatmap(count_wide_pa_selected,
         cluster_rows = FALSE,  # Keep predefined genome order
         cluster_cols = TRUE,  # Cluster samples
         annotation_row = annotation_df,  # Genome order annotations
         annotation_col = sample_annotation_df,  # Sample location annotations
         annotation_colors = annotation_colors,  # Define color coding for both annotations
         color = colorRampPalette(c("white", "black"))(50),  # White (absent) to black (present)
         fontsize_row = 8,  # Adjust genome font size
         fontsize_col = 8,  # Adjust sample font size
         show_rownames = TRUE,  # Show genome names
         show_colnames = FALSE)  # Show sample names

##export the figure as Suplementary_Figure_6


######################### Import Data core microbiome and phylogenetic distance###################

######################### Data for figures 3C, 3D, S7, S8###################


master_genome_groupings <- read.table(file = "~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Main_Figures/Figure_4/Data/Master_Bacterial_Groupings.txt", header = T, sep = "\t") 
instrain_compare <- read.table(file = "~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Main_Figures/Figure_4/Data/instrain_compare_all_genomeWide_compare.tsv", header = T, sep = "\t")
rp16_tree <- read.tree("~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Main_Figures/Figure_4/Data/unicom_2021_All_Genomes_wOG_iqtree.treefile")

# ### Defining Functions For Analysis
# 

#   
## (3) sample summary statistics - expects samples in columns and featurs in rows
sample.summary <- function(mat_raw, metaD = NULL, std.asv.names = NULL){
  
  # calculate summary statistics
  Sample_Summaries <- data.frame(Sample = colnames(mat_raw),
                                 Total_Cts = colSums(mat_raw),
                                 Mean_Cts = round(apply(mat_raw, 2, function(x) mean(x)), digits = 0),
                                 Mean_Cts_no0 = round(apply(mat_raw, 2, function(x) mean(x[x>0])), digits = 0),
                                 Min_Cts = apply(mat_raw, 2, function(x) min(x)),
                                 Max_Cts = apply(mat_raw, 2, function(x) max(x)),
                                 Gm_Mean = apply(mat_raw, 2, function(x) gm.mean(x + 1)), # adding 1 so that we can deal with zeros 
                                 Gm_Mean_no0 = apply(mat_raw, 2, function(x) gm.mean(x[x>0])),
                                 Has_Cts = colSums(mat_raw > 0),
                                 SD = apply(mat_raw, 2, FUN = sd),
                                 CV = apply(mat_raw, 2, function(x) sd(x)/mean(x)),
                                 Max_Feature = colnames(t(mat_raw))[max.col(t(mat_raw))],
                                 Max_Feature_Frac = apply(mat_raw, 2, function(x) max(x)) / colSums(mat_raw))
  
  # add total counts w/o standards if std asvs provided
  if(is.null(std.asv.names) == F){
    mat_no_stds <- mat_raw[rownames(mat_raw) %notin% std.asv.names,]
    Sample_Summaries$Cts_NoStd <- colSums(mat_no_stds)
    Sample_Summaries$Frac_NoStd <- Sample_Summaries$Cts_NoStd / Sample_Summaries$Total_Cts
  }
  
  # remove row names
  rownames(Sample_Summaries) <- NULL
  
  # merge in metadata based on "Sample"
  if(is.null(metaD) == F){
    Sample_Summaries <- merge(Sample_Summaries, metaD, by = "Sample")
  }
  
  # return output
  return(Sample_Summaries)
  
}

## NOT IN Operator for Arg Parseing
'%notin%' <- Negate('%in%')




####
### Direct Visualization and Assessment of Genome Abundances Across Culture Samples 
####
## Create vector of cyano genomes to filter out
cyano_genomes <- c("unicom2021_StrawCreek_S_6803_BW_A_idba_maxbin.001", #6803
                   "unicom2021_StrawCreek_S_A7120_W_A_spades_metabat.3", #A7120
                   "unicom2021_Dbay_Mouth_L0902_BW_C_idba_metabat.8", #L0902
                   "unicom2021_Eel_Fox_3055_W_B_spades_maxbin.002") #7942 and 3055

## Create vector of source samples that should be removed sometimes
source_samples <- subset(metadata, Sample_Type == "Source")$Sample

## Create matrix of normalized (fractional) abundance data for cultures only + bad removed

# create new variable for wide abundance matrix with only cultures and bad removed
abund_wide_mat_cult <- abund_wide

# make into a matrix (counts)
genome_names <- abund_wide_mat_cult$Genome
abund_wide_mat_cult <- as.matrix(abund_wide_mat_cult[,-1])
rownames(abund_wide_mat_cult) <- genome_names

# remove bad and source samples
abund_wide_mat_cult <- abund_wide_mat_cult[,which(colnames(abund_wide_mat_cult) %notin% sample_drop)]
abund_wide_mat_cult <- abund_wide_mat_cult[,which(colnames(abund_wide_mat_cult) %notin% source_samples)]

# Perform sum total normalization on remaining samples
abund_wide_mat_cult <- vegan::decostand(abund_wide_mat_cult, method = "total", MARGIN = 2) 

## Create metadata set for cultures only + bad removed
cult_metadata <- subset(metadata, Sample %notin% sample_drop)
cult_metadata <- subset(cult_metadata, Sample %notin% source_samples)


### Calculate Presence Absence Summary for Genomes Across Strains and Locations (using relative abundance wide matrix as input)

## remake relative abundance wide df
abund_wide_df_cult <- data.frame(Genome = genome_names,
                                 abund_wide_mat_cult)

## Melt relative abundance matrix for easy summaries
abund_long_df_cult <- pivot_longer(abund_wide_df_cult,
                                   cols = starts_with("unicom"),
                                   names_to = "Sample",
                                   values_to = "Fraction")

## Merge in metadata
abund_long_df_cult <- merge(abund_long_df_cult, cult_metadata[,c(1,10,15,17)])


## Calculate sums of presence absence
genome_presabs_summaries <- data.table(abund_long_df_cult %>%
                                         group_by(Genome) %>%
                                         summarise(total_obs = sum(Fraction > 0), # total number of samples a genome was observed in
                                                   mean_perc_non0 = mean(Fraction[Fraction > 0]), # mean of read Fractions that are > 0 across all observed samples
                                                   uniq_loc = length(unique(Location[Fraction > 0])), # number of unique locations a genome was observed in (any number = a hit)
                                                   uniq_strain = length(unique(Strain[Fraction > 0])), # number of unique species a genome was observed with (any number = a hit)
                                                   DiscBay_obs = sum(Fraction[Location == "Disc_Bay"] > 0), # number of times a genome was observed in DiscBay samples 
                                                   EelRiver_obs = sum(Fraction[Location == "Eel_River"] > 0), # number of times a genome was observed in Eel River samples 
                                                   StrawCreek_obs = sum(Fraction[Location == "Straw_Creek"] > 0), # number of times a genome was observed in Straw_Creek samples 
                                                   X3055_obs = sum(Fraction[Strain == "3055"] > 0), # number of times a genome was observed with a 3055 sample 
                                                   X7942_obs = sum(Fraction[Strain == "7942"] > 0), # number of times a genome was observed with a 7942 sample 
                                                   X6803_obs = sum(Fraction[Strain == "6803"] > 0), # number of times a genome was observed with a 6803 sample 
                                                   A7120_obs = sum(Fraction[Strain == "A7120"] > 0), # number of times a genome was observed with a A7120 sample 
                                                   L0902_obs = sum(Fraction[Strain == "L0902"] > 0), # number of times a genome was observed with a L0902 sample 
                                                   uniq_loc_3 = sum(c(DiscBay_obs >= 3, EelRiver_obs >= 3, StrawCreek_obs >= 3)), # number of unique sites where a site is only counted if a genome was observed in >= 3 samples from that site
                                                   uniq_strain_3 = sum(c(X3055_obs >= 3, X7942_obs >= 3, X6803_obs >= 3, A7120_obs >= 3, L0902_obs >= 3)), # number of unique host strains where a strain is only counted if a genome was observed in >= 3 samples with that strain
                                                   core_mb = ifelse((uniq_loc_3 >= 2 & uniq_strain_3 >= 4), TRUE, FALSE), # select genomes that have reproducible presence with >= 4 strains and have been seen in cultures from >= 2 source inoculum locations...represent core microbiome
                                                   DiscBay_obs_3 = ifelse(sum(Fraction[Location == "Disc_Bay"] > 0) >= 3,TRUE,FALSE), # is genome observed in >=3 samples from this site 
                                                   EelRiver_obs_3 = ifelse(sum(Fraction[Location == "Eel_River"] > 0) >= 3,TRUE,FALSE), # is genome observed in >=3 samples from this site
                                                   StrawCreek_obs_3 = ifelse(sum(Fraction[Location == "Straw_Creek"] > 0) >= 3,TRUE,FALSE), # is genome observed in >=3 samples from this site
                                                   X3055_obs_3 = ifelse(sum(Fraction[Strain == "3055"] > 0) >= 3,TRUE,FALSE), # is genome observed in >=3 samples with this host
                                                   X7942_obs_3 = ifelse(sum(Fraction[Strain == "7942"] > 0) >= 3,TRUE,FALSE), # is genome observed in >=3 samples with this host
                                                   X6803_obs_3 = ifelse(sum(Fraction[Strain == "6803"] > 0) >= 3,TRUE,FALSE), # is genome observed in >=3 samples with this host
                                                   A7120_obs_3 = ifelse(sum(Fraction[Strain == "A7120"] > 0) >= 3,TRUE,FALSE), # is genome observed in >=3 samples with this host
                                                   L0902_obs_3 = ifelse(sum(Fraction[Strain == "L0902"] > 0) >= 3,TRUE,FALSE))) # is genome observed in >=3 samples with this host


## Filter out Host Cyano Genomes and Genomes Never Observed in Cultures
genome_presabs_summaries <- subset(genome_presabs_summaries, Genome %notin% cyano_genomes) # filter out cyano
genome_presabs_summaries_cultonly <- subset(genome_presabs_summaries, total_obs >= 1) # filter out genomes never observed in cultures 

## Add in genome metadata
genome_presabs_summaries <- merge(genome_presabs_summaries, genome_info, by = "Genome")
genome_presabs_summaries_cultonly <- merge(genome_presabs_summaries_cultonly, genome_info, by = "Genome")

## Save Data Table on Presence Absence Summary  
write.table(genome_presabs_summaries, paste0(table_output_dir,"PresAbs_Summaries_All.txt"), quote = F, sep = "\t", row.names = F)
write.table(genome_presabs_summaries_cultonly, paste0(table_output_dir,"PresAbs_Summaries_CultOnly.txt"), quote = F, sep = "\t", row.names = F)

# 
### Create Upset Plots for Strain and Overlaps

## Load upset package ***NOTE: Upset Plots need boolean input columns 
library(ComplexUpset)


########################################################  Figure 3 D ############################################################################

# make upset plot - Only Strains involved in Core and Specific Microbiomes

# create metadata data frame - info on genomes
genome_md <- subset(genome_info, Genome %notin% cyano_genomes)

# create data frame of intersections - info on if genome observed with strain >= 3 times
genome_presabs_strain <- genome_presabs_summaries_cultonly[,c("Genome", "X3055_obs_3", "X7942_obs_3", "X6803_obs_3", "A7120_obs_3", "L0902_obs_3")]
colnames(genome_presabs_strain) <- c("Genome", "3055", "7942", "6803", "A7120", "L0902")

# merge tables + create vector for columns of interest
genome_presabs_strain_upset <- merge(genome_md, genome_presabs_strain, by = "Genome")

# filter to Core and Specific
core_genomes <- subset(master_genome_groupings, Assoc_Class == "Core")$Genome
genome_presabs_strain_upset_core_specific <- subset(genome_presabs_strain_upset, Genome %in% core_genomes)

# create vector for columns of interest 
strains <- c("3055", "7942", "6803", "A7120", "L0902")

# Build Upset - Just Core - Color by Order
ComplexUpset::upset(genome_presabs_strain_upset_core_specific, strains,
                    name='Host-Strain Instersection',
                    width_ratio = 0.2,
                    min_size = 2,
                    #group_by = 'sets',
                    sort_intersections_by=c('degree','cardinality'),
                    sort_intersections = 'descending',
                    sort_sets = 'descending',
                    annotations = list(
                      'Phylum' = (
                        ggplot(mapping = aes(fill = Order)) +
                          geom_bar(position = "stack", color = "black") + #position=fill if we want realative
                          #scale_fill_manual(values = tax_cols_culture) +
                          scale_fill_brewer(palette = "Set3") +
                          ylab("Fractional Number of Orders"))),
                    queries = list(
                      upset_query(
                        intersect=c('6803', 'L0902','A7120','3055'),
                        color='steelblue',
                        fill='steelblue',
                        only_components=c('intersections_matrix', 'Intersection size')),
                      upset_query(
                        intersect=c('6803', 'L0902','A7120','3055','7942'),
                        color='steelblue',
                        fill='steelblue',
                        only_components=c('intersections_matrix', 'Intersection size'))))

# Save Plot
ggsave(paste0(plot_output_dir,"Figure_3D_Core_microbiome.pdf"), width = 6, height = 8)



########################################################  Figure 3 C ############################################################################


### Plot Relative Abundance of Core, Specific, and Aux Microbiomes

## Merge in Class of Genome (e.g. Core, Specific, Aux)

# subset data to just genome name and class
genome_abundance_class <- master_genome_groupings[,c("Genome", "Assoc_Class")]

# merge into abund_long_df_cult
abund_long_df_cult_class <- merge(abund_long_df_cult, genome_abundance_class, by = "Genome")

# remove genomes classified as Source or NA
abund_long_df_cult_class <- subset(abund_long_df_cult_class, Assoc_Class != "Source" & is.na(Assoc_Class) == F)

# sum up relative abundances by Strain + Assoc_Class
abund_long_df_cult_class_plot <- data.table(abund_long_df_cult_class %>% 
                                              group_by(Strain, Assoc_Class, Location) %>%
                                              summarise(Fraction = sum(Fraction)))


## Plot stacked bars - FIGURE 3C - Mean Community Compositions of Core, Specific, Aux Species - Without Cyanobacteria
ggplot(subset(abund_long_df_cult_class_plot, Assoc_Class != "Host_Cyano"), aes(x = Strain, y = Fraction, fill = reorder(Assoc_Class, -Fraction))) +
  geom_bar(stat = "identity", position = "fill", color = "black") + 
  #geom_hline(aes(yintercept = 0.50)) +
  labs(fill = "Assoc_Class") +
  scale_fill_manual(values = c("grey70", "steelblue")) +
  xlab(NULL) +
  ylab("Fraction of Coverage") +
  facet_wrap( .~ Location, scales = "free_x") +
  theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))

ggsave(paste0(plot_output_dir,"Figure_3D_Core_and_auxiliary_fractions.pdf"), width = 6, height = 4, units = "in", dpi = 300)



### Summarize fratctions of Host, Core, Specfic, and Aux

## Load R companion Library
library(rcompanion)

## Calculate Sums of Fractions 
assoc_class_sums <- groupwiseSum(Fraction ~ Assoc_Class + Sample, data = abund_long_df_cult_class)

## Calculate Groupwise Means
assoc_class_means <- groupwiseMean(Sum ~ Assoc_Class, data = assoc_class_sums, boot = T, bca = T)

## Calculate Core Microbiome range and SD
sd(subset(assoc_class_sums, Assoc_Class == "Core")$Sum)
range(subset(assoc_class_sums, Assoc_Class == "Core")$Sum)



########################################################  Figure S7 ############################################################################

####
### Statistical Analysis of InStrain Data 
####


### Load Required Packages
library(betareg)
library(ggallin)


### Generate Location and Strain Columns for inStrain Compare Data

## make new data frame to hold outpits
instrain_compare_anno <- data.frame(instrain_compare,
                                    loc1 = NA,
                                    loc2 = NA,
                                    strain1 = NA,
                                    strain2 = NA)

## run for loop to populate annotations
for (i in 1:nrow(instrain_compare_anno)){
  
  # logicals for loc1
  if(grepl("Dbay", instrain_compare_anno$name1[i]) == T){
    instrain_compare_anno[i,11] <- "Disc_Bay"
  }
  
  if(grepl("Eel", instrain_compare_anno$name1[i]) == T){
    instrain_compare_anno[i,11] <- "Eel_River"
  }
  
  if(grepl("Straw", instrain_compare_anno$name1[i]) == T){
    instrain_compare_anno[i,11] <- "Straw_Creek"
  }
  
  # logicals for loc2
  if(grepl("Dbay", instrain_compare_anno$name2[i]) == T){
    instrain_compare_anno[i,12] <- "Disc_Bay"
  }
  
  if(grepl("Eel", instrain_compare_anno$name2[i]) == T){
    instrain_compare_anno[i,12] <- "Eel_River"
  }
  
  if(grepl("Straw", instrain_compare_anno$name2[i]) == T){
    instrain_compare_anno[i,12] <- "Straw_Creek"
  }
  
  # logicals for strain1
  if(grepl("3055", instrain_compare_anno$name1[i]) == T){
    instrain_compare_anno[i,13] <- "3055"
  }
  
  if(grepl("7942", instrain_compare_anno$name1[i]) == T){
    instrain_compare_anno[i,13] <- "7942"
  }
  
  if(grepl("6803", instrain_compare_anno$name1[i]) == T){
    instrain_compare_anno[i,13] <- "6803"
  }
  
  if(grepl("A7120", instrain_compare_anno$name1[i]) == T){
    instrain_compare_anno[i,13] <- "A7120"
  }
  
  if(grepl("L0902", instrain_compare_anno$name1[i]) == T){
    instrain_compare_anno[i,13] <- "L0902"
  }
  
  # logicals for strain2
  if(grepl("3055", instrain_compare_anno$name2[i]) == T){
    instrain_compare_anno[i,14] <- "3055"
  }
  
  if(grepl("7942", instrain_compare_anno$name2[i]) == T){
    instrain_compare_anno[i,14] <- "7942"
  }
  
  if(grepl("6803", instrain_compare_anno$name2[i]) == T){
    instrain_compare_anno[i,14] <- "6803"
  }
  
  if(grepl("A7120", instrain_compare_anno$name2[i]) == T){
    instrain_compare_anno[i,14] <- "A7120"
  }
  
  if(grepl("L0902", instrain_compare_anno$name2[i]) == T){
    instrain_compare_anno[i,14] <- "L0902"
  }
  
}

## Create columns for same/diff strain and location
instrain_compare_anno$loc_same <- ifelse(instrain_compare_anno$loc1 == instrain_compare_anno$loc2, "YES", "NO")
instrain_compare_anno$strain_same <- ifelse(instrain_compare_anno$strain1 == instrain_compare_anno$strain2, "YES", "NO")

## Create a fixed popANI varaible that will allow beta regression
instrain_compare_anno$popANI_fix <- instrain_compare_anno$popANI - .0000001

## Merge in genome info 
instrain_compare_anno <- merge(instrain_compare_anno, master_genome_groupings, by.x = "genome", by.y = "Genome")


### Summarize InStrain Comparison Data and Generate Summary Statistics 

## Summarize data on comparisons for each genome
instrain_compare_anno_summary <- data.table(instrain_compare_anno %>% 
                                              group_by(genome) %>%
                                              summarise(comparisons = n(),
                                                        same_popANI = sum(popANI >= 0.9999900),
                                                        loc_same_compare = sum(loc_same == "YES"),
                                                        loc_same_same_popANI = sum((popANI[loc_same == "YES"] >= 0.9999900)),
                                                        loc_diff_compare = sum(loc_same == "NO"),
                                                        loc_diff_same_popANI = sum((popANI[loc_same == "NO"] >= 0.9999900)),
                                                        mean_popANI = mean(popANI),
                                                        sd_popANI = sd(popANI),
                                                        mean_perccomp = mean(percent_compared),
                                                        sd_percomp = sd(percent_compared),
                                                        Assoc_Class = unique(Assoc_Class)))

## Merge in genome_metadata
instrain_compare_anno_summary <- merge(instrain_compare_anno_summary, genome_info, by.x = "genome", by.y = "Genome")

## Create filtered version without cyano
instrain_compare_anno_summary_filt <- subset(instrain_compare_anno_summary, genome %notin% cyano_genomes)

write.table(instrain_compare_anno_summary_filt, "~/Berkeley_Postdoc/Unicom_analyses/Figures/figures_all/instrain_summary.txt", quote = F, sep = "\t", row.names = T)

## Get numbers of core and aux species
assoc_class_counts <- plyr::count(instrain_compare_anno_summary_filt$Assoc_Class)

## Get comparison counts of core, specific, aux
assoc_class_comparisons <- aggregate(comparisons ~ Assoc_Class, data = instrain_compare_anno_summary_filt, FUN = sum)

## Calculate fractions of identical comparisons

# fraction of identical strains from same source
sum(instrain_compare_anno_summary_filt$loc_same_same_popANI) / sum(instrain_compare_anno_summary_filt$loc_same_compare)

# fraction of identical strains from different sources
sum(instrain_compare_anno_summary_filt$loc_diff_same_popANI) / sum(instrain_compare_anno_summary_filt$loc_diff_compare)


### Run Beta Regression Modeling

## Filter out cyanobacterial genomes
instrain_compare_anno_filt <- subset(instrain_compare_anno, genome %notin% cyano_genomes)

## Model outcomes with beta regression

# build beta regression model
instrain_betareg_mod <- betareg(popANI_fix ~ loc_same + strain_same, data = instrain_compare_anno_filt)

# check model + diagnostic plots
summary(instrain_betareg_mod)
plot(instrain_betareg_mod)
lmtest::lrtest(instrain_betareg_mod) # test significance of model

# check significance of terms 
emmeans::joint_tests(instrain_betareg_mod) # this method is preferred as emmeans explicitly works with beta regression

# get estimated marginal means for loc_same
instrain_betareg_mod_em <- emmeans::emmeans(instrain_betareg_mod, ~ loc_same + strain_same)
instrain_betareg_mod_em

# generate compact letter display for statistically different values between comparisons
emmeans::contrast(instrain_betareg_mod_em, method = "pairwise")
instrain_betareg_mod_em_cld <- multcomp::cld(instrain_betareg_mod_em, Letters = letters, , reversed = T, )
print(instrain_betareg_mod_em_cld)

## plot marginal means - w/o strain

# put estimated marginal means data into a data.frame
marginal_means_plot <- data.frame(instrain_betareg_mod_em)

# add combined variables for loc and strain
marginal_means_plot$comb_var <- paste0(marginal_means_plot$loc_same," - ",marginal_means_plot$strain_same)
instrain_compare_anno_filt$comb_var <- paste0(instrain_compare_anno_filt$loc_same," - ",instrain_compare_anno_filt$strain_same)

ggplot() +
  geom_jitter(data = instrain_compare_anno_filt, aes(x = comb_var, y = popANI), width = 0.1, height = 0, shape = 21, size = 2.5, alpha = 0.2, fill = "grey50") +
  geom_errorbar(data = marginal_means_plot, aes(x = comb_var, ymin=asymp.LCL, ymax=asymp.UCL), width = 0.1, size = 1) +
  geom_point(data = marginal_means_plot, aes(x = comb_var, y = emmean, fill = strain_same), shape = 22, size = 5, alpha = 0.7) +
  geom_hline(aes(yintercept = 0.99999), linetype = 2, color = "firebrick") +
  scale_fill_brewer(palette = "Set1") +
  xlab("Same Source Inoculum") +
  ylab("popANI") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8)) 

ggsave(paste0(plot_output_dir,"Suplementary_Figure_8.pdf"), width = 5, height = 4, units = "in", dpi = 300)



### Summarize fratctions of Host, Core, Specfic, and Aux

## Load R companion Library
library(rcompanion)

## Calculate Sums of Fractions 
assoc_class_sums <- groupwiseSum(Fraction ~ Assoc_Class + Sample, data = abund_long_df_cult_class)

## Calculate Groupwise Means
assoc_class_means <- groupwiseMean(Sum ~ Assoc_Class, data = assoc_class_sums, boot = T, bca = T)

## Calculate Core Microbiome range and SD
sd(subset(assoc_class_sums, Assoc_Class == "Core")$Sum)
range(subset(assoc_class_sums, Assoc_Class == "Core")$Sum)


########################################################  Figure S7 ############################################################################


####
### Statistical Analysis of Phylogenetic Distance Between Core, Aux, and Source 
####  


### Generate all pairwise distances of phylogenetic tree branch tips
rp16_tree_tipdist <- cophenetic(rp16_tree)


### Generate all pariwise comparisons of genome sets: Core, Aux, and Source 

## Filter all genomes to those in tree and remove cyanobacteria
rp16_tree_genomes_nocyano <- subset(master_genome_groupings, Genome %in% rp16_tree[["tip.label"]] & Phylum != "p__Cyanobacteria")

## Create bacterial groupings 

# core
rp16_tree_genomes_nocyano_core <- subset(rp16_tree_genomes_nocyano, Assoc_Class == "Core")$Genome

# aux
rp16_tree_genomes_nocyano_aux <- subset(rp16_tree_genomes_nocyano, Assoc_Class == "Aux")$Genome

# source
rp16_tree_genomes_nocyano_source <- subset(rp16_tree_genomes_nocyano, Assoc_Class == "Source")$Genome

## Generate data.frames of all pairs for each set

# core 
rp16_tree_core_pairs <- t(combn(rp16_tree_genomes_nocyano_core, 2))
rp16_tree_core_pairs <- data.frame(rp16_tree_core_pairs, Assoc_Class = "Core")

# aux 
rp16_tree_aux_pairs <- t(combn(rp16_tree_genomes_nocyano_aux, 2))
rp16_tree_aux_pairs <- data.frame(rp16_tree_aux_pairs, Assoc_Class = "Aux")

# source  
rp16_tree_source_pairs <- t(combn(rp16_tree_genomes_nocyano_source, 2))
rp16_tree_source_pairs <- data.frame(rp16_tree_source_pairs, Assoc_Class = "Source")

## Combine data into master data.frame and rename columns
rp16_tree_all_pairs <- rbind(rp16_tree_core_pairs, rp16_tree_aux_pairs, rp16_tree_source_pairs)
colnames(rp16_tree_all_pairs)[1:2] <- c("Genome1", "Genome2")


### Lookup All Pairwise Phylogenetc Distances into Master Table

## Add empty column to accept lookup values
rp16_tree_all_pairs$branch_dist <- 0

## Run for loop to add values 
for (i in 1:nrow(rp16_tree_all_pairs)){
  
  # get genomes for iteration
  tmp_genome1 <- rp16_tree_all_pairs[i,1]
  tmp_genome2 <- rp16_tree_all_pairs[i,2]
  
  # get branch length distance
  tmp_dist <- rp16_tree_tipdist[tmp_genome1, tmp_genome2]
  
  # populate output data.frame
  rp16_tree_all_pairs[i,4] <- tmp_dist
  
}

### Calculate Summary Statistics + Statistical Differences

## Summary Stats - Need bootstrapped percentile confidence inervals for non-normal data (BCa not working)
rp16_tree_all_pairs_means <- rcompanion::groupwiseMean(branch_dist ~ Assoc_Class, data = rp16_tree_all_pairs, R = 5000, boot = T, percentile = T)

## Pairwise Wilcox Test - Data does not appear to be normal
pairwise.wilcox.test(rp16_tree_all_pairs$branch_dist, rp16_tree_all_pairs$Assoc_Class, p.adjust.method = "fdr")



ggplot(rp16_tree_all_pairs) +
  #geom_jitter(data = rp16_tree_all_pairs[sample(nrow(rp16_tree_all_pairs), 10000),], aes(x = reorder(Assoc_Class, branch_dist), y = branch_dist) ,width = 0.2, alpha = 0.1) +
  geom_violin(data = rp16_tree_all_pairs, aes(x = reorder(Assoc_Class, branch_dist), y = branch_dist, fill = Assoc_Class), alpha = 0.5) +
  geom_errorbar(data = rp16_tree_all_pairs_means, aes(x = Assoc_Class, ymin = Percentile.lower, ymax = Percentile.upper), width = 0.1,) +
  geom_point(data = rp16_tree_all_pairs_means, aes(x = Assoc_Class, y = Boot.mean), shape = 22, size = 4, fill = "grey") +
  scale_fill_manual(values = c("grey70", "steelblue","#984EA3")) +
  xlab("Genome Group") +
  ylab("Branch Distance") +
  theme_bw()

ggsave(paste0(plot_output_dir,"Suplementary_Figure_7.pdf"), width = 3, height = 4, units = "in", dpi = 300)








