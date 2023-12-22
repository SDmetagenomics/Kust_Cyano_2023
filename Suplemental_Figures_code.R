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



##
#################################################### Import  Data 1. ####################################################
####

####singleM analysis data####


###
### Import  Data
####
#######NEED DATA FROM SPENCER########
#####Data used for parts of Figure S1, Figure S4 and Figure S5


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

ggsave(paste0(output,"Plot1_Markers_per_Sample.pdf"), width = 20, height = 6)

## Plot Results - Total Counts
ggplot(marker_summary, aes(x = marker, y = total)) +
  geom_boxplot(outlier.shape = NA, fill = "steelblue3") +
  geom_hline(aes(yintercept = median(marker_summary$total)), color = "firebrick", linetype = 2) +
  geom_jitter(shape = 21, fill = "grey50", alpha = 0.5, width = 0.2, height = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10() +
  ggtitle("Total Counts of Marker per Sample")

ggsave(paste0(output,"Plot2_Marker_Counts_per_Sample.pdf"), width = 20, height = 6)



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




########################################################  Figure S1 A ############################################################################

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
otutab_robust_phy_counts$Phylum_small <- sub("Patescibacteria", "CPR",  otutab_robust_phy_counts$Phylum_small)

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
              "Proteobacteria" = "#B7B180",
              "Verrucomicrobiota" = "#B4FF00",
              "CPR" = "#79C5D1",
              "Other" = "grey")



# Convert 'Phylum_small' to a factor with levels in the desired order
desired_order = c("Planctomycetota", "Verrucomicrobiota", "Actinobacteriota", "Cyanobacteria",
                  "Proteobacteria", "Bacteroidota", "CPR")  # specify all levels in order
otutab_robust_phy_counts$Phylum_small <- factor(otutab_robust_phy_counts$Phylum_small, levels = desired_order)

# rpL6
ggplot(subset(otutab_robust_phy_counts, marker == "S3.39.ribosomal_protein_L6_rplF"), aes(x = variable, y = counts, fill = Phylum_small)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phy_cols) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 9, angle = 60, hjust = 1)) +

  ylab("Relative_Abundance") +
  xlab(NULL)

ggsave(paste0(plot_output_dir,"Plot11_Phylum_Relative_Abundance_Source_singlem.pdf"), width = 10, height = 10)


##
#################################################### Import  Data 2. ####################################################
####

####genome resolved metagenomics data####


###
### Import  Data
####

#####Data used for parts of Figure S1, Figure S5 


###
### Import All Data + Define Analysis Functions
####


### Input Data
abund_long <- read.table(file = "~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Main_Figures/Figure_3/Data/Fig3_H-K_Data/coverM_Full_Stats_Filt.txt", header = T, sep = "\t") 
counts_wide <- read.table(file = "~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Main_Figures/Figure_3/Data/Fig3_H-K_Data/coverM_Counts_Wide_Filt.txt", header = T, sep = "\t")
abund_wide <- read.table(file = "~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Main_Figures/Figure_3/Data/Fig3_H-K_Data/coverM_Coverage_Wide_Filt.txt", header = T, sep = "\t") 
genome_info <- read.table(file = "~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Main_Figures/Figure_3/Data/Fig3_H-K_Data/Final_Genomes_Summary.tsv", header = T, sep = "\t") 
metadata <- read.table(file = "~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Main_Figures/Figure_3/Data/Fig3_H-K_Data/All_Full_Sample_Metadata_wCov.txt", header = T, sep = "\t") 


### Output Directories
plot_output_dir <- "~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Main_Figures/Figure_3/Fig3_H-K_Code/plots/"
table_output_dir <- "~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Main_Figures/Figure_3/Fig3_H-K_Code/stats/"


### Check sample names are same in all metadata
abund_samp <- unique(abund_long$Sample)
metadata_samp <- unique(metadata$Sample)
abund_samp %in% metadata_samp


### Check all genome names are same in abundance and genome_info
abund_genome <- unique(abund_long$Genome)
genome_info_genome <- unique(genome_info$Genome)
sum(abund_genome %in% genome_info_genome) # should be 537; extra genome in abund_genome is "Unmapped"


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



########################################################  Figure S1 B ############################################################################


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

## Generate vector of colors for phyla based on detected phyla above
tax_cols <- c(p__Acidobacteriota = "#5EFFA1",
              p__Actinobacteriota = "#2444FF",
              p__Armatimonadota = "#7082B2",
              p__Bacteroidota = "#FDE4A5",
              p__Chlamydiota = "#F00000",
              p__Cyanobacteria = "#35B779FF",
              p__Firmicutes = "#687a6a",
              p__Gemmatimonadota = "#EA15D6",
              p__Planctomycetota = "#4DB6FF",
              p__Proteobacteria = "#B7B180",
              p__Spirochaetota = "#F36579",
              p__Verrucomicrobiota = "#B4FF00")

## Plot relative abundance stacked bar for source samples
ggplot(abund_long_plot_filt_source, aes(x = Sample, y = Coverage_Filt, fill = reorder(Phylum, -Coverage_Filt))) +
  geom_bar(stat = "identity", position = "fill") + 
  #geom_hline(aes(yintercept = 0.50)) +
  scale_fill_manual(values = tax_cols) +
  labs(fill = "Phylum") +
  xlab(NULL) +
  ylab("Fraction of Filtered Coverage") +
  theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))

ggsave(paste0(plot_output_dir,"Plot14_Phylum_Relative_Abundance_Source.pdf"), width = 10, height = 10)



###
#################################################### Import  Data 3. ####################################################
####

####16S rRNA data

#####Data used for parts of Figure S3



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




####
### Alpha Diversity Analysis
####


### Rarefy data to normalize between samples (needs samples in rows and ASVs in columns)

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


### Perform Linear Mixed Effects Modeling for Time Series + Repeated Measures Data on Lagged Beta Diversity
## Add factorized data

# w/Cyano
asvmat_filt_rare_alpha_plot$Total_Days_Fac <- as.factor(asvmat_filt_rare_alpha_plot$Total_Days)

# w/o Cyano
asvmat_filt_rare_nocyano_alpha_plot$Total_Days_Fac <- as.factor(asvmat_filt_rare_nocyano_alpha_plot$Total_Days)


#####

########################################################  Figure S3 A ############################################################################

#####

## Shannon Diversity - w/Cyano

# run LME with factorized time points
alpha_shan_mod <- lmer(diversity_shannon ~ General_Site + Strain + Passage_Rate + Total_Days_Fac + (1|Tube_No), data = asvmat_filt_rare_alpha_plot)
summary(alpha_shan_mod)
plot(alpha_shan_mod)

# assess term significance and model fit 
anova(alpha_shan_mod)
write.table(data.frame(anova(alpha_shan_mod)), file = "lme_stats/alpha_shan_mod_anova.txt", sep = "\t", quote = F)
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
write.table(data.frame(alpha_shan_em_cld), file = "lme_stats/alpha_shan_em_cld.txt", sep = "\t", quote = F, row.names = F)


## Shannon Diversity - w/Cyano
ggplot() +
  geom_jitter(data = asvmat_filt_rare_alpha_plot, aes(x = Total_Days_Fac, y = diversity_shannon), fill = "grey90",height = 0, width = 0.05, shape = 21, alpha = 0.5, size = 3) +
  geom_errorbar(data = alpha_shan_em_df, aes(x = Total_Days_Fac, ymin = lower.CL, ymax = upper.CL), width = 0.2,) +
  geom_point(data = alpha_shan_em_df, aes(x = Total_Days_Fac, y = emmean), shape = 22, size = 3, fill = "firebrick4") +
  geom_line(data = alpha_shan_em_df, aes(x = as.numeric(Total_Days_Fac), y = emmean), linetype = 2, color = "firebrick4") +
  #geom_smooth() +
  xlab("Time From Inoculation (d)") +
  ylab("Shannon Diversity") +
  theme_bw() 

ggsave(paste0(plot_out,"Plot3_Shannon_wPF.pdf"), width = 5, height = 3, units = "in")  


########################################################  Figure S3 B ############################################################################


## Shannon Diversity - w/o Cyano

# run LME with factorized time points
alpha_shan_nocyano_mod <- lmer(diversity_shannon ~ General_Site + Strain + Passage_Rate + Total_Days_Fac + (1|Tube_No), data = asvmat_filt_rare_nocyano_alpha_plot)
summary(alpha_shan_nocyano_mod)
plot(alpha_shan_nocyano_mod)

# assess term significance and model fit 
anova(alpha_shan_nocyano_mod)
write.table(data.frame(anova(alpha_shan_nocyano_mod)), file = "lme_stats/alpha_shan_nocyano_mod_anova.txt", sep = "\t", quote = F)
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
write.table(data.frame(alpha_shan_nocyano_em_cld), file = "~/Berkeley_Postdoc/Unicom_analyses/Figures/figures_all/alpha_shan_nocyano_em_cld.txt", sep = "\t", quote = F, row.names = F)

## Shannon Diversity - w/o Cyano
ggplot() +
  geom_jitter(data = asvmat_filt_rare_nocyano_alpha_plot, aes(x = Total_Days_Fac, y = diversity_shannon), fill = "grey90",height = 0, width = 0.05, shape = 21, alpha = 0.5, size = 3) +
  geom_errorbar(data = alpha_shan_nocyano_em_df, aes(x = Total_Days_Fac, ymin = lower.CL, ymax = upper.CL), width = 0.2,) +
  geom_point(data = alpha_shan_nocyano_em_df, aes(x = Total_Days_Fac, y = emmean), shape = 22, size = 3, fill = "firebrick4") +
  geom_line(data = alpha_shan_nocyano_em_df, aes(x = as.numeric(Total_Days_Fac), y = emmean), linetype = 2, color = "firebrick4") +
  #geom_smooth() +
  xlab("Time From Inoculation (d)") +
  ylab("Shannon Diversity") +
  theme_bw() 

ggsave(paste0(plot_out,"Plot4_Shannon_nocyano_wPF.pdf"), width = 5, height = 3, units = "in") 


########################################################  Figure S3 C ############################################################################


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



### Plot Data

## Richness - w/o Cyano
ggplot() +
  geom_jitter(data = asvmat_filt_rare_nocyano_alpha_plot, aes(x = Total_Days_Fac, y = observed), fill = "grey90",height = 0, width = 0.1, shape = 21, alpha = 0.5, size = 3) +
  geom_errorbar(data = alpha_obs_nocyano_em_df, aes(x = Total_Days_Fac, ymin = lower.CL, ymax = upper.CL), width = 0.2,) +
  geom_point(data = alpha_obs_nocyano_em_df, aes(x = Total_Days_Fac, y = emmean), shape = 22, size = 3, fill = "steelblue4") +
  geom_line(data = alpha_obs_nocyano_em_df, aes(x = as.numeric(Total_Days_Fac), y = emmean), linetype = 2, color = "steelblue4") +
  #geom_smooth() +
  xlab("Time From Inoculation (d)") +
  ylab("Richness") +
  theme_bw() 

ggsave(paste0(plot_out,"Plot2_Richness_nocyano_wPF.pdf"), width = 5, height = 3, units = "in")



########################################################  Figure S3 D ############################################################################

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


## w/o Cyano
asvmat_filt_nocyano_no0_alr <- alr.tfm(asvmat_filt_nocyano_no0, features = "rows", stds = pred_stand$ID)


### Calculate Achinson Distance


## w/o Cyano
asvmat_filt_nocyano_no0_alr_dist <- vegan::vegdist(t(asvmat_filt_nocyano_no0_alr), method = "euclidian")



### Generate Table of Lagged Beta-Diversity Comparisons - w/o CyAno

## Create parsable matrix out of distance object
asvmat_filt_nocyano_no0_alr_dist_mat <- as.matrix(asvmat_filt_nocyano_no0_alr_dist)

## Get tube numbers
tube_numbers <- unique(metadata_filt$Tube_No)

## Time points for evaluation
#tp_eval <- c(0,14,28,42,56,70,84)
tp_eval <- c(0,14,28,42,56,70,84,91,105,119)

## Create df to hold output
beta.long.nocyano.out <- data.frame()

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
      tmp_dist <- asvmat_filt_nocyano_no0_alr_dist_mat[tmp_sample1,tmp_sample2]
      
      tmp_out <- data.frame(Sample1 = tmp_sample1,
                            Sample2 = tmp_sample2,
                            Tube_No = tmp_tube_no,
                            Comparison = tmp_compare,
                            Ach_Dist = tmp_dist)
      
      # bind data to output df 
      beta.long.nocyano.out <- rbind(beta.long.nocyano.out, tmp_out)
      
    } else {
      next
    }
  }
}

## Merge in metadata to beta.long.nocyano.out
beta.long.nocyano.out <- merge(beta.long.nocyano.out, good_samples[,-c(7,8)], by = "Tube_No", all.x = T)
beta.long.nocyano.out$Tube_No_Char <- as.character(beta.long.nocyano.out$Tube_No)
beta.long.nocyano.out$Comparison_Fac <- as.factor(beta.long.nocyano.out$Comparison)


### Perform Linear Mixed Effects Modeling for Time Series + Repeated Measures Data on Lagged Beta Diversity


## w/o Cyano 

# run LME with factorized time points
beta_nocyano_mod <- lmer(Ach_Dist ~ General_Site + Strain + Passage_Rate + Comparison_Fac + (1|Tube_No), data = beta.long.nocyano.out)
summary(beta_nocyano_mod)
plot(beta_nocyano_mod)

# assess term significance and model fit 
anova(beta_nocyano_mod)
write.table(data.frame(anova(beta_nocyano_mod)), file = "lme_stats/beta_nocyano_mod_anova.txt", sep = "\t", quote = F)
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
# write.table(data.frame(beta_nocyano_mod_em_cld), file = "~/Berkeley_Postdoc/Unicom_analyses/Figures/figures_all/beta_nocyano_mod_em_cld.txt", sep = "\t", quote = F, row.names = F)



## w/o Cayno

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

ggsave(paste0(plot_out,"Plot6_deltaBeta_nocyano_wPF.pdf"), width = 5, height = 3, units = "in")   


### Perform PERMANOVA analysis on beta-diversity as a function of metadata variables 


## w/o Cyano
adonis_nocyano_out <- vegan::adonis2(asvmat_filt_nocyano_no0_alr_dist ~ General_Site + Strain + Passage_Rate + as.factor(Total_Days), data = metadata_filt, by = "margin", permutations = 999, parallel = 8)
adonis_nocyano_out
write.table(data.frame(adonis_nocyano_out), file = "lme_stats/adonis_nocyano_out.txt", sep = "\t", quote = F)



######

########################################################  Figure S3 E ############################################################################

#####


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
  
  ## Phylum Colors 
  
  
  
  ## w/ Cyano
      
    # Order by General Site
      
      # summarize data
      area_ord_site <- data.frame(asvtab_filt_area_long %>%
                                    group_by(General_Site, Ord_Small, Total_Days) %>%
                                    summarize(Total = sum(Counts)))
      #area plot 
      ggplot(area_ord_site, aes(x = Total_Days, y = Total, group = reorder(Ord_Small, -Total))) +
        geom_area(aes(fill = Ord_Small), stat = "identity", position = "fill", linewidth = 0.2, colour = "black") +
        #geom_vline(aes(xintercept = 84)) +
        facet_wrap(.~General_Site,nrow = 3) +
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

      ggsave(paste0(plot_out,"Plot8_Area_Order_Site.pdf"), width = 7, height = 10, units = "in")


      
      
      
      
      
#####
      ########################################################  Figure S3 F, G ############################################################################
      
#####
      ### Load Library
      library(Maaslin2)
      
      
      ### Parse Data Inputs for Maaslin2
      
      ## Subset Metadata of Interest 
      
      # subset metadata
      maaslin2_metadata <- as.data.frame(subset(metadata_filt, Total_Days %in% c(0,14,28,42,56,70,84,91,105,119)))
      
      # create rows named by sample
      rownames(maaslin2_metadata) <- maaslin2_metadata$Sample
      
      # keep only columns we want
      maaslin2_metadata <- maaslin2_metadata[,c("Strain","Tube_No","General_Site", "Passage_Rate", "Total_Days")]
      
      # fix stuff in data
      maaslin2_metadata$General_Site <- as.factor(sub(" ", "_", maaslin2_metadata$General_Site))
      maaslin2_metadata$Strain <- as.factor(maaslin2_metadata$Strain)
      maaslin2_metadata$Passage_Rate <- as.factor(maaslin2_metadata$Passage_Rate)
      maaslin2_metadata$Total_Days_Fac <- factor(maaslin2_metadata$Total_Days, levels = c("84","14","28","42","56","70","0","91","105","119"))
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
      sum(good_asv) # 1682 out of 2129 ASVs
      
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
      
      # 0d --> 84d Comparison
      maaslin2_fit_TD84ref_0d <- subset(maaslin2_fit_TD84ref$results, metadata == "Total_Days_Fac" & value == "0" & qval <= 0.05)
      
      # 84d --> 119d Comparison
      maaslin2_fit_TD84ref_119d <- subset(maaslin2_fit_TD84ref$results, metadata == "Total_Days_Fac" & value == "119" & qval <= 0.05)
      
      ## Add direction of coefficient relative to 84d time point
      
      # 0d --> 84d Comparison
      maaslin2_fit_TD84ref_0d$Coeff_dir <- ifelse(maaslin2_fit_TD84ref_0d$coef < 0, "Increase", "Decrease")
      
      # 84d --> 119d Comparison
      maaslin2_fit_TD84ref_119d$Coeff_dir <- ifelse(maaslin2_fit_TD84ref_119d$coef > 0, "Increase", "Decrease")
      
      ## Merge in taxonomic information
      
      # 0d --> 84d Comparison
      maaslin2_fit_TD84ref_0d <- merge(maaslin2_fit_TD84ref_0d, taxtab, by.x = "feature", by.y = "ASV_ID", all.x = T)
      
      # 84d --> 119d Comparison
      maaslin2_fit_TD84ref_119d <- merge(maaslin2_fit_TD84ref_119d, taxtab, by.x = "feature", by.y = "ASV_ID", all.x = T)
      
      ## Add in column for streamlined Order taxonomy (Top10)
      
      # 0d --> 84d Comparison
      maaslin2_fit_TD84ref_0d$Ord_Small <- forcats::fct_lump_n(maaslin2_fit_TD84ref_0d$Order, 10, other_level = "Other", ties.method = "first")
      
      # 84d --> 119d Comparison
      maaslin2_fit_TD84ref_119d$Ord_Small <- forcats::fct_lump_n(maaslin2_fit_TD84ref_119d$Order, 10, other_level = "Other", ties.method = "first")
      
      
      ### Plot Results as Volcano Scatter
      
      ## 0d --> 84d Comparison
      
      # Set Colors
      
      # Plot 
      ggplot(maaslin2_fit_TD84ref_0d, aes(x = (-1 * coef), y = -log10(qval), fill = Ord_Small)) +
        geom_point(shape = 21, size = 3) +
        xlim(c(-6,6)) +
        scale_y_log10() +
        scale_fill_brewer(palette = "Set3") +
        xlab("Coefficent (84 d / 0 d)") +
        
        theme(panel.background = element_blank(),
              panel.border = element_rect(fill = NA),
              legend.position = "top",
              strip.background = element_rect(colour = "black") ,
              legend.key.size = unit(0.3, "cm") ,
              axis.text.x = element_text(colour = "black")) 
      
      
      ggsave(paste0(plot_out,"Plot11_maaslin_84dv0d.pdf"), width = 6, height = 4, units = "in")
      
      
      ## 84d --> 119d Comparison
      
      # Set Colors
      
      # Plot
      ggplot(maaslin2_fit_TD84ref_119d, aes(x = coef, y = -log10(qval), fill = Ord_Small)) +
        geom_point(shape = 21, size = 3) +
        xlim(c(-4,4)) +
        scale_y_log10() +
        scale_fill_brewer(palette = "Set3") +
        xlab("Coefficent (119 d / 84 d)") +
        
        theme(panel.background = element_blank(),
              panel.border = element_rect(fill = NA),
              legend.position = "top",
              strip.background = element_rect(colour = "black") ,
              legend.key.size = unit(0.3, "cm") ,
              axis.text.x = element_text(colour = "black")) 
      
      ggsave(paste0(plot_out,"Plot12_maaslin_119dv84d.pdf"), width = 6, height = 4, units = "in")

      
      
      
########################################################  Figure S4  ############################################################################

####################################################  Data 1. ####################################################
      
####singleM analysis data####
  
####
### Alpha Diversity Analysis
####


### Filter out bad samples (based on metagenomic)
    
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


### Calculate a rarefaction curve from the data 
  
  ## Use vegan to calculate and output df (needs row = samp / col = species input)
  rarecurve_dat <- vegan::rarecurve(t(otumat_L6_filtered), step = 50, tidy = T)
  
  ## Merge in metadata with only site names
  rarecurve_dat <- merge(rarecurve_dat, metadata, by.x = "Site", by.y = "Sample")
  
  ## Plot rarefaction curve in ggplot
  ggplot(rarecurve_dat) +
    geom_line(aes(x = Sample, y = Species, group = Site, color = Location), size = 0.5) +
    geom_vline(aes(xintercept = 965), linetype = 2, color = "red") +
    scale_color_brewer(palette = "Set1") +
    xlab("Sample Size") +
    ylab("Unique OTUs") +
    facet_wrap(.~ Sample_Type + Location, scales = "free_y")
  
  ggsave(paste0(output,"Plot10_Rare_Curve_L6_freey.pdf"), width = 20, height = 8)
  
  
  ## Plot rarefaction curve in ggplot
  ggplot(rarecurve_dat) +
    geom_line(aes(x = Sample, y = Species, group = Site, color = Location), size = 0.5) +
    geom_vline(aes(xintercept = 965), linetype = 2, color = "red") +
    scale_color_brewer(palette = "Set1") +
    xlab("Sample Size") +
    ylab("Unique OTUs") +
    facet_wrap(.~ Sample_Type + Location)
  
  ggsave(paste0(output,"Plot11_Rare_Curve_L6_samescale.pdf"), width = 20, height = 8)
  
  
  
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
plot(density(otumat_L6_filtered_rare_alpha_plot$veg_rare)) # unimodal with extreme values
plot(density(otumat_L6_filtered_rare_alpha_plot$chao1)) # unimodal with extreme values
plot(density(otumat_L6_filtered_rare_alpha_plot$diversity_gini_simpson)) # bimodal
plot(density(otumat_L6_filtered_rare_alpha_plot$diversity_inverse_simpson)) # unimodal with extreme values
plot(density(otumat_L6_filtered_rare_alpha_plot$diversity_shannon)) # unimodal with extreme values

########################################################  Figure S4 A, B ############################################################################

### Analysis of Source vs Culture Samples


## Summarize group means

# Richness/Observed Diversity
sourcecult_rich_summary <- rcompanion::groupwiseMean(observed ~ Sample_Type, data = otumat_L6_filtered_rare_alpha_plot, bca = T)

# Shannon Diversity
sourcecult_shan_summary <- rcompanion::groupwiseMean(diversity_shannon ~ Sample_Type, data = otumat_L6_filtered_rare_alpha_plot, bca = T)  


## Plot Diversity Metrics

# Observed Diversity
ggplot() +
  geom_jitter(data = otumat_L6_filtered_rare_alpha_plot, aes(x = Sample_Type, y = observed), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = sourcecult_rich_summary, aes(x = Sample_Type, ymin = Bca.lower, ymax = Bca.upper), width = 0.1) +
  geom_point(data = sourcecult_rich_summary, aes(x = Sample_Type, y = Mean, fill = Sample_Type), shape = 22, size = 5) +
  xlab(NULL) +
  ylab("Richness") +
  ylim(c(0,600)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0(output,"Plot14_Alpha_SourceCult_Richness.pdf"), width = 4, height = 4)

# Shannon Diversity - more influenced by rare species
ggplot() +
  geom_jitter(data = otumat_L6_filtered_rare_alpha_plot, aes(x = Sample_Type, y = diversity_shannon), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
  geom_errorbar(data = sourcecult_shan_summary, aes(x = Sample_Type, ymin = Bca.lower, ymax = Bca.upper), width = 0.1) +
  geom_point(data = sourcecult_shan_summary, aes(x = Sample_Type, y = Mean, fill = Sample_Type), shape = 22, size = 5) +
  xlab(NULL) +
  ylab("Shannon") +
  ylim(c(0,6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0(output,"Plot15_Alpha_SourceCult_Shannon.pdf"), width = 4, height = 4)

## Statistical Testing of Diversity Differences Between Source and Culture

# Observed Diversity
t.test(observed ~ Sample_Type, data = otumat_L6_filtered_rare_alpha_plot)

# Shannon Diversity
t.test(diversity_shannon ~ Sample_Type, data = otumat_L6_filtered_rare_alpha_plot)

########################################################  Figure S4 C, D ############################################################################


### Analysis of Source Samples Individually (Relative to Each Other)

  ## Subset Data to Source Only
  otumat_L6_filtered_rare_alpha_source_plot <- subset(otumat_L6_filtered_rare_alpha_plot, Sample_Type == "Source")
  
  ## Summarize group means
  
    # Richness/Observed Diversity
    source_rich_summary <- rcompanion::groupwiseMean(observed ~ Location, data = otumat_L6_filtered_rare_alpha_source_plot, bca = T)
    
    # Shannon Diversity
    source_shan_summary <- rcompanion::groupwiseMean(diversity_shannon ~ Location, data = otumat_L6_filtered_rare_alpha_source_plot, bca = T)
  
  ## Plot Diversity Metrics
  
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
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    ggsave(paste0(output,"Plot12_Alpha_Source_Richness.pdf"), width = 4, height = 4)
    
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
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    
    ggsave(paste0(output,"Plot13_Alpha_Source_Shannon.pdf"), width = 4, height = 4)
  
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


########################################################  Figure S4 E, F ############################################################################
    

### Analysis of Culture Samples Individually (Relative to Each Other)

  ## Subset data to cultures only
  otumat_L6_filtered_rare_alpha_culture_plot <- subset(otumat_L6_filtered_rare_alpha_plot, Sample_Type == "Culture")

  ## Richness: 
  
    # Statistically evaluate impacting features: Set up model to test global significance 
    mod_cult_rich <- lm(observed ~ as.factor(Location) + as.factor(Strain) + as.factor(Passage_Rate), data = otumat_L6_filtered_rare_alpha_culture_plot)
    summary(mod_cult_rich) # summarize model
    car::Anova(mod_cult_rich, type="III") # use type III (no interactions) accounts for each variance individually while controlling for other variables  
    
    # Calculate factor importance 
    calc.relimp(mod_cult_rich, diff = T)

    # Perform multiple comparisons on significant main effects with emmeans
    em_out_rich <- emmeans::emmeans(mod_cult_rich, pairwise ~ Location)
    em_out_rich_em <- as.data.frame(em_out_rich$emmeans)
    em_out_rich
  
    
    
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
    
    
    # Plot Richness / Observed Diversity
    ggplot() +
      geom_jitter(data = otumat_L6_filtered_rare_alpha_culture_plot, aes(x = Location, y = observed), fill = "grey",height = 0, width = 0.05, shape = 21, alpha = 0.7, size = 3) +
      geom_errorbar(data = em_out_rich_em, aes(x = Location, ymin = lower.CL, ymax = upper.CL), width = 0.1) +
      geom_point(data = em_out_rich_em, aes(x = Location, y = emmean, fill = Location), shape = 22, size = 5) +
      xlab(NULL) +
      ylab("Richness") +
      ylim(c(0,60)) +
      scale_fill_brewer(palette = "Set1") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    ggsave(paste0(output,"Plot16_Alpha_Cult_Richness.pdf"), width = 4, height = 4)
    
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
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    ggsave(paste0(output,"Plot17_Alpha_Cult_Shannon.pdf"), width = 4, height = 4)

     
    
########################################################  Figure S5 A ############################################################################
    
####################################################  Data 1. ####################################################
    
####singleM analysis data####
    
    ### Create Data filtered for Cyanobacteria
    
    otu_df_L6_filtered_noCyano <- subset(otu_df_L6_filtered, Phylum != "Cyanobacteria")
    
    
    ### Create Matricies out of filtered data
    
    ## L6_filtered_noCyano
    otu_names_tmp <- otu_df_L6_filtered_noCyano$OTU_ID
    otumat_L6_filtered_noCyano <- as.matrix(otu_df_L6_filtered_noCyano[,4:124]) #adjust based on numeric columns in otumat
    rownames(otumat_L6_filtered_noCyano) <- otu_names_tmp
    rm(otu_names_tmp)
    
    ### For all non-cyanobacterial genomes (PCA Fig 5a)
    
    
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
    
    ## plot PC1 and PC2 by location
    ggplot(otumat_L6_filtered_noCyano_zcompzero_pca_plot, aes(x = PC1, y = PC2, fill = Location, shape = Sample_Type)) +
      geom_point(size = 4, alpha = 0.8) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") +
      scale_shape_manual(values = c(21,22)) +
      #stat_ellipse(aes(group = Location, color = Location), linetype = 2) +
      xlab("PC1 (90.3 %)") +
      ylab("PC2 (0.9 %)") +
      theme_bw()
    
    ggsave(paste0(output,"Plot19_PCA_Source_v_Culture_NoCyano.pdf"), width = 5, height = 4)    

    

########################################################  Figure S5 B, C, D ############################################################################
    
####################################################  Data 2. ########################################################
    
####genome resolved metagenomics data####

    
    
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
    #do the same thing for the good samples and source data
    #load the count table and metadata for source 

    ########################################################  Figure S5 B ############################################################################
    
    ### Perform Zero Imputation on Data
    
    ## identify genomes with >=2 positive values and only keep these, i adjusted to keep more samples, particualryl source
    good_genomes <- apply(count_wide_mat_nocyano, 1, function(x) sum(x > 0)) >= 2
    length(good_genomes[which(good_genomes == F)]) # how many genomes to be removed
    nrow(count_wide_mat_nocyano) - length(good_genomes[which(good_genomes == F)]) # how many genomes will remain
    
    ## remove genomes with < 2 positive values in any sample
    count_wide_mat_all_zcompzero <- count_wide_mat_nocyano[good_genomes,]
    
    # impute remaining zeros using zCompositions package
    count_wide_mat_all_zcompzero <- zCompositions::cmultRepl(t(count_wide_mat_all_zcompzero), output = "p-counts", z.warning = 1) # Zcompositions needs samples in rows and features in columns
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
    
    ## plot PC1 and PC2 by location

    ggplot(count_wide_mat_all_zcompzero_pca_plot, aes(x = PC1, y = PC2, fill = Location, shape = Sample_Type)) +
      geom_point(size = 4, alpha =0.8) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") +
      scale_shape_manual(values = c(21,22)) +
      scale_x_reverse() +  # Reverse the x-axis
      #stat_ellipse(aes(group = Location, color = Location), linetype = 2) +
      xlab("PC1 (51.8 %)") +
      ylab("PC2 (6.4 %)") +
      theme_bw()
    
    # save plot 
    ggsave(paste0("~/Berkeley_Postdoc/Unicom_analyses/Figures/figures_all/PCA_GENOME_ALL.pdf"), width = 5, height = 4)
    
    
    
    
    
########################################################  Figure S5  C, D ############################################################################
    ####genome resolved metagenomics data####   
    
    ## Create metadata set withouth source samples using already clean data
    metadata_cult <- subset(all_metadata, Sample_Type == "Culture")
    #create vector with only culture samples
    sample_names <- metadata_cult$Sample
    #use the vectro to get only cultures count data
    colnames(count_wide_mat_all) <- as.character(colnames(count_wide_mat_all))
    count_wide_mat_cult <- count_wide_mat_all[, colnames(count_wide_mat_all) %in% sample_names]
    
    ### Perform Zero Imputation on Data
    
    ## identify genomes with >=2 positive values and only keep these, i adjusted to keep more samples, particualryl source
    good_genomes_cult <- apply(count_wide_mat_cult, 1, function(x) sum(x > 0)) >= 2
    length(good_genomes_cult[which(good_genomes_cult == F)]) # how many genomes to be removed
    nrow(count_wide_mat_cult) - length(good_genomes_cult[which(good_genomes_cult == F)]) # how many genomes will remain
    
    ## remove genomes with < 2 positive values in any sample
    count_wide_mat_cult_zcompzero <- count_wide_mat_cult[good_genomes_cult,]
    
    # impute remaining zeros using zCompositions package
    count_wide_mat_cult_zcompzero <- zCompositions::cmultRepl(t(count_wide_mat_cult_zcompzero), output = "p-counts", z.warning = 1) # Zcompositions needs samples in rows and features in columns
    count_wide_mat_cult_zcompzero <- t(count_wide_mat_cult_zcompzero)
    
    
    ### Perform CLR Transformation (Also lite normalization through CLR)
    count_wide_mat_cult_zcompzero_clr <- clr.tfm(count_wide_mat_cult_zcompzero, features = "rows")
    
    
    ### Estimate impacts of factors on beta-diversity
    
    ## Calculate distance matrix
    count_wide_mat_cult_zcompzero_clr_ach <- vegdist(t(count_wide_mat_cult_zcompzero_clr), method = "euc") 
 
    ## Run Adonis
    adonis2(count_wide_mat_cult_zcompzero_clr_ach ~ Strain + Location, data = metadata_cult)

    ### Visualize beta-diversity using PCA
    
    ## run PCA 
    count_wide_mat_cult_zcompzero_pca <- prcomp(count_wide_mat_cult_zcompzero_clr, center = T, scale = T) # center and scale so value magnitude is negated 
    plot(count_wide_mat_cult_zcompzero_pca) # plot scree 
    summary(count_wide_mat_cult_zcompzero_pca)
    
    ## merge in metadata 
    count_wide_mat_cult_zcompzero_pca_plot <- data.frame(Sample = rownames(count_wide_mat_cult_zcompzero_pca$rotation), count_wide_mat_cult_zcompzero_pca$rotation)
    count_wide_mat_cult_zcompzero_pca_plot <- merge(count_wide_mat_cult_zcompzero_pca_plot, metadata_cult, by = "Sample")
    
   
    
    # plot PC1 and PC2 by Strain
    ggplot(count_wide_mat_cult_zcompzero_pca_plot, aes(x = PC1, y = PC2, fill = Strain)) +
      geom_point(shape = 21, size = 4, alpha = 0.8) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") 
      #stat_ellipse(aes(group = Strain, color = Strain), linetype = 2) +
  
    # save plot 
    #ggsave(paste0("~/Berkeley_Postdoc/Unicom_analyses/Figures/figures_all/PCA_genome_cult.pdf"), width = 5, height = 4)
    
    
    # plot PC1 and PC2 by Strain
    ggplot(count_wide_mat_cult_zcompzero_pca_plot, aes(x = PC1, y = PC2, fill = Location)) +
      geom_point(shape = 21, size = 4, alpha = 0.8) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") 
    
      #stat_ellipse(aes(group = Strain, color = Strain), linetype = 2) +
   
    
    
    
    ### Visualize beta-diversity using UMAP
    
    ## Load UMAP library
    library(umap)
    
    ## set seed
    set.seed(1234)
    
    ## run umap
    count_wide_mat_cult_zcompzero_umap <- umap(t(count_wide_mat_cult_zcompzero_clr))
    
    ## create umap component df 
    count_wide_mat_cult_zcompzero_umap_df <- count_wide_mat_cult_zcompzero_umap$layout
    count_wide_mat_cult_zcompzero_umap_df <- data.frame(Sample = rownames(count_wide_mat_cult_zcompzero_umap_df), count_wide_mat_cult_zcompzero_umap_df)
    
    ## merge metadata for plotting
    count_wide_mat_cult_zcompzero_umap_plot <- merge(count_wide_mat_cult_zcompzero_umap_df, all_metadata, by = "Sample")
    
    ## plot umap for Location
    ggplot(count_wide_mat_cult_zcompzero_umap_plot, aes(x = X1, y = X2, fill = Location)) +
      geom_point(shape = 21, size = 4, alpha = 0.8) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") +
      #geom_text(aes(label = BioRep), check_overlap = F) +
      #stat_ellipse(aes(group = tax_fam, color = tax_fam), linetype = 2) +
      theme_bw() +
      xlab("UMAP1") +
      ylab("UMAP2")
    
    
    ## plot umap for Location
    ggplot(count_wide_mat_cult_zcompzero_umap_plot, aes(x = X1, y = X2, fill = Strain)) +
      geom_point(shape = 21, size = 4, alpha = 0.8) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") +
      #geom_text(aes(label = BioRep), check_overlap = F) +
      #stat_ellipse(aes(group = tax_fam, color = tax_fam), linetype = 2) +
      theme_bw() +
      xlab("UMAP1") +
      ylab("UMAP2")
    
    ## save plot
    plot_output_dir <- "~/Berkeley_Postdoc/Unicom_analyses/Figures/Codes/Plots/"
    
    
 