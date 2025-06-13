library(data.table)
library(dplyr)
library(ggplot2)



##### Load Data #####
  
  #### Load Master Feature Table
  cyano_ko_master <- fread("Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Manuscript_Submissions/ISME_3/New_Sup_Tables/Table_S24.txt")

  #### Make Long Format Master Table
  cyano_ko_master_long <-  reshape2::melt(cyano_ko_master)
  colnames(cyano_ko_master_long)[5] <- "Species"

  
  
##### Get Summary Statistics on FunCats #####

  #### Get Unique FunCat Counts by Species
  cyano_ko_unique_funcat_summary_byspecies <- cyano_ko_master_long %>%
    group_by(Species, FunCat) %>%
    summarise(FunCat_Total = sum(value),
              FunCat_Uniq = sum(value[Unique == TRUE]),
              FunCat_Uniq_Frac = FunCat_Uniq / FunCat_Total)
    
  
  #### Get Unique FunCat Counts and Variation Stats
  cyano_ko_unique_funcat_summary_all <- cyano_ko_unique_funcat_summary_byspecies %>%
    group_by(FunCat) %>%
    summarise(Total = sum(FunCat_Total),
              Uniq = sum(FunCat_Uniq),
              Frac_Uniq = Uniq / Total,
              Mean_Frac_Uniq = mean(FunCat_Uniq_Frac, na.rm = T),
              SD_Frac_Uniq = sd(FunCat_Uniq_Frac, na.rm = T),
              Spec_Obs_Total = sum(FunCat_Total > 1),
              Spec_Obs_Uniq = sum(FunCat_Uniq > 1))
  
  

##### Identify Specific Functional Categories Over or Underrepresented In Unique Set #####

  #### Test Functional Categories by Genome
    
    ### Set Up Testing Vectors 
    funcat <- unique(cyano_ko_master_long$FunCat)
    species <- unique(cyano_ko_master_long$Species)
    
    ### Create Output Dataframe
    funcat_unique_enrichment_byspecies <- data.frame()
    
    ### testing for loop
    for(i in 1:length(species)){
      
      for(j in 1:length(funcat)){
      
        ## fisher test values 
        a <- sum(subset(cyano_ko_master_long, Species == species[i] & FunCat == funcat[j] & Unique == T & value == 1)$value)  # unique genes in category X in genome Y
        b <- sum(subset(cyano_ko_master_long, Species != species[i] & FunCat == funcat[j] & Unique == T & value == 1)$value)  # unique genes in category X not in genome Y
        c <- sum(subset(cyano_ko_master_long, Species == species[i] & FunCat == funcat[j] & Unique == F & value == 1)$value)  # non-unique genes in category X in genome Y
        d <- sum(subset(cyano_ko_master_long, Species != species[i] & FunCat == funcat[j] & Unique == F & value == 1)$value)  # non-unique genes in category X not in genome Y
      
        ## create fisher test matrix
        tmp_mat <- matrix(c(a, b, c, d),
                          nrow = 2, 
                          byrow = TRUE,
                          dimnames = list(
                            "Gene Status" = c("Unique", "Non-Unique"),
                            "Category" = c("In X", "Not in X")))
      
        ## run fisher test
        tmp_fish <- fisher.test(tmp_mat)
      
        ## create results dataframe
        tmp_res <- data.frame(Species = species[i],
                              FunCat = funcat[j],
                              OR = log(unname(tmp_fish$estimate)),
                              p_value = tmp_fish$p.value)
      
        ## merge to output dataframe
        funcat_unique_enrichment_byspecies <- rbind(funcat_unique_enrichment_byspecies, tmp_res)
      
      }
    }
    
    ### Create Final Output
      
      ## FDR correct p-values
      funcat_unique_enrichment_byspecies$fdr <- p.adjust(funcat_unique_enrichment_byspecies$p_value, method = "fdr")
      funcat_unique_enrichment_byspecies$enriched <- ifelse(funcat_unique_enrichment_byspecies$fdr <= 0.05, T, F)
      
      ## merge tables
      cyano_ko_unique_funcat_summary_byspecies <- merge(cyano_ko_unique_funcat_summary_byspecies, funcat_unique_enrichment_byspecies, by = c("Species", "FunCat"))
    
    
  #### Test Functional Categories Overall
    
    ### Set Up Testing Vectors 
    funcat <- unique(cyano_ko_master_long$FunCat)
    
    ### Create Output Dataframe
    funcat_unique_enrichment <- data.frame()
    
    ### testing for loop
    for(i in 1:length(funcat)){
      
      ## fisher test values 
      a <- sum(subset(cyano_ko_master_long, FunCat == funcat[i] & Unique == T & value == 1)$value) # unique genes in category X
      b <- sum(subset(cyano_ko_master_long, FunCat != funcat[i] & Unique == T & value == 1)$value) # unique genes not in category X
      c <- sum(subset(cyano_ko_master_long, FunCat == funcat[i] & Unique == F & value == 1)$value) # non-unique genes in category X
      d <- sum(subset(cyano_ko_master_long, FunCat != funcat[i] & Unique == F & value == 1)$value) # non-unique genes not in category X
      
      ## create fisher test matrix
      tmp_mat <- matrix(c(a, b, c, d),
                        nrow = 2, 
                        byrow = TRUE,
                        dimnames = list(
                          "Gene Status" = c("Unique", "Non-Unique"),
                          "Category" = c("In X", "Not in X")))
      
      ## run fisher test
      tmp_fish <- fisher.test(tmp_mat)
      
      ## create results dataframe
      tmp_res <- data.frame(FunCat = funcat[i],
                            log_OR = log(unname(tmp_fish$estimate)),
                            p_value = tmp_fish$p.value)
      
      ## merge to output dataframe
      funcat_unique_enrichment <- rbind(funcat_unique_enrichment, tmp_res)
      
    }
    
    ### Create Final Output
    
      ## FDR correct p-values
      funcat_unique_enrichment$fdr <- p.adjust(funcat_unique_enrichment$p_value, method = "fdr")
      funcat_unique_enrichment$enriched <- ifelse(funcat_unique_enrichment$fdr <= 0.05, T, F)
      
      ## merge tables
      cyano_ko_unique_funcat_summary_all <- merge(cyano_ko_unique_funcat_summary_all, funcat_unique_enrichment, by = "FunCat")
      
      ## write table
      fwrite(cyano_ko_unique_funcat_summary_all, "~/Dropbox/Diamond_Lab_Files/Manuscripts/2023_Unicom2021_Cyano/Manuscript_Submissions/ISME_3/New_Sup_Tables/Table_S25.txt", quote = F, sep = "\t")
  
  
##### Plot Important Data for Figures ##### 
  
  #### Figure SX_A: How many unique KEGG functions are found per genome
  cyano_unique_by_species <- subset(cyano_ko_master_long, Unique == T & value == 1)
  cyano_unique_by_species <- plyr::count(cyano_unique_by_species$Species)
  colnames(cyano_unique_by_species) <- c("Species", "FunCat_Uniq")
          
  ggplot(cyano_unique_by_species, aes(x = reorder(Species, -FunCat_Uniq), y = FunCat_Uniq, fill = Species)) +
    geom_bar(stat = "identity", color = "black") +
    xlab("Species") +
    ylab("Unique KEGG Functions") +
    scale_fill_brewer(palette = "Set3") +
    theme_bw()
  #ggsave("~/Downloads/FigSX_A.pdf", width = 5, height = 3)
  
  #### Figure SX_B: Species Breakdown of each Functional Category
  ggplot(cyano_ko_unique_funcat_summary_all, aes(x = reorder(FunCat, -Uniq), y = Uniq)) +
    geom_bar(stat = "identity", fill = "steelblue4", color = "black") +
    xlab("Functional Category") +
    ylab("Unique KEGG Functions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #ggsave("~/Downloads/FigSX_B.pdf", width = 7, height = 5)
      
      
  #### Figure SX_C: Enrichment of unique functions by category across species
  ggplot(cyano_ko_unique_funcat_summary_all, aes(x = reorder(FunCat,-log_OR), y = log_OR, fill = enriched)) +
    geom_bar(stat = "identity", color = "black") +
    xlab(NULL) +
    ylab("Log (Odds Ratio)") +
    scale_fill_manual(values = c("grey", "firebrick3")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #ggsave("~/Downloads/FigSX_C.pdf", width = 7, height = 5)
  
  
  #### Figure SX_C: Species Breakdown of each Functional Category
  ggplot(subset(cyano_ko_unique_funcat_summary_byspecies, FunCat %in% c("Phosphonate_Metabolism", "Nitrogen_Metabolism", "Other")), aes(x = reorder(FunCat, -FunCat_Uniq), y = FunCat_Uniq, fill = Species)) +
    geom_bar(stat = "identity", position = "fill", color = "black") +
    xlab("Functional Category") +
    ylab("Category Fraction by Species") +
    scale_fill_brewer(palette = "Set3") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #ggsave("~/Downloads/FigSX_D.pdf", width = 4, height = 4)
  
  
    

    
    
    
    