library(dplyr)
library(networkD3)
library(tidyverse)
library(dplyr)
# # Install if not already installed
# if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
#   install.packages("ComplexUpset")
# }
# if (!requireNamespace("ggplot2", quietly = TRUE)) {
#   install.packages("ggplot2")
# }

# Load libraries
library(ComplexUpset)
library(ggplot2)
library(viridis)


#load the gathered information form literature and this study culutres
metadata_literature <- read.table('~/Berkeley_Postdoc/UniCom_analyses/literature/metadata_all_filtered_clust_core.txt', header=TRUE, sep="\t", stringsAsFactors = FALSE)

#create a directory for saving figures
plot_out <- "~/Berkeley_Postdoc/UniCom_analyses/literature/figures/"


# Creating the binary presence/absence matrix with studies as columns 
binary_matrix <- metadata_literature %>%
  distinct(Study, secondary_cluster, phylum, order) %>%  # Ensure we are working with unique entries for each study-phylum combination
  mutate(presence = 1) %>%  # Indicate presence
  pivot_wider(names_from = Study, values_from = presence, values_fill = list(presence = 0))
str(binary_matrix)

##ONLY OVERLAPS BETWEEN ALL STUDIES
# Filter rows where the sum of columns from the fourth to the last is greater than 1
filtered_data <- binary_matrix %>%
   filter(rowSums(select(., 4:ncol(.)) == 1, na.rm = TRUE) > 1)

# View the filtered metadata_literature
print(filtered_data)
# View the filtered metadata_literature
print(metadata_literature)

# Ensure the dataset is correctly formatted
filtered_data$secondary_cluster <- as.factor(filtered_data$secondary_cluster)
filtered_data$phylum <- as.factor(filtered_data$phylum)
filtered_data$order <- as.factor(filtered_data$order)
# 
unique<-unique(binary_matrix$order)
print (unique)
ord_colors_all <- c("o__Acetobacterales" = "#8CCFC3", "o__Bacillales" = "#687a6a",
                    "o__Chitinophagales"= "#FAF5B5","o__Cytophagales" = "#BDB8D7","o__Rhizobiales" = "#F7CCE0",
                    "o__Sphingomonadales" = "#BA80B5","o__Fimbriimonadales" = "#E3F4B7","o__Caulobacterales" = "#B8E3BF",
                    "o__Flavobacteriales" = "#CEAABD","o__Gemmatimonadales" = "#E59395","o__Burkholderiales" = "#99A6BE",
                    "o__UBA7662" = "#95B1BF","o__Actinomycetales"  = "#C4B294","o__Ferrovibrionales" = "#41B369",
                    "o__M30B66" = "#E6C164","o__CACIAM-22H2"= "#F681D2", "o__Pseudomonadales"="#C9D066", 
                    "o__Xanthomonadales"="#9B6A68",
                    "o__Rhodobacterales" = "#D8D8D7",
                    "o__GCA-2729495" = "#FB8072",
                    "o__NS11-12g" = "#FDB462",
                    "o__Leptospirales" = "#80B1D3"
                    
)

# # Generate the UpSet plot order ovelaps between studies species level
# ComplexUpset::upset(
#   filtered_data,
#   #sets = study,
#   c("Bai_2015", "Berthold_2023", "Bouma-Gregson_2019", "Cornet_2018", "Levy_2018", "Li_2018", "Li_2023", "Perez-Carrascal_2021", "this_study", "Zhao_2023", "Zorz_2019"), # directly using the study columns
#   width_ratio=0.2,
#   min_size=1,
#   sort_intersections_by=c('degree', 'cardinality'),
#   sort_intersections='descending',
#   sort_sets='descending',
#   #sort_sets='input',
#   #sort_sets=FALSE,
#   annotations=list(
#     'Phylum Distribution' = (
#       ggplot() +
#         geom_bar(aes(fill=order), position="stack", color="black", size = 0.25) + #or "fill if i want realtive abundance
#         scale_fill_manual(values = ord_colors_all) +
#         theme_bw() +
#         theme(legend.position = "right",
#               legend.key.size = unit(0.35, "cm")) +
#         ylab("Number of Orders")+
#         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
#     )
#   ))
# # 
# ggsave(paste0(plot_out,"upset_species_overlaps_only.pdf"), width = 8, height = 8)


###show all the data overlaps and uniques

# Generate a fallback palette with unique colors
remaining_orders <- setdiff(unique(binary_matrix$order), names(ord_colors_all))
library(viridis)
fallback_palette <- viridis(length(remaining_orders))


# Combine specified colors with fallback colors
additional_colors <- setNames(fallback_palette, remaining_orders)
all_colors <- c(ord_colors_all, additional_colors)

# Ensure factor levels for 'order' are set in the correct order
binary_matrix$order <- factor(binary_matrix$order, levels = names(all_colors))

# Generate the UpSet plot for order overlaps between studies at the species level
upset_plot <- ComplexUpset::upset(
  binary_matrix,
  intersect = c("Bai_2015", "Berthold_2023", "Bouma-Gregson_2019", "Cornet_2018", "Levy_2018", "Li_2018", "Li_2023", "Perez-Carrascal_2021", "this_study", "Zhao_2023", "Zorz_2019"), # directly using the study columns
  width_ratio = 0.2,
  min_size = 1,
  sort_intersections_by = c('degree', 'cardinality'),
  sort_intersections = 'descending',
  sort_sets = 'descending',
  annotations = list(
    'Order Distribution' = (
      ggplot() +
        geom_bar(aes(fill = order), position = "fill", color = "black", size = 0.25) + # or "fill" if you want relative abundance
        scale_fill_manual(values = ord_colors_all) +
        theme_bw() +
        theme(
          legend.position = "right",
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.box.spacing = unit(0.45, "cm"),
          legend.spacing.x = unit(0.45, "cm"),
          legend.spacing.y = unit(0.45, "cm"),
          legend.key.height = unit(0.1, "cm"),
          legend.key.width = unit(0.2, "cm")
        ) +
        guides(fill = guide_legend(ncol = 3)) +  # Control the number of columns in the legend
        ylab("Number of Orders") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    )
  )
)

print(upset_plot)


ggsave(paste0(plot_out,"Suplementary_Figure_9A_Species_overlap_between_all_studies.pdf"), width = 10, height = 10)

write.table(binary_matrix, "~/Berkeley_Postdoc/UniCom_analyses/literature/figures/binary_matrix_overlaps.txt", quote = F, sep = "\t", row.names = T)



##################################### order overlaps with our study at species level##########################

# Creating the binary presence/absence matrix with studies as columns and phyla as rows
binary_matrix <- metadata_literature %>%
  distinct(Study, secondary_cluster, phylum, order) %>%  # Ensure we are working with unique entries for each study-phylum combination
  mutate(presence = 1) %>%  # Indicate presence
  pivot_wider(names_from = Study, values_from = presence, values_fill = list(presence = 0))
str(binary_matrix)

# # Ensure the dataset is correctly formatted

filtered_data_study <- binary_matrix %>%
  # Filter rows where the sum of columns from the third to the last is greater than 1
  filter(rowSums(.[, 4:ncol(.)] == 1, na.rm = TRUE) > 1)

# View the filtered metadata_literature
print(filtered_data_study)


# View the filtered metadata_literature
print(metadata_literature)

# Ensure the dataset is correctly formatted
filtered_data_study$secondary_cluster <- as.factor(filtered_data_study$secondary_cluster)

filtered_data_study$phylum <- as.factor(filtered_data_study$phylum)
filtered_data_study$order <- as.factor(filtered_data_study$order)


our_study <- subset(filtered_data_study, this_study == 1)
write.table(our_study, "~/Berkeley_Postdoc/UniCom_analyses/literature/figures/overlapping_genomes_with_our_study.txt", quote = F, sep = "\t", row.names = T)


# Calculate the total number of overlapping genomes for all sets
our_study$total_overlap <- rowSums(our_study[, 4:ncol(our_study)] == 1)

# Calculate and print the number of total overlapping genomes
total_overlapping_genomes <- sum(our_study$total_overlap > 0)

columns_of_interest <- c( "Bouma-Gregson_2019",  "Li_2018", 
                         "Li_2023", "Perez-Carrascal_2021", "this_study", "Zhao_2023")

# Initialize a matrix to store the count of overlaps
overlap_matrix <- matrix(0, nrow = length(columns_of_interest), ncol = length(columns_of_interest))
rownames(overlap_matrix) <- columns_of_interest
colnames(overlap_matrix) <- columns_of_interest

# Calculate overlaps for each pair of columns
for (i in 1:length(columns_of_interest)) {
  for (j in 1:length(columns_of_interest)) {
    if (i != j) {
      overlap_matrix[i, j] <- sum(our_study[[columns_of_interest[i]]] == 1 & our_study[[columns_of_interest[j]]] == 1)
    }
  }
}


# Convert the matrix to a dataframe for better readability
overlap_df <- as.data.frame(overlap_matrix)


##manuallly oad coloring for orders
ord_colors <- c("o__Acetobacterales" = "#8CCFC3",
                "o__Chitinophagales"= "#FAF5B5","o__Cytophagales" = "#BDB8D7","o__Rhizobiales" = "#F7CCE0",
                "o__Sphingomonadales" = "#BA80B5","o__Fimbriimonadales" = "#E3F4B7","o__Caulobacterales" = "#B8E3BF",
                "o__Flavobacteriales" = "#CEAABD","o__Gemmatimonadales" = "#E59395","o__Burkholderiales" = "#99A6BE",
                "o__UBA7662" = "#95B1BF","o__Actinomycetales"  = "#C4B294","o__Ferrovibrionales" = "#41B369",
                "o__M30B66" = "#E6C164","o__CACIAM-22H2"= "#F681D2")


#make a vector of study names:
  study_order <- c("Bai_2015", "Berthold_2023", "Bouma-Gregson_2019", "Cornet_2018", "Levy_2018", "Li_2018", "Li_2023", 
                   "Perez-Carrascal_2021", "this_study", "Zhao_2023", "Zorz_2019")



  # Generate the UpSet plot order overlaps with our study at species level
ComplexUpset::upset(
  our_study,
  study_order,
  width_ratio = 0.2,
  min_size = 1,
  sort_intersections_by = c('degree', 'cardinality'),
  #sort_intersections = 'descending', # Ensure intersections are sorted by degree/cardinality
  sort_sets='descending',
  annotations = list(
    'Phylum Distribution' = (
      ggplot() +
        geom_bar(aes(fill = order), position = "stack", color = "black", size = 0.25) +
        scale_fill_manual(values = ord_colors) +
        theme_bw() +
        theme(legend.position = "right",
              legend.key.size = unit(0.35, "cm")) + 
        
        ylab("Fractional Number of Orders")
      
    )
  )
)

ggsave(paste0(plot_out,"Figure_3E_Species_overlap_this_study_vs_oothers.pdf"), width = 8, height = 8)


######studies and bacteria found in them #######

#ORDER level

binary_matrix_study_or <- metadata_literature %>%
  distinct(order, Study) %>%  # Ensure we are working with unique entries
  mutate(presence = 1) %>%  # Indicate presence
  pivot_wider(names_from = order, values_from = presence, values_fill = list(presence = 0))

# Filter columns where 'this_study' row has value 1
binary_matrix_study_or <- binary_matrix_study_or %>%
  dplyr::select(Study, which(binary_matrix_study_or[binary_matrix_study_or$Study == "this_study", -1] == 1) + 1)

# Summarize the metadata_literature
order_sums <- colSums(binary_matrix_study_or[, -1])

# Convert the summarized metadata_literature to a metadata_literature frame
order_summary <- data.frame(
  Order = names(order_sums),
  Count = as.vector(order_sums)
)


# Convert the metadata_literature to long format
data_long <- binary_matrix_study_or %>%
  pivot_longer(cols = -Study, names_to = "Order", values_to = "Presence") %>%
  filter(Presence == 1)


###arrange by the number of studies detected
# Count the number of studies for each order
order_counts <- data_long %>%
  group_by(Order) %>%
  summarise(Count = n())

# Reorder the Order factor based on the count
data_long <- data_long %>%
  mutate(Order = factor(Order, levels = order_counts$Order[order(order_counts$Count, decreasing = TRUE)]))
# 

library(RColorBrewer)
# Create a color palette with more colors
colors <- brewer.pal(n = 9, name = "Set1")
# Extend the palette if needed
colors <- colorRampPalette(colors)(length(unique(data_long$Study)))


# Create the segmented bar chart
ggplot(data_long, aes(x = Order, fill = Study)) +
  geom_bar(position = "stack") +
  theme_minimal() +
  labs(title = "Number of Studies Each Order is Found In",
       x = "Order",
       y = "Count of Studies",
       fill = "Study") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = colors)

ggsave(paste0(plot_out,"Suplementary_Figure_9B_This_study_orders_in_other_studies.pdf"), width = 8, height = 8)


#### GENUS level

binary_matrix_study_gen <- metadata_literature %>%
  distinct(genus, Study) %>%  # Ensure we are working with unique entries
  mutate(presence = 1) %>%  # Indicate presence
  pivot_wider(names_from = genus, values_from = presence, values_fill = list(presence = 0))

# Filter columns where 'this_study' row has value 1
binary_matrix_study_gen <- binary_matrix_study_gen %>%
  dplyr::select(Study, which(binary_matrix_study_gen[binary_matrix_study_gen$Study == "this_study", -1] == 1) + 1)

# Summarize the metadata_literature
gen_sums <- colSums(binary_matrix_study_gen[, -1])

# Convert the summarized metadata_literature to a metadata_literature frame
gen_summary <- data.frame(
  Genus = names(gen_sums),
  Count = as.vector(gen_sums)
)

# Convert the long format
data_long_gen <- binary_matrix_study_gen %>%
  pivot_longer(cols = -Study, names_to = "Genus", values_to = "Presence") %>%
  filter(Presence == 1)


###arrange by the number of studies detected
# Count the number of studies for each order
gen_counts <- data_long_gen %>%
  group_by(Genus) %>%
  summarise(Count = n())

# Reorder the Genus factor based on the count
data_long_gen <- data_long_gen %>%
  mutate(Genus = factor(Genus, levels = gen_counts$Genus[order(gen_counts$Count, decreasing = TRUE)]))

# Create the segmented bar chart
ggplot(data_long_gen, aes(x = Genus, fill = Study)) +
  geom_bar(position = "stack") +
  theme_minimal() +
  labs(title = "Number of Studies Each Genus is Found In",
       x = "Genus",
       y = "Count of Studies",
       fill = "Study") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = colors)

ggsave(paste0(plot_out,"Suplementary_Figure_9C_This_study_genera_in_other_studies.pdf"), width = 15, height = 10)



