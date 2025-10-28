#####################################################################
###################  Neutrality Test Script ######################### 
################### Miguel Campos January 2025 ######################
#####################################################################

# Load libraries
library(PopGenome)
library(ggplot2)
library(readr)
library(reshape2)

# Load the data from FASTA files for different species
dataHD <- readData("./Bhybridum_D/Aligned", format = "fasta")
dataHS <- readData("./Bhybridum_S/Aligned", format = "fasta")
dataD <- readData("./Bdistachyon/Aligned", format = "fasta")
dataS <- readData("./Bstacei/Aligned", format = "fasta")

# Perform neutrality statistics for each dataset
resultsHD = neutrality.stats(dataHD, detail=TRUE)
resultsHS = neutrality.stats(dataHS, detail=TRUE)
resultsD = neutrality.stats(dataD, detail=TRUE)
resultsS = neutrality.stats(dataS, detail=TRUE)

# Create data frames for the neutrality results
neutrality_resultsHD <- data.frame(gene=resultsHD@region.names, n_segregating_sites=resultsHD@n.segregating.sites, Tajima_D=resultsHD@Tajima.D)
neutrality_resultsHS <- data.frame(gene=resultsHS@region.names, n_segregating_sites=resultsHS@n.segregating.sites, Tajima_D=resultsHS@Tajima.D)
neutrality_resultsD <- data.frame(gene=resultsD@region.names, n_segregating_sites=resultsD@n.segregating.sites, Tajima_D=resultsD@Tajima.D)
neutrality_resultsS <- data.frame(gene=resultsS@region.names, n_segregating_sites=resultsS@n.segregating.sites, Tajima_D=resultsS@Tajima.D)

# View the results
print(neutrality_resultsHD)
print(neutrality_resultsHS)
print(neutrality_resultsD)
print(neutrality_resultsS)

# Read the manually created CSV files with tab-separated values
setwd("./Neutrality")

TajimaD <- read_delim("TajimaD.csv", delim="\t", col_names=TRUE)

# Ensure the data is in the correct format (melt the data)
TajimaD_melted <- melt(TajimaD, id.vars = "Gene", variable.name = "Species", value.name = "TajimaD")

# Create the plot for Tajima D
TajimaD_Plot <- ggplot(TajimaD_melted, aes(x=TajimaD, fill=Species, color=Species)) +
  geom_density(alpha=0.4) +
  scale_fill_manual(values=c("blue", "purple", "pink", "darkred")) +
  scale_color_manual(values=c("blue", "purple", "pink", "darkred")) +
  theme_minimal() +
  labs(title="Tajima D",
       x="Tajima D",
       y="Density",
       fill="Species",
       color="Species") +
  theme(legend.position="top")

# Display the plots
print(TajimaD_Plot)
