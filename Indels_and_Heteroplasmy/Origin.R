######################################################################
###################  Origin assignment Script ########################
################### Miguel Campos January 2025 #######################
######################################################################

# Load libraries
library(Biostrings)
library(dplyr)

setwd("./Origin")

# Load parental sequences
father_seq <- readDNAStringSet("Bdistachyon.fasta")
mother_seq <- readDNAStringSet("Bstacei.fasta")

######################
#### B. hybridum S ###
######################

# Read the table from a CSV file
HybridumS <- read.csv("Heteroplasmy_HybridumS.csv", sep=";")

# Function to extract alleles from sequences at specific positions
extract_allele <- function(seq, position) {
  position <- as.integer(position)  # Ensure the position is an integer
  as.character(subseq(seq, start=position, end=position))
}

# Extract alleles at heteroplasmic positions
father_alleles <- sapply(HybridumS$Minimum, extract_allele, seq=father_seq)
mother_alleles <- sapply(HybridumS$Minimum, extract_allele, seq=mother_seq)

# Function to check allele matches
check_maternity <- function(change, father_allele, mother_allele) {
  minor_allele <- strsplit(change, " -> ")[[1]][2]
  match_father <- minor_allele == father_allele
  match_mother <- minor_allele == mother_allele
  return(data.frame(match_father, match_mother))
}

# Apply the function to each row of the table
results <- mapply(check_maternity, HybridumS$Change, father_alleles, mother_alleles, SIMPLIFY = FALSE)
results_df <- do.call(rbind, results)
results_df <- cbind(HybridumS, results_df)

write.csv(results_df, file='Putative_BhybridumS.csv')

# Read the table from a CSV file
HybS <- read.csv("Putative_BhybridumS.csv", sep=";")

# Group by Sequence.Name and count the matches
results_summary <- HybS %>%
  group_by(Sequence.Name) %>%
  summarize(
    match_father_only = sum(match_father & !match_mother),
    match_mother_only = sum(!match_father & match_mother),
    match_both = sum(match_father & match_mother),
    match_neither = sum(!match_father & !match_mother)
  )

# Display the results
print(results_summary)


######################
#### B. hybridum D ###
######################

# Read the table from a CSV file
HybridumD <- read.csv("Heteroplasmy_HybridumD.csv", sep=";")

# Function to extract alleles from sequences at specific positions
extract_allele <- function(seq, position) {
  position <- as.integer(position)  # Ensure the position is an integer
  as.character(subseq(seq, start=position, end=position))
}

# Extract alleles at heteroplasmic positions
father_alleles <- sapply(HybridumD$Minimum, extract_allele, seq=father_seq)
mother_alleles <- sapply(HybridumD$Minimum, extract_allele, seq=mother_seq)

# Function to check allele matches
check_maternity <- function(change, father_allele, mother_allele) {
  minor_allele <- strsplit(change, " -> ")[[1]][2]
  match_father <- minor_allele == father_allele
  match_mother <- minor_allele == mother_allele
  return(data.frame(match_father, match_mother))
}

# Apply the function to each row of the table
results <- mapply(check_maternity, HybridumD$Change, father_alleles, mother_alleles, SIMPLIFY = FALSE)
results_df <- do.call(rbind, results)
results_df <- cbind(HybridumD, results_df)

write.csv(results_df, file='Putative_BhybridumD.csv')

# Read the table from a CSV file
HybD <- read.csv("Putative_BhybridumD.csv", sep=";")

# Group by Sequence.Name and count the matches
results_summary <- HybD %>%
  group_by(Sequence.Name) %>%
  summarize(
    match_father_only = sum(match_father & !match_mother),
    match_mother_only = sum(!match_father & match_mother),
    match_both = sum(match_father & match_mother),
    match_neither = sum(!match_father & !match_mother)
  )

# Display the results
print(results_summary)
