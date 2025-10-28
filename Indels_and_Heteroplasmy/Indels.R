######################################################################
###################  Indels detection Script #########################
################### Miguel Campos January 2025 #######################
######################################################################

# Load libraries
library(Biostrings)
library(DECIPHER)
library(dplyr)

# Load the alignment from the FASTA file
setwd("./Indels")
alignment <- readDNAStringSet("Plastomes.fasta")

# Transform the alignment to a DECIPHER object
alignment <- AlignSeqs(DNAStringSet(alignment))

# Transform the alignment to a matrix and add column names
alignment_matrix <- as.matrix(alignment)
rownames(alignment_matrix) <- names(alignment)

# Function to detect gap islands and determine if they are insertions or deletions
detect_gap_islands2 <- function(alignment_matrix) {
  gap_islands_list <- list()
  ref_seq <- alignment_matrix[1, ]  # Use the first sequence as the reference
  
  for (i in 2:nrow(alignment_matrix)) {
    query_seq <- alignment_matrix[i, ]
    gap_positions <- which(ref_seq != query_seq & (ref_seq == "-" | query_seq == "-"))
    
    if (length(gap_positions) > 0) {
      gaps <- data.frame(Position = gap_positions, Sequence = query_seq[gap_positions])
      gaps$Group <- cumsum(c(1, diff(gaps$Position) != 1))  # Identify continuous groups
      
      gap_islands <- gaps %>%
        group_by(Group) %>%
        summarize(
          Start = min(Position),
          End = max(Position),
          Length = n(),
          Type = ifelse(all(ref_seq[min(Position):max(Position)] == "-"), "Deletion", "Insertion")
        )
      
      gap_islands_list[[rownames(alignment_matrix)[i]]] <- gap_islands
    }
  }
  return(gap_islands_list)
}

# Detect gap islands
gap_islands <- detect_gap_islands2(alignment_matrix)

# Create a data frame with the results
gap_islands_df <- do.call(rbind, lapply(names(gap_islands), function(seq_name) {
  data.frame(Sequence = seq_name, gap_islands[[seq_name]], stringsAsFactors = FALSE)
}))

# View the results
head(gap_islands_df)

# Save the results to a CSV file
write.csv(gap_islands_df, "gap_islands2.csv", row.names = FALSE)

# Function to identify shared gap islands
identify_shared_gaps2 <- function(gap_islands_df) {
  gap_islands_grouped <- gap_islands_df %>%
    group_by(Start, End, Type) %>%
    summarize(
      Species_Count = n_distinct(Sequence),
      Species = paste(unique(Sequence), collapse = ", "),
      .groups = 'drop'
    )
  
  return(gap_islands_grouped)
}

# Identify shared gap islands
shared_gaps <- identify_shared_gaps2(gap_islands_df)

# Classify the gap islands
classified_gaps <- shared_gaps %>%
  mutate(Classification = case_when(
    Species_Count == 1 ~ "Species Specific",
    Species_Count == 2 ~ "Shared by 2 Species",
    Species_Count == 3 ~ "Shared by 3 Species"
  ))

# View the results
head(classified_gaps)

# Guardar los resultados en un archivo CSV
write.csv(classified_gaps, "classified_gap_islands2.csv", row.names = FALSE)
