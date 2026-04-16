# Required packages
require(data.table)
require(Biostrings)
require(tidyverse)

# Read input files
classification <- fread("data/full_ITS/eukaryome_ITS.classification")
fasta <- readDNAStringSet("data/full_ITS/eukaryome_ITS.fasta")
fungal_traits <- fread("data/fungal_traits.txt")

# Get a list of families with ecm genera
ecm_families <- classification %>%
  filter(genus %in% fungal_traits$genus[fungal_traits[["primary_lifestyle"]] == "ectomycorrhizal"]) %>%
  pull(family) %>%
  unique()

# Filter classification to ecm
ecm_classification <- classification %>%
  filter(
    family %in% ecm_families,
    # Filter to identified sequences
    species != "unclassified"
  )

# Check unique families
ecm_classification %>%
  group_by(family) %>%
  summarise(n = n()) %>%
  print(n = Inf)

# Check unique genera
ecm_classification %>%
  group_by(genus) %>%
  summarise(n = n()) %>%
  print(n = Inf)

# Filter the fasta to ecm
ecm_fasta <- fasta[names(fasta) %in% ecm_classification$id]

# Save the filtered files
writeXStringSet(ecm_fasta, "data/full_ITS/ecm_family.fasta")
fwrite(ecm_classification, "data/full_ITS/ecm_family.classification")