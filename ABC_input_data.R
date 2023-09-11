### Create ABC input data files
## MURI MMARINeDNA 
## August 2023

library(tidyverse)
library(phyloseq)
library(microViz)

wd <- ""
dada2_output <- ""
sample_data <- "new_prey_meta_7.26.22.csv"
project_name<- "MURI_0304"

setwd(wd)

### Make ps object from dada2 output ---------------------------------------------

# load dada2 output
load(dada2_output)

# get sample metadata
samdf <- read.csv(sample_data) 

# create master phyloseq object
ps.raw <- phyloseq(otu_table(cleaned.seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(joined_old_new_taxa))

# shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.raw))
names(dna) <- taxa_names(ps.raw)
ps.raw <- merge_phyloseq(ps.raw, dna)
taxa_names(ps.raw) <- paste0("ASV", seq(ntaxa(ps.raw)))

### Create dataframe for ABC ---------------------------------------------------

# convert ps back to dataframes
readcount.table <- as.data.frame(cleaned.seqtab.nochim)
taxon.table <- as.data.frame(joined_old_new_taxa) %>% rownames_to_column(var = "ASV")
metadata.table <- samdf %>% rownames_to_column(var = "Sample")
reference.table <- as.data.frame(refseq(ps.raw)) %>% rownames_to_column("ASVname") %>% 
  rename("DNAseq" = x)

# merge dataframes into long table, keep ASVs without species ID
readcount.table.long <- readcount.table %>% 
  rownames_to_column(var = "Sample") %>% 
  pivot_longer(-Sample, names_to = "ASV", values_to = "read.count") %>% 
  left_join(metadata.table, by = c("Sample" = "Sample")) %>% 
  select(Sample, read.count, ASV) %>% 
  left_join(taxon.table, by = c("ASV" = "ASV")) %>% 
  select(Sample, Species, ASV) %>% 
  mutate(Species = case_when(is.na(Species) ~ ASV, TRUE ~ Species)) %>% 
  select(-ASV) %>% 
  group_by(Sample, Species) %>% 
  mutate(sp.read.count = sum(read.count)) %>% 
  slice_head() %>% 
  select(-read.count)

# pivot wider by species, remove ASVs without species assignment
readcount.table.species <- readcount.table.long %>% 
  filter(!grepl("ASV",Species)) %>% 
  pivot_wider(names_from = "Species", values_from = sp.read.count) %>% 
  mutate("bio_rep" = 1, .after = Sample)

# pivot wider again, store ASVs without species assignment
readcount.table.asv <- readcount.table.long %>% 
  filter(grepl("ASV",Species)) %>% 
  pivot_wider(names_from = "Species", values_from = sp.read.count) %>% 
  mutate("bio_rep" = 1, .after = Sample)

# Export for amplification bias correction
write.csv(readcount.table.species, file = paste0(project_name, "_ABC_readcount_species_", format(Sys.Date(), "%m%d%Y"), ".csv"), row.names = FALSE)
write.csv(readcount.table.asv, file = paste0(project_name, "_ABC_readcount_ASV_", format(Sys.Date(), "%m%d%Y"), ".csv"), row.names = FALSE)
save(reference.table, file = paste0(project_name, "_ASVreference.Rdata"))
