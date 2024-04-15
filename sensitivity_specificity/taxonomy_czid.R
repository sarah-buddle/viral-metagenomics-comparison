## Make species-level taxonomy database for CZID ##

library(tidyverse)

input <- read.csv(input)

fungal_phyla <- c("Basidiomycota", "Ascomycota", "Chytridiomycota", "Neocallimastigomycota",
                  "Blastocladiomycota", "Zoopagomycota", "Mucoromycota", "Glomeromycota", "Opisthosporidia")


findType <- function(superkingdom_name, phylum_name, species_name) {
  
  if (superkingdom_name == "Archaea") {
    
    type <- "Archaea"
    
  } else if (superkingdom_name == "Viroids") {
    
    type <- "Viroid"
    
  } else if (superkingdom_name == "Bacteria") {
    
    type <- "Bacteria"
    
  } else if (superkingdom_name == "Viruses") {
    
    type <-  "Virus"
    
  } else if (superkingdom_name == "Eukaryota") {
    
    if (species_name == "Homo sapiens") {
      
      type <- "Human"
      
    } else if (phylum_name %in% fungal_phyla) {
      
      type <- "Fungi"
      
    } else {
      
      type <- "Other eukaryote"
      
    }
  } else {
    
      type <- NA
      
      print(species_name)
    
  }
  
  return(type)
  
}

findType <- Vectorize(findType)

findRank <- function(genus_taxid, species_taxid) {
  
  if (genus_taxid > 0) {
    
    if (species_taxid > 0) {
      
      rank <- "species"
      
    } else {
      
      rank <- "genus"
      
    }
    
  } else {
    
    rank <- NA
    
  }
} 

findRank <- Vectorize(findRank)


output <- input %>% 
  dplyr::filter(version_end == "2021-01-22") %>% 
  dplyr::mutate(type = findType(superkingdom_name, phylum_name, species_name)) %>% 
  dplyr::mutate(rank = findRank(genus_taxid, species_taxid)) %>% 
  dplyr::select(c(tax_name, taxid, type, rank, species_name, species_taxid)) %>% 
  dplyr::rename(name = tax_name, species = species_name)

write.csv(output, output_filepath, row.names = FALSE, quote = FALSE)


