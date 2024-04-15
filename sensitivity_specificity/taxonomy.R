## Make species-level taxonomy database for all tools other than CZID ##

library(dplyr)
library(taxonomizr)

# Prepare database 
# prepareDatabase(getAccession = FALSE)

# Read in existing taxonomy file
taxonomy_csv <- ""

taxonomy_input <- read.csv(taxonomy_csv) %>% 
  distinct() %>% 
  filter(!is.na(taxid))

# New taxids to be added

input <- read.delim(classifier_output, sep = "\t", header = TRUE, row.names = NULL)

# Find any taxids that aren't alreaddy in taxonomy.csv from input file
new_taxids_mixed <- setdiff(input$taxid, taxonomy_input$taxid) %>% 
  unique()

# Sometimes import doesn't work correctly - just use numeric taxids
new_taxids <- as.numeric(grep('^-?[0-9.]+$', new_taxids_mixed, val = T))

# Function to identify type(e.g. bacteria), rank, species name and ID from taxid
taxonomy <- function(taxid) {
  
  rawTax <- taxonomizr::getRawTaxonomy(taxid)
  
  superkingdom <- NA
  species <- NA
  type <- NA
  name <- NA
  rank <- NA
  species_taxid <- NA
  
  # Determine whether bacteria, virus, fungi etc
  
  if (!is.null(rawTax[[1]])) {
    
    if ("superkingdom" %in% names(rawTax[[1]])) {
      superkingdom <- rawTax[[1]][["superkingdom"]]
    }
    
    if ("species" %in% names(rawTax[[1]])) {
      species <- rawTax[[1]][["species"]]
    }
    
    if (!is.na(superkingdom) & superkingdom == "Archaea") {
      
      type <- "Archaea"
      
    } else if (!is.na(superkingdom) & superkingdom == "Viruses") {
      
      type <- "Virus"
      
    } else if (!is.na(superkingdom) & superkingdom == "Bacteria") {
      
      type <- "Bacteria"
      
    } else if (!is.na(superkingdom) & superkingdom == "Eukaryota") {
      
      if ("Fungi" %in% rawTax[[1]]) {
        type <- "Fungi"
        
      } else if (!is.na(species) & species == "Homo sapiens") {
        type <- "Human"
        
      } else {
        
        type <- "Other eukaryote"
        
      }
      
    } else {
      
      type <- NA
    }
  
  # Find taxonomic rak e.g. species, genus  
    ranksRawTax <- names(rawTax[[1]][1])
    
    if (is.na(ranksRawTax)) {
      
      rank <- "unclassified"
      
    } else if (ranksRawTax == "no rank"){
      
      if ("strain" %in% names(rawTax[[1]]) | "species" %in% names(rawTax[[1]])){
        
        rank <- "strain"
        
      } else {
        
        rank <- ranksRawTax
      }
      
    } else {
      
      rank <- ranksRawTax
      
    }
    
    # Find species taxid
    if (!is.na(species)) {
    
      species_taxid <- getId(species)
      
    }
    
    # Find species name
    name <- rawTax[[1]][[1]]
    
  }
  
  if (is.na(species)) {
    
    species <- name
    
  }
    
    # Output
    
    output <- c(name, taxid, type, rank, species, species_taxid)
    
    names(output) <- c("name", "taxid", "type", "rank", "species", "species_taxid")
  
    return(output)
  
}

# Need to add something here to remove commas
taxonomy_output <- lapply(new_taxids, taxonomy) %>%
  do.call(rbind, .) %>%
  as.data.frame(.) %>%
  dplyr::filter(!is.na(taxid)) %>% 
  dplyr::distinct()

write.table(taxonomy_output, taxonomy_csv, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
