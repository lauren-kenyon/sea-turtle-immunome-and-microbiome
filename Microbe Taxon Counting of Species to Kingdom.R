library(stringr)

# Purpose: Track the lineage of microbial hits, get info on genus, phylum, etc. counts
#
# Therefore the % in this output are of the total microbial hits for that species and will add to > 100%
# since if 90% of the hits are classified as 'Bacteria', it will report the taxa underneath it. For instance,
# 80% Proteobacteria, 5% Chloroflexi, and 5% Firmicutes could make up the total 90%.
#
# NOTE: this script requires an intermediate file from https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
# - Use output file from this script as input for the Taxonomy name/is status report
# - This taxon report supplies the names for the taxon ids for the final output

### Get taxon report ###
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\"

#Saved from NCBI with lineage information
tax_report_name <- "Cc_95_lineages"
fp2 <- paste(data_drive_with_tax_report, tax_report_name, ".txt", sep = "")
taxon <- read.delim(fp2, header=TRUE, sep="\t")

taxon <- taxon[c(TRUE,FALSE)] #Get rid of even "|" columns

#remove last 2 extran. lines
length <- nrow(taxon)
taxon <- taxon[-c(length-1, length),]

  
is_microbe <- function (lineage){
  my_tax_categories <- c("2","2157","10239")
  sep_lin <- strsplit(lineage, "\\s+")
  microbe <- FALSE
  
  for(id in my_tax_categories){
    if(is.element(id, sep_lin[[1]])){
      microbe <- TRUE
      break
    }
  }
  return(microbe)
}


lineages = taxon$lineage
bool <- sapply(lineages,is_microbe)

microbe_taxa <- taxon[bool,]
lineages = microbe_taxa$lineage

separate_lineage <- function (lineage){
  sep_lin <- strsplit(lineage, "\\s+")
  return (sep_lin[1])
}

ls_lineages <- sapply(lineages,separate_lineage)

all_taxa_level_IDS <- unique(unlist(ls_lineages, use.names=FALSE))


tax_counts <- vector(mode="numeric", length=length(all_taxa_level_IDS))
names(tax_counts) <- all_taxa_level_IDS


### Count taxonomy occurrences  ###
for(i in 1:length(ls_lineages)){
  lineage_list <- ls_lineages[[i]]
  for(id in lineage_list){
    tax_counts[id] <- tax_counts[id]+1
  }
}

### Write to excel ###

#Need names for taxa, return to NCBI taxonomy
fp_r <- paste(data_drive_with_tax_report, tax_report_name, "_all_taxa_levels_ID.xls", sep = "") 
write.table(all_taxa_level_IDS, fp_r, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t")


####### STOP AND GET TAX REPORT FILE #########
#name as Species_blastpercent_lineages_all_taxa.txt


fp_names <- paste(data_drive_with_tax_report, tax_report_name, "_all_taxa.txt", sep = "") 
all_taxa <- read.delim(fp_names, header=TRUE, sep="\t")

#remove last 2 extran. lines
length <- nrow(all_taxa)
all_taxa <- all_taxa[-c(length-1, length),]
#Get names



fp_1 <- paste(data_drive_with_tax_report, tax_report_name, "_all_taxa_counts.xls", sep = "")
ls <- list(tax_counts)
percent_of_blasts <- ls[[1]]/nrow(microbe_taxa)

taxon.table <- data.frame(
  tax_names = all_taxa$taxname,
  tax_ids = all_taxa_level_IDS,
  counts = tax_counts,
  percent = percent_of_blasts
)

write.table(taxon.table, fp_1, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t")
