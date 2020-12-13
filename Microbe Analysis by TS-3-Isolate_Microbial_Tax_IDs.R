library(stringr)


### Get taxon report ###
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\GT_quants\\"

#Saved from NCBI with lineage information
fp2 <- paste(data_drive_with_tax_report, "tax_report_for_GT_cds.txt", sep = "")
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




taxnames <- microbe_taxa$taxname
all_names <- unique(taxnames)
taxids <- unique(microbe_taxa$taxid)




### Write to excel ###
fp_1 <- paste(data_drive_with_tax_report, "microbes_in_GT_cds.xls", sep = "")


taxon.table <- data.frame(
  tax_names = all_names,
  tax_ids = taxids
)

write.table(taxon.table, fp_1, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t")
