library(stringr)


### Get taxon report ###
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\"

#Saved from NCBI with lineage information
tax_report_name <- "Ei_95_lineages"
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




taxnames <- microbe_taxa$taxname
all_names <- unique(taxnames)
taxids <- unique(microbe_taxa$taxid)


tax_counts <- vector(mode="numeric", length=length(all_names))
names(tax_counts) <- all_names

### Count taxonomy occurrences  ###
for(i in 1:nrow(microbe_taxa)){
  row <- microbe_taxa [i,]
  tax_name <- row[["taxname"]]
  
  tax_counts[tax_name] <- tax_counts[tax_name] + 1
}

### Write to excel ###
fp_1 <- paste(data_drive_with_tax_report, tax_report_name, "_counts.xls", sep = "")
ls <- list(tax_counts)
percent_of_blasts <- ls[[1]]/nrow(microbe_taxa)

taxon.table <- data.frame(
  tax_names = names(ls[[1]]),
  tax_ids = taxids,
  counts = ls[[1]],
  percent = percent_of_blasts
)

write.table(taxon.table, fp_1, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t")
