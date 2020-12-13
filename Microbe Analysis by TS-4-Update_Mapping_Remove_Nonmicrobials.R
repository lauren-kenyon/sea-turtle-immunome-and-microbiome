library(stringr)


### Get taxon report ###
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\GT_quants\\"

fp_Map <- paste(data_drive_with_tax_report, "GT_transcriptID_taxonIDs_map.xls", sep = "")
fp_Microbes <- paste(data_drive_with_tax_report, "microbes_in_GT_cds.xls", sep = "")


map <- read.delim(fp_Map, header=TRUE, sep="\t")

microbes <- read.delim(fp_Microbes, header=TRUE, sep="\t")
microbe_IDs <- microbes$tax_ids


for(i in 1:nrow(map)){
  tax_ids <- map[i, 2]
  if(grepl(";",tax_ids)){
    t<- strsplit(map[i,2],";")[[1]]
    a <- t[is.element(t,microbe_IDs)]
    map[i, 2] <- paste(a, collapse=" ")
  }
}

new_map <- data.frame(
  tax_names = map$transcript_IDs,
  tax_ids = map$taxon_ids
)


fp_transcript_with_ID <- paste(data_drive_with_tax_report, "GT_microbe_tax_ID_map.xls", sep = "") 
write.table(new_map, fp_transcript_with_ID, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t" )


