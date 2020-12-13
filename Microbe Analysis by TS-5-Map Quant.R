library(stringr)


### Get taxon report ###
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\GT_quants\\"

fp_transcript_map <- paste(data_drive_with_tax_report, "GT_microbe_tax_ID_map.xls", sep = "") 
map <- read.delim(fp_transcript_map, header=TRUE, sep="\t")

fp_all_taxon <- paste(data_drive_with_tax_report, "microbes_in_GT_cds.xls", sep = "") 
microbes <- read.delim(fp_all_taxon, header=TRUE, sep="\t")

microbe_tax_ids <- sapply(microbes$tax_ids, toString)

quant_mappings <- data.frame(
  tax_names = microbes$tax_names,
  tax_ids = microbe_tax_ids)
  
#For each quant
for(k in 1:44){
  fp_quant <- paste(data_drive_with_tax_report, "quant", k, ".sf", sep = "")
  quant <- read.delim(fp_quant, header=TRUE, sep="\t")
  
  
  
  transcript_names <- map$tax_names
  quant_names <- quant$Name
  
  v_num_reads <- vector(mode="numeric", length=length(microbe_tax_ids))
  v_TPM <- vector(mode="numeric", length=length(microbe_tax_ids))
  
  error <- 0
  error2 <- 0
  for(i in 1:length(quant_names)){ #For every mapped transcript in salmon quant file
    quant_n <- quant_names[i]
    if(is.element(quant_n, transcript_names)){
      index <- which(quant_n == transcript_names)[[1]] #Find which row it maps back to in orginal .cds
      tax_ids <- map[index, 2] #Get associated taxon id(s)
      all_tax_ids <- strsplit(tax_ids, "\\s+")[[1]]
      
      for(n in 1:length(all_tax_ids)){ #For each taxon ID -> update count
        taxon_id <- all_tax_ids[n]
        if(is.element(taxon_id,microbe_tax_ids)){
          index_microbe <- which(taxon_id == microbe_tax_ids)[[1]]
          
          v_num_reads[index_microbe] <- v_num_reads[index_microbe] + quant[i,]$NumReads
          v_TPM[index_microbe] <- v_TPM[index_microbe] + quant[i,]$TPM
        }else{
          error2 <- error2 + 1
        }
      }
    }else{a
      error <- error + 1
    }
  }
  
  
  quant_mappings <- cbind(quant_mappings,v_num_reads)
  quant_mappings <- cbind(quant_mappings,v_TPM)
}

fp_quant <- paste(data_drive_with_tax_report, "quant_mappings.xls", sep = "")
write.table(quant_mappings, fp_quant, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t")