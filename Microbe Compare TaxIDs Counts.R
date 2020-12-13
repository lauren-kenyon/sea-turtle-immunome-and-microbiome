### Import Data ###

#problem with multiple tax IDs
tax_report_name <- "_95_lineages"
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\"


rename_ids <- function(taxon_table){
  tax_ids <- taxon_table$tax_ids
  if(length(which(tax_ids == 1760)) >0){
    row_num <- which(tax_ids == 1760)
    taxon_table[row_num, 1] <- "Actinobacteria high GC Gram+"
  }
  if(length(which(tax_ids == 203693)) >0){
    row_num <- which(tax_ids == 203693)
    taxon_table[row_num, 1] <- "Nitrospira [class]"
    
  }
  if(length(which(tax_ids == 1234)) >0){
    row_num <- which(tax_ids == 1234)
    taxon_table[row_num, 1] <- "Nitrospira [Genus]"
  }

  return(taxon_table)
}





GT <- read.delim(paste(data_drive_with_tax_report, "GT", tax_report_name, "_all_taxa_counts.xls", sep = ""), header=TRUE, sep="\t")
GT <- rename_ids(GT)
rownames(GT) <- GT$tax_names

Dc <- read.delim(paste(data_drive_with_tax_report, "Dc", tax_report_name, "_all_taxa_counts.xls", sep = ""), header=TRUE, sep="\t")
Dc <- rename_ids(Dc)
rownames(Dc) <- Dc$tax_names

Cc <- read.delim(paste(data_drive_with_tax_report, "Cc", tax_report_name, "_all_taxa_counts.xls", sep = ""), header=TRUE, sep="\t")
Cc <- rename_ids(Cc)
rownames(Cc) <- Cc$tax_names

Ei <- read.delim(paste(data_drive_with_tax_report, "Ei", tax_report_name, "_all_taxa_counts.xls", sep = ""), header=TRUE, sep="\t")
Ei <- rename_ids(Ei)
rownames(Ei) <- Ei$tax_names





#List of all microbes
full_names <- unique(c(GT$tax_names, Dc$tax_names, Cc$tax_names, Ei$tax_names))
full_taxids <- unique(c(GT$tax_ids, Dc$tax_ids, Cc$tax_ids, Ei$tax_ids))

### Get counts between species ###

species <- list(Cc, GT, Dc, Ei)
counts_and_percents <- vector(mode="list", length=length(species)*2)
index <- 1

for(i in 1:length(species)){
  s <- species[[i]]
  counts <- vector(mode="numeric", length=length(full_names))
  percents <- vector(mode="numeric", length=length(full_names))
  
  for(n in 1:length(full_names)){
    taxname <- full_names[n]
    
    if(is.na(s[taxname,][1])){ #If taxname is not in this species' list
      counts[n] <- 0
      percents[n] <- 0.0
    }else{
      counts[n] <- s[taxname,]$counts
      percents[n] <- s[taxname,]$percent
    }
  }
  counts_and_percents[[index]] <- counts
  index <- index +1
  counts_and_percents[[index]] <- percents
  index <- index +1
}


### Construct output ###
taxon.table <- data.frame(
  tax_names = full_names,
  tax_ids = full_taxids,
  Cc_counts = counts_and_percents[[1]],
  Cc_percent = counts_and_percents[[2]],
  Cm_counts = counts_and_percents[[3]],
  Cm_percent = counts_and_percents[[4]],
  Dc_counts = counts_and_percents[[5]],
  Dc_percent = counts_and_percents[[6]],
  Ei_counts = counts_and_percents[[7]],
  Ei_percent = counts_and_percents[[8]]
)


fp_r <- paste(data_drive_with_tax_report, "all_taxa_lineages_with_counts", tax_report_name, "_counts.xls", sep = "")
write.table(taxon.table, fp_r, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t")


