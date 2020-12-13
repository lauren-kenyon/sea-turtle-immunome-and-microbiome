library(stringr)
library(myTAI)

help(package="myTAI")

### Get diamond blast tophits ###
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\"
species <- "Ei"
blast_file <-"_nr_tophits.m8" #Change for input

fp <- paste(data_drive_with_tax_report, species, blast_file, sep = "") 
top_hits <- read.delim(fp, header=FALSE, sep="\t")


#Label output
colnames(top_hits) <- c("taxids", "sskingdoms", "skingdoms", "sphylums", "pident", "qseqid", "sseqid", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

### Subset for blast hits with high percent identity ###
before <- nrow(top_hits)

perc <- 95.0
top_hits <- top_hits[top_hits$pident >= perc,]
after <- nrow(top_hits)

removed <- before-after
removed_perc <- (100*removed/before)

sprintf("Removed %i blast hits under %.2f%%", (before-after), perc)
sprintf("Removed %.2f%% of hits", removed_perc)

### Isolate taxids ####
my_tax_ids <- top_hits$taxids

#Separate hits with multiple taxids (separted with ;)
singular_tax_ids <- my_tax_ids[!grepl(";",my_tax_ids)]
mult_tax_ids <- my_tax_ids[grepl(";",my_tax_ids)]

separate_tax_ids<- vector(mode="character", length=length(mult_tax_ids)*5) #Allocate enough memory, add many empty spots
index <- 1
for(i in 1:length(mult_tax_ids)){
  tax_ids <- strsplit(mult_tax_ids[i],";")[[1]]
  
  for(t in tax_ids){
    separate_tax_ids[index] <- t
    index <- index + 1
  }
}

separate_tax_ids <- separate_tax_ids[separate_tax_ids!=""] #Remove empty spots from list

all_tax_ids <- append(singular_tax_ids, separate_tax_ids)
#Check
length(all_tax_ids) #Should match error file with # blast hits

fp_r <- paste(data_drive_with_tax_report, species, "__tax_ids_over_pident", perc, ".xls",sep = "") 

write.table(all_tax_ids, fp_r, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t" )

