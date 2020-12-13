## Get fasta transcript IDs to keep

library(stringr)
library(myTAI)


### Get diamond blast tophits ###
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\"
species <- "GT"
blast_file <-"_nr_tophits.m8" #Change for input

fp <- paste(data_drive_with_tax_report, species, blast_file, sep = "") 
top_hits <- read.delim(fp, header=FALSE, sep="\t")

fp_fasta_output <- fp <- paste(data_drive_with_tax_report, "trimmed_Cm_longest_orfs_cds.fasta", sep = "") 

#Label output
colnames(top_hits) <- c("taxids", "sskingdoms", "skingdoms", "sphylums", "pident", "qseqid", "sseqid", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

### Subset for blast hits with high percent identity ###
before <- nrow(top_hits)

perc <- 95.0
top_hits <- top_hits[top_hits$pident >= perc,]
after <- nrow(top_hits)

my_transcripts <- str_sort(top_hits$qseqid)



fp <- paste(data_drive_with_tax_report, "longest_orfs.cds", sep = "") 
fasta <- read.delim(fp, header=FALSE, sep=">")

fasta_seqs <- fasta[[1]]
fasta_seqs <- fasta_seqs[fasta_seqs!=""]

fasta_transcripts <- fasta[[2]]
fasta_transcripts <- fasta_transcripts[fasta_transcripts!=""]


#isolate_symbol <- function(fasta_line){
#  return(strsplit(fasta_line, "[.]")[[1]][1])
#}

#my_transcripts <- sapply(my_transcripts, isolate_symbol) 




keep <- vector(mode="logical", length=length(fasta_seqs))

for(i in 1:length(fasta_seqs)){
  fasta_line <- fasta_transcripts[i]
  
  for(n in 1:length(my_transcripts)){
    current_transcript <- my_transcripts[n]
    
    if(grepl(current_transcript, fasta_line, fixed=TRUE)){
      keep[i] <- TRUE
      break
    }
  }
}


put_carrot_in <- function(fasta_line){
  return (paste(">",fasta_line,sep=""))
}


keep_fasta <- fasta_transcripts[keep]
keep_fasta <- sapply(keep_fasta, put_carrot_in)

keep_seqs <- fasta_seqs[keep]

my_subsetted_fasta <- c(rbind(keep_fasta,keep_seqs))


fileConn <- file(fp_fasta_output)    
writeLines(my_subsetted_fasta, fileConn)    
close(fileConn)





