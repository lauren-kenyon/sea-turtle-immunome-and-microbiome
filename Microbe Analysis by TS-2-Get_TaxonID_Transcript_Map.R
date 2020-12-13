## Get fasta transcript IDs to keep

library(stringr)
library(myTAI)


### Get diamond blast tophits over 95% identity ###
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\"
blast_file <-"GT_nr_tophits.m8" #Change for input

fp <- paste(data_drive_with_tax_report, blast_file, sep = "") 
top_hits <- read.delim(fp, header=FALSE, sep="\t")

#Label output
colnames(top_hits) <- c("taxids", "sskingdoms", "skingdoms", "sphylums", "pident", "qseqid", "sseqid", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

### Subset for blast hits with high percent identity ###
before <- nrow(top_hits)

perc <- 95.0
top_hits <- top_hits[top_hits$pident >= perc,]
after <- nrow(top_hits)

blast_hits<- str_sort(top_hits$qseqid)

### Get longest orfs cds (what'll be mapped to) ###
fp <- paste(data_drive_with_tax_report, "GT_quants\\trimmed_Cm_longest_orfs_cds.fasta", sep = "") 
fasta <- read.delim(fp, header=FALSE, sep=">")

fasta_seqs <- fasta[[1]]
fasta_seqs <- fasta_seqs[fasta_seqs!=""]

fasta_transcripts <- fasta[[2]]
fasta_transcripts <- fasta_transcripts[fasta_transcripts!=""]


isolate_symbol <- function(fasta_line){
  return(strsplit(fasta_line, "[ ]")[[1]][1])
}

trimmed_gene_symbols <- sapply(fasta_transcripts, isolate_symbol) #list of all possible gene symbols from trimmed cds


### Get associated taxon IDs to gene names/IDs from longest orf cds ###
#
taxon_IDs <- vector(mode="character", length=length(trimmed_gene_symbols))


all_tax_IDs <- top_hits$taxids
qseqs <- top_hits$qseqid

for(i in 1:length(trimmed_gene_symbols)){
  for(r in 1:length(qseqs)){ #Find same 
    qseq = qseqs[r]
    if(trimmed_gene_symbols[i] == qseq){
      taxon_IDs[i] <- all_tax_IDs[r]
      break
    }
  }
}



### ^ have list of gene IDs from cds with associated taxon ids

### now to isolate taxon ids
### Isolate taxids ####
my_tax_ids <- taxon_IDs

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
uniq_tax_ids <- unique(all_tax_ids)

fp_IDs <- paste(data_drive_with_tax_report, "GT_trimmed_cds_taxon_IDs.xls", sep = "") 
write.table(uniq_tax_ids, fp_IDs, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t" )

transcript_ID_map <- data.frame(
  transcript_IDs = trimmed_gene_symbols,
  taxon_ids = taxon_IDs
)



fp_transcript_with_ID <- paste(data_drive_with_tax_report, "GT_transcriptID_taxonIDs_map.xls", sep = "") 
write.table(transcript_ID_map, fp_transcript_with_ID, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t" )


