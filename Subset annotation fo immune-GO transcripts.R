library(GO.db)
library(biomaRt)
immune_genes <- as.list(GOBPOFFSPRING[["GO:0002376"]])

#Filepath for annotations & output
data_drive <- "C:\\Users\\laure\\Documents\\Thesis\\Annotations\\"
annot_name <- "GT_trinotate_annotation_report.xls" #found in data_drive
fp <- paste(data_drive, annot_name, sep = "")
fp_for_immunome <- paste(data_drive, "immunome_", annot_name, sep = "")

#get annotation, as created by Trinotate pipeline
mydata <- read.delim(fp, header=TRUE, sep="\t")


#isolating pertinent data (transcript IDs and GO terms)
newDF <- mydata[,c(1,2,13,14)]
#empty vector to contain row #s of newDF that have a immune-related GO term
immune_gene_indexes <- vector(mode = "numeric", length =20000)

index <- 1

for (i in 1:nrow(newDF)){
  row <- newDF[i,]
  blast <- row[3]
  pfam <- row[4] 
  
  if(grepl("GO", blast, fixed=TRUE) || grepl("GO", pfam, fixed=TRUE)){
    for (gene in immune_genes) {
      b <- grepl(gene, blast, fixed=TRUE)
      p <- grepl(gene, pfam, fixed=TRUE)
      if(p || b){
        immune_gene_indexes[index] <- i
        index <- index + 1
        break
      }
    }
  }
  
}

immunome <- mydata[immune_gene_indexes,]
write.table(immunome, fp_for_immunome, row.names = FALSE, quote = FALSE, sep = "\t")