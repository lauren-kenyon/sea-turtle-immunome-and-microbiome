library(GO.db)
library(biomaRt)
library(stringr)
#This script creates a unique, sorted list of the immune GO terms present in the provided annotation


immune_GO_terms <- as.list(GOBPOFFSPRING[["GO:0002376"]])

#Filepath for immunome
data_drive <- "C:\\Users\\laure\\Documents\\Thesis\\Annotations\\"

immunome_annot_name <- "immunome_GT_trinotate_annotation_report.xls" #found in data_drive
fp <- paste(data_drive, immunome_annot_name, sep = "")

#get immunome, as created by isolate GO terms.R
mydata <- read.delim(fp, header=TRUE, sep="\t")

#isolating pertinent data (transcript IDs and GO terms)
newDF <- mydata[,c(1,2,13,14)]

#empty vector to contain all the immunome's (mydata) GO terms
all_annot_GO_terms <- vector(mode = "character", length =20000)

index <- 1
for (i in 1:nrow(newDF)){
  row <- newDF[i,]
  
  #isolate GO term columns
  blast <- row[3]
  pfam <- row[4] 
  
  #Get list of GO terms from each column
  GO_terms_blast <- regmatches(blast, gregexpr("GO:\\d+",blast))
  GO_terms_pfam <- regmatches(pfam, gregexpr("GO:\\d+",pfam))

  list_GO_bl <- GO_terms_blast[[1]] 
  list_GO_pf <- GO_terms_pfam[[1]]
  
  #Add all the GO terms found to the grand list of GO terms
  for (term in list_GO_bl) {
    all_annot_GO_terms[index] <- term
    index <- index + 1
  }
  for (term in list_GO_pf) {
    all_annot_GO_terms[index] <- term
    index <- index + 1
  }
}

#Get rid of repeat GO terms
uniq_annotation_GO_terms <- unique(all_annot_GO_terms)

#Subset our annotation's GO terms for only the ones under the immune category
immunome_GO_terms<-intersect(uniq_annotation_GO_terms, immune_GO_terms)
sorted_GO <-str_sort(immunome_GO_terms)

GO_fp_for_list <- paste(data_drive, "GO_list_of_", immunome_annot_name, sep = "")
write.table(sorted_GO, GO_fp_for_list, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")

                    

