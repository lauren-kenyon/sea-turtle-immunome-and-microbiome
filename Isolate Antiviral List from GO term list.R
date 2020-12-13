library(GO.db)
library(biomaRt)
library(stringr)
immune_genes <- as.list(GOBPOFFSPRING[["GO:0002376"]])

#Filepath for annotations & output
data_drive <- "C:\\Users\\laure\\Documents\\Thesis\\Annotations\\"

immunome_annot_name <- "immunome_GT_trinotate_annotation_report.xls" #found in data_drive

GO_fp_for_list <- paste(data_drive, "\\GO_lists\\GO_list_of_", immunome_annot_name, sep = "")
antiviral_fp_GO_list <- paste(data_drive, "\\GO_lists\\antiviral_GO_list_of_", immunome_annot_name, sep = "")

#input file is just a column ofimmune GO terms
mydata <- read.delim(GO_fp_for_list , header=FALSE, sep="\t")

data_drive_with_GO <- "C:\\Users\\laure\\Documents\\Thesis\\"
fp_antiv <- paste(data_drive_with_GO, "My_GO_slim_antiviral.csv", sep = "")
antiv_table <- read.csv(fp_cat, header=TRUE, sep=",")
GO_antiv <- antiv_table[,1]

my_antiviral_GO_terms <- c()

for(GO_term in GO_antiv){
  my_antiviral_GO_terms <- append(my_antiviral_GO_terms, as.list(GOBPOFFSPRING[[GO_term]]))
}


my_antiviral_GO_terms <- str_sort(unique(my_antiviral_GO_terms))
my_antiviral_GO_terms <- my_antiviral_GO_terms[my_antiviral_GO_terms != "NA"]

#empty vector to contain all the immunome's (mydata) GO terms
all_antiviral_GO_terms <- vector(mode = "character", length =100)

index <- 1


GT_GO_terms <- mydata[,1]
for(i in 1:length(GT_GO_terms)){
  term <- GT_GO_terms[i]
  if(is.element(term, my_antiviral_GO_terms)){
    all_antiviral_GO_terms[index] <- term
    index <- index + 1
  }
}

#Remove empty elements
all_antiviral_GO_terms <- unique(all_antiviral_GO_terms[all_antiviral_GO_terms != ""])

#Get GO descriptions from a GO term list
#input: GO term list

get_GO_table <- function(immunome){
  uniq_GO_table <-matrix(nrow=length(immunome), ncol=3)
  colnames(uniq_GO_table) <- c("GO ID", "Term Description", "Definition")
  for(i in 1:length(immunome)){
    GO_ID <- immunome[i]
    GO_data <- GOTERM[[GO_ID]]
    uniq_GO_table[i, 1] = GO_data@GOID
    uniq_GO_table[i, 2] = GO_data@Term
    uniq_GO_table[i, 3] = GO_data@Definition
  }
  return(uniq_GO_table)
}

my_go_table <- get_GO_table(all_antiviral_GO_terms) #includes descriptions
only_the_GOs <- my_go_table[,1] #Good for UPsetR

write.table(only_the_GOs, antiviral_fp_GO_list, col.names=FALSE, row.names = FALSE, quote = FALSE, sep = "\t")