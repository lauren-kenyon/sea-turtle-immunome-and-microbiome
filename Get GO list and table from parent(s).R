library(GO.db)
library(biomaRt)
library(stringr)

#Input file with GO terms
data_drive <- "C:\\Users\\laure\\Documents\\Thesis\\WD\\"
#table_name <- "oxid_burst_go_terms.csv"
#GO_fp_for_list <- paste(data_drive, table_name, sep = "")

#immune_GOs <- read.csv(GO_fp_for_list, header=TRUE)
#unsorted_GOs <- immune_GOs[,1]
#immune_GOs <- str_sort(unsorted_GOs)

immune_GOs <- c("GO:0002376")

#Output file
fp_1<- paste(data_drive, "all_go_terms_GO_table.xls", sep = "")

all_my_terms_and_offspring <- c()


#Add offspring
for(i in 1:length(immune_GOs)){
  immune_list <- c(GOBPOFFSPRING[[immune_GOs[i]]])
  immune_list <- c(immune_GOs[i], immune_list)
  all_my_terms_and_offspring <- c(all_my_terms_and_offspring, immune_list)
}

all_my_terms_and_offspring <- str_sort(unique(all_my_terms_and_offspring))
all_my_terms_and_offspring <- all_my_terms_and_offspring[!is.na(all_my_terms_and_offspring)] #Some terms have no offspring -> NA was added to list

get_GO_table <- function(immunome){
  uniq_GO_table <-matrix(nrow=length(immunome), ncol=3)
  colnames(uniq_GO_table) <- c("GO ID", "Term Description", "Definition")
  for(i in 1:length(immunome)){
    GO_ID <- immunome[i]
    GO_data < GOTERM[[GO_ID]-]
    uniq_GO_table[i, 1] = GO_data@GOID
    uniq_GO_table[i, 2] = GO_data@Term
    uniq_GO_table[i, 3] = GO_data@Definition
  }
  return(uniq_GO_table)
}


### Workspace ###
my_GO_table <- get_GO_table(all_my_terms_and_offspring)
write.table(my_GO_table, fp_1, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")