library(UpSetR)
library(GO.db)
library(biomaRt)

data_drive_with_GO <- "C:\\Users\\laure\\Documents\\Thesis\\Annotations\\GO_list_of_immunome_"

annot_name <- "_trinotate_annotation_report.xls" #found in data_drive

#Get data
Cc_fp <- paste(data_drive_with_GO, "Cc", annot_name, sep = "")
Dc_fp <- paste(data_drive_with_GO, "Dc", annot_name, sep = "")
Ei_fp <- paste(data_drive_with_GO, "Ei", annot_name, sep = "")
GT_fp <- paste(data_drive_with_GO, "GT", annot_name, sep = "")

Cc_Go_list <- read.csv(Cc_fp, header=FALSE, sep="\t")
Dc_Go_list <- read.csv(Dc_fp, header=FALSE, sep="\t")
Ei_Go_list <- read.csv(Ei_fp, header=FALSE, sep="\t")
GT_Go_list <- read.csv(GT_fp, header=FALSE, sep="\t")

Cc_ls <- Cc_Go_list[,1]
Dc_ls <- Dc_Go_list[,1]
Ei_ls <- Ei_Go_list[,1]
GT_ls <- GT_Go_list[,1]

#UpsetR

listInput <- list("C. carretta" = Cc_ls, "D. coriacea" = Dc_ls, "E. imbricata" = Ei_ls, "C. mydas" = GT_ls)


upset(fromList(listInput), order.by = "freq", point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Number of Shared Immunome GO terms", sets.x.label = "Immunome GO terms", 
      #text.scale = c(1.3, 1.3, 1.3, 1, 1.5, 2),
      text.scale = c(1.3, 1.3, 1, 1, 2, 2),
      #main.bar.color = "darkred", sets.bar.color = "firebrick3", matrix.color = "red")
      queries = list(list(query = intersects, params = list("D. coriacea", "C. mydas", "C. carretta", "E. imbricata"), color = "red", active = T), list(query = intersects, params = list("C. mydas", "C. carretta", "E. imbricata"), color = "blue", active = T)))


others <- union (Cc_ls, Ei_ls)
uniq_Dc_best <- setdiff(Dc_ls, others)


#find all gO terms unique to immunome1 compared to other immunomes 
#input: GO term list and a  vectors of GO term lists

uniq_to <- function(immunome1, other_immunomes) {
  others <- unlist(other_immunomes[1])
  
  i <- 2
  while(i <= length(other_immunomes)) {
    others <- union(others, unlist(other_immunomes[i]))
    i <- i + 1
  }
  uniq_to_immunome1 <- setdiff(immunome1, others)
  return (uniq_to_immunome1)
}


#find all gO terms uniquely shared in allimmunomes of interest compared to other immunomes 
#input: two lists of GO term lists of any size

uniq_shared <- function(immunomes_of_interest, other_immunomes){
  grand_immunome_list <- unlist(immunomes_of_interest[1])
  i <- 2
  
  while(i <= length(immunomes_of_interest)) {
    grand_immunome_list <- intersect(grand_immunome_list, unlist(immunomes_of_interest[i]))
    i <- i + 1
                                     
  }
  others <- unlist(other_immunomes[1])
  
  i <- 2
  while(i <= length(other_immunomes)) {
    others <- union(others, unlist(other_immunomes[i]))
    i <- i + 1
  }
  uniq_shared_to_immunomes_of_interest <- setdiff(grand_immunome_list, others)
  
  return(uniq_shared_to_immunomes_of_interest)
}

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


###### Workspace/Scrap
#list(Ei_ls, Cc_ls, Dc_ls, GT_ls)

uniq_Dc <- uniq_to (Dc_ls, list(Ei_ls, Cc_ls, GT_ls))
dc <- get_GO_table(uniq_Dc)   
uniq_Cc <- uniq_to (Cc_ls, list(Ei_ls, Dc_ls, GT_ls))
cc <- get_GO_table(uniq_Cc)   
uniq_Ei <- uniq_to (Ei_ls, list(Cc_ls, Dc_ls, GT_ls))
ei <- get_GO_table(uniq_Ei)   
uniq_GT <- uniq_to (GT_ls, list(Ei_ls, Cc_ls, Dc_ls))
gt <- get_GO_table(uniq_GT)   




uniq_s <- uniq_shared(list(Dc_ls, Cc_ls, GT_ls, Ei_ls), list())

#a <- get_GO_table(uniq_Cc_Ei)
c <- get_GO_table(uniq_s)   
    
fp_for_table <- "C:\\Users\\laure\\Documents\\Thesis\\core_GO_term_set.xls"
write.table(c, fp_for_table, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
