library(GO.db)
library(biomaRt)
library(stringr)
library(scales)
library(DataCombine)
#input; CHANGE 

data_drive1 <- "C:\\Users\\laure\\Documents\\Thesis\\WD\\"
fp_1<- paste(data_drive1, "all_go_terms_GO_table_counts.xls", sep = "")

#getimmune GO term list
table_name <- "all_go_terms_GO_table.xls"
GO_fp_for_list <- paste(data_drive1, table_name, sep = "")

immune_GOs <- read.delim(GO_fp_for_list, header=TRUE, sep="\t")
unsorted_GOs <- immune_GOs[,1]
immune_GOs <- str_sort(unsorted_GOs)

#Data from annotation stats; used for %
num_trans_in_annot <- c(132146,	489355, 347717, 280711)
num_trans_annotated <- c(69723,	188047, 150987, 117483)
num_immune_trans_annotated <- c(5888, 12555, 9603, 8706)


species <- c("Cc", "Cm", "Dc", "Ei")


names(num_trans_in_annot) <- species
names(num_trans_annotated) <- species
names(num_immune_trans_annotated) <- species

#Store all the data between different species for output table
GO_term_counts <- vector(mode="list", length=length(species))
GO_term_trans_percentages <- vector(mode="list", length=length(species))
GO_term_immmune_trans_percentages <- vector(mode="list", length=length(species))

names(GO_term_counts) <- species
names(GO_term_trans_percentages) <- species
names(GO_term_immmune_trans_percentages) <- species


for (my_species in species){

  
  #Filepath for annotations & output
  data_drive <- "C:\\Users\\laure\\Documents\\Thesis\\"
  immunome_annot_name <- paste("immunome_", my_species, "_trinotate_annotation_report.xls", sep="") #found in data_drive
  fp <- paste(data_drive, "Annotations\\", immunome_annot_name, sep = "")
  fp_output <- paste(data_drive, "Genes_under_GO_terms_for_", immunome_annot_name, sep = "")
  
  #get immunome, as created by Trinotate pipeline
  immunome <- read.delim(fp, header=TRUE, sep="\t")
  

  

  #Isolating pertinent data (transcript IDs and GO terms)
  newDF <- immunome[,c(1,2,13,14)]
  
  rows_with_immune_GO_term <- vector(mode="list", length=length(immune_GOs)) #Each element corresponds with a list of rows (in immunome) for a specific GO term
  
  ### Get Counts ###
  
  for (i in 1:nrow(newDF)){ #For each immune transcript in immunome there may be multiple GO terms
    
    row <- newDF[i,]
    blast <- row[3]
    pfam <- row[4] 
    
    for(t in 1:length(immune_GOs)) { #Use whole list of immune GO terms and make a list of which rows in immunome have it
      term <- immune_GOs[t]
      b <- grepl(term, blast, fixed=TRUE)
      p <- grepl(term, pfam, fixed=TRUE)
      if(p || b){ #If specific GO term is found, add row # to the end of its row list
        row_list <- rows_with_immune_GO_term[[t]]
        
        if(is.null(row_list[1])){ #Don't want NULL as 1st element
          rows_with_immune_GO_term[[t]] <- c(i)
        }else{
          rows_with_immune_GO_term[[t]] <- append(row_list, i)
        }
      }
    }
  }
  
  #### Collapse on isoforms ###
  before_collapse <- rows_with_immune_GO_term #save the compleete list of rows
  
  for (i in 1:length(immune_GOs)){
    
    GO_term_row <- rows_with_immune_GO_term[[i]] #Get row list for specific GO term
    
    first_row_num <- GO_term_row[1]
    last_transcript_id <- newDF[first_row_num,2] #Get first transcript ID, save as last
    current_transcript_id <- ""
   
    remove <- c() #isoform indexes to remove
    
    if(length(GO_term_row) > 1){
      for (t in 2:length(GO_term_row)){ #Already got transcript #1
        row_num <- GO_term_row[t]
        current_transcript_id = newDF[row_num, 2]
        
        if(current_transcript_id == last_transcript_id){
          remove <- append(remove, t)
        }
        last_transcript_id <- current_transcript_id #Update last transcript ID
      }
      if(length(remove) != 0){ #Only if there's something to remove
        GO_term_row <- GO_term_row[-remove]
        rows_with_immune_GO_term[[i]] <- GO_term_row
      }
    }
  }
  
  
  ### Get counts ###
  rows_counts_no_iso <- vector(mode="numeric", length=length(immune_GOs))
  rows_counts_with_all_iso<- vector(mode="numeric", length=length(immune_GOs))
  
  for(i in 1:length(immune_GOs)){
    rows_counts_with_all_iso[i] <- length(before_collapse[[i]])
    rows_counts_no_iso[i] <- length(rows_with_immune_GO_term[[i]])
  }
  
  names(rows_counts_with_all_iso) <- immune_GOs
  names(rows_counts_no_iso) <- immune_GOs
  
  #diff <- rows_counts_with_all_iso - rows_counts_no_iso #Get # of isoforms removed for each GO term


  percent_of_annotated <- rows_counts_no_iso/num_trans_annotated[my_species]
  percent_of_immune_annotated <- rows_counts_no_iso/num_immune_trans_annotated[my_species]
  
  GO_term_counts[[my_species]] <- rows_counts_no_iso
  GO_term_trans_percentages[[my_species]]
  GO_term_immmune_trans_percentages[[my_species]] <-percent_of_immune_annotated
  
}
  
GO_term_counts

GO.data <- data.frame(
  GO_terms = immune_GOs,
  Cc_num = GO_term_counts[["Cc"]],
  Cm_num = GO_term_counts[["Cm"]],
  Dc_num = GO_term_counts[["Dc"]],
  Ei_num = GO_term_counts[["Ei"]],
  Cc_p = label_percent(accuracy = 0.001)(GO_term_immmune_trans_percentages[["Cc"]]),
  Cm_p = label_percent(accuracy = 0.001)(GO_term_immmune_trans_percentages[["Cm"]]),
  Dc_p = label_percent(accuracy = 0.001)(GO_term_immmune_trans_percentages[["Dc"]]),
  Ei_p = label_percent(accuracy = 0.001)(GO_term_immmune_trans_percentages[["Ei"]])
)



broad_cat <- c("", "Counts of transcripts mapped to GO term","","","","Percent of immune transcripts annotated","","","")
colnames(GO.data) <- broad_cat
my_cols <- c("GO term ID", species, species)

GO.data <- InsertRow(GO.data, NewRow = my_cols, RowNum = 1)

write.table(GO.data, fp_1, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t" )

