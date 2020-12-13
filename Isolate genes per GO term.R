library(GO.db)
library(biomaRt)

#input; CHANGE 
my_species <- "GT"


#Data from annotation stats; used for %
num_trans_in_annot <- c(280711,	132146,	347717,	489355)
num_trans_annotated <- c(117483,	69723,	150987,	188047)
species <- c("Ei", "Cc", "Dc", "GT")
names(num_trans_in_annot) <- species
names(num_trans_annotated) <- species


#Filepath for annotations & output
data_drive <- "C:\\Users\\laure\\Documents\\Thesis\\Annotations\\"
immunome_annot_name <- paste("immunome_", my_species, "_trinotate_annotation_report.xls", sep="") #found in data_drive
fp <- paste(data_drive, immunome_annot_name, sep = "")
fp_output <- paste(data_drive, "Genes_under_GO_terms_for_", immunome_annot_name, sep = "")

#get immunome, as created by Trinotate pipeline
immunome <- read.delim(fp, header=TRUE, sep="\t")

#getimmune GO term list
GO_fp_for_list <- paste(data_drive, "GO_lists\\GO_list_of_", immunome_annot_name, sep = "")
immune_GOs <- read.delim(GO_fp_for_list, header=TRUE, sep="\t")
immune_GOs <- immune_GOs[,1]

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

diff <- rows_counts_with_all_iso - rows_counts_no_iso #Get # of isoforms removed for each GO term

########################### Scrap
### Get Data ###

#Get # of transcripts by GO term
rows_counts_no_iso["GO:2000914"]
  #NA means no GO term associated
rows_counts_no_iso["GO:2000974"]


subsetted_immunome_per_GO <- vector(mode="list", length=length(immune_GOs))
for(i in 1:length(immune_GOs)){
  subsetted_immunome_per_GO[[i]] <- list(immunome[rows_counts_no_iso[i],])
}
names(subsetted_immunome_per_GO) <- immune_GOs

#Replace with desired GO
write.table(subsetted_immunome_per_GO["GO:0001777"], fp_output , row.names = FALSE, quote = FALSE, sep = "\t")
write.table(subsetted_immunome_per_GO[[5]], fp_output , row.names = FALSE, quote = FALSE, sep = "\t")



#TODO: 
#With subsetted immunome per GO term get isolated list of gene symbols
#Find a way to print all subsets in 1 excel file
#Add statistics
