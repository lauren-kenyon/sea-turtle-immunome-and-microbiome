library(biomaRt)
library(stringr)

### Filepath for annotations & immunome outputs ###
#Edit filepaths if naming doesn't follow the same convection 

data_drive <- "C:\\Users\\laure\\Documents\\Thesis\\Annotations\\"

#Should only have to adjust annot_name based on naming of outputs in other scripts
annot_name <- "GT_trinotate_annotation_report.xls" #found in data_drive
annot_fp <- paste(data_drive, annot_name, sep = "")

immunome_name <- paste("immunome_", annot_name, sep = "")
immunome_fp <- paste(data_drive, immunome_name, sep = "")

GO_list_immunome_name <- paste("GO_list_of_", immunome_name, sep = "")
GO_list_immunome_fp <- paste(data_drive, GO_list_immunome_name, sep = "")


#### Get data ####
my_annot <- read.delim(annot_fp, header=TRUE, sep="\t")
my_GO_list_immunome <- read.delim(GO_list_immunome_fp, header=TRUE, sep="\t")
my_immunome <-read.delim(immunome_fp, header=TRUE, sep="\t")


#### Gather transcript count data ####
#Some transcripts and their isoforms are listed multiple times
all_transcripts <-  my_annot[,2]
num_transcripts <- length(unique(all_transcripts)) #repeated isoforms


num_unique_immunume_GO_terms <- nrow(my_GO_list_immunome)
all_immune_transcripts <- my_immunome[,2]
num_immune_transcripts <- length(unique(all_immune_transcripts)) #repeated isoforms


#Keep track of # of pfam hits, blast hits, etc. by column in annotation .xls
col_counts<- vector(mode = "numeric", length = 14)
col_bool<- vector(mode = "logical", length = 14)

#### Go through initial annotation and get stats ####
num_GO_blast_annotated <- 0
num_GO_pfam_annotated <- 0
transcripts_with_GO_hits <- 0 #GO pfam and/or GO blastp
num_annotated <- 0 #num of transcripts annotated

last_transcript_id <- my_annot[1,2] #first transcript ID, used to check for redundant transcript IDs
current_transcript_id <- ""

annotated <- FALSE #Used to check if transcript has any annotations
annot_length <- nrow(my_annot)

for (i in 1:annot_length){
  row <- my_annot[i,]
  current_transcript_id = row[2]

  

  for(c in 3:14){
      if(row[c] != "."){
        col_bool[c] <- TRUE
        annotated <- TRUE
      }
  }
  if((current_transcript_id != last_transcript_id) || i == annot_length){ #TODO: check for string comparison
    if(col_bool[13] || col_bool[14]){ #GO blast or pfam hits
        transcripts_with_GO_hits <- transcripts_with_GO_hits + 1
    }
    if(annotated){
        num_annotated <- num_annotated + 1
    }
    for(k in 3:14){
      if(col_bool[k]){
        col_counts[k] <- col_counts[k]+ 1
        col_bool[k] <- FALSE #Reset
      }
    }
    annotated <- FALSE #Reset 
  }
  last_transcript_id = row[2]
}




result_row_names<- c(
                     "Total transcripts with blastp hits", 
                     "Total transcripts with pfam hits",
                     "Total transcripts with signalP hits", 
                     "Total protein sequences with Eggnog terms",
                     "Total protein sequences matching Kegg pathways",
                     "Total transcripts with GO blastp hits", 
                     "Total transcripts with GO pfam hits", 
                     "Total transcripts with GO blastp and/or GO pfam hits",
                     "Total transcripts with immune GO hits",
                     "Number of unique immune GO terms in annotation",
                     "Number of transcripts in annotation",
                     "Number of transcripts annotated")
values <- c(
            col_counts[7],
            col_counts[8],
            col_counts[9], 
            col_counts[11],
            col_counts[12],
            col_counts[13],
            col_counts[14], 
            transcripts_with_GO_hits,
            num_immune_transcripts,
            num_unique_immunume_GO_terms,
            num_transcripts,
            num_annotated)

annot_results <-matrix(data=values, nrow=12, ncol=1)
rownames(annot_results) <- result_row_names

#look at annot_results variable for stats


