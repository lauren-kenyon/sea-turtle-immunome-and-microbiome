library(GO.db)
library(biomaRt)
library(stringr)

### Get GO term list ###
data_drive_with_GO <- "C:\\Users\\laure\\Documents\\Thesis\\"
#GO_list_name <- "GO_list_of_immunome_Ei_trinotate_annotation_report.xls" #found in data_drive
#fp <- paste(data_drive_with_GO, "Annotations\\GO_lists\\",  GO_list_name, sep = "")

fp <- paste(data_drive_with_GO, "core_GO_term_list.csv", sep = "")

GO_table <- read.csv(fp, header=FALSE, sep="\t")
GO_list <- GO_table[[1]]

### Define GO slims ###
fp_cat <- paste(data_drive_with_GO, "My_GO_slim_cat_for_core.csv", sep = "")
cat_table <- read.csv(fp_cat, header=TRUE, sep=",")
cat <- cat_table[,2]
GO_cat <- cat_table[,1]

#Format GO category labels, "\n" included in excel
for(c in 1:length(cat)){
 cat[c] <-  gsub("([\\])n","\n", cat[c])
}

#Get all descendants 
immune_offspring <- list()

num_in_cat <- vector(mode="numeric", length=length(GO_cat))

for(i in 1:length(GO_cat)){
  immune_list <- list(GOBPOFFSPRING[[GO_cat[i]]])
  immune_offspring <- append(immune_offspring, immune_list)
  num_in_cat[i] <- length(immune_list[[1]])
}

#All possible offspring for GO category
names(immune_offspring) <- GO_cat

observed_offspring <- vector(mode="list", length=length(GO_cat))
names(observed_offspring) <- GO_cat


### Get category counts ###
count_vec <- vector(mode="numeric", length=length(GO_cat))

for(i in 1:length(GO_list)){
  GO_term <- GO_list[i]
  for(c in 1:length(GO_cat)){
    if(is.element(GO_term, immune_offspring[[c]])){
      count_vec[c] = count_vec[c]+1
      observed_offspring[[c]] <- append(observed_offspring[[c]], GO_term)
    }
  }
}


### Plot

num_absent <- num_in_cat - count_vec
max <- max(count_vec)

at <- seq(from = 0, to = 120, by = 10)
par(mar=c(4,14,4,4)+.1, mgp=c(8,0,3,0))
GO_slim <- barplot(count_vec, main="GO slim categories and counts of GO terms", 
        col=c("darkblue"),
        horiz=TRUE, las = 1,axes=FALSE,
        xlim=c(0,120),
        names.arg=cat, cex.names = 0.7)
axis(side = 1, at = at)
mtext("Number of GO terms in category", side = 1, line = 2, cex = 0.9)
mtext("Immune GO categories", side = 2, line = 12, cex = 0.9)
text(x = count_vec, y = GO_slim, label = count_vec, pos = 4, cex = 0.8, col = "black")
text(x = (count_vec + 10), y = GO_slim, label = num_absent, pos = 4, cex = 0.8, col = "grey")
#Rerun once plot size is large enough
legend("topright", 
       legend=c("Observed", "Absent"),
       col=c("black", "grey"),
       pch=c(19),
       cex=.9,
       pt.cex=.9,
       text.col="black")



       