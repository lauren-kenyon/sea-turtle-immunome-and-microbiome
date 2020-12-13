
library(RColorBrewer)
library(ggplot2)
library(stringr)
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\WD\\"

fp_counts <- paste(data_drive_with_tax_report, "Table_4_.txt", sep="")
my_data <- read.delim(fp_counts, header=TRUE, sep="\t")

fp_GO_term_list <- paste(data_drive_with_tax_report, "tumor_table.txt", sep="")
GO_terms <- read.delim(fp_GO_term_list, header=TRUE, sep="\t")

GO_IDs <- GO_terms[,1]
#GO_IDs <- c("GO:0002845")#, "GO:0002859")
#GO_IDs <- c("GO:0002651", "GO:0002663", "GO:0002666")
GO_IDs <- c("GO:0002230", "GO:0050689") 
#GO_IDs <- c( "GO:0039532", "GO:0140374")
            #), "GO:0071346")
#GO_IDs <- c("GO:0045630", "GO:2000553", "GO:0002830","GO:0002829","GO:0001812", "GO:0002827", "GO:0002826", "GO:2000321", "GO:0045627")
#"GO:0001868", "GO:0001867" elctin pathway
GO_IDs <- str_sort(GO_IDs)
species <- c("Cc", "Cm", "Dc", "Ei")
percent_sum <- vector(mode="numeric", length=4)
names(percent_sum) <- species

percents <- my_data[, c(7,8,9,10)]
row.names(percents) <- my_data[,2]

to_subset <- is.element(my_data$GO.term.ID, GO_IDs)
subsetted <- my_data[to_subset, c(2,7,8,9,10)]


for(i in 1:length(GO_IDs)){
  ID <- GO_IDs[i]
  row <- percents[ID,]
  
  for(s in species){
    percent_sum[s] <- percent_sum[s] + as.numeric(sub("%","",row[s]))/100
  }
}



percent_sum <- percent_sum *100


add <- max(percent_sum)*.15
# Simple Bar Plot
coul <- brewer.pal(5, "Set2") 
par(mar=c(4,5,7,1)+.1)
BP <- barplot(percent_sum, main="Percent of Immune Transcripts Mapping to GO:0050900\n",
        xlab="Marine Turtle Species",
        ylab="Percent of Species'\nImmune Transcripts Mapped",
        ylim = c(0,max(percent_sum) + add),
        col=coul 
        )
text(x = BP, y = percent_sum + add, label = percent_sum, pos = 1, cex = 0.8, col = "black")



num_IDs <- length(GO_IDs)
specie <- c(rep("Cc" , num_IDs) , rep("Cm" , num_IDs) , rep("Dc" , num_IDs) , rep("Ei" , num_IDs) )
condition <- rep(subsetted[,1] , 4) #4 species

convert_to_num <- function(percent){
  return(as.numeric(sub("%","",percent)))
}

### Visualize ###
Cc_p <- sapply(subsetted["Cc"],convert_to_num)
Cm_p <- sapply(subsetted["Cm"],convert_to_num)
Dc_p <- sapply(subsetted["Dc"],convert_to_num)
Ei_p <- sapply(subsetted["Ei"],convert_to_num)

#as.numeric(sub("%","",row[s]))/100
# condition <- rep(c("Negative Regulation", "Positive Regulation"), 4)
# num_IDs <- 2
# specie <- c(rep("Cc" , num_IDs) , rep("Cm" , num_IDs) , rep("Dc" , num_IDs) , rep("Ei" , num_IDs) )
# value <- c(0.93,0,0.66,0.1,0.70,0.02, 0.60, 0.05)

value <- c(Cc_p, Cm_p, Dc_p, Ei_p )
data <- data.frame(specie,condition,value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity") + 
  ggtitle("Regulation of Defense Response to Virus GO Term Percentages") + 
  xlab("Marine Turtle Species") + ylab("Percent") + 
  labs(fill = "GO IDs") + scale_fill_manual(values=c("palegreen3", "indianred2"))
# + scale_fill_manual(values=c("indianred2","palegreen3", "brown2","mediumseagreen"))#, "brown", "springgreen4" ))
