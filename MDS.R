
library(GO.db)
library(biomaRt)
library(stringr)
immune_genes <- as.list(GOBPOFFSPRING[["GO:0002376"]])

my_GO_immune_matrix <- matrix(nrow=4, ncol=length(immune_genes))
my_immunomes <- c("Cc", "GT", "Dc", "Ei")
colnames(my_GO_immune_matrix) = immune_genes
rownames(my_GO_immune_matrix) = my_immunomes



#Filepath for GO list folder
data_drive <- "C:\\Users\\laure\\Documents\\Thesis\\Annotations\\GO_lists\\"


#Ex Go list name:" GO_list_of_immunome_GT_trinotate_annotation_report.xls" #found in data_drive
immunome_annot_name <- "_trinotate_annotation_report.xls"

num_immunomes <- length(my_immunomes)
#Comparing GO term prescense for all immune GO terms
for (r in 1:num_immunomes){
  GO_fp_for_list <- paste(data_drive, "GO_list_of_immunome_", my_immunomes[r], immunome_annot_name, sep = "")
  
  mydata <- read.delim(GO_fp_for_list , header=FALSE, sep="\t")
  GO_terms <- mydata[,1]
  
  for(c in 1:length(immune_genes)){
    if(is.element(immune_genes[c], GO_terms)){
      my_GO_immune_matrix[r,c] = 1 #If immunome has particular GO term -> gets 1
    }else{
      my_GO_immune_matrix[r,c] = 0 #Otherwise -> matrix gets 0
    }

  }
}

#(aka asymmetric binary): The vectors are regarded as binary bits, so non-zero elements are 'on' and zero elements are 'off'. The distance is the proportion of bits in which only one is on amongst those in which at least one is on.
d <- dist(my_GO_immune_matrix, method="binary")


mds1 = cmdscale(d, k = 2)

# plot
plot(mds1[,1], mds1[,2], type = "n", xlab = "", ylab = "", axes = FALSE,
     main = "MDS Plot for Immune GO Term Prescence between Immunomes")
text(mds1[,1], mds1[,2], c("C.carretta","C. mydas","D. coriacea","E. imbricata"), cex = 1, xpd = TRUE, col = c("darkgreen"))




