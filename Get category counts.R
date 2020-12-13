library(GO.db)
library(biomaRt)

### Get GO term list ###
data_drive_with_GO <- "C:\\Users\\laure\\Documents\\Thesis\\"
GO_list_name <- "core_GO_term_set.xls" #found in data_drive
fp <- paste(data_drive_with_GO, GO_list_name, sep = "")
GO_table <- read.csv(fp, header=TRUE, sep="\t")
GO_list <- GO_table[[1]]



### Define GO categories to slim ###
all_immune_genes <- as.list(GOBPOFFSPRING[["GO:0002376"]])

immune_categories <- c("GO:0002684", "GO:0002683","GO:0002418", "GO:0002507", 
                       "GO:0002437" , "GO:0045087","GO:0002250", "GO:0045321", "GO:0051607", "GO:0045917", "GO:0045916")
names(immune_categories) <- c("positive regulation of\nimmune response", "negative regulation of\nimmune response", "immune response\nto tumor cell", "tolerance\ninduction", 
                              "inflammatory response\nto antigenic stimulus", "innate\nimmune response", "adaptive\nimmune response", "leukocyte\nactivation", "antiviral\nresponse", 
                              "positive regulation of\ncomplement activation", "negative regulation of\ncomplement activation")

immune_offspring <- vector("list", length = length(immune_categories))

#for(i in 1:length(immune_categories)){
#  immune_offspring[i] <- as.list(GOBPOFFSPRING[[immune_categories[i]]])
#}

pos_imm_ls <- as.list(GOBPOFFSPRING[[immune_categories[1]]])
neg_imm_ls <- as.list(GOBPOFFSPRING[[immune_categories[2]]])
tumor_ls <- as.list(GOBPOFFSPRING[[immune_categories[3]]])
tolerance_ls <- as.list(GOBPOFFSPRING[[immune_categories[4]]])
inflam_ls <- as.list(GOBPOFFSPRING[[immune_categories[5]]])
innate_ls <- as.list(GOBPOFFSPRING[[immune_categories[6]]])
adapt_ls <- as.list(GOBPOFFSPRING[[immune_categories[7]]])
leuko_ls <- as.list(GOBPOFFSPRING[[immune_categories[8]]])
antiv_ls <- as.list(GOBPOFFSPRING[[immune_categories[9]]])
pcompl_ls <- as.list(GOBPOFFSPRING[[immune_categories[10]]])
ncompl_ls <- as.list(GOBPOFFSPRING[[immune_categories[11]]])

immune_offspring <- list(pos_imm_ls, neg_imm_ls, tumor_ls, tolerance_ls, inflam_ls, innate_ls,adapt_ls, leuko_ls, antiv_ls, pcompl_ls, ncompl_ls)


### Get category counts ###
count_vec <- vector(mode="numeric", length=length(immune_categories))

for(i in 1:length(GO_list)){
  GO_term <- GO_list[i]
  for(c in 1:length(immune_categories)){
    if(is.element(GO_term, immune_offspring[[c]])){
      count_vec[c] = count_vec[c]+1
    }
  }
}

my_labels <- names(immune_categories)
barplot(count_vec, main="GO slim categories and counts of GO terms", 
        ylab="Number of GO terms GO category", col=c("darkblue","red", "purple"),
        xlab="Immune GO categories",
        names.arg=my_labels, cex.names = 0.9)


count_vec[9]
