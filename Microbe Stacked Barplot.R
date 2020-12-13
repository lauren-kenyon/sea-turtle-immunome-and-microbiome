library(ggplot2)
library(stringr)
### Get Data ###

data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\"
table <- "taxa_counts.txt"
fp <- paste(data_drive_with_tax_report,table, sep="")


my_data <- read.delim(fp, header=TRUE, sep="\t")
rownames(my_data) <- my_data$tax_names

### Specify Microbes ###

microbes <- my_data[,1]# c("Acinetobacter baumannii", "Klebsiella pneumoniae")
#microbes <- c("Avian leukosis virus", "Bovine viral diarrhea virus", "Avian carcinoma Mill Hill virus 2","Escherichia virus phiX174", "Pestivirus", "Moloney murine sarcoma virus", "Y73 sarcoma virus","Escherichia virus phiX174")
m <- c("Rhodoferax koreense")
microbes <-  str_sort(m)

num_microbes <- length(microbes)
specie <- c(rep("Cc" , num_microbes) , rep("Cm" , num_microbes) , rep("Dc" , num_microbes) , rep("Ei" , num_microbes) )
condition <- rep(microbes , 4) #4 species

### Get Microbe Counts ###
microbe_counts <- vector(mode="list", length=4)
names(microbe_counts) <- c("Cc", "Cm", "Dc", "Ei")

#col_num <- 4 #3 for counts, 4 for percent
col_num <- 4
for(i in 1:4){
  ls <- c()
  for(m in 1:length(microbes)){
    #ls <- append(ls, my_data[microbe, col_num])
    ls <- append(ls,as.numeric(sub("%","",my_data[microbes[m], col_num])))
  }
  microbe_counts[[i]] <- ls
  col_num <- col_num + 2#2
}


### Visualize ###
value <- c(microbe_counts[[1]],microbe_counts[[2]], microbe_counts[[3]], microbe_counts[[4]] )
data <- data.frame(specie,condition,value)

stack_mode <- !TRUE
legend <- "Species"
title <- "Percent of Microbial Hits to E. coli vs. E.coli's Phage "
y_lab <- "Percent of Microbial Hits"

if(!stack_mode){
  title <- paste("Relative ",title, sep="")
}

if(stack_mode){
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity") + 
  ggtitle(title) + 
  xlab("Marine Turtle Species") + ylab(y_lab)+ labs(fill = legend) 
}else{
  ggplot(data, aes(fill=condition, y=value, x=specie)) + 
    geom_bar(position="fill", stat="identity") + 
    ggtitle(title) + 
    xlab("Marine Turtle Species") + ylab(y_lab)+ labs(fill = legend)  + scale_y_continuous(labels = scales::percent_format(suffix=""))
}


