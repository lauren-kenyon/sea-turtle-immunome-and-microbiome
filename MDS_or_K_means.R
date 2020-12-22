
library(biomaRt)
library(stringr)
library(magrittr)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(fpc)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

### Adjust parameters for graph output ###
plot_MDS <- !FALSE

plot_kmeans <- !FALSE
num_clusters <- 9
kmeans_elbow_test <- TRUE #edit kmeans number of clusters to match result
kmeans_silhoutette_test <- TRUE

color_code_by_tumor <- TRUE #False is by location

by_region <- !TRUE        #collapse location by overall region e.g. California
where <- "topleft"       #where the legend should be
general_title <- "of Enterobacter hormaechei and\nRhodoferax spp. Microbial Mappings"
#,\nand Rhodoferax sp. OTU1

MDS_title <- paste("MDS Plot of ", general_title, sep="" )
k_title <- paste ("K-means Clustering of MDS Plot of\n", general_title, sep="" )

##########################################



### Get data 
data_drive <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\GT_quants\\"
### Get quant mappings ###

fp_quant <- paste(data_drive, "quant_mappings.xls", sep = "") 
#fp_quant <- paste(data_drive, "viral_quant_mappings.txt", sep = "") 
quants <- read.delim(fp_quant, header=TRUE, sep="\t")

### Get metadata
fp_chmy_meta <- paste(data_drive, "chmy_metadata_JAS_7.29.19.csv", sep = "") 
meta <- read.delim(fp_chmy_meta, header=TRUE, sep=",")
 
row.names(meta) <- meta$samp_ID
sorted_meta <- meta[str_sort(meta$samp_ID) , ] #Sort by sample number as the quant file is sorted by sample #
samps <- sorted_meta$samp_ID

tumor_sc <- sorted_meta$Tumor.Score
if(color_code_by_tumor){
  location <- sorted_meta$Sampling_Location
  location <- sapply(location, trimws)
  Hawaii <- c("Kailua","Kiholo","Kapoho")
  loc_bool <- is.element(location, Hawaii)
  
  tumor_sc <- tumor_sc[loc_bool]
  
  #Color code by tumor score
  colors <- c("black", rgb(255,165,38, maxColorValue = 255), rgb(255,114,0, maxColorValue = 255), rgb(204,0,0, maxColorValue = 255))
  names(colors) <- c("0","1","2", "3") #tumor scores
  
  pt_colors <- vector(mode="character", length=length(tumor_sc)) #assign each pt a color
  
  for(i in 1:length(tumor_sc)){
    t <- toString(tumor_sc[i])
    pt_colors[i] <- colors[t] #choose color based on tumor score
  }
}else{
  #Color code by sampling location
  location <- sorted_meta$Sampling_Location
  location <- sapply(location, trimws)
  colors <- c()
  
  if(!by_region){
    colors <- c("darkolivegreen3", "green3", "chartreuse4","#A6CEE3", "#1F78B4", "#E31A1C", "#FDBF6F", "#FF7F00" )
  }else{
    colors <- c("forestgreen", "forestgreen", "forestgreen", "blue", "blue", "red", "red", "red" )
  }
  
  names(colors) <- c("Kailua","Kiholo","Kapoho","CNMI-TI", "CNMI-SA", "SBNWR","SGR", "SDB")
  pt_colors <- vector(mode="character", length=length(location))
  
  for(i in 1:length(location)){
    t <- trimws(location[i])
    pt_colors[i] <- colors[t]
  }
}


#Looking at specific microbes of interest
keep_microbes <- c("Enterobacter hormaechei","Rhodoferax sp. BAB1", "Rhodoferax koreense","Rhodoferax sp. OTU1")


keep <- vector(mode="logical", length=nrow(quants))
for(i in 1:length(quant_names)){
  keep[i] <- is.element(quant_names[i], keep_microbes)
}
g <- quants$tax_names[keep]
TPMS <- quants[keep,]



TPMS <- TPMS[,c(FALSE,TRUE)] #for TPM
#TPMS <- quants[,c(TRUE,FALSE)] #for counts
TPMS <- TPMS[,-1] #remove taxid column

#TPMS <- TPMS[-c(147:177),] #none in these rows if looking at all microbe comparasion
microbe_hits_genes <- quants$tax_names

x <- 44
cnames <- vector(mode="character", length=x)

for(i in 1:x){
  n <- ""
  if(i < 10){
    n <- "0"
  }
  cnames[i] <- paste("samp_", n, i, sep = "")
}

 #Name quant mappings correctly

colnames(TPMS) <- cnames
TPMS_subsetted <- TPMS[samps] #Remove quants that aren't in metadata excel


if(color_code_by_tumor){
  TPMS_subsetted <- TPMS_subsetted[loc_bool]
}

samples <- t(TPMS_subsetted) #Change orientation of data (flip rows and columns)
  
d <- dist(samples, method="euclidean")

fit <- cmdscale(d, eig=TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]



if(plot_MDS){
  
  plot(x, y, pch = 19, xlim = range(x), col = pt_colors, cex = 1.1)
  title(main= MDS_title)
  #text(x,y, cex = 0.5,labels=rownames(samples))
  if(color_code_by_tumor){
    legend(where,
           legend = c("0 - no visible tumors", "1", "2", "3 - heavily tumored"),
           title = "Tumor Score",
           col = colors,
           pch = c(19),
           bty = "n",
           pt.cex = 1.2,
           cex = 0.8,
           text.col = "black",
           horiz = F ,
           inset = c(0.1, 0.1))
  }else{
    if(by_region){
      colors <- unique(colors)
      names(colors) <- c("Hawaii", "Mariana Island", "California")
      
    }
    legend(where, 
           legend = names(colors), 
           title = "Sampling Populations",
           col = colors, 
           pch = c(19), 
           bty = "n", 
           pt.cex = 1.2, 
           cex = 0.8, 
           text.col = "black", 
           horiz = F , 
           inset = c(0.1, 0.1))
  }

}

if(kmeans_silhoutette_test){
  fviz_nbclust(samples, kmeans, method = "silhouette")
}
if(kmeans_elbow_test){
  # # # Determine number of clusters: 4
  wss <- (nrow(samples)-1)*sum(apply(samples,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(samples,
                                       centers=i)$withinss)
  plot(1:15, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}

if(plot_kmeans){
  mds <- samples %>%
    dist() %>%
    cmdscale() %>%
    as_tibble()
  colnames(mds) <- c("Dim.1", "Dim.2")
  # Plot MDS
  ggscatter(mds, x = "Dim.1", y = "Dim.2",
            label = rownames(samples),
            size = 1,
            repel = TRUE)
  
  # K-means clustering
  set.seed(123) #same results with and without
  clust <- kmeans(mds, num_clusters ,nstart=2500)$cluster %>%
    as.factor()
  mds <- mds %>%
    mutate(groups = clust)
  # Plot and color by groups
  ggscatter(mds, x = "Dim.1", y = "Dim.2", xlim = range(x) + c(0, 0),
            label = tumor_sc,
            #label = sorted_meta[rownames(samples),8],
            title = k_title,
            color = "groups",
            palette = "jco",
            size = 1.7,
            fill = 0.1,
            ellipse = TRUE,
            ellipse.type = "convex",
            repel = TRUE)
}






