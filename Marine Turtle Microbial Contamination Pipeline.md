# Marine Turtle Microbial Contamination Pipeline
Literature Review in box: Sea turtle Immunomes -Kenyon Honors Thesis/Thesis Files

Scripts in box and below:  Sea turtle Immunomes -Kenyon Honors Thesis/Scripts/Microbial Contamination 

Questions:
-	What pathogenic bacteria are present in the blood?
    -	Group hits by taxon [DONE]
-	Is there an association between microbial communities and tumor score?
    -	Map individuals to their designated transcriptome and get a shared list of microbial transcripts (instead of making new transcriptomes for individual GT)
    -	Group by tumor score
-	Is there an association between microbial communities and sample locations?
    -	Compare taxon % between groups of turtles from different sampling locations

Addressing:
- Underlying blood infection?
- Microbial communities & disease
- Environmental niches


Cmd line diamond params:
https://github.com/bbuchfink/diamond_docs/blob/master/3%20Command%20line%20options.MD


## 1. Download non-redundant database and supporting files:
NR:
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

Mapping file:
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
- gunzip

Nodes.dmp (required for taxonomy features):
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip.
- unzip for names.dmp and nodes.dmp

## 2. Make diamond database
- This will create a binary DIAMOND database file with the specified name (nr.dmnd).
```
#BSUB -q long
#BSUB -W 20:00
#BSUB -R rusage[mem=5000]
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -e make_nr.err
#BSUB -oo make_nr.log

module load diamond/2.0.4
diamond  makedb --in nr.gz -d nr --taxonmap ./prot.accession2taxid --taxonnodes ./nodes.dmp --taxonnames ./names.dmp

```
## 3. Transdecoder to convert .fasta to .pep
    - extract candidate LORF
    - input for blastp
```
#!/bin/bash
#BSUB -q long
#BSUB -W 100:00
#BSUB -n 20
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6400]"
#BSUB -oo GT_to_pep.log
#BSUB -e GT_to_pep.err


module load transdecoder/5.5.0

assembly_dir=/project/uma_lisa_komoroske/Full_Turtle_RNAseq_Spring2019/green_turtle_assembly/Trinity_miniGT_assembly/Trinity_HIGT_small_assembly
cd $assembly_dir
TransDecoder.LongOrfs -t mini_chmy_Trinity.fasta

```


## 4. Blastp for top hit
```
#!/bin/bash
#use diamond to search uniprot database with input protein sequences

#BSUB -q short
#BSUB -W 3:00
#BSUB -R rusage[mem=6000]
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -e GT_nr.err
#BSUB -o GT_nr.log


module load diamond/2.0.4
pep=/project/uma_lisa_komoroske/Full_Turtle_RNAseq_Spring2019/green_turtle_assembly/Trinity_miniGT_assembly/Trinity_HIGT_small_assembly/mini_chmy_Trinity.fasta.transdecoder_dir/longest_orfs.pep


#need a downloaded protein dabase file in fasta format AND a file of translated DNA reads

#Just getting taxon data + regular output
diamond blastp -d ./nr.dmnd -q $pep -o GT_nr_tophits.m8 --max-target-seqs 1 --max-hsps 1 --taxonlist 2,2157,10239 -f 6 staxids sskingdoms skingdoms sphylums pident qseqid sseqid length mismatch gapopen qstart qend sstart send evalue bitscore


```
### Output:

![](https://i.imgur.com/1vg8znw.png)

![](https://i.imgur.com/3tzfUBI.png)


### Quick stats:
![](https://i.imgur.com/iRRbdJ0.png)



## 5. Count Species Hits
a. Isolate taxon IDs above a percent identity
Scripts: [located in Scripts/Microbial Contamination-related Scripts]
1. 1_Get_TaxonID_Transcript_Map.R
2. 2_Isolate_Microbial_Tax_IDs.R
3. 3_Update_Mapping_Remove_Nonmicrobials.R


b. Map tax ID to organism with 
https://www.ncbi.nlm.nih.gov/taxonomy > Name/ID Status
direct link: https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi 

1. Upload tax id output from script above
a. Please enter a list of names or tax_id(s): > Click "Choose file" > select output
2. Select full taxid lineage
3. Click save in file
a. tax_report.txt will download
b. rename to [species]_95_lineages.txt

![](https://i.imgur.com/GvkaVaI.png)


c. Use tax_report to get counts of microbe hits

Note: staxids will report multiple taxids when the gene maps to other species (even if you specified a particular taxon list)

In the script below, I remove extra taxon ids that aren't in the bacteria, archae, or virus lineage


[script below]
```
library(stringr)


### Get taxon report ###
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\"

#Saved from NCBI with lineage information
tax_report_name <- "Ei_95_lineages"
fp2 <- paste(data_drive_with_tax_report, tax_report_name, ".txt", sep = "")
taxon <- read.delim(fp2, header=TRUE, sep="\t")

taxon <- taxon[c(TRUE,FALSE)] #Get rid of even "|" columns

#remove last 2 extran. lines
length <- nrow(taxon)
taxon <- taxon[-c(length-1, length),]

  
is_microbe <- function (lineage){
  my_tax_categories <- c("2","2157","10239")
  sep_lin <- strsplit(lineage, "\\s+")
  microbe <- FALSE
  
  for(id in my_tax_categories){
    if(is.element(id, sep_lin[[1]])){
      microbe <- TRUE
      break
    }
  }
  return(microbe)
}


lineages = taxon$lineage
bool <- sapply(lineages,is_microbe)

microbe_taxa <- taxon[bool,]




taxnames <- microbe_taxa$taxname
all_names <- unique(taxnames)
taxids <- unique(microbe_taxa$taxid)


tax_counts <- vector(mode="numeric", length=length(all_names))
names(tax_counts) <- all_names

### Count taxonomy occurrences  ###
for(i in 1:nrow(microbe_taxa)){
  row <- microbe_taxa [i,]
  tax_name <- row[["taxname"]]
  
  tax_counts[tax_name] <- tax_counts[tax_name] + 1
}

### Write to excel ###
fp_1 <- paste(data_drive_with_tax_report, tax_report_name, "_counts.xls", sep = "")
ls <- list(tax_counts)
percent_of_blasts <- ls[[1]]/nrow(microbe_taxa)

taxon.table <- data.frame(
  tax_names = names(ls[[1]]),
  tax_ids = taxids,
  counts = ls[[1]],
  percent = percent_of_blasts
)

write.table(taxon.table, fp_1, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t")

```

d. Compare species
```
### Import Data ###


tax_report_name <- "_95_lineages"
data_drive_with_tax_report <- "C:\\Users\\laure\\Documents\\Thesis\\Microbe\\"

GT <- read.delim(paste(data_drive_with_tax_report, "GT", tax_report_name, "_counts.xls", sep = ""), header=TRUE, sep="\t")
rownames(GT) <- GT$tax_names

Dc <- read.delim(paste(data_drive_with_tax_report, "Dc", tax_report_name, "_counts.xls", sep = ""), header=TRUE, sep="\t")
rownames(Dc) <- Dc$tax_names

Cc <- read.delim(paste(data_drive_with_tax_report, "Cc", tax_report_name, "_counts.xls", sep = ""), header=TRUE, sep="\t")
rownames(Cc) <- Cc$tax_names

Ei <- read.delim(paste(data_drive_with_tax_report, "Ei", tax_report_name, "_counts.xls", sep = ""), header=TRUE, sep="\t")
rownames(Ei) <- Ei$tax_names

#List of all microbes
full_names <- unique(c(GT$tax_names, Dc$tax_names, Cc$tax_names, Ei$tax_names))
full_taxids <- unique(c(GT$tax_ids, Dc$tax_ids, Cc$tax_ids, Ei$tax_ids))

### Get counts between species ###

species <- list(Cc, GT, Dc, Ei)
counts_and_percents <- vector(mode="list", length=length(species)*2)
index <- 1

for(i in 1:length(species)){
  s <- species[[i]]
  counts <- vector(mode="numeric", length=length(full_names))
  percents <- vector(mode="numeric", length=length(full_names))
  
  for(n in 1:length(full_names)){
    taxname <- full_names[n]
    
    if(is.na(s[taxname,][1])){ #If taxname is not in this species' list
      counts[n] <- 0
      percents[n] <- 0.0
    }else{
      counts[n] <- s[taxname,]$counts
      percents[n] <- s[taxname,]$percent
    }
  }
  counts_and_percents[[index]] <- counts
  index <- index +1
  counts_and_percents[[index]] <- percents
  index <- index +1
}


### Construct output ###
taxon.table <- data.frame(
  tax_names = full_names,
  tax_ids = full_taxids,
  #Cc_counts = counts_and_percents[[1]],
  Cc_percent = counts_and_percents[[2]],
  #Cm_counts = counts_and_percents[[3]],
  Cm_percent = counts_and_percents[[4]],
  #Dc_counts = counts_and_percents[[5]],
  Dc_percent = counts_and_percents[[6]],
  #Ei_counts = counts_and_percents[[7]],
  Ei_percent = counts_and_percents[[8]]
)


fp_r <- paste(data_drive_with_tax_report, "all", tax_report_name, "_counts.xls", sep = "")
write.table(taxon.table, fp_r, row.names = FALSE, col.names <- FALSE, quote = FALSE, sep = "\t")
```

Then:
1. Open in Excel
2. Right click in percentage => Sort => Sort largest to smallest






## 6. Visualize

with
/Scripts/Microbial Contamination-related Scripts/Microbe Stacked Barplot.R

![](https://i.imgur.com/sKvTQDp.png)

