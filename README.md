# Sea turtle immunome analysis

## Pipeline
- Outline to create annotation in Trinotate Pipeline.md
## Scripts
- Subset annotation for immune-GO transcripts
  - Create immunome by removing transcripts without atleast one immune GO term association
  
# Sea turtle microbial analysis
## Pipeline
- Outlined in Marine Turtle Microbial Contamination Pipeline.md

### Microbial analysis on each species' transriptome

#### Scripts
- Microbial Analysis by TS [Scripts 2-4]
- Step 6 of pipeline uses: Microbe Stacked Barplot.R

#### Results
- Microbial Mappings to Assemblies [with/without lineage].xlsx
  - Microbial (archaea, bacteria, or virus) mapping rates (with over 95% sequence identity) from each assembly 
  
### Microbial analysis by tumor score
- Microbial Analysis by TS [Scripts 1-5]
  - Microbial analysis by tumor score (TS) involves 5 scripts to take the unfiltered Cm assembly, filter it to contain only microbrial mappings, then use Salmon to quanitfy individual 
  green turtle mappings in transcripts per million (TPM)
- MDS_or_Kmeans 
  - Visualize results of Microbial Analysis by TS between individuals using an MDS plot


