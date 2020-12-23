# Sea turtle immunome analysis

## Pipeline
- Outline to create annotation in Trinotate Pipeline.md
## Scripts
- Subset annotation for immune-GO transcripts
  - Create immunome by removing transcripts without atleast one immune GO term association
- Annotation and immunome stats_updated.R
  - Get summary of annotation/immunome findings
- GO_slim_figure_and_counts.R
  - Compare presence of GO terms by category and make visual bar graph
- Get GO list and table from parent(s).R
  - Get list of GO terms under a parent GO term and a table with their ID and description
- Get Percent Sum between Species.R
  - Create visual bar graph showing percent of immune transcripts associated with GO term(s)
- Isolate genes per GO term.R
  - Find annotated genes from immunome based on desired GO term
- Isolate list of immunome GO term.R
  - Get list of all immunome GO terms found in annotation
 - UpsetR with GO tables.R
  - Create upsetR plot based on shared GO terms between species and create a table of GO term IDs/descriptions to see what the bins in the upsetR plot hold
## Results
- All Immune GO term Counts and Percentages.xls
  - Species immunomes compared by GO term with relative counts and percentages for comparison
 
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

# Important links/hackmds

Overall project hackmd: https://hackmd.io/hkPu-fbuRMu8zp7v9NhQ-A?edit

Trinotate pipeline: https://hackmd.io/v5CYlWM-RBaA0A5TlySsjg

Creating immunome: https://hackmd.io/76f_GZKnQ6Gl9uAslL5xzA?both


## Microbial-related pipelines

Microbial characterization of transcriptomes pipeline: https://hackmd.io/SH-2SWL1RmWvE0JXHJ8igg

Microbial analysis by tumor score: https://hackmd.io/ATlxYyZSSvewGez9wwPw6g
