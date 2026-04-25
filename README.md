# H3K9me2-Schizophrenia-Analysis
Investigating whether differential H3K9me2 modifications correlate with the expression of nearby genes in postmortem tissue samples of individuals diagnosed with Schizophrenia.
## Overview
Despite decades of research, common treatments for Schizophrenia remain largely unsuccesfull, pushing researchers to investigae the disorder's genetic & epigenetic basis. Schizophrenia (SCZ) is a severe, chronic psychiatric disorder dirven in part by epigenetic dysregulation in the brain. H3K9me2, a repressive histone marker, is understudied yet suspected to alter gene expression in SCZ. This project integrates ChIP-se & RNA-seq datasets to postmortem prefrontal cortex tissue to determine if differential H3K9me2 modifications correlate with the expression of nearby genes. <br> 

**Key findings:**
- H3K9me2 is negatively correlated with gene expression, consistent with its role as a repressive mark
- SCZ samples show significantly more H3K9me2 markers per gene than controls, suggesting elevated repression
- However, direction of H3K9me2 change at individual genes does not predict direction of expression change
- This suggests H3K9me2 operates at the chromatin domain level rather than as a precise per-gene switch

## Data 
| Dataset | Source | Samples | Region |
|---------|--------|---------|--------|
| ChIP-seq (H3K9me2) | [GSE215991](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215991) / [PRJNA891668](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA891668) | 15 SCZ, 15 controls | PFC, Brodmann area 9 |
| RNA-seq | [PRJNA695206](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA695206) | 19 SCZ, 18 controls | PFC, Brodmann area 10 |

## Pipeline Overview
**Preprocessing**
- Quality control (FastQC, MultiQC)
- Read trimming (Trimmomatic)
- Alignment to hg38 (Bowtie2)
- BAM to BED conversion

**ChIP-seq analysis**
- Differential peak calling (DiffReps) <br>
- Volcano plot of differential H3K9me2 peaks <br>

**RNA-seq analysis**
- Read quantification (featureCounts) <br>
- TMM normalization to log2 CPM (edgeR) <br>
- Differential expression analysis (LIMMA) <br>

**Integration**
- Map ChIP-seq peaks to nearest genes <br>
- Define regulatory windows around transcription start sites (1kb, 5kb, 10kb upstream) <br>
- Identify peak-associated differentially expressed genes <br>

**Downstream analysis**
- Pearson correlation between H3K9me2 signal and gene expression <br>
- Per-sample correlation comparison (SCZ vs control) <br>
- Gene Ontology / KEGG enrichment (clusterProfiler) <br>
- Co-modification network analysis <br>
## Repository Structure 

H3K9me2-Schizophrenia-Analysis/
```
├── README.md  
├── .gitignore  
├── scripts/  
│   ├── preprocessing/      # Trimming, alignment, QC  
│   ├── chip-seq/           # Peak calling, volcano plots  
│   ├── rna-seq/            # Quantification, normalization, DEG  
│   ├── integration/        # Peak-gene mapping, regulatory windows  
│   ├── analysis/           # Correlation, scatter plots, box plots  
│   └── network/            # Co-modification networks  
├── annotations/            # Reference annotation files  
├── poster/                 # Project poster  
```
## Dependancies 

**Command-line tools**
- Trimmomatic
- Bowtie2
- SAMtools
- BEDTools
- DiffReps
- Subread (featureCounts)
- MultiQC

**R packages**
- edgeR
- LIMMA
- clusterProfiler
- org.Hs.eg.db
- enrichplot
- tidyverse
- ggpubr
- ggplot2

**Python packages**
- pandas
- pybedtools
- networkx
- matplotlib

**Reference genome**
- hg38 (GENCODE v38 annotation)

## References
1. Frontiers in Psychiatry (2022). *Differential H3K9me2 heterochromatin levels and concordant mRNA expression in postmortem brain tissue of individuals with schizophrenia, bipolar, and controls.* [Link](https://www.frontiersin.org/journals/psychiatry/articles/10.3389/fpsyt.2022.1006109/full)
2. Schizophrenia Bulletin Open (2021). *Concordance of Immune-Related Markers in Lymphocytes and Prefrontal Cortex in Schizophrenia.* [Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC7865130/)

## Author
Isabella Castellano <br>
Mentor: Xander Krohannon, MS
