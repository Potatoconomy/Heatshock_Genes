# Heatshock_Genes
Returns Heat Shock Response and Molecular Chaperone genes from a DESeq2 Results Table.

# Usage

Simply load in your log-fold change shrinkage (LFC) results from DESeq2 (or other DE analysis). Best is if your rownames have been converted to HGNC symbols, otherwise basic modifications to the code need to be made. Feedback is greatly appreciated. Keep an eye out for false results by checking gene descriptions in the dataframe output.

2 User inputs are required. They are clearly notated in the code.

# Version

### Dependencies
 
DESeq2_1.26.0 
hash_2.2.6.1  
conflicted_1.0.4.9000

### My R Session

R version 3.6.1 (2019-07-05)  
Platform: x86_64-conda_cos6-linux-gnu (64-bit)  
Running under: Ubuntu 20.04.1 LTS  
