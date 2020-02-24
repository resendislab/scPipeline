##### This file is the flow control of all pipeline
##### see pipeline.md for more information.

########### Needed Folders
### it Creates the folder to store numerical results.
  if (file.exists("./results/data") == FALSE){
    dir.create("./results/data", recursive = TRUE)
  }

### it Creates the folder to store plots.
if (file.exists("./results/plots") == FALSE){
  dir.create("./results/plots", recursive = TRUE)
}

### it Creates the folder to store processed data.
  if (file.exists("./data/processed") == FALSE){
    dir.create("./data/processed", recursive = TRUE)
  }

########### Sample Processing
### Raw data should be stored at ./data/raw/
### There is a count matrix for every line from every sample

########### Knee Filter
### -i ddseeker output files
### -o Filtered count matrix
  source("./Src/Data_Conversion.R")


########### Gene Over-dispersion
### -i Filtered count matrix
### -o Over-dispersion matrix
  source("./Src/Gene_Dispersion.R")

########### Gene number selection
### -i Over-dispersion matrix
### -o Variance plots to select the number of considered genes
  source("./analysis/Variance.R")

########### Clustering
## NOTE: Dimensionality reduction step is embedded in the clustering script
######  Kmeans
### -i Over-dispersion matrix
### -o optimal number of groups
###    uMAP representation of clustering
###    Data re-grouped according to kmeans clustering algorithm
###    proportions plots
  source("./analysis/Clustering_kmeans_umap.R")

######  Expectation-Maximization
### -i Over-dispersion matrix
### -o optimal number of groups
###    uMAP representation of clustering
###    Data re-grouped accoding to EM clustering algorithm
###    proportions plots
  source("./analysis/Clustering_EM_umap.R")

########### Data re-assignation
### -i raw count matrix
### -o Count matrix rearranged given clustering methods
  source("./analysis/Grouping.R")

########## Enrichment Analysis
## NOTE: As a output of the previous step, one matrix file
##       is generated that contains information of all clusters.
##       the matrix is needed to use GSEA tool

########### Differential Gene expression analysis
## NOTE: Grouping.R script also saves the matrix for comparisons
##       by pairs needed for DGEA
### -i Count matrix re-arranged by comparison pairs
### -o Data file with log2 fold change and the statistical
###    parameters associated for each comparison.
  source("./analysis/Diferential_Gene_Expression.R")

######  Genetic Signatures
  ## p-values
### -i Data from DGEA
### -o Genes that over the threshold of log2 fold change >= 2.0 and a
###    p-value >= 0.01
  source("./analysis/p_values genes.R")
 
  ## Gene list intersection
### -i Filtered genes by log2 fold change and p-values for
###    a cluster given all possible comparisons
### -o Global differentially expressed genes for every cluster
  source("./analysis/Gene_list_Intersection.R")
 
 ## Heatmap
### -i genetic signatures for every cluster
###    raw count matrix  
### -o heatmap plot of global differentially expressed genes
  source("./analysis/Gene_expression_heatmap.R")
