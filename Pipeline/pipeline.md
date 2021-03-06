<p align="center">
  <img src="https://user-images.githubusercontent.com/60892768/74993425-88d39900-5410-11ea-8643-b701551d0472.png">
</p>

# Single Cell Pipeline

The full pipeline starts with the fastq files and ends with the differentially expressed genes and the interaction maps according pathway enrichment. The **Pipeline.R** file contains the pipeline as R code, which it is structured by the following modules:

* [Sample Processing](#samples-processing)
  - [Data Availability](#data-availability)
* [Knee filter](#knee-filter)
* [Gene Over-dispersion](#gene-over-dispersion)
* [Gene number selection](#gene-number-selection)
* [Dimensionality Reduction](#Dimensionality-Reduction)
* [Clustering](#Clustering)
  - [Kmeans](#kmeans)
  - [EM](#expectation-maximization-algorithm)
* [Enrichment Analysis](#Enrichment-Analysis)
  - [Interaction Map](#Interaction-Map)
* [Differential Gene Expression Analysis](#Differential-Gene-Expression-Analysis)
  - [Genetic Signatures](#Genetic-signatures)

## Samples Processing

Raw files were processed through ddSeeker pipeline ([Rogmanoli et al. 2018](https://link.springer.com/epdf/10.1186/s12864-018-5249-x?author_access_token=5GkMGb1JtnR887KgJLyh-m_BpE1tBhCbnbw3BuzI2ROCtPkWTeG4740r3l1fSwcikEVJgHZm9jmSpiShpuk2FV8ae_KUm1O2Kb8nf4xLIHEJfbwJ3tasCADHjVdZ23iWPECc69GpWSYTzS6lSpBaJA%3D%3D)) using **STAR** as aligner. As output, BAM files were created for every sample. Count matrix were computed by **DigitalExpression** function of the [Drop-Seq tools](https://github.com/broadinstitute/Drop-seq) with a value of NUM_CORE_BARCODES=5000.

There is a count matrix for every of the 4 lines. All files can be merged into one.

### Data Availability
The data used and discussed have been deposited in NCBI's Gene Expression Omnibus and are accessible through the accession number [GSE145633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145633).

## Knee filter
Since we used droplet based technology for the high throughput sequencing, data must be filtered by taking out droplets with no RNA from the biological samples, so, only informative cells remain. To implement knee filter we plotted the cumulative UMI sum   for all cells. The first inflection point was set as threshold.  

<p align="center">
  <img src="https://user-images.githubusercontent.com/60892768/74989579-6177ce80-5406-11ea-8e20-866bc38cdefd.png">
</p>

As output file, the count matrix of cells leftward of the threshold is saved. All samples are saved into one big matrix

## Gene Over-dispersion
Instead of using all genes into the analysis, we selected the most informative ones. To do so, first taking hand of ![SCDE error models](https://hms-dbmi.github.io/scde/index.html) overdispersion matrix is computed of every gene and every sample.
The basic concept is if a gene is affected by a process or condition it will be reflected in their variance. A non affected gene will have the same value across data. An example is shown below.

<p align="center">
  <img src="https://user-images.githubusercontent.com/60892768/74991777-98e97980-540c-11ea-90ca-a90fe7ff38e1.png">
</p>

Overdispersion matrix for all data is saved in csv file

## Gene number selection
The number of genes used for further analysis has a relevant consideration. If only genes with greater variance were taken (e.g. top 500) gene expression analysis and pathway enrichment will give information about most significant processes. However, you will discard subtle processes with relevant information about the phenomena. Therefore, we proposed to discard genes with a variation less than 15% of the total dynamic range.

<p align="center">
  <img width="450" height="450" src="https://user-images.githubusercontent.com/60892768/74992683-7193ac00-540e-11ea-948d-91af7130f1b1.png">
</p>

## Dimensionality Reduction
Since data is in a multidimensional space it is hard to associate similar characteristics among data. Therefore, a dimension reduction is a proper analysis to visualize data properties in a bi-dimensional space. Particularly for our data set, uMAP method works better due its objective function definition (R package **umap**).

## Clustering
Instead of compared data according the sample labels, we mixed all data and performed clustering methods to re-group data by similarities in their gene expression profile. Two clustering methods were used to validate results and reduce possible induced bias by each one of them.

Clustering was performed in two spaces. Firstly, in the uMAP bidimensional space. and secondly, in the multidimensional space of the distance of every cell given the variance. This was done to evaluate the analysis robustness and set possible considerations.

### Kmeans
It is a supervised method based on distance of a centroid to the data, for kmeans implementation the number of cluster must be predefined. To select the optimal number of clusters, an screening was performed fot several number of groups, from 2 to 18. In each one, Sum of Squared estimate of Errors (SSE) was computed. So, the optimal clusters number was assigned when the elbow plot indicates the maximum change in SSE. This was done by computing the second derivative of SSE to find the optimal point. Kmeans method was implemented through the R package **kmeans**

<p align="center">
  <img width="450" height="450" src="https://user-images.githubusercontent.com/60892768/75059691-1eb80400-54a3-11ea-85b0-deb5973e2a44.png">
</p>

The projection of both inputs spaces clustered by kmeans are as follows. As can be seen, despite cluster assignation, uMAP input space has a more refined assignation. This consideration must be taken carefully,because multidimensional space clustering was not performed in the uMAP space but it is plotted in the uMAP space.

<p align="center">
  <img width="900" height="450" src="https://user-images.githubusercontent.com/60892768/75061267-fe3d7900-54a5-11ea-9651-d80937df830a.png">
</p>

### Expectation-Maximization Algorithm
Expectation-Maximization Algorithm is non supervised method, that computes the maximum likelihood of the probability distribution.

EM algorithm was implemented through the R package **mclust**. mclust gives you the optimal clusters number by computing the associated error. EM was used with the same two inputs as kmeans. An example is presented below, as can be seen EM fractures the space into more groups, however as further discussed, there is no difference between analysis. An important consideration is that based on the sample size, fewer groups are better because is easier to be associated with a particular process.

<p align="center">
  <img width="500" height="450" src="https://user-images.githubusercontent.com/60892768/75064692-0947d780-54ad-11ea-9542-d76210a13f57.png">
</p>

As a result of clustering methods, data is re-grouped with their respective cluster assignation in a csv file where columns sets for cluster name with the number of sample and rows sets for the selected genes.

## Enrichment Analysis
Since enrichment analysis performed by [GSEA tool](https://www.gsea-msigdb.org/gsea/index.jsp) can be done comparing every cluster with the rest of them, only one file is needed containing all clusters data. Moreover, there is a need to convert the count matrix into a specific formatted files according [GSEA file formats](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats).

### Interaction Map
As a result from gsea analysis, we got the enriched pathways according the expression profile. A filter is needed to keep those pathways with statistical relevance (FDR ≤ 0.05 and a p-value ≤ 0.01). The passing filter pathways were used to build interaction maps to give insight about general processes. This map was computed using Cytoscape **Enrichment Map** app. Final map was refined with custom algorithms. For example, Jak Stat group is confirmed for several enriched pathways that has shared genes.

<p align="center">
  <img width="450" height="450" src="https://user-images.githubusercontent.com/60892768/75078790-3bb4fd00-54cc-11ea-8501-86410c7cb78e.png">
 </p>

## Differential Gene Expression Analysis
Differential gene expression analysis (DGEA) is done by pairs, so, one file is needed per possible comparison of two given the number of clusters (eg. If there are 3 clusters, there 3 possible combinations, AvsB, AvsC and BvsC). DGEA was performed by the **SCDE** R package. Error models were fitted to split drop-out events(technical errors) from biological data.

### Genetic Signatures

Differentially expressed genes were selected with a threshold of |Log<sub>2</sub>(Fold Change)| ≥ 4 and p-value ≤ 0.01. So, there is a list of differentially expressed genes for each comparison. To get globally expressed genes for each cluster, we constructed intersection lists, they contained genes that are overexpressed in all comparasissons related to a specific cluster. For example, a gene **GeneX** either must be overexpressed for Cluster A in AvsB and AvsC comparisons to be selected as global(in the case of 3 clusters). This process can represented as a restructuration of a Venn diagram as shown below.

<p align="center">
  <img width=900" height="450" src="https://user-images.githubusercontent.com/60892768/75173471-8ca73a00-56f4-11ea-9a60-4142bcda8bfe.png">
 </p>

Finally, globally differentially expressed genes are visualized in a expression heatmap. Upper bar represents Clusters.

<p align="center">
  <img width="450" height="450" src="https://user-images.githubusercontent.com/60892768/75082905-8ab65e80-54db-11ea-9365-2e7bd3da90dc.png">
 </p>
