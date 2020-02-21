<p align="center">
  <img src="https://user-images.githubusercontent.com/60892768/74993425-88d39900-5410-11ea-8643-b701551d0472.png">
</p>

# Bicycle-Gummy

The pipeline described was developed by Aarón Vázquez-Jiménez in the [Human Systems Biology Group](https://resendislab.github.io/) at INMEGEN.

## Atributtion
This pipeline was used to study functionally heterogeneity in the MCF7 cell line using Multicelluas Tumor Spheroids as tumor model. If you have any questions contact us at avazquezj(at)inmegen.gob.mx and oresendis (at) inmegen.gob.mx. 

## Getting Started

The complete pipeline starts with the fastq files. It ends with the files related to the differential gene expression analysis and the ones necceary to use the GSEA tool. The complete [pipeline](Pipeline/pipeline.md) is formed by the following steps:

* [Sample Processing](Pipeline/pipeline.md#samples\processing)
  - [Data Availabity](Pipeline/pipeline.md#samples\data\availability)
* [Knee filter](Pipeline/pipeline.md#knee\filter)
* [Gene Over-dispersion](Pipeline/pipeline.md#gene\over-disperssion)
* [Gene number selection](Pipeline/pipeline.md#gene\number\selection)
* [Dimensionality Reduction](Pipeline/pipeline.md#Dimensionality\Reduction)
* [Clustering](Pipeline/pipeline.md#Clustering)
  - [Kmeans](Pipeline/pipeline.md#kmeans)
  - [EM](Pipeline/pipeline.md#expectation-maxinization-algorithm)
* [Enrichment Analysis](Pipeline/pipeline.md#Enrichment\Analysis)
