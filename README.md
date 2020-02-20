# Bicycle-Gummy

The pipeline described was developed by Aarón Vázquez-Jiménez in the [Human Systems Biology Group](https://resendislab.github.io/) at INMEGEN.

## Atributtion
This pipeline was used to study functionally heterogeneity in the MCF7 cell line using Multicelluas Tumor Spheroids as tumor model. If you have any questions contact us at avazquezj(at)inmegen.gob.mx and oresendis (at) inmegen.gob.mx. 

## Getting Started

The complete pipeline is formed by the following steps:

### Samples Processing

Raw files were processed througth ddSeeker pipeline ([Rogmanoli et al. 2018](https://link.springer.com/epdf/10.1186/s12864-018-5249-x?author_access_token=5GkMGb1JtnR887KgJLyh-m_BpE1tBhCbnbw3BuzI2ROCtPkWTeG4740r3l1fSwcikEVJgHZm9jmSpiShpuk2FV8ae_KUm1O2Kb8nf4xLIHEJfbwJ3tasCADHjVdZ23iWPECc69GpWSYTzS6lSpBaJA%3D%3D)) using **STAR** as aligner. As output, BAM files for every sample were created. Count matrix were generated by **DigitalExpression** function of the [Drop-Seq tools](https://github.com/broadinstitute/Drop-seq) with a value of NUM_CORE_BARCODES=5000.

A count matrix is generated for every of the 4 lines. All files can be merged into one.

#### Data Availabity 
The data used and discussed have been deposited in NCBI's Gene Expression Omnibus and are accessible through the accession number [GSE145633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145633).

### knee filter
Since we used droplet based tecnology for the hight throuput sequencing, data must be filtered by taking out droplets with no RNA from the biological samples, so, only informative cells remain. To implement knee filter we ploted the cumulaitve UMI sum   for all cells. The point where the firts inflection point occurs set the threshold. 

![Knee filter](https://github.com/Aaron-Vazquez/Bicycle-Gummy/tree/master/images/Knee_plot.png)
