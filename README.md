# Code and data for "A gene regulatory network reveals features and regulators of the root response to elevated CO2 in Arabidopsis", O. Cassan et al., 2023.


This repository contains data and Rmarkdown scripts to reproduce the figures and analyses from "A gene regulatory network reveals features and regulators of the root response to elevated CO2 in Arabidopsis", O. Cassan et al., 2023, published in the New Pyhtologist : https://doi.org/10.1111/nph.18788.


## Data

The folder `Data` contains phenotypic observations and expression counts under the combinatorial design in csv format. `Data/Candidates` contains the phenptypes of the candidate genes identified in GRN inference.

## Phenotypic analysis

Contains the code for biomass, N content, Fe content, and N absorption plots under the combinatorial design (`Phenotype_anamysis.Rmd`).
It also contains the code for the validation of candidate genes (genotype x CO2 interaction tests and figure in `Candidate_regulators.Rmd`).

## Transcriptomic analysis

Contains the scripts for : 

+ Transcriptome normalisation, PCA, and differential expression analyses (`PCA_DEA.Rmd`)

+ Heatmaps for nitrate and iron nutrition genes (`Hetamaps.Rmd`)

+ Co-expression clustering with [Coseq](https://www.bioconductor.org/packages/release/bioc/html/coseq.html) (`Coexpression_clustering.Rmd`)

+ Gene Regulatory Network inference with [DIANE](https://oceanecsn.github.io/DIANE/) (`GRN_inference.Rmd`)

+ Gene Regulatory Network exploration, validation, and visualisation (`GRN_analysis.Rmd`)


### Dependencies


Many packages are used in those scripts, and can by installed from the CRAN in the usual way using `install.packages()`, except for DIANE, that should be installed as follows :


```
library(remotes)
install_github("OceaneCsn/DIANE")
```