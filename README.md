# Code and data for "A gene regulatory network reveals features and regulators of the root response to elevated CO2 in Arabidopsis", O. Cassan et al., 202X.


This repository contains data and Rmarkdown scripts to reproduce the figures and analyses from "A gene regulatory network reveals features and regulators of the root response to elevated CO2 in Arabidopsis", O. Cassan et al., 202X.


## Data

The folder `Data` contains phenotypic observations and expression counts under the combinatorial design in csv format. `Data/Candidates` contains the phenptypes of the candidate genes identified in GRN inference.

## Phenotypic analysis

Contains the code for biomass, N content, Fe content, and N absorption plots under the combinatorial design.
It also contains the code for the validation of candidate genes (genotype x CO2 interaction tests, and figure).

## Transcriptomic analysis

Contains the scripts for : 

+ Transcriptome normalisation, PCA, and differential expression analyses

+ Co-expression clustering with [Coseq](https://www.bioconductor.org/packages/release/bioc/html/coseq.html)

+ Gene Regulatory Network inference with [DIANE](https://oceanecsn.github.io/DIANE/)

Many packages are used in those scripts, and can by installed from the CRAN in the usual way using `install.packages()`, except for DIANE, that should be installed as follows :

### Dependencies

```
library(remotes)
install_github("OceaneCsn/DIANE")
```