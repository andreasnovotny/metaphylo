# metaphylo R package V 0.1.0

This is an R package with the purpose of sharing data and data analysis.

- Author: Andreas Novotny
- Contact: andreas.novotny@su.se
- source: https://github.com/andreasnovotny/metaphylo

**Plankton Diet Metabarcoding Project, Stockholm University 2018**


## Installation

The package is in a private repository at GitHub.
Install with devtools:

```
library(devtools)

install_github("andreasnovotny/metaphylo", auth_token = "cc00cf689d3a412aa3e3082d19f4133a286452ed", build_vignettes=TRUE)

library("metaphylo")
```

## Dependencies

Before installation, make shore that these packages are installed:
(dada2, phyloseq, DESeq2, ggplot2, bipartite, magrittr)

## Instructions

Instructions for using the package are in the package vignette:

```
browseVignettes("metaphylo")
```

