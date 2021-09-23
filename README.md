# Tollervey-iCLIP

This repo contains the code for the re-analysis of the human brain iCLIP data from the Tollervey et al., 2011 [paper](https://www.nature.com/articles/nn.2778).
The data is preprocessed, mapped to the GRCh38 genome and cross-links are identified with the [iCount](https://icount.readthedocs.io/en/latest/) pipeline.
iCount is run using [snakemake](https://doi.org/10.12688/f1000research.29032.1) and specfically a setup adapted from the [ARMOR workflow](https://doi.org/10.1534/g3.119.400185)([github](https://github.com/csoneson/ARMOR)).
Downstream analysis is performed with R (see the Rmd/ directory).
