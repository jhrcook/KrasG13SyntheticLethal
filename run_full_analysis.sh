#!/bin/bash

module load gcc R

# prepare data
Rscript subscripts/data_preparation.R

# linear model
Rscript subscripts/linear_model.R

# render README
Rscript -e 'rmarkdown::render("README.Rmd")'