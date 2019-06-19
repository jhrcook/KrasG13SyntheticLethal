#!/bin/bash

module load gcc R

# prepare data
Rscript subscripts/data_preparation.R

# linear model
Rscript subscripts/linear_model.R

# predicting KRAS mutation by genetic dependencies
Rscript subscripts/predict_rasallele.R

# check overlap with Alex (UCSF)
Rscript subscripts/overlap_alex.R

# functional annotation of genetic dependency hits
Rscript subscripts/hit_annotation.R

# render README
Rscript -e 'rmarkdown::render("README.Rmd")'

