#!/bin/bash

R -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc/');rmarkdown::render('supplementary.Rmd',output_file='Bacci_et_al_2020_SM.html')"
