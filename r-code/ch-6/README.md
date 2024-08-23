# /r-code/ch-6

Code to reproduce the analysis and figures of Chap. 6.

## sup-material-ch6.qmd and sup-material-ch6.html

The Quarto document to reproduce the analysis and figures in Chap. 6. The `.qmd` file contains the Quarto source code with impeded R code and commentary and the `.html` file is the rendered report so that the reader can be easily followed along.

## elppd.R

An R script that runs elppd calculations in Sec. 6.6 and outputs Table 6.1.

## pit.R and ppc-obj.rds

An R script that runs the posterior predictive checks using the test quantity $Min(z)$ from Sec. 6.1. The code saves the output as `ppc-obj.rds` which is then read in when the main quarto document (`sup-material-ch6.qmd`) is run.

## ft-cdf.R and ft-\*-obs-\*.rds

`ft-cdf.R` is an R script that generates the failure time distributions from the first seven, eight, and nine observations from both models and saves them as the `.rds` objects `ft-\*-obs-\*.rds` to be read in by `sup-material-ch6.qmd` to generate the failure time plots in Figs. 6.13 amd 6.14.

## compiled-stan-model-\*.rds and stan-data.rds

The defined stan model and data to avoid the need to compile the models and prep the data in `elppd.R`, `pit.R`, and `ft-cdf.R`.

## belt-cm-example.PNG and animation_\*.gif

Supplementary figures that are called by `sup-material-ch6.qmd`.
