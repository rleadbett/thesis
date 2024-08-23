# /r-code/ch-2

Code to reproduce the analysis and figures of Chap. 2.

## background-figures.R

An R script that generates Figs. 2.1--2.4 in Sec. 2.1 of the Thesis.

## left-trunc-sim-study-\*.R and LT_experiment_results_df\*.rds

The R scripts to run the simulation study in Sec. 2.5 of the chapter (for different values of beta) and their output R objects. These scripts are computationally intensive and were run on Pawsey's cloud compute service Nimbus.

## sup-material-ch2.qmd and sup-material-ch2.html

The Quarto document to reproduce the main analysis in Sec. 2.3--2.5. The `.qmd` file contains the Quarto source code with impeded R code and commentary and the `.html` file is the rendered report so that the reader can be easily followed along.
