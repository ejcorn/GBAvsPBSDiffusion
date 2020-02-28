# GBAvsPBSDiffusion

Code to reproduce all analysis in Henderson et al. 2019 ("Glucocerebrosidase activity modulates neuronal susceptibility to pathological αsynuclein insult"). In these analyses, we use a linear network diffusion model to explain the spread of α-synuclein pathology through the brain over time after injection of misfolded α-synuclein into the caudoputamen in mice that were coinjected with either vehicle (DPBS group) or glucocerebrocidase inhibitor (CBE group).

## Requirements:
  - R 3.3.3 or later. Requisite packages are listed in code/misc/packages.R

## Directory structure

Master branch contains 2 major folders:
  - code/process contains code that deals with loading and processing data
  - code/misc contains various helper functions for analysis and plotting that are called in other scripts
  - code/diffmodel contains code that uses linear diffusion models to predict spread of protein through structural connectome.
  - code/PBSvsCBE contains code that compares time constant values between PBS and CBE-injected mice
  - data/ contains csv and xlsx files with 1) experimentally obtained pathology data and 2) parcellated Snca expression data and connectome data from Allen Brain Institute. 
  
Using gene expression and linear dynamics of spread along the connectome, we attempt to predict the spatial distribution of experimentally observed pathology at 1 month post injection, and test whether injection of a glucocerebrocidase inhibitor impacts spreading dynamics.

## Input specification

The file ‘pipeline.R’ is located in the main directory. This file will coordinate the sequential execution of all scripts within the code/ folder, generating all the figures in the paper and more. Custom specification of the following inputs at the top of ‘pipeline.R’ is required:
  - basedir:  path to the main directory containing the 'code' and 'data' folders 
  - opdir: name of output directory that contains all results, which will be housed in basedir. the line params$opdir <- paste(...) can be replaced with any string that ends with a '/', and that string will serve as the name of the output directory.
  - grps: character vector containing the name of groups in data file to test. For our data set, these were 'DPBS' and 'CBE'
  - c.min: minimum value for time constant sweep (0 is fine)
  - c.max: maximum value for time constant sweep (10 is fine)
  - c.n: number of time constants to test, linearly spaced between c.min and c.max (100 is fine)

The script moveresults.sh compiles the pipeline output for a specific output folder name and creates a new folder with the figures that ultimately were included in the manuscript, along with relevant stats output as a .log file.

## Questions, suggestions, comments?

Please contact Eli Cornblath (Eli ~`DOT`~ Cornblath ~`AT`~ pennmedicine.upenn.edu) with any questions regarding network analysis and code, and contact Mike Henderson (hendm ~`AT`~ pennmedicine.upenn.edu) with any questions regarding experiments and data.
