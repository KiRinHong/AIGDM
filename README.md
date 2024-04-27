## File Description

This repository is used for depositing the analysis scripts, derived data/results and corresponding figures and tables for the paper "A marginal regression model for longitudinal compositional counts with application to microbiome data". 

Here are brief descriptions for each folder, detailed explanations are commented within code and described in the paper:

* a folder named **Analysis**, which contains the real data analysis R scripts. 
  - 0.prepareData.R
  - 1.runModel.R
  - 2.generateRslt.R
  - 3.realdataVis.R
  - \*.utility.R

* a folder named **Data**, which contains raw data and derived data.

* a folder named **Figs**, which contains the figures in the main text and supplementary information.

* a folder named **Simulation**, which contains the simulation shell scripts and R scripts.
  - **ProposedModel**: Type1 and Power simulations.
  - realpara_gdmNlnm.R and realpara_utility.R: Generate the para.Rdata, which is used as a basis in the simulation.
  - generateTxt.R: Generate the input.txt, which is different scenario used in simulation.

## R Package

The proposed method **AIGDM** in the paper is implemented in the Analysis/1a.AIGDM_utility.R, will be incoporated in the R pacakge [miLineage](https://github.com/tangzheng1/tanglab/blob/gh-pages/software/miLineage_v3.1.zip) in the future.

## Contact

* Qilin (Kirin) Hong - qhong8@wisc.edu
