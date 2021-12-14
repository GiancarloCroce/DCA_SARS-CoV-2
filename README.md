# Epistatic models predict mutable sites in SARS-CoV-2 proteins and epitopes 

GitHub repository accompanying the paper [Epistatic models predict mutable sites in SARS-CoV-2 proteins and epitopes](https://www.biorxiv.org/content/10.1101/2021.12.11.472202v1)(Juan Rodriguez-Rivas, Giancarlo Croce, Maureen Muscat, Martin Weigt, biorXiv 2021). Doi: [https://doi.org/10.1101/2021.12.11.472202](https://doi.org/10.1101/2021.12.11.472202)

Run the Jupyter-notebook ```dca_sarscov2.ipynb``` to reproduce key results from the paper and guide the data analysis.  You can run the notebook directly on Google Colab at [this link](https://colab.research.google.com/github/GiancarloCroce/DCA_SARS-CoV-2/blob/main/dca_sarscov2.ipynb) without downloading the GitHub repository.

# Data

## Mutability score

Direct Coupling Analysis (DCA) mutability score for each site of the SARS-CoV-2 proteome included in a PFAM domain: ```./data/data_dca_whole_proteome.csv```.

## Mutability score and IEDB response frequency (RBD domain)

We combine our **DCA-Mutability Score** predictions with the **IEDB-Response Frequency ** to identify sites that are predicted to be mutable and are shared by multiple positively responding epitopes ```./data/IEDB_updated_data/```
The section "Predicting immunologically relevant mutable sites" of the ```dca_sarscov2.ipynb``` notebook, explains how to analyse and interpret the data.

## Beyond the Spike protein 

A key advantage of our data-driven modeling approach is the possibility to obtain predictions for all the protein domains in the SARS-CoV-2 proteome. Run the Jupyter-notebook to extend the DCA predictions to all 39 protein domains covering 81% of the entire proteome (8037 out of 9748 positions), and combine them with immunological IEDB data.
Raw data are available at  ```./data/data_dca_proteome.csv``` and ```./data/IEDB_updated_data/```. 
