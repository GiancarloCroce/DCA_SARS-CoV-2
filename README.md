# Epistatic prediction of mutable sites in SARS-CoV-2 proteins and epitopes

Jupyter notebooks ```dca_sarscov2.ipynb``` to reproduce the results of [papername](link_to_paper).
You can run the notebook Google Colab at [this link](https://colab.research.google.com/github/GiancarloCroce/DCA_SARS-CoV-2/blob/main/dca_sarscov2.ipynb) without downloading the GitHub repository.

# Data

## Mutability score

Direct Coupling Analysis (DCA) mutability score for each site of the SARS-CoV-2 proteome included in a PFAM domain: ```./data/data_dca_whole_proteome.csv```.

## Mutability score and IEDB response frequency (RBD domain)

For the RDB domain of the spike protein we combine our predictions with the [IEDB Response Frequency](https://www.iedb.org/immunomebrowser.php?cookie_id=638356&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC2&source_antigen_name=Spike+glycoprotein)
to identify sites that are predicted to be mutable and are shared by multiple positively responding epitopes ```./data/data_dca_iedb_RDB_domain.csv ```

See also the section "Predicting immunologically relevant mutable sites" in the ```dca_sarscov2.ipynb``` notebook.

## Protein domains

Direct Coupling Analysis (DCA) and Indepented (IND) models are inferred for all protein domains in ```./data/data_meff.csv ```
Meff is the effective number of sequences available for a protein domain. 
As a thumb rule: more sequence data -> better models -> increased predictive power.

