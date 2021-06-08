# Epistatic model predicts SARS-CoV-2 mutability of proteins and epitopes

Jupyter notebooks ```paper_dca_sarscov2.ipynb``` to reproduce the results of [papername](link_to_paper).

# Data

## Mutability score

Direct Coupling Analysis (DCA) mutability score for each site of the SARS-CoV-2 proteome included in a PFAM domain: ```./data/data_dca_whole_proteome.csv```.

## Mutability score and IEDB response frequency (RBD domain)

For the RDB domain of the spike protein we combine our predictions with the [IEDB Response Frequency](https://www.iedb.org/immunomebrowser.php?cookie_id=638356&source_organism=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCBITaxon_2697049&source_organism_name=SARS-CoV2&source_antigen=http%3A%2F%2Fwww.uniprot.org%2Funiprot%2FP0DTC2&source_antigen_name=Spike+glycoprotein)
to identify sites that are predicted to be mutable and are shared by multiple positively responding epitopes ```./data/data_dca_iedb_RDB_domain.csv ```

See also the section "Predicting immunologically relevant mutable sites" in the ```paper_dca_sarscov2.ipynb``` notebook.
