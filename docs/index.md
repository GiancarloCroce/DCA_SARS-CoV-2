---
title: DCA for SARS-CoV-2 
---

**Epistatic model predicts SARS-CoV-2 mutability of proteins and epitopes**

For each site of the SARS-CoV-2 Wuhan-Hu-1 proteome (Accession [NC_045512](https://www.genome.jp/dbget-bin/www_bget?refseq:NC_045512)) included in a [PFAM domain](http://pfam.xfam.org/), we introduce a  [Direct Coupling Analysis](https://en.wikipedia.org/wiki/Direct_coupling_analysis) **mutability score to predict mutable and constrained sites**.

The predictions are validated using the genomes of SARS-CoV-2 strains from the [GISAID](https://www.gisaid.org/) database.

Paper: [papername](link_to_paper).
![](pipeline.png)

The results of our analysis are available in the ```./data ``` folder of the [Github page](https://github.com/GiancarloCroce/DCA_SARS-CoV-2/)

The data structure is:
```
protein  domain	      position_protein  position_domain  aa_Wuhan-Hu-1  mutability_score(IND)  mutability_score(DCA)  observed_mut_Mar2021  observed_mut_Nov2020  observed_mut_Jul2020
Spike  	 bCoV_S1_RBD  349               1.0              S              -3.7992                -7.0714                0.0                   0.0                   0.0
Spike  	 bCoV_S1_RBD  350               2.0              V              -5.4405                -6.8487                1.0                   0.0                   0.0
Spike  	 bCoV_S1_RBD  351               3.0              Y              -4.9535                -7.4424                6.0                   1.0                   1.0
Spike    bCoV_S1_RBD  352               4.0              A              -2.7878                -6.6612                17.0                  3.0                   2.0
```

For sites in the [RDB domain](http://pfam.xfam.org/family/bCoV_S1_RBD), we also include the [IEDB](https://www.iedb.org/) site response frequency and the corresponding 95% confidence interval upper/lowerbound:
```
protein  domain	      position_protein	... 	IEDB_upperbound  IEDB_lowerbound  IEDB_response_frequency
Spike  	 bCoV_S1_RBD  349             	... 	0.10             0.05             0.073563
Spike  	 bCoV_S1_RBD  350             	... 	0.10             0.05             0.069284
Spike  	 bCoV_S1_RBD  351             	... 	0.11             0.06             0.083333
Spike    bCoV_S1_RBD  352             	... 	0.10             0.05             0.069212
```

We provide the  Jupyter-notebook ```dca_sarscov2.ipynb``` to guide data analysis and reproduce the paper results.
