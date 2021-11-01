---
title: DCA for SARS-CoV-2 
---

**Epistatic models predict mutable sites in SARS-CoV-2 proteins and epitopes**

For each site of the SARS-CoV-2 Wuhan-Hu-1 proteome (Accession [NC_045512](https://www.genome.jp/dbget-bin/www_bget?refseq:NC_045512)) included in a [PFAM domain](http://pfam.xfam.org/), we introduce a  [Direct Coupling Analysis](https://en.wikipedia.org/wiki/Direct_coupling_analysis) **mutability score to predict mutable and constrained sites**.

We validate our mutability predictions with the mutations observed in SARS-CoV-2 proteomes deposited in the [GISAID](https://www.gisaid.org/) database.

Paper: [papername](link_to_paper).

Run the Jupyter-notebook ```dca_sarscov2.ipynb``` on Google Colab at [this link](https://colab.research.google.com/github/GiancarloCroce/DCA_SARS-CoV-2/blob/main/dca_sarscov2.ipynb) to reproduce key results from the paper and guide the data analysis. Colab allows you to execute the Python code through your browser (no need to clone the GitHub repository on your local machine). 

![](pipeline2.png)

The results of our analysis are also available in the ```./data ``` folder on the [Github page](https://github.com/GiancarloCroce/DCA_SARS-CoV-2/)

The data structure is:
```
protein  domain	      position_protein  position_domain  aa_Wuhan-Hu-1  mutability_score(IND) 	mutability_score(DCA)  observed_mut_May2021  observed_mut_Dec2020  observed_mut_Jul2020
Spike  	 bCoV_S1_RBD  349               1.0              S              -1.3818			-1.2046			0.0                   0.0                   0.0
Spike  	 bCoV_S1_RBD  350               2.0              V              -1.9788			-1.1667			3.0                   0.0                   0.0
Spike  	 bCoV_S1_RBD  351               3.0              Y              -1.8017			-1.2678			8.0                   1.0                   1.0
Spike    bCoV_S1_RBD  352               4.0              A              -1.0140			-1.1347			23.0                  3.0                   0.0
```

For sites in the [RDB domain](http://pfam.xfam.org/family/bCoV_S1_RBD), we also include the [IEDB](https://www.iedb.org/) site response frequency and the corresponding 95% confidence interval upper/lowerbound:
```
protein  domain	      position_protein	... 	IEDB_upperbound  IEDB_lowerbound  IEDB_response_frequency
Spike  	 bCoV_S1_RBD  349             	... 	0.10             0.05             0.073563
Spike  	 bCoV_S1_RBD  350             	... 	0.10             0.05             0.069284
Spike  	 bCoV_S1_RBD  351             	... 	0.12             0.07             0.093023
Spike    bCoV_S1_RBD  352             	... 	0.10             0.05             0.069212
```

In [papername](link_to_paper) we limit our analysis to the case of SARS-CoV-2. However, our approach requires only a single reference genome to identify distant homologs and make predictions. It can potentially be extended to any virus, as long as sufficient sequence data are available to train reliable models. Code for training sequence-based models to predict mutability scores is available [here](https://github.com/juan-rodriguez-rivas/covmut). It computes both independent and DCA epistatic models, with the latter providing a better prediction of the mutability in the vast majority of cases.
