# Homoplasy
R code for the Relative Homoplasy Index, a metric for quantifying homoplasy in discrete phylogenetic datasets

Relative Homoplasy Index: A New Cross-comparable Metric for Quantifying Homoplasy in Discrete Character Datasets
https://doi.org/10.5061/dryad.34tmpg4rp

This dataset includes all nexus, phylip and .R files to reproduce analysis included in 'Relative Homoplasy Index: A New Cross-comparable Metric for Quantifying Homoplasy in Discrete Character Datasets'.
Description of the data and file structure
Six published datasets are included:
Passeriformes a & b (Steell et al 2023)
Neornithes (Field et al 2020)
Telluraves (Ksepka et al 2019)
Avialae a & b (Benito et al 2022)
Each dataset has a .R file to prepare for further analysis (e.g., 'Passeriformes_a_prep.R') which can be used with the appropriate data files (either .nex or .phy for trees, and .nex for matrices).

The following new R functions are included:
RHI() - A function to calculate the relative homoplasy index from a given matrix and tree.
subset.trees() - A function to randomly subset tips in phylogenetic trees n times.
subset.matrix() - A function to subset matrices based on a list of subset trees, to be used after subset.trees().
statef() - A function to give the character state frequencies in a given matrix.
contrast.matrix() - A function to create a contrast matrix for use with the 'phangorn' package tools such as the 'phyDat' object. This function is in development and currently can only create contrast matrices if there are no polymorphic states present in the dataset.

Each analysis (Analyses 1-6 in the manuscript) has accompanying annotated R code incorporating the datasets and new functions.
Sharing/Access information
Links to other publicly accessible locations of the data:
	•	GitHub (including code updates): http://github.com/LizzySteell/Homoplasy
Data was derived from the following sources:
	•	Steell EM, Nguyen JM, Benson RB, Field DJ. Comparative anatomy of the passerine carpometacarpus helps illuminate the early fossil record of crown Passeriformes. Journal of Anatomy. 2023 Mar;242(3):495-509.
	•	Field DJ, Benito J, Chen A, Jagt JW, Ksepka DT. Late Cretaceous neornithine from Europe illuminates the origins of crown birds. Nature. 2020 Mar 19;579(7799):397-401.
	•	Benito J, Kuo PC, Widrig KE, Jagt JW, Field DJ. Cretaceous ornithurine supports a neognathous crown bird ancestor. Nature. 2022 Dec 1;612(7938):100-5.
	•	Ksepka DT, Grande L, Mayr G. Oldest finch-beaked birds reveal parallel ecological radiations in the earliest evolution of passerines. Current Biology. 2019 Feb 18;29(4):657-63.
Code/Software
All code was run in R version 4.2.3 (2023-03-15) -- "Shortstop Beagle".

Code for this project is being curated on GitHub (http://github.com/LizzySteell/Homoplasy)

The following packages and their dependencies were used in the development of this code:
ape v5.7-1
Paradis, E., Claude, J. and Strimmer, K. (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20, 289–290.
Paradis, E. and Schliep, K. (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics, 35, 526–528.

phytools v1.5-1
Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223.

phangorn v2.11.1
Schliep KP. phangorn: phylogenetic analysis in R. Bioinformatics. 2011 Feb 15;27(4):592-3.
Schliep K, Paradis E, de Oliveira Martins L, Potts A, White TW, Stachniss C, Kendall M. Package ‘phangorn’. Available at< Available at https://cran. r-project. org/web/packages/phangorn/phangorn. pdf>. Accessed on 01th August. 2019 Jun 19.

TreeTools v1.8.0
DOI: 10.5281/zenodo.3522726
Smith MR.

TreeSearch v1.2.0
Smith MR. TreeSearch: morphological phylogenetic analysis in R. bioRxiv. 2021 Nov 10:2021-11.

geiger v2.0.10
LJ Harmon, J Weir, C Brock, RE Glor, W Challenger, G Hunt, R FitzJohn, MW Pennell, GJ Slater, JW Brown, J Uyeda, and JM Eastman

dispRity v1.7.0
Guillerme T. dispRity: a modular R package for measuring disparity. Methods in Ecology and Evolution. 2018 Jul;9(7):1755-63.

fitdistrplus v1.1-8
Marie-Laure Delignette-Muller and Christophe Dutang.

