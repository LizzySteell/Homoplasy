# Homoplasy
R code for the Relative Homoplasy Index, a method for quantifying homoplasy in discrete phylogenetic datasets

Preprint published at BioRxiv: doi: https://doi.org/10.1101/2023.10.10.561677

Relative Homoplasy Index: A New Cross-comparable Index for Quantifying Homoplasy in Discrete Character Datasets
https://doi.org/10.5061/dryad.34tmpg4rp

This dataset includes all nexus, phylip and .R files to reproduce the analyses included in 'Relative Homoplasy Index: A New Cross-comparable Metric for Quantifying Homoplasy in Discrete Character Datasets'.

Description of the data and file structure
CHARACTER-TAXON MATRICES
Four published phylogenetic discrete character matrices are included in nexus file format:
- "Passeriformes_matrix.nex" for datasets Passeriformes a & b (Steell et al 2023)
"Neornithes_matrix.nex" for dataset Neornithes (Field et al 2020)
"Telluraves_matrix.nex" for dataset Telluraves (Ksepka et al 2019)
"Avialae_matrix.nex" for datasets Avialae a & b (Benito et al 2022)
These files can be opened using a text editor application, or opened in a phylogenetics software programme such as Mesquite. Each nexus file consists of:
A block of taxon labels (except for Avialae) and the dimensions of the dataset (NTAX is number of taxa, NCHAR is number of characters).
A matrix with taxon labels at the start of each row, and a character state scoring for each character. Character is a categorical morphological character, where '?' indicates an ambiguity and '-' indicates a non-applicable scoring (e.g., for a dependent character). Some matrices contain polymorphic states (indicated by either '{01}' or '(01)' and can contain any number of combination of states present for that character).
Some matrices contain a list of ordered characters in a MrBayes format. They may also contain MrBayes phylogenetic analysis code.

PHYLOGENETIC TREE FILES
Each phylogenetic matrix is accompanied by one or two phylogenetic tree files. Each file is in either phylip or nexus format, and depicts the phylogenetic relationships between taxa through parentheses. These files can be opened in a text editor.
The following tree files are included:
"Passeriformes_a_tree.phy" for dataset Passeriformes a.
"Passeriformes_b_tree.phy" for dataset Passeriformes b.
"Neornithes_tree.nex" for dataset Neornithes.
"Telluraves_tree.phy" for dataset Telluraves.
"Avialae_a_tree.nex" for dataset Avialae a.
"Avialae_b_tree.nex" for dataset Avialae b.

PREPARATORY R SCRIPTS FOR EACH DATASET
Each matrix and dataset-specific tree file should be read into R using the appropriate accompanying dataset-specific R file, which standardises the datasets. These R scripts are as follows:
"Passeriformes_a_prep.R" for dataset Passeriformes a.
"Passeriformes_b_prep.R" for dataset Passeriformes b.
"Neornithes_prep.R" for dataset Neornithes.
"Telluraves_prep.R" for dataset Telluraves.
"Avialae_a_prep.R" for dataset Avialae a.
"Avialae_b_prep.R" for dataset Avialae b.

R FUNCTIONS
The following R script files describe code for new R functions that are used within the manuscript:
"RHI().R" - A function to calculate the relative homoplasy index (RHI) from a given matrix and tree.
"subset.trees().R" - A function to randomly subset tips in phylogenetic trees n times.
"subset.matrix().R" - A function to subset matrices based on a list of subset trees, to be used after subset.trees().
"statef().R" - A function to give the character state frequencies in a given matrix.
"contrast.matrix().R" - A function to create a contrast matrix for use with the 'phangorn' package tools such as the 'phyDat' object. This function is in development and currently can only create contrast matrices if there are no polymorphic states present in the dataset.

R SCRIPTS FOR ANALYSES 1-6
Each analysis (Analyses 1-6 in the manuscript) has accompanying annotated R code incorporating the datasets and new functions, in the following files:
"Analyses_1.R" Code for inflating the transition rate for characters in simulated matrices.
"Analyses_2.R" Code for subtree sampling in simulated matrices.
"Analyses_3.R" Code for subsetting characters in simulated matrices.
"Analyses_4.R" Code for calculating relative homoplasy index and other homoplasy indexes for empirical datasets.
"Analyses_5.R" Code for subsetting empirical trees.
"Analyses_6.R" Code for subsetting empirical characters.

Sharing/Access information
Links to other publicly accessible locations of the data:

GitHub (including code updates): http://github.com/LizzySteell/Homoplasy
Data was derived from the following sources:

Steell EM, Nguyen JM, Benson RB, Field DJ. Comparative anatomy of the passerine carpometacarpus helps illuminate the early fossil record of crown Passeriformes. Journal of Anatomy. 2023 Mar;242(3):495-509.
Field DJ, Benito J, Chen A, Jagt JW, Ksepka DT. Late Cretaceous neornithine from Europe illuminates the origins of crown birds. Nature. 2020 Mar 19;579(7799):397-401.
Benito J, Kuo PC, Widrig KE, Jagt JW, Field DJ. Cretaceous ornithurine supports a neognathous crown bird ancestor. Nature. 2022 Dec 1;612(7938):100-5.
Ksepka DT, Grande L, Mayr G. Oldest finch-beaked birds reveal parallel ecological radiations in the earliest evolution of passerines. Current Biology. 2019 Feb 18;29(4):657-63.
Code/Software
All code was run in R version 4.2.3 (2023-03-15) -- "Shortstop Beagle".

Code for this project is being curated on GitHub (http://github.com/LizzySteell/Homoplasy)

The following packages and their dependencies were used in the development of this code:
ape v5.7-1
Paradis, E., Claude, J. and Strimmer, K. (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20, 289–290.
Paradis, E. and Schliep, K. (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics, 35, 526–528.

phytools v1.5-1
Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223.

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

Methods
Data used in this manuscript was downloaded from published papers with available nexus files. Data from the following papers were used in analyses:

Datasets 'Passeriformes a' and 'Passeriformes b': Steell EM, Nguyen JM, Benson RB, Field DJ. Comparative anatomy of the passerine carpometacarpus helps illuminate the early fossil record of crown Passeriformes. Journal of Anatomy. 2023 Mar;242(3):495-509.

Datasets 'Avialae a' and 'Avialae b': Benito J, Kuo PC, Widrig KE, Jagt JW, Field DJ. Cretaceous ornithurine supports a neognathous crown bird ancestor. Nature. 2022 Dec 1;612(7938):100-5.

Dataset 'Neoaves': Field DJ, Benito J, Chen A, Jagt JW, Ksepka DT. Late Cretaceous neornithine from Europe illuminates the origins of crown birds. Nature. 2020 Mar 19;579(7799):397-401.

Dataset 'Telluraves': Ksepka DT, Grande L, Mayr G. Oldest finch-beaked birds reveal parallel ecological radiations in the earliest evolution of passerines. Current Biology. 2019 Feb 18;29(4):657-63.

Other datasets used in this paper were generated through simulations using the phylogenetic topologies from the above six datasets.

All analyses were carried out in the R programming language with R version 4.2.3 (2023-03-15) -- "Shortstop Beagle".

Funding
Natural Environment Research Council, Award: NE/S007164/1

Natural Environment Research Council, Award: MR/S032177/1

Swedish Research Council Starting Grant Within Natural and Engineering Sciences, Award: ÄR-NT 2020-03515
