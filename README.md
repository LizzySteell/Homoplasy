# Homoplasy

**R code for the 'Revealing patterns of homoplasy in discrete phylogenetic datasets with a cross-comparable index' in _Zoological Journal of the Linnean Society_ (2025)**

Article DOI: doi.org/10.1093/zoolinnean/zlaf024

Authors: Elizabeth M. Steell (University of Cambridge/Girton College, Cambridge); Allison Y. Hsiang (University of Stockholm); Daniel J. Field (University of Cambridge).

This dataset includes all Nexus, Phylip and .R files to reproduce the analyses included in the above article. The custom functions are regularly maintained and may differ from the R code provided in the Supporting Information associated with the published article.

**DATA FILES**

_Empirical tree files_

These Newark format PHYLIP files represent the topologies derived from the empirical studies and were used in Analysis 1. Some are rooted and some are not fully bifurcating.

•	Telluraves (Ksepka et al., 2019): ‘Telluraves_tree.phy’.

•	Passeriformes a (Steell et al., 2023): ‘Passeriformes_a_tree.phy’.

•	Passeriformes b (Steell et al., 2023): ‘Passeriformes_b_tree.phy’.

•	Neornithes (Field et al., 2020): ‘Neornithes_tree.nex’.

•	Avialae a (Benito et al., 2022): ‘Avialae_a_tree.nex’.

•	Avialae b (Benito et al., 2022): ‘Avialae_b_tree.nex’.

_Dichotomous trees for simulations_

These Newark format PHYLIP files represent fully dichotomous trees where nodes were randomly resolved when populating branch lengths using the accelerated transformation algorithm (see extended methods above). Branch lengths represent morphological change. These trees were used in all simulations for this study (Analyses 2-5). 

•	Telluraves (Ksepka et al., 2019): ‘Telluraves_tree_branchlengths_dichot.phy’.

•	Passeriformes a (Steell et al., 2023): ‘Passeriformes_tree_branchlengths_dichot.phy’.

•	Neornithes (Field et al., 2020): ‘Neornithes_tree_branchlengths_dichot.phy’.

•	Avialae a (Benito et al., 2022): ‘Avialae_tree_branchlengths_dichot.phy’.

_Empirical character taxon matrices_

NEXUS files of empirical character taxon matrices used in Analysis 1.

•	Telluraves (Ksepka et al., 2019): ‘Telluraves_matrix.nex’.

•	Passeriformes a (Steell et al., 2023): ‘Passeriformes_matrix_a.nex’.

•	Passeriformes b (Steell et al., 2023): ‘Passeriformes_matrix_b.nex’.

•	Neornithes (Field et al., 2020): ‘Neornithes_matrix.nex’.

•	Avialae a and b (Benito et al., 2022): ‘Avialae_matrix.nex’.

**R SCRIPT FILES**

•	Functions: Script containing all custom functions necessary to carry out below analyses; ‘Functions_SuppInfo.R’.

•	Analysis 1: Comparison of indexes for empirical and random matrices; ‘Analysis_1_v2.R’.

•	Analysis 2: Inflating transition rate for simulated matrices; ‘Analysis_2_v1.R’.

•	Analysis 3: Varying number of taxa for simulated matrices at slow, medium and fast transition rate categories; ‘Analysis_3_u-6.R’, ‘Analysis_3_u-4.R’, ‘Analysis_3_u-2.R’.

•	Analysis 4: Varying number of characters for simulated matrices at slow, medium and fast transition rate categories; ‘Analysis_4_u-6.R’, ‘Analysis_4_u-4.R’, ‘Analysis_4_u-2’.

•	Analysis 5: Phylomorphospace patterns at different homoplasy levels; ‘Analysis_5_u-6.R’, ‘Analysis_5_u-4.R’, ‘Analysis_5_u-2.R’, ‘Analysis_5_u-0’.

**FUNCTION USAGE NOTES**

_HER()_

Function to calculate homoplasy excess ratio (Archie, 1989). Should only be used for topologies that are inferred via maximum parsimony from the input matrix with deviations from this method (e.g., applying partial or complete topological constraints or using implied weighting).

Arguments:

‘data’ - A matrix of morphological data in the phangorn package ‘phydat’ (Schliep, 2011) format. Can use custom function ‘morph.phydat’, see below.

‘tree’ - A phylogenetic tree of class ‘phylo’ (Paradis et al., 2004).

‘n’ - Number of permutations of the matrix.

‘maxit’ - Maximum number of iterations in parsimony ratchet, as in ‘pratchet’ function in phangorn (Schliep, 2011). Default is 500.

‘minit’ - Minimum number of iterations in parsimony ratchet, as in ‘pratchet’ function in phangorn (Schliep, 2011). Default is 50.

‘k’ - Maximum number of rounds before ratchet is stopped when there is no improvement, as in ‘pratchet’ function in phangorn (Schliep, 2011). Default is 10.

_HSR()_

Function to calculate homoplasy slope ratio (Meier et al., 1991). Includes two options: to calculate the HSR as in Meier et al. (1991) or to calculate the modified homoplasy slope ratio (HSRm) as described in the present manuscript. HSRm includes polymorphic states and missing data into the matrix randomisation step.

Arguments:

‘data’ - A matrix of morphological data in the phangorn package ‘phydat’ (Schliep, 2011) format. Can use custom function ‘morph.phydat’, see below.

‘tree’ - A phylogenetic tree of class ‘phylo’ (Paradis et al., 2004).

‘n’ - Number of permutations of the matrix.

‘modified’ – To state whether HSR should be calculated from the original method (Meier et al. 1991; modified = FALSE) or the modified method (modified = TRUE; default).

_RHI()_

Function to calculate the relative homoplasy index described here. 

Arguments:

‘data’ - A matrix of morphological data in the phangorn package ‘phydat’ (Schliep, 2011) format. Can use custom function ‘morph.phydat’, see below.

‘tree’ - A phylogenetic tree of class ‘phylo’ (Paradis et al., 2004).

‘n’ - Number of randomisations of the tip order in the tree.

‘ord’ – A vector containing the numbers for the specific ordered characters in the dataset. Default = NULL.

‘cost’ – A cost matrix specifying the cost per character state change (transition). See custom function ‘cost.matrix’ below. Must be specified if ‘ord’ is specified. Default = NULL.

_RHI.char()_

Function to calculate the per character relative homoplasy using the RHI formula. Currently only available for unordered characters (specified ordered multistate characters are in development).

Arguments:

‘data’ - A matrix of morphological data in the phangorn package ‘phydat’ (Schliep, 2011) format. Can use custom function ‘morph.phydat’, see below.

‘tree’ - A phylogenetic tree of class ‘phylo’ (Paradis et al., 2004).

‘n’ - Number of randomisations of the tip order in the tree.

_RHI.multi()_

Function to calculate RHI for a distribution of trees with the same taxa and the same matrix. Ordered characters are allowed.

Arguments:

‘data’ - A matrix of morphological data in the phangorn package ‘phydat’ (Schliep, 2011) format. Can use custom function ‘morph.phydat’, see below.

‘tree’ - A list of phylogenetic trees of class ‘MultiPhylo’ (Paradis et al., 2004).

‘n’ - Number of randomisations of the tip order in the tree.

‘ord’ – A vector containing the numbers for the specific ordered characters in the dataset. Default = NULL.

‘cost’ – A cost matrix specifying the cost per character state change (transition). See custom function ‘cost.matrix’ below. Must be specified if ‘ord’ is specified. Default = NULL.

_contrast.matrix()_

Function to generate a dataset-specific contrast matrix. This is required to generate a ‘phydat’ object required for homoplasy index calculations (Schliep, 2011). Contrast matrices specify which character scorings should behave as character states or combinations of states (e.g., ambiguities or polymorphisms). Here, polymorphisms and uncertainties are treated the same, and ‘gaps’ (‘-‘) and ambiguities (‘?’) are treated the same. A new contrast matrix should be generated for each phylogenetic character matrix. This function is implemented in custom function ‘morph.phydat’ detailed below.

Arguments:

‘data’ – A character-taxon matrix of class ‘matrix array’. Note that NAs should not be present.

_cost.matrix()_

Function to generate a dataset-specific cost matrix for use with ordered characters.

Arguments:

‘data’ – A character-taxon matrix of class ‘matrix array’. Note that NAs should not be present.

_get.states()_

Function that pastes the character states within a matrix. It ignores non-state scorings (polymorphisms and ambiguities).

Arguments:

‘data’ – A character-taxon matrix of class ‘matrix array’. Note that NAs should not be present.

_morph.phydat()_

Function to generate a ‘phydat’ object (Schliep, 2011) to specifically handle morphological data as opposed to other data types more typically used in the phangorn package (e.g., nucleotide or amino acid data). Makes other phangorn functions compatible with morphological data that includes polymorphisms. Necessary for use with homoplasy index functions in phangorn and detailed here.

Arguments:

‘matrix’ – A character-taxon matrix of class ‘matrix array’. Note that NAs should not be present.

_parenthesize()_

Internal function to place round braces around polymorphic character scorings.

Arguments:

‘x’ – Character string in which to enclose with brackets.

_poly.to.amb()_

Function to replace all polymorphisms and uncertainties with ambiguities (‘?’).

Arguments:

‘matrix’ – A character-taxon matrix of class ‘matrix array’. Note that NAs should not be present.

_statef()_

Function to paste the relative character state frequencies of a matrix. Excludes non-state scorings (polymorphisms and ambiguities).

Arguments:

‘data’ – A character-taxon matrix of class ‘matrix array’. Note that NAs should not be present.

_subset.matrix()_

Function to subset matrices based on a list of subset trees. This function subsets a matrix by removing specific taxa, but does not remove characters.

Arguments:

‘trees’ – A list of subset trees (object class ‘phylo’) from the ‘subset.trees’ function detailed below.

‘matrix’ – A character taxon matrix to be subsetted.

_uncert.to.poly()_

Internal function to change uncertainties denoted with curly braces to polymorphisms (denoted by round braces).

Arguments:

‘matrix’ – A character-taxon matrix of class ‘matrix array’. Note that NAs should not be present.

**SOURCES AND REFERENCES**

Datasets were derived from the following published articles:

Steell EM, Nguyen JM, Benson RB, Field DJ. Comparative anatomy of the passerine carpometacarpus helps illuminate the early fossil record of crown Passeriformes. Journal of Anatomy. 2023 Mar;242(3):495-509.

Field DJ, Benito J, Chen A, Jagt JW, Ksepka DT. Late Cretaceous neornithine from Europe illuminates the origins of crown birds. Nature. 2020 Mar 19;579(7799):397-401.

Benito J, Kuo PC, Widrig KE, Jagt JW, Field DJ. Cretaceous ornithurine supports a neognathous crown bird ancestor. Nature. 2022 Dec 1;612(7938):100-5.

Ksepka DT, Grande L, Mayr G. Oldest finch-beaked birds reveal parallel ecological radiations in the earliest evolution of passerines. Current Biology. 2019 Feb 18;29(4):657-63.

_Code/Software_

All code was run in R version 4.2.3 (2023-03-15) -- "Shortstop Beagle".

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

tictoc v2014.1
Izrailev S. tictoc: Functions for timing R scripts, as well as implementations of Stack and List structures.
