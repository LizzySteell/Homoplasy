##### 2024-02-28

### Relative homoplasy index

### Analysis 5: phylomorphospace of simulated matrices across different rate categories

### Four rate categories: -6, -4, -2, 0

### This script is for rate category: -2

# Packages

packages <- c("ape", "phangorn", "StatMatch", "TreeTools", "phytools", "janitor", "stringr", "geiger", "tictoc", "dispRity", "pracma")

# lapply(packages, install.packages)

lapply(packages, library, character.only=TRUE)

# Don't set seed for these analyses as each matrix needs to be random

# Import data - trees with resolved polytomies and branch lengths are saved from Analysis 2

av_tree <- read.tree("Avialae_tree_branchlengths_dichot.phy")
neo_tree <- read.tree("Neornithes_tree_branchlengths_dichot.phy")
pas_tree <- read.tree("Passeriformes_tree_branchlengths_dichot.phy")
tel_tree <- read.tree("Telluraves_tree_branchlengths_dichot.phy")

# Eliminate zero-length branches by adding 0.01 to each branch length in each tree

av_tree$edge.length <- av_tree$edge.length + 0.01
neo_tree$edge.length <- neo_tree$edge.length + 0.01
pas_tree$edge.length <- pas_tree$edge.length + 0.01
tel_tree$edge.length <- tel_tree$edge.length + 0.01

### Avialae

# tree has 85 taxa, dataset has 282 characters
av_char <- 282

# Simulate 1 matrix per rate category, with the empirical topology (fully resoloved), but 0.9 binary chars and 0.1 3-state chars.

av_matrix <- sim.morpho(av_tree, characters=av_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-2, sdlog=0.1), invariant=FALSE)

# Generate Gower distance matrices per simulated morphological matrix

av_dist <- gower.dist(av_matrix)

# Principal coordinate analyses

av_pco <- pcoa(av_dist)
av_pco$values

av_1 <- av_pco$values[,3][1]
av_2 <- av_pco$values[,3][2]
av_3 <- av_pco$values[,3][3]
av_4 <- av_pco$values[,3][4]
av_5 <- av_pco$values[,3][5]

pdf("Analysis5_meanlog-2_Avialae.pdf")

phylomorphospace(av_tree, av_pco$vectors[,1:2], label="off", node.size=c(0.1,1), pch=16, xlab=av_pco$values[,3][1], ylab=av_pco$values[,3][2])

dev.off()

# Calculate homoplasy metrics for each one
av_phydat <- morph.phydat(av_matrix)

rhi1 <- RHI(av_phydat, av_tree, 100)
# her1 <- HER(av_phydat, av_tree, 100)
ci1 <- CI(av_tree, av_phydat)
ri1 <- RI(av_tree, av_phydat)
hsrm1 <- HSR(av_phydat, av_tree, 100)
hsr1 <- HSR(av_phydat, av_tree, 100, noise=FALSE)

### Neornithes ###

# tree has 39 taxa, dataset has 295 characters
neo_char <- 295

# Simulate 1 matrix per rate category, with the empirical topology (fully resoloved), but 0.9 binary chars and 0.1 3-state chars.

neo_matrix <- sim.morpho(neo_tree, characters=neo_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-2, sdlog=0.1), invariant=FALSE)

# Generate Gower distance matrices per simulated morphological matrix

neo_dist <- gower.dist(neo_matrix)

# Principal coordinate analyses

neo_pco <- pcoa(neo_dist)
neo_pco$values

neo_1 <- neo_pco$values[,3][1]
neo_2 <- neo_pco$values[,3][2]
neo_3 <- neo_pco$values[,3][3]
neo_4 <- neo_pco$values[,3][4]
neo_5 <- neo_pco$values[,3][5]

pdf("Analysis5_meanlog-2_Neornithes.pdf")

phylomorphospace(neo_tree, neo_pco$vectors[,1:2], label="off", node.size=c(0.1,1), pch=16, xlab=neo_pco$values[,3][1], ylab=neo_pco$values[,3][2])

dev.off()

# Calculate homoplasy metrics for each one
neo_phydat <- morph.phydat(neo_matrix)

rhi2 <- RHI(neo_phydat, neo_tree, 100)
# her2 <- HER(neo_phydat, neo_tree, 100)
ci2 <- CI(neo_tree, neo_phydat)
ri2 <- RI(neo_tree, neo_phydat)
hsrm2 <- HSR(neo_phydat, neo_tree, 100)
hsr2 <- HSR(neo_phydat, neo_tree, 100, noise=FALSE)

### Passeriformes ###


pas_char <- 49

# Simulate 1 matrix per rate category, with the empirical topology (fully resoloved), but 0.9 binary chars and 0.1 3-state chars.

pas_matrix <- sim.morpho(pas_tree, characters=pas_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-2, sdlog=0.1), invariant=FALSE)

# Generate Gower distance matrices per simulated morphological matrix

pas_dist <- gower.dist(pas_matrix)

# Principal coordinate analyses

pas_pco <- pcoa(pas_dist)
pas_pco$values

pas_1 <- pas_pco$values[,3][1]
pas_2 <- pas_pco$values[,3][2]
pas_3 <- pas_pco$values[,3][3]
pas_4 <- pas_pco$values[,3][4]
pas_5 <- pas_pco$values[,3][5]

pdf("Analysis5_meanlog-2_Passeriformes.pdf")

phylomorphospace(pas_tree, pas_pco$vectors[,1:2], label="off", node.size=c(0.1,1), pch=16, xlab=pas_pco$values[,3][1], ylab=pas_pco$values[,3][2])

dev.off()

# Calculate homoplasy metrics for each one
pas_phydat <- morph.phydat(pas_matrix)

rhi3 <- RHI(pas_phydat, pas_tree, 100)
# her3 <- HER(pas_phydat, pas_tree, 100)
ci3 <- CI(pas_tree, pas_phydat)
ri3 <- RI(pas_tree, pas_phydat)
hsrm3 <- HSR(pas_phydat, pas_tree, 100)
hsr3 <- HSR(pas_phydat, pas_tree, 100, noise=FALSE)


### Telluraves ###


tel_char <- 146

# Simulate 1 matrix per rate category, with the empirical topology (fully resoloved), but 0.9 binary chars and 0.1 3-state chars.

tel_matrix <- sim.morpho(tel_tree, characters=tel_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-2, sdlog=0.1), invariant=FALSE)

# Generate Gower distance matrices per simulated morphological matrix

tel_dist <- gower.dist(tel_matrix)

# Principal coordinate analyses

tel_pco <- pcoa(tel_dist)
tel_pco$values

tel_1 <- tel_pco$values[,3][1]
tel_2 <- tel_pco$values[,3][2]
tel_3 <- tel_pco$values[,3][3]
tel_4 <- tel_pco$values[,3][4]
tel_5 <- tel_pco$values[,3][5]

pdf("Analysis5_meanlog-2_Telluraves.pdf")

phylomorphospace(tel_tree, tel_pco$vectors[,1:2], label="off", node.size=c(0.1,1), pch=16, xlab=tel_pco$values[,3][1], ylab=tel_pco$values[,3][2])

dev.off()

# Calculate homoplasy metrics for each one
tel_phydat <- morph.phydat(tel_matrix)

rhi4 <- RHI(tel_phydat, tel_tree, 100)
# her4 <- HER(tel_phydat, tel_tree, 100)
ci4 <- CI(tel_tree, tel_phydat)
ri4 <- RI(tel_tree, tel_phydat)
hsrm4 <- HSR(tel_phydat, tel_tree, 100)
hsr4 <- HSR(tel_phydat, tel_tree, 100, noise=FALSE)
