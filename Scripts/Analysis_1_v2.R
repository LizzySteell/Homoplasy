# 2024-02-05

### Relative Homoplasy Index

### Analysis 1a: Comparison between methods for empirical datasets and trees

# Packages

packages <- c("ape", "phangorn", "StatMatch", "TreeTools", "phytools", "janitor", "stringr", "geiger", "tictoc")

lapply(packages, install.packages)

lapply(packages, library, character.only=TRUE)

set.seed(100)

# Import data

Avial_matrix <- ReadCharacters("Avialae_matrix.nex")
Avial_a <- read.tree("Avialae_a_tree.nex")
Avial_b <- read.tree("Avialae_b_tree.nex")

Neorn_matrix <- ReadCharacters("Neornithes_matrix.nex")
Neorn_tree <- read.tree("Neornithes_tree.nex")

Passe_matrix_a <- ReadCharacters("Passeriformes_matrix_a.nex")
Passe_matrix_b <- ReadCharacters("Passeriformes_matrix_b.nex")
Passe_a <- read.tree("Passeriformes_a_tree.phy")
Passe_b <- read.tree("Passeriformes_b_tree.phy")

Tell_matrix <- ReadCharacters("Telluraves_matrix.nex")
Tell_tree <- read.tree("Telluraves_tree.phy")

# Remove non-varying characters and mismatching taxa

#### Avialae

for(i in 1:ncol(Avial_matrix)){
	if (sd(Avial_matrix[,i], na.rm=TRUE) == 0){
		print(i)
	}
} #43, 101, 194 need to be removed

Avial_matrix <- Avial_matrix[,-c(43, 101, 194)]

# Name check

name.check(Avial_a, Avial_matrix)
name.check(Avial_b, Avial_matrix)
# Data not tree - "Apsaravis_ukhaana"
Avial_matrix <- Avial_matrix[row.names(Avial_matrix) != "Apsaravis_ukhaana",]

#### Neornithes

for(i in 1:ncol(Neorn_matrix)){
	if (sd(Neorn_matrix[,i], na.rm=TRUE) == 0){
		print(i)
	}
} #201, 282 need to be removed

Neorn_matrix <- Neorn_matrix[,-c(201, 282)]

# Name check

name.check(Neorn_tree, Neorn_matrix) # OK

#### Passeriformes

# Name check

name.check(Passe_a, Passe_matrix_a) # OK
name.check(Passe_b, Passe_matrix_b) # OK

#### Telluraves

for(i in 1:ncol(Tell_matrix)){
	if (sd(Tell_matrix[,i], na.rm=TRUE) == 0){
		print(i)
	}
}

# Name check

name.check(Tell_tree, Tell_matrix) # OK

########

# Make phydats of each matrix

av_phydat <- morph.phydat(Avial_matrix)
neorn_phydat <- morph.phydat(Neorn_matrix)
passe_a_phydat <- morph.phydat(Passe_matrix_a)
passe_b_phydat <- morph.phydat(Passe_matrix_b)
tell_phydat <- morph.phydat(Tell_matrix)

# Calculate RHI, L (tree length) and Lmin (theoretical minimum tree length), H (extra steps for whole matrix), ES (extra steps per character)

### Avialae
tic()
rhi_av_a <- RHI(av_phydat, Avial_a, 100) # RHI 0.4107649
toc()
av_a_L <- rhi_av_a[4] #1407
av_a_Lmin <- rhi_av_a[5] #392
av_a_H <- av_a_L[[1]]-av_a_Lmin[[1]] #1015

tic()
av_a_ci <- CI(Avial_a, av_phydat)
toc()
tic()
av_a_ri <- RI(Avial_a, av_phydat)
toc()

tic()
av_a_hsr <- HSR(av_phydat, Avial_a, 100, noise=FALSE) #0.1138785
toc()
tic()
av_a_hsrm <- HSR(av_phydat, Avial_a, 100) #0.4069665
toc()

tic()
av_a_her <- HER(av_phydat, Avial_a, 100)
toc()

# Avial b

tic()
rhi_av_b <- RHI(av_phydat, Avial_b, 1000) # RHI 0.4696853
toc()

tic()
av_b_L <- rhi_av_b[4] #1616
toc()
av_b_Lmin <- rhi_av_b[5] #392
av_b_H <- av_b_L[[1]]-av_b_Lmin[[1]] #1224

tic()
av_b_ci <- CI(Avial_b, av_phydat)
toc()
tic()
av_b_ri <- RI(Avial_b, av_phydat)
toc()

tic()
av_b_hsr <- HSR(av_phydat, Avial_b, 100) #0.4614826
toc()
tic()
av_b_hsrm <- HSR(av_phydat, Avial_b, 100, noise=FALSE) #0.1256474
toc()

tic()
av_b_her <- HER(av_phydat, Abial_b, 100)
toc()

### Neornithes
tic()
rhi_neorn <- RHI(neorn_phydat, Neorn_tree, 100) # RHI 0.5039225
toc()
neorn_L <- rhi_neorn[4] #1498
neorn_Lmin <- rhi_neorn[5] #406
neorn_H <- neorn_L[[1]]-neorn_Lmin[[1]] #1092

tic()
neorn_ci <- CI(Neorn_tree, neorn_phydat)
toc()
tic()
neorn_ri <- RI(Neorn_tree, neorn_phydat)
toc()

tic()
neorn_hsr <- HSR(neorn_phydat, Neorn_tree, 100) #0.5034043
toc()
tic()
neorn_hsrm <- HSR(neorn_phydat, Neorn_tree, 100, noise=FALSE) #0.2758129
toc()

tic()
norn_her <- HER(neorn_phydat, Neorn_tree, 100)
toc()

### Passeriformes

# dataset a

tic()
rhi_passe_a <- RHI(passe_a_phydat, Passe_a, 100) # RHI 0.6177632
toc() #0.328 sec elapsed
passe_a_L <- rhi_passe_a[4] #995
passe_a_Lmin <- rhi_passe_a[5] #56
passe_a_H <- passe_a_L[[1]]-passe_a_Lmin[[1]] #939
passe_a_ES <- passe_a_H/dim(Passe_matrix_a)[2] #19.16327

tic()
ci_passe_a <- CI(Passe_a, passe_a_phydat) #0.05628141
toc()
tic()
ri_passe_a <- RI(Passe_a, passe_a_phydat) #0.5063091
toc()

tic()
HSR(passe_a_phydat, Passe_a, 100) #0.6195172
toc()
tic()
HSR(passe_a_phydat, Passe_a, 100, noise=FALSE) #0.3901193
toc()

# dataset b

tic()
rhi_passe_b <- RHI(passe_b_phydat, Passe_b, 1000)
toc()
 
passe_b_L <- rhi_passe_b[4]
passe_b_Lmin <- rhi_passe_b[5]
passe_b_H <- passe_b_L[[1]]-passe_b_Lmin[[1]] 
passe_b_ES <- passe_b_H/dim(Passe_matrix_b)[2] 

tic()
ci_passe_b <- CI(Passe_b, passe_b_phydat)
toc()
tic()
ri_passe_b <- RI(Passe_b, passe_b_phydat)
toc()

tic()
HSR(passe_b_phydat, Passe_b, 100) 
toc()
tic()
HSR(passe_b_phydat, Passe_b, 100, noise=FALSE) 
toc()


# Telluraves

tic()
CI(Tell_tree, tell_phydat)
toc()

tic()
RI(Tell_tree, tell_phydat)
toc()

tic()
RHI(tell_phydat, Tell_tree, 1000)
toc()

tic()
HSR(tell_phydat, Tell_tree, 100) #0.4436284
toc()

tic()
HSR(tell_phydat, Tell_tree, 100, noise=FALSE) #0.1636464
toc()

### HERs 

tic()
HER(av_phydat, Avial_a, 100)
toc()

tic()
HER(av_phydat, Avial_b, 100)
toc()

tic()
HER(neorn_phydat, Neorn_tree, 100)
toc()

tic()
HER(passe_b_phydat, Passe_b, 100) 
toc()

tic()
HER(passe_a_phydat, Passe_a, 100) 
toc()

tic()
HER(tell_phydat, Tell_tree, 100)
toc()

# Analysis 1b: generate random matrices for each dataset with same numbers of characters and taxa, and same proportions of binary and multistate characters. Then infer MP trees for each random matrix using pratchet(). Then repeat analyses (CI, RI, HER, HSR, HSRm, RHI)

# function to populate matrices based on empirical states

populate <- function(matrix, states){	
for(k in 1:length(states)){
	matrix[,k] <- sample(states[[k]], size=dim(matrix)[1], replace=TRUE)
}
return(matrix)
}

# Avialae 

Avial_matrixR <- matrix(nrow=dim(Avial_matrix)[1], ncol=dim(Avial_matrix)[2])

nchar <- dim(Avial_matrix)[2]

	states <- apply(Avial_matrix, 2, get.states, simplify=FALSE)

Avial_matrixR <- populate(Avial_matrixR, states)

rows <- c(Avial_a$tip.label)

row.names(Avial_matrixR) <- rows

table(Avial_matrixR)
		
# get parsimony topology from matrix

Avial_phydatR <- morph.phydat(Avial_matrixR)
Avial_Rtree <- pratchet(Avial_phydatR)
plot(Avial_Rtree, cex=0.5)

# Neornithes

# remove uncertanties 

Neorn_matrix_ <- uncert.to.poly(Neorn_matrix)

Neorn_matrixR <- matrix(nrow=dim(Neorn_matrix_)[1], ncol=dim(Neorn_matrix_)[2])

nchar <- dim(Neorn_matrix_)[2]

	states <- apply(Neorn_matrix_, 2, get.states, simplify=FALSE)

Neorn_matrixR <- populate(Neorn_matrixR, states)

rows <- c(Neorn_tree$tip.label)

row.names(Neorn_matrixR) <- rows

table(Neorn_matrixR)
		
# get parsimony topology from matrix

Neorn_phydatR <- morph.phydat(Neorn_matrixR)
Neorn_Rtree <- pratchet(Neorn_phydatR)
plot(Neorn_Rtree, cex=0.5)

# Passeriformes (only do matrix a)

Passe_matrixR <- matrix(nrow=dim(Passe_matrix_a)[1], ncol=dim(Passe_matrix_a)[2])

nchar <- dim(Passe_matrix_a)[2]

	states <- apply(Passe_matrix_a, 2, get.states, simplify=FALSE)

Passe_matrixR <- populate(Passe_matrixR, states)

rows <- c(Passe_a$tip.label)

row.names(Passe_matrixR) <- rows

table(Passe_matrixR)
		
# get parsimony topology from matrix

Passe_phydatR <- morph.phydat(Passe_matrixR)
Passe_Rtree <- pratchet(Passe_phydatR)
plot(Passe_Rtree, cex=0.5)

# Telluraves

# remove uncertainites from matrix

Tell_matrix_ <- uncert.to.poly(Tell_matrix)

Tell_matrixR <- matrix(nrow=dim(Tell_matrix_)[1], ncol=dim(Tell_matrix_)[2])

nchar <- dim(Tell_matrix_)[2]

	states <- apply(Tell_matrix_, 2, get.states, simplify=FALSE)

Tell_matrixR <- populate(Tell_matrixR, states)

rows <- c(Tell_tree$tip.label)

row.names(Tell_matrixR) <- rows

table(Tell_matrixR)
		
# get parsimony topology from matrix

Tell_phydatR <- morph.phydat(Tell_matrixR)
Tell_Rtree <- pratchet(Tell_phydatR)
plot(Tell_Rtree, cex=0.5)


# Homplasy calculations for randomly generated matrices and their resulting minimum length Maximum Parsimony trees

### Avialae
RHI(Avial_phydatR, Avial_Rtree, 100)

CI(Avial_Rtree, Avial_phydatR)

RI(Avial_Rtree, Avial_phydatR)

HSR(Avial_phydatR, Avial_Rtree, 100, noise=FALSE)
HSR(Avial_phydatR, Avial_Rtree, 100)

### Neornithes
RHI(Neorn_phydatR, Neorn_Rtree, 100)

CI(Neorn_Rtree, Neorn_phydatR)

RI(Neorn_Rtree, Neorn_phydatR)

HSR(Neorn_phydatR, Neorn_Rtree, 100, noise=FALSE)
HSR(Neorn_phydatR, Neorn_Rtree, 100)

### Passeriformes
RHI(Passe_phydatR, Passe_Rtree, 100)

CI(Passe_Rtree, Passe_phydatR)

RI(Passe_Rtree, Passe_phydatR)

HSR(Passe_phydatR, Passe_Rtree, 100, noise=FALSE)
HSR(Passe_phydatR, Passe_Rtree, 100)

### Telluraves
RHI(Tell_phydatR, Tell_Rtree, 100)

CI(Tell_Rtree, Tell_phydatR)

RI(Tell_Rtree, Tell_phydatR)

HSR(Tell_phydatR, Tell_Rtree, 100, noise=FALSE)
HSR(Tell_phydatR, Tell_Rtree, 100)


### HERs 

a <- HER(Avial_phydatR, Avial_Rtree, 100)

b <- HER(Neorn_phydatR, Neorn_Rtree, 100)

c <- HER(Passe_phydatR, Passe_Rtree, 100)

d <- HER(Tell_phydatR, Tell_Rtree, 100)


