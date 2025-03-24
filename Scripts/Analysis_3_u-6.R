# 2024-02-20

### Relative Homoplasy Index

### Analysis 3: Decreasing taxon number by subsetting trees

### Three rate categories: -6, -4 and -2 meanlog. Run for all four simulated datasets.
# This simulation method generates whole matrices using the sim.morpho() feature from the 'DispRity' package.

#### This script generates matrices with meanlog = -6.

# Packages

packages <- c("ape", "phangorn", "StatMatch", "TreeTools", "phytools", "janitor", "stringr", "geiger", "tictoc", "dispRity", "pracma")

# lapply(packages, install.packages)

lapply(packages, library, character.only=TRUE)

set.seed(100)

# Import data - trees with resolved polytomies and branch lengths are saved from Analysis 2

av_tree <- read.tree("Avialae_tree_branchlengths_dichot.phy")
neo_tree <- read.tree("Neornithes_tree_branchlengths_dichot.phy")
pas_tree <- read.tree("Passeriformes_tree_branchlengths_dichot.phy")
tel_tree <- read.tree("Telluraves_tree_branchlengths_dichot.phy")

#################

### Avialae ###

# tree has 85 taxa, dataset has 282 characters
av_char <- 282
av_t <- 85

# Subset from 10 taxa to 82 taxa by 8
# Generate 100 subsetted trees per taxon number
# 900 trees in total, 900 matrices in total (per rate category)

av_n_subsets <- seq(10, 74, by=8)

# Subset trees:

av_sub_trees <- c()
for(i in 1:length(av_n_subsets)){
	av_sub_trees[[i]] <- subset.trees(100, av_tree, av_n_subsets[i])
}

#Simulate 1 matrix per subsetted tree per rate category, with the same number of characters as empirical matrix, but 0.9 binary chars and 0.1 3-state chars.

av_matrices <- c()
for(i in 1:length(av_n_subsets)){
	av_matrices[[i]] <- lapply(av_sub_trees[[i]], sim.morpho, characters=av_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-6, sdlog=0.1), invariant=FALSE)
	print(i)
	}

# Generate 100 matrices at that the transition rate for the full set of taxa
av_matrices_ <- replicate(100, sim.morpho(av_tree, characters=av_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-6, sdlog=0.1), invariant=FALSE))

# Convert to list

sim <- list()
	for(i in 1:dim(av_matrices_)[3]){
		sim[[i]] <- av_matrices_[,,i]
		}
		
av_phydats_ <- lapply(sim, morph.phydat)

# Convert matrices to phydat objects
av_phydats <- c()
for(i in 1:length(av_matrices)){
	av_phydats[[i]] <- lapply(av_matrices[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
av_cis <- c()
for(j in 1:length(av_sub_trees)){
	av_cis[[j]] <- mapply(CI, av_sub_trees[[j]], av_phydats[[j]])
}

#Make numeric

for(i in 1:length(av_cis)){
	av_cis[[i]] <- as.numeric(av_cis[[i]])
}

#Retention index
av_ris <- c()
for(j in 1:length(av_sub_trees)){
	av_ris[[j]] <- mapply(RI, av_sub_trees[[j]], av_phydats[[j]])
}

#Make numeric

for(i in 1:length(av_ris)){
	av_ris[[i]] <- as.numeric(av_ris[[i]])
}

# Relative homoplasy index
# n = 100
tic()
av_rhis <- c()
for(j in 1:length(av_sub_trees)){
	av_rhis[[j]] <- mapply(RHI, av_phydats[[j]], av_sub_trees[[j]], n=100)
}
toc()

#Extract RHI values
av_rhi_vals <- c()

for(j in 1:length(av_rhis)){
	
		av_rhi_vals[[j]] <- av_rhis[[j]][2,]
	}

# Get medians and quantiles

rhi_med <- as.numeric(lapply(av_rhi_vals, quantile, prob=0.5))
rhi_05 <- as.numeric(lapply(av_rhi_vals, quantile, prob=0.05))
rhi_95 <- as.numeric(lapply(av_rhi_vals, quantile, prob=0.95))

ri_med <- as.numeric(lapply(av_ris, quantile, prob=0.5))
ri_05 <- as.numeric(lapply(av_ris, quantile, prob=0.05))
ri_95 <- as.numeric(lapply(av_ris, quantile, prob=0.95))

ci_med <- as.numeric(lapply(av_cis, quantile, prob=0.5))
ci_05 <- as.numeric(lapply(av_cis, quantile, prob=0.05))
ci_95 <- as.numeric(lapply(av_cis, quantile, prob=0.95))

# For full topology

ci_ <- c()
for(x in 1:length(av_phydats_)){
	ci_[[x]] <- CI(av_tree, av_phydats_[[x]])
}
ci_ <- unlist(lapply(ci_, as.numeric))
ci_ <- quantile(ci_, prob=0.5)[[1]]
ci_05_ <- quantile(ci_, prob=0.05)[[1]]
ci_95_ <- quantile(ci_, prob=0.95)[[1]]

ri_ <- c()
for(x in 1:length(av_phydats_)){
	ri_[[x]] <- RI(av_tree, av_phydats_[[x]])
}
ri_ <- unlist(lapply(ri_, as.numeric))
ri_ <- quantile(ri_, prob=0.5)[[1]]
ri_05_ <- quantile(ri_, prob=0.05)[[1]]
ri_95_ <- quantile(ri_, prob=0.95)[[1]]

rhi_ <- lapply(av_phydats_, RHI, av_tree, 100)
rhi__ <- c()
for(j in 1:length(rhi_)){
		rhi__[[j]] <- rhi_[[j]][2]
	}
rhi_ <- quantile(as.numeric(unlist(rhi__)), prob=0.5)[[1]]
rhi_05_ <- quantile(as.numeric(unlist(rhi__)), prob=0.05)[[1]]
rhi_95_ <- quantile(as.numeric(unlist(rhi__)), prob=0.95)[[1]]

# Plot

pdf(file="Analysis3_meanlog-6_Avialae_2024-02-20.pdf")
plot(xlim=c(10,av_t), ri_med, ylim=c(0,1.1), col="white")
lines(av_n_subsets, 1-ri_med, col="blue")
lines(av_n_subsets, 1-ri_05, col="light blue")
lines(av_n_subsets, 1-ri_95, col="light blue")
lines(av_n_subsets, 1-ci_med, col="red")
lines(av_n_subsets, 1-ci_05, col="pink")
lines(av_n_subsets, 1-ci_95, col="pink")
lines(av_n_subsets, rhi_med, col="black")
lines(av_n_subsets, rhi_05, col="gray")
lines(av_n_subsets, rhi_95, col="gray")
points(av_t, 1-ci_, col="red")
points(av_t, 1-ci_05_, col="pink")
points(av_t, 1-ci_95_, col="pink")
points(av_t, 1-ri_, col="blue")
points(av_t, 1-ri_05_, col="light blue")
points(av_t, 1-ri_95_, col="light blue")
points(av_t, rhi_, col="black")
points(av_t, rhi_05_, col="gray")
points(av_t, rhi_95_, col="gray")
dev.off()

# Y axis range differences:

a1 <- max(1-ci_med) - min(1-ci_med)

b1 <- max(1-ri_med) - min(1-ri_med)

c1 <- max(1-rhi_med) - min(1-rhi_med)

#####

### Neornithes ###

# tree has 39 taxa, dataset has 295 characters
neo_char <- 295
neo_t <- 39

# Subset from 10 taxa to 37 taxa by 3
# Generate 100 subsetted trees per taxon number
# 900 trees in total, 900 matrices in total (per rate category)

neo_n_subsets <- seq(10, 37, by=3)

# Subset trees:

neo_sub_trees <- c()
for(i in 1:length(neo_n_subsets)){
	neo_sub_trees[[i]] <- subset.trees(100, neo_tree, neo_n_subsets[i])
}

#Simulate 1 matrix per subsetted tree per rate category, with the same number of characters as empirical matrix, but 0.9 binary chars and 0.1 3-state chars.

neo_matrices <- c()
for(i in 1:length(neo_n_subsets)){
	neo_matrices[[i]] <- lapply(neo_sub_trees[[i]], sim.morpho, characters=neo_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-6, sdlog=0.1), invariant=FALSE)
	print(i)
	}
	
# Generate 100 matrices at that the transition rate for the full set of taxa
neo_matrices_ <- replicate(100, sim.morpho(neo_tree, characters=neo_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-6, sdlog=0.1), invariant=FALSE))

# Convert to list

sim <- list()
	for(i in 1:dim(neo_matrices_)[3]){
		sim[[i]] <- neo_matrices_[,,i]
		}
		
neo_phydats_ <- lapply(sim, morph.phydat)

# Convert matrices to phydat objects
neo_phydats <- c()
for(i in 1:length(neo_matrices)){
	neo_phydats[[i]] <- lapply(neo_matrices[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
neo_cis <- c()
for(j in 1:length(neo_sub_trees)){
	neo_cis[[j]] <- mapply(CI, neo_sub_trees[[j]], neo_phydats[[j]])
}

#Make numeric

for(i in 1:length(neo_cis)){
	neo_cis[[i]] <- as.numeric(neo_cis[[i]])
}

#Retention index
neo_ris <- c()
for(j in 1:length(neo_sub_trees)){
	neo_ris[[j]] <- mapply(RI, neo_sub_trees[[j]], neo_phydats[[j]])
}

#Make numeric

for(i in 1:length(neo_ris)){
	neo_ris[[i]] <- as.numeric(neo_ris[[i]])
}

# Relative homoplasy index
# n = 100
tic()
neo_rhis <- c()
for(j in 1:length(neo_sub_trees)){
	neo_rhis[[j]] <- mapply(RHI, neo_phydats[[j]], neo_sub_trees[[j]], n=100)
}
toc()

#Extract RHI values
neo_rhi_vals <- c()

for(j in 1:length(neo_rhis)){
	
		neo_rhi_vals[[j]] <- neo_rhis[[j]][2,]
	}

# Get medians and quantiles

rhi_med <- as.numeric(lapply(neo_rhi_vals, quantile, prob=0.5))
rhi_05 <- as.numeric(lapply(neo_rhi_vals, quantile, prob=0.05))
rhi_95 <- as.numeric(lapply(neo_rhi_vals, quantile, prob=0.95))

ri_med <- as.numeric(lapply(neo_ris, quantile, prob=0.5))
ri_05 <- as.numeric(lapply(neo_ris, quantile, prob=0.05))
ri_95 <- as.numeric(lapply(neo_ris, quantile, prob=0.95))

ci_med <- as.numeric(lapply(neo_cis, quantile, prob=0.5))
ci_05 <- as.numeric(lapply(neo_cis, quantile, prob=0.05))
ci_95 <- as.numeric(lapply(neo_cis, quantile, prob=0.95))

# For full topology

ci_ <- c()
for(x in 1:length(neo_phydats_)){
	ci_[[x]] <- CI(neo_tree, neo_phydats_[[x]])
}
ci_ <- unlist(lapply(ci_, as.numeric))
ci_ <- quantile(ci_, prob=0.5)[[1]]
ci_05_ <- quantile(ci_, prob=0.05)[[1]]
ci_95_ <- quantile(ci_, prob=0.95)[[1]]

ri_ <- c()
for(x in 1:length(neo_phydats_)){
	ri_[[x]] <- RI(neo_tree, neo_phydats_[[x]])
}
ri_ <- unlist(lapply(ri_, as.numeric))
ri_ <- quantile(ri_, prob=0.5)[[1]]
ri_05_ <- quantile(ri_, prob=0.05)[[1]]
ri_95_ <- quantile(ri_, prob=0.95)[[1]]

rhi_ <- lapply(neo_phydats_, RHI, neo_tree, 100)
rhi__ <- c()
for(j in 1:length(rhi_)){
		rhi__[[j]] <- rhi_[[j]][2]
	}
rhi_ <- quantile(as.numeric(unlist(rhi__)), prob=0.5)[[1]]
rhi_05_ <- quantile(as.numeric(unlist(rhi__)), prob=0.05)[[1]]
rhi_95_ <- quantile(as.numeric(unlist(rhi__)), prob=0.95)[[1]]

# Plot

pdf(file="Analysis3_meanlog-6_Neornithes_2024-02-20.pdf")
plot(xlim=c(10,neo_t), ri_med, ylim=c(0,1.1), col="white")
lines(neo_n_subsets, 1-ri_med, col="blue")
lines(neo_n_subsets, 1-ri_05, col="light blue")
lines(neo_n_subsets, 1-ri_95, col="light blue")
lines(neo_n_subsets, 1-ci_med, col="red")
lines(neo_n_subsets, 1-ci_05, col="pink")
lines(neo_n_subsets, 1-ci_95, col="pink")
lines(neo_n_subsets, rhi_med, col="black")
lines(neo_n_subsets, rhi_05, col="gray")
lines(neo_n_subsets, rhi_95, col="gray")
points(neo_t, 1-ci_, col="red")
points(neo_t, 1-ci_05_, col="pink")
points(neo_t, 1-ci_95_, col="pink")
points(neo_t, 1-ri_, col="blue")
points(neo_t, 1-ri_05_, col="light blue")
points(neo_t, 1-ri_95_, col="light blue")
points(neo_t, rhi_, col="black")
points(neo_t, rhi_05_, col="gray")
points(neo_t, rhi_95_, col="gray")
dev.off()

# Y axis range differences:

a2 <- max(1-ci_med) - min(1-ci_med)

b2 <- max(1-ri_med) - min(1-ri_med)

c2 <- max(1-rhi_med) - min(1-rhi_med)

#####

### Passeriformes ###

# tree has 143 taxa, dataset has 49 characters
pas_char <- 49
pas_t <- 143

# Subset from 10 taxa to 135 taxa by 15
# Generate 100 subsetted trees per taxon number
# 900 trees in total, 900 matrices in total (per rate category)

pas_n_subsets <- seq(10, 135, by=15)

# Subset trees:

pas_sub_trees <- c()
for(i in 1:length(pas_n_subsets)){
	pas_sub_trees[[i]] <- subset.trees(100, pas_tree, pas_n_subsets[i])
}

#Simulate 1 matrix per subsetted tree per rate category, with the same number of characters as empirical matrix, but 0.9 binary chars and 0.1 3-state chars.

pas_matrices <- c()
for(i in 1:length(pas_n_subsets)){
	pas_matrices[[i]] <- lapply(pas_sub_trees[[i]], sim.morpho, characters=pas_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-6, sdlog=0.1), invariant=FALSE)
	print(i)
	}
	
# Generate 100 matrices at that the transition rate for the full set of taxa
pas_matrices_ <- replicate(100, sim.morpho(pas_tree, characters=pas_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-6, sdlog=0.1), invariant=FALSE))

# Convert to list

sim <- list()
	for(i in 1:dim(pas_matrices_)[3]){
		sim[[i]] <- pas_matrices_[,,i]
		}
		
pas_phydats_ <- lapply(sim, morph.phydat)

# Convert matrices to phydat objects
pas_phydats <- c()
for(i in 1:length(pas_matrices)){
	pas_phydats[[i]] <- lapply(pas_matrices[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
pas_cis <- c()
for(j in 1:length(pas_sub_trees)){
	pas_cis[[j]] <- mapply(CI, pas_sub_trees[[j]], pas_phydats[[j]])
}

#Make numeric

for(i in 1:length(pas_cis)){
	pas_cis[[i]] <- as.numeric(pas_cis[[i]])
}

#Retention index
pas_ris <- c()
for(j in 1:length(pas_sub_trees)){
	pas_ris[[j]] <- mapply(RI, pas_sub_trees[[j]], pas_phydats[[j]])
}

#Make numeric

for(i in 1:length(pas_ris)){
	pas_ris[[i]] <- as.numeric(pas_ris[[i]])
}

# Relative homoplasy index
# n = 100
tic()
pas_rhis <- c()
for(j in 1:length(pas_sub_trees)){
	pas_rhis[[j]] <- mapply(RHI, pas_phydats[[j]], pas_sub_trees[[j]], n=100)
}
toc()

#Extract RHI values
pas_rhi_vals <- c()

for(j in 1:length(pas_rhis)){
	
		pas_rhi_vals[[j]] <- pas_rhis[[j]][2,]
	}

# Get medians and quantiles

rhi_med <- as.numeric(lapply(pas_rhi_vals, quantile, prob=0.5))
rhi_05 <- as.numeric(lapply(pas_rhi_vals, quantile, prob=0.05))
rhi_95 <- as.numeric(lapply(pas_rhi_vals, quantile, prob=0.95))

ri_med <- as.numeric(lapply(pas_ris, quantile, prob=0.5))
ri_05 <- as.numeric(lapply(pas_ris, quantile, prob=0.05))
ri_95 <- as.numeric(lapply(pas_ris, quantile, prob=0.95))

ci_med <- as.numeric(lapply(pas_cis, quantile, prob=0.5))
ci_05 <- as.numeric(lapply(pas_cis, quantile, prob=0.05))
ci_95 <- as.numeric(lapply(pas_cis, quantile, prob=0.95))

# For full topology

ci_ <- c()
for(x in 1:length(pas_phydats_)){
	ci_[[x]] <- CI(pas_tree, pas_phydats_[[x]])
}
ci_ <- unlist(lapply(ci_, as.numeric))
ci_ <- quantile(ci_, prob=0.5)[[1]]
ci_05_ <- quantile(ci_, prob=0.05)[[1]]
ci_95_ <- quantile(ci_, prob=0.95)[[1]]

ri_ <- c()
for(x in 1:length(pas_phydats_)){
	ri_[[x]] <- RI(pas_tree, pas_phydats_[[x]])
}
ri_ <- unlist(lapply(ri_, as.numeric))
ri_ <- quantile(ri_, prob=0.5)[[1]]
ri_05_ <- quantile(ri_, prob=0.05)[[1]]
ri_95_ <- quantile(ri_, prob=0.95)[[1]]

rhi_ <- lapply(pas_phydats_, RHI, pas_tree, 100)
rhi__ <- c()
for(j in 1:length(rhi_)){
		rhi__[[j]] <- rhi_[[j]][2]
	}
rhi_ <- quantile(as.numeric(unlist(rhi__)), prob=0.5)[[1]]
rhi_05_ <- quantile(as.numeric(unlist(rhi__)), prob=0.05)[[1]]
rhi_95_ <- quantile(as.numeric(unlist(rhi__)), prob=0.95)[[1]]

# Plot

pdf(file="Analysis3_meanlog-6_Passeriformes_2024-02-20.pdf")
plot(xlim=c(10,pas_t), ri_med, ylim=c(0,1.1), col="white")
lines(pas_n_subsets, 1-ri_med, col="blue")
lines(pas_n_subsets, 1-ri_05, col="light blue")
lines(pas_n_subsets, 1-ri_95, col="light blue")
lines(pas_n_subsets, 1-ci_med, col="red")
lines(pas_n_subsets, 1-ci_05, col="pink")
lines(pas_n_subsets, 1-ci_95, col="pink")
lines(pas_n_subsets, rhi_med, col="black")
lines(pas_n_subsets, rhi_05, col="gray")
lines(pas_n_subsets, rhi_95, col="gray")
points(pas_t, 1-ci_, col="red")
points(pas_t, 1-ci_05_, col="pink")
points(pas_t, 1-ci_95_, col="pink")
points(pas_t, 1-ri_, col="blue")
points(pas_t, 1-ri_05_, col="light blue")
points(pas_t, 1-ri_95_, col="light blue")
points(pas_t, rhi_, col="black")
points(pas_t, rhi_05_, col="gray")
points(pas_t, rhi_95_, col="gray")
dev.off()

# Y axis range differences:

a3 <- max(1-ci_med) - min(1-ci_med)

b3 <- max(1-ri_med) - min(1-ri_med)

c3 <- max(1-rhi_med) - min(1-rhi_med)

#####

### Telluraves ###

# tree has 62 taxa, dataset has 146 characters
tel_char <- 146
tel_t <- 62

# Subset from 10 taxa to 55 taxa by 5
# Generate 100 subsetted trees per taxon number
# 1000 trees in total, 1000 matrices in total (per rate category)

tel_n_subsets <- seq(10, 55, by=5)

# Subset trees:

tel_sub_trees <- c()
for(i in 1:length(tel_n_subsets)){
	tel_sub_trees[[i]] <- subset.trees(100, tel_tree, tel_n_subsets[i])
}

#Simulate 1 matrix per subsetted tree per rate category, with the same number of characters as empirical matrix, but 0.9 binary chars and 0.1 3-state chars.

tel_matrices <- c()
for(i in 1:length(tel_n_subsets)){
	tel_matrices[[i]] <- lapply(tel_sub_trees[[i]], sim.morpho, characters=tel_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-6, sdlog=0.1), invariant=FALSE)
	print(i)
	}
	
# Generate 100 matrices at that the transition rate for the full set of taxa
tel_matrices_ <- replicate(100, sim.morpho(tel_tree, characters=tel_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-6, sdlog=0.1), invariant=FALSE))

# Convert to list

sim <- list()
	for(i in 1:dim(tel_matrices_)[3]){
		sim[[i]] <- tel_matrices_[,,i]
		}
		
tel_phydats_ <- lapply(sim, morph.phydat)

# Convert matrices to phydat objects
tel_phydats <- c()
for(i in 1:length(tel_matrices)){
	tel_phydats[[i]] <- lapply(tel_matrices[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
tel_cis <- c()
for(j in 1:length(tel_sub_trees)){
	tel_cis[[j]] <- mapply(CI, tel_sub_trees[[j]], tel_phydats[[j]])
}

#Make numeric

for(i in 1:length(tel_cis)){
	tel_cis[[i]] <- as.numeric(tel_cis[[i]])
}

#Retention index
tel_ris <- c()
for(j in 1:length(tel_sub_trees)){
	tel_ris[[j]] <- mapply(RI, tel_sub_trees[[j]], tel_phydats[[j]])
}

#Make numeric

for(i in 1:length(tel_ris)){
	tel_ris[[i]] <- as.numeric(tel_ris[[i]])
}

# Relative homoplasy index
# n = 100
tic()
tel_rhis <- c()
for(j in 1:length(tel_sub_trees)){
	tel_rhis[[j]] <- mapply(RHI, tel_phydats[[j]], tel_sub_trees[[j]], n=100)
}
toc()

#Extract RHI values
tel_rhi_vals <- c()

for(j in 1:length(tel_rhis)){
	
		tel_rhi_vals[[j]] <- tel_rhis[[j]][2,]
	}

# Get medians and quantiles

rhi_med <- as.numeric(lapply(tel_rhi_vals, quantile, prob=0.5))
rhi_05 <- as.numeric(lapply(tel_rhi_vals, quantile, prob=0.05))
rhi_95 <- as.numeric(lapply(tel_rhi_vals, quantile, prob=0.95))

ri_med <- as.numeric(lapply(tel_ris, quantile, prob=0.5))
ri_05 <- as.numeric(lapply(tel_ris, quantile, prob=0.05))
ri_95 <- as.numeric(lapply(tel_ris, quantile, prob=0.95))

ci_med <- as.numeric(lapply(tel_cis, quantile, prob=0.5))
ci_05 <- as.numeric(lapply(tel_cis, quantile, prob=0.05))
ci_95 <- as.numeric(lapply(tel_cis, quantile, prob=0.95))

# For full topology

ci_ <- c()
for(x in 1:length(tel_phydats_)){
	ci_[[x]] <- CI(tel_tree, tel_phydats_[[x]])
}
ci_ <- unlist(lapply(ci_, as.numeric))
ci_ <- quantile(ci_, prob=0.5)[[1]]
ci_05_ <- quantile(ci_, prob=0.05)[[1]]
ci_95_ <- quantile(ci_, prob=0.95)[[1]]

ri_ <- c()
for(x in 1:length(tel_phydats_)){
	ri_[[x]] <- RI(tel_tree, tel_phydats_[[x]])
}
ri_ <- unlist(lapply(ri_, as.numeric))
ri_ <- quantile(ri_, prob=0.5)[[1]]
ri_05_ <- quantile(ri_, prob=0.05)[[1]]
ri_95_ <- quantile(ri_, prob=0.95)[[1]]

rhi_ <- lapply(tel_phydats_, RHI, tel_tree, 100)
rhi__ <- c()
for(j in 1:length(rhi_)){
		rhi__[[j]] <- rhi_[[j]][2]
	}
rhi_ <- quantile(as.numeric(unlist(rhi__)), prob=0.5)[[1]]
rhi_05_ <- quantile(as.numeric(unlist(rhi__)), prob=0.05)[[1]]
rhi_95_ <- quantile(as.numeric(unlist(rhi__)), prob=0.95)[[1]]

# Plot

pdf(file="Analysis3_meanlog-6_Telluraves_2024-02-20.pdf")
plot(xlim=c(10,tel_t), ri_med, ylim=c(0,1.1), col="white")
lines(tel_n_subsets, 1-ri_med, col="blue")
lines(tel_n_subsets, 1-ri_05, col="light blue")
lines(tel_n_subsets, 1-ri_95, col="light blue")
lines(tel_n_subsets, 1-ci_med, col="red")
lines(tel_n_subsets, 1-ci_05, col="pink")
lines(tel_n_subsets, 1-ci_95, col="pink")
lines(tel_n_subsets, rhi_med, col="black")
lines(tel_n_subsets, rhi_05, col="gray")
lines(tel_n_subsets, rhi_95, col="gray")
points(tel_t, 1-ci_, col="red")
points(tel_t, 1-ci_05_, col="pink")
points(tel_t, 1-ci_95_, col="pink")
points(tel_t, 1-ri_, col="blue")
points(tel_t, 1-ri_05_, col="light blue")
points(tel_t, 1-ri_95_, col="light blue")
points(tel_t, rhi_, col="black")
points(tel_t, rhi_05_, col="gray")
points(tel_t, rhi_95_, col="gray")
dev.off()

# Y axis range differences:

a4 <- max(1-ci_med) - min(1-ci_med)

b4 <- max(1-ri_med) - min(1-ri_med)

c4 <- max(1-rhi_med) - min(1-rhi_med)

#####



