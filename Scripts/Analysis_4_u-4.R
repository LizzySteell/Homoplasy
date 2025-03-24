# 2024-02-26

### Relative Homoplasy Index

### Analysis 4: Decreasing character number by simulating matrices with decreasing characters

### Three rate categories: -6, -4 and -2 meanlog. Run for all four simulated datasets.
# This simulation method generates whole matrices using the sim.morpho() feature from the 'DispRity' package.

#### This script generates matrices with meanlog = -4.

# Packages

packages <- c("ape", "phangorn", "StatMatch", "TreeTools", "phytools", "janitor",  "stringr", "geiger", "tictoc", "dispRity", "pracma")

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

# Subset from 10 chars to 280 chars by 30
# Generate 100 matrices per character number category
# 10 character categories, 1000 matrices in total (per rate category)

av_chars <- seq(10, 280, by=30)

#Simulate 100 matrices per character number category per rate category, with the empirical topology (fully resoloved), but 0.9 binary chars and 0.1 3-state chars.

av_matrices <- c()
for(i in 1:length(av_chars)){
	av_matrices[[i]] <- replicate(100, sim.morpho(av_tree, characters=av_chars[[i]], states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-4, sdlog=0.1), invariant=FALSE), simplify=FALSE)
	print(i)
	}

# Generate 100 matrices at that the transition rate for the full set of characters
av_matrices_ <- replicate(100, sim.morpho(av_tree, characters=av_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-4, sdlog=0.1), invariant=FALSE))

# Convert to list - the full set of chars

sim <- list()
	for(i in 1:dim(av_matrices_)[3]){
		sim[[i]] <- av_matrices_[,,i]
		}
		
av_phydats_ <- lapply(sim, morph.phydat)

# Convert decreasing chars matrices to phydat objects
av_phydats <- c()
for(i in 1:length(av_matrices)){
	av_phydats[[i]] <- lapply(av_matrices[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
av_cis <- c()
av_cis_ <- c()
for(i in 1:length(av_phydats)){
	for(j in 1:length(av_phydats[[i]])){
	av_cis_[[j]]	 <- CI(av_tree, av_phydats[[i]][[j]])
	av_cis[[i]] <- av_cis_
	}
		}	

#Make numeric

for(i in 1:length(av_cis)){
	av_cis[[i]] <- as.numeric(av_cis[[i]])
}

#Retention index
av_ris <- c()
av_ris_ <- c()
for(i in 1:length(av_phydats)){
	for(j in 1:length(av_phydats[[i]])){
	av_ris_[[j]]	 <- RI(av_tree, av_phydats[[i]][[j]])
	av_ris[[i]] <- av_ris_
	}
		}

#Make numeric

for(i in 1:length(av_ris)){
	av_ris[[i]] <- as.numeric(av_ris[[i]])
}

# Relative homoplasy index
# n = 100

av_rhis <- c()
for(j in 1:length(av_phydats)){
	av_rhis[[j]] <- lapply(av_phydats[[j]], RHI, av_tree, n=100)
}


#Extract RHI values
av_rhi_vals <- c()
av_rhi_vals_ <- c()
for(i in 1:length(av_rhis)){
	for(j in 1:length(av_rhis[[i]])){
		av_rhi_vals_[[j]] <- av_rhis[[i]][[j]][[2]]
		av_rhi_vals[[i]] <- av_rhi_vals_
		}
	}	

for(i in 1:length(av_rhi_vals)){
	av_rhi_vals[[i]] <- as.numeric(av_rhi_vals[[i]])
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

pdf(file="Analysis4_meanlog-4_Avialae_2024-02-26.pdf")
plot(xlim=c(10,av_char), ri_med, ylim=c(0,1.1), col="white")
lines(av_chars, 1-ri_med, col="blue")
lines(av_chars, 1-ri_05, col="light blue")
lines(av_chars, 1-ri_95, col="light blue")
lines(av_chars, 1-ci_med, col="red")
lines(av_chars, 1-ci_05, col="pink")
lines(av_chars, 1-ci_95, col="pink")
lines(av_chars, rhi_med, col="black")
lines(av_chars, rhi_05, col="gray")
lines(av_chars, rhi_95, col="gray")
points(av_char, 1-ci_, col="red")
points(av_char, 1-ci_05_, col="pink")
points(av_char, 1-ci_95_, col="pink")
points(av_char, 1-ri_, col="blue")
points(av_char, 1-ri_05_, col="light blue")
points(av_char, 1-ri_95_, col="light blue")
points(av_char, rhi_, col="black")
points(av_char, rhi_05_, col="gray")
points(av_char, rhi_95_, col="gray")
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

# Subset from 10 chars to 280 chars by 30
# Generate 100 matrices per character number category
# 10 character categories, 1000 matrices in total (per rate category)

neo_chars <- seq(10, 280, by=30)

#Simulate 100 matrices per character number category per rate category, with the empirical topology (fully resoloved), but 0.9 binary chars and 0.1 3-state chars.

neo_matrices <- c()
for(i in 1:length(neo_chars)){
	neo_matrices[[i]] <- replicate(100, sim.morpho(neo_tree, characters=neo_chars[[i]], states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-4, sdlog=0.1), invariant=FALSE), simplify=FALSE)
	print(i)
	}

# Generate 100 matrices at that the transition rate for the full set of characters
neo_matrices_ <- replicate(100, sim.morpho(neo_tree, characters=neo_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-4, sdlog=0.1), invariant=FALSE))

# Convert to list - the full set of chars

sim <- list()
	for(i in 1:dim(neo_matrices_)[3]){
		sim[[i]] <- neo_matrices_[,,i]
		}
		
neo_phydats_ <- lapply(sim, morph.phydat)

# Convert decreasing chars matrices to phydat objects
neo_phydats <- c()
for(i in 1:length(neo_matrices)){
	neo_phydats[[i]] <- lapply(neo_matrices[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
neo_cis <- c()
neo_cis_ <- c()
for(i in 1:length(neo_phydats)){
	for(j in 1:length(neo_phydats[[i]])){
	neo_cis_[[j]]	 <- CI(neo_tree, neo_phydats[[i]][[j]])
	neo_cis[[i]] <- neo_cis_
	}
		}	

#Make numeric

for(i in 1:length(neo_cis)){
	neo_cis[[i]] <- as.numeric(neo_cis[[i]])
}

#Retention index
neo_ris <- c()
neo_ris_ <- c()
for(i in 1:length(neo_phydats)){
	for(j in 1:length(neo_phydats[[i]])){
	neo_ris_[[j]]	 <- RI(neo_tree, neo_phydats[[i]][[j]])
	neo_ris[[i]] <- neo_ris_
	}
		}

#Make numeric

for(i in 1:length(neo_ris)){
	neo_ris[[i]] <- as.numeric(neo_ris[[i]])
}

# Relative homoplasy index
# n = 100

neo_rhis <- c()
for(j in 1:length(neo_phydats)){
	neo_rhis[[j]] <- lapply(neo_phydats[[j]], RHI, neo_tree, n=100)
}


#Extract RHI values
neo_rhi_vals <- c()
neo_rhi_vals_ <- c()
for(i in 1:length(neo_rhis)){
	for(j in 1:length(neo_rhis[[i]])){
		neo_rhi_vals_[[j]] <- neo_rhis[[i]][[j]][[2]]
		neo_rhi_vals[[i]] <- neo_rhi_vals_
		}
	}	

for(i in 1:length(neo_rhi_vals)){
	neo_rhi_vals[[i]] <- as.numeric(neo_rhi_vals[[i]])
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

pdf(file="Analysis4_meanlog-4_Neornithes_2024-02-26.pdf")
plot(xlim=c(10,neo_char), ri_med, ylim=c(0,1.1), col="white")
lines(neo_chars, 1-ri_med, col="blue")
lines(neo_chars, 1-ri_05, col="light blue")
lines(neo_chars, 1-ri_95, col="light blue")
lines(neo_chars, 1-ci_med, col="red")
lines(neo_chars, 1-ci_05, col="pink")
lines(neo_chars, 1-ci_95, col="pink")
lines(neo_chars, rhi_med, col="black")
lines(neo_chars, rhi_05, col="gray")
lines(neo_chars, rhi_95, col="gray")
points(neo_char, 1-ci_, col="red")
points(neo_char, 1-ci_05_, col="pink")
points(neo_char, 1-ci_95_, col="pink")
points(neo_char, 1-ri_, col="blue")
points(neo_char, 1-ri_05_, col="light blue")
points(neo_char, 1-ri_95_, col="light blue")
points(neo_char, rhi_, col="black")
points(neo_char, rhi_05_, col="gray")
points(neo_char, rhi_95_, col="gray")
dev.off()

# Y axis range differences:

a2 <- max(1-ci_med) - min(1-ci_med)

b2 <- max(1-ri_med) - min(1-ri_med)

c2 <- max(1-rhi_med) - min(1-rhi_med)



### Passeriformes ###

# tree has 143 taxa, dataset has 49 characters
pas_char <- 49
pas_t <- 143

# Subset from 10 chars to 45 chars by 5
# Generate 100 matrices per character number category
# 8 character categories, 800 matrices in total (per rate category)

pas_chars <- seq(10, 45, by=5)

#Simulate 100 matrices per character number category per rate category, with the empirical topology (fully resoloved), but 0.9 binary chars and 0.1 3-state chars.

pas_matrices <- c()
for(i in 1:length(pas_chars)){
	pas_matrices[[i]] <- replicate(100, sim.morpho(pas_tree, characters=pas_chars[[i]], states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-4, sdlog=0.1), invariant=FALSE), simplify=FALSE)
	print(i)
	}

# Generate 100 matrices at that the transition rate for the full set of characters
pas_matrices_ <- replicate(100, sim.morpho(pas_tree, characters=pas_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-4, sdlog=0.1), invariant=FALSE))

# Convert to list - the full set of chars

sim <- list()
	for(i in 1:dim(pas_matrices_)[3]){
		sim[[i]] <- pas_matrices_[,,i]
		}
		
pas_phydats_ <- lapply(sim, morph.phydat)

# Convert decreasing chars matrices to phydat objects
pas_phydats <- c()
for(i in 1:length(pas_matrices)){
	pas_phydats[[i]] <- lapply(pas_matrices[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
pas_cis <- c()
pas_cis_ <- c()
for(i in 1:length(pas_phydats)){
	for(j in 1:length(pas_phydats[[i]])){
	pas_cis_[[j]]	 <- CI(pas_tree, pas_phydats[[i]][[j]])
	pas_cis[[i]] <- pas_cis_
	}
		}	

#Make numeric

for(i in 1:length(pas_cis)){
	pas_cis[[i]] <- as.numeric(pas_cis[[i]])
}

#Retention index
pas_ris <- c()
pas_ris_ <- c()
for(i in 1:length(pas_phydats)){
	for(j in 1:length(pas_phydats[[i]])){
	pas_ris_[[j]]	 <- RI(pas_tree, pas_phydats[[i]][[j]])
	pas_ris[[i]] <- pas_ris_
	}
		}

#Make numeric

for(i in 1:length(pas_ris)){
	pas_ris[[i]] <- as.numeric(pas_ris[[i]])
}

# Relative homoplasy index
# n = 100

pas_rhis <- c()
for(j in 1:length(pas_phydats)){
	pas_rhis[[j]] <- lapply(pas_phydats[[j]], RHI, pas_tree, n=100)
}


#Extract RHI values
pas_rhi_vals <- c()
pas_rhi_vals_ <- c()
for(i in 1:length(pas_rhis)){
	for(j in 1:length(pas_rhis[[i]])){
		pas_rhi_vals_[[j]] <- pas_rhis[[i]][[j]][[2]]
		pas_rhi_vals[[i]] <- pas_rhi_vals_
		}
	}	

for(i in 1:length(pas_rhi_vals)){
	pas_rhi_vals[[i]] <- as.numeric(pas_rhi_vals[[i]])
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

pdf(file="Analysis4_meanlog-4_Passeriformes_2024-02-26.pdf")
plot(xlim=c(10,pas_char), ri_med, ylim=c(0,1.1), col="white")
lines(pas_chars, 1-ri_med, col="blue")
lines(pas_chars, 1-ri_05, col="light blue")
lines(pas_chars, 1-ri_95, col="light blue")
lines(pas_chars, 1-ci_med, col="red")
lines(pas_chars, 1-ci_05, col="pink")
lines(pas_chars, 1-ci_95, col="pink")
lines(pas_chars, rhi_med, col="black")
lines(pas_chars, rhi_05, col="gray")
lines(pas_chars, rhi_95, col="gray")
points(pas_char, 1-ci_, col="red")
points(pas_char, 1-ci_05_, col="pink")
points(pas_char, 1-ci_95_, col="pink")
points(pas_char, 1-ri_, col="blue")
points(pas_char, 1-ri_05_, col="light blue")
points(pas_char, 1-ri_95_, col="light blue")
points(pas_char, rhi_, col="black")
points(pas_char, rhi_05_, col="gray")
points(pas_char, rhi_95_, col="gray")
dev.off()

# Y axis range differences:

a3 <- max(1-ci_med) - min(1-ci_med)

b3 <- max(1-ri_med) - min(1-ri_med)

c3 <- max(1-rhi_med) - min(1-rhi_med)

### Telluraves ###

# tree has 62 taxa, dataset has 146 characters
tel_char <- 146
tel_t <- 62

# Subset from 10 chars to 130 chars by 15
# Generate 100 matrices per character number category
# 9 character categories, 900 matrices in total (per rate category)

tel_chars <- seq(10, 130, by=15)

#Simulate 100 matrices per character number category per rate category, with the empirical topology (fully resoloved), but 0.9 binary chars and 0.1 3-state chars.

tel_matrices <- c()
for(i in 1:length(tel_chars)){
	tel_matrices[[i]] <- replicate(100, sim.morpho(tel_tree, characters=tel_chars[[i]], states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-4, sdlog=0.1), invariant=FALSE), simplify=FALSE)
	print(i)
	}

# Generate 100 matrices at that the transition rate for the full set of characters
tel_matrices_ <- replicate(100, sim.morpho(tel_tree, characters=tel_char, states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog=-4, sdlog=0.1), invariant=FALSE))

# Convert to list - the full set of chars

sim <- list()
	for(i in 1:dim(tel_matrices_)[3]){
		sim[[i]] <- tel_matrices_[,,i]
		}
		
tel_phydats_ <- lapply(sim, morph.phydat)

# Convert decreasing chars matrices to phydat objects
tel_phydats <- c()
for(i in 1:length(tel_matrices)){
	tel_phydats[[i]] <- lapply(tel_matrices[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
tel_cis <- c()
tel_cis_ <- c()
for(i in 1:length(tel_phydats)){
	for(j in 1:length(tel_phydats[[i]])){
	tel_cis_[[j]]	 <- CI(tel_tree, tel_phydats[[i]][[j]])
	tel_cis[[i]] <- tel_cis_
	}
		}	

#Make numeric

for(i in 1:length(tel_cis)){
	tel_cis[[i]] <- as.numeric(tel_cis[[i]])
}

#Retention index
tel_ris <- c()
tel_ris_ <- c()
for(i in 1:length(tel_phydats)){
	for(j in 1:length(tel_phydats[[i]])){
	tel_ris_[[j]]	 <- RI(tel_tree, tel_phydats[[i]][[j]])
	tel_ris[[i]] <- tel_ris_
	}
		}

#Make numeric

for(i in 1:length(tel_ris)){
	tel_ris[[i]] <- as.numeric(tel_ris[[i]])
}

# Relative homoplasy index
# n = 100

tel_rhis <- c()
for(j in 1:length(tel_phydats)){
	tel_rhis[[j]] <- lapply(tel_phydats[[j]], RHI, tel_tree, n=100)
}


#Extract RHI values
tel_rhi_vals <- c()
tel_rhi_vals_ <- c()
for(i in 1:length(tel_rhis)){
	for(j in 1:length(tel_rhis[[i]])){
		tel_rhi_vals_[[j]] <- tel_rhis[[i]][[j]][[2]]
		tel_rhi_vals[[i]] <- tel_rhi_vals_
		}
	}	

for(i in 1:length(tel_rhi_vals)){
	tel_rhi_vals[[i]] <- as.numeric(tel_rhi_vals[[i]])
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

pdf(file="Analysis4_meanlog-4_Telluraves_2024-02-26.pdf")
plot(xlim=c(10,tel_char), ri_med, ylim=c(0,1.1), col="white")
lines(tel_chars, 1-ri_med, col="blue")
lines(tel_chars, 1-ri_05, col="light blue")
lines(tel_chars, 1-ri_95, col="light blue")
lines(tel_chars, 1-ci_med, col="red")
lines(tel_chars, 1-ci_05, col="pink")
lines(tel_chars, 1-ci_95, col="pink")
lines(tel_chars, rhi_med, col="black")
lines(tel_chars, rhi_05, col="gray")
lines(tel_chars, rhi_95, col="gray")
points(tel_char, 1-ci_, col="red")
points(tel_char, 1-ci_05_, col="pink")
points(tel_char, 1-ci_95_, col="pink")
points(tel_char, 1-ri_, col="blue")
points(tel_char, 1-ri_05_, col="light blue")
points(tel_char, 1-ri_95_, col="light blue")
points(tel_char, rhi_, col="black")
points(tel_char, rhi_05_, col="gray")
points(tel_char, rhi_95_, col="gray")
dev.off()

# Y axis range differences:

a4 <- max(1-ci_med) - min(1-ci_med)

b4 <- max(1-ri_med) - min(1-ri_med)

c4 <- max(1-rhi_med) - min(1-rhi_med)