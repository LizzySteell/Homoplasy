# 2024-02-19

### Relative Homoplasy Index

### Analysis 2: Inflating transition rate for simulated matrices

# Packages

packages <- c("ape", "phangorn", "StatMatch", "TreeTools", "phytools", "janitor", "stringr", "geiger", "tictoc", "dispRity", "pracma")

# lapply(packages, install.packages)

lapply(packages, library, character.only=TRUE)

set.seed(100)

# Import data

Avial_matrix <- ReadCharacters("Avialae_matrix.nex")
Avial_a <- read.tree("Avialae_a_tree.nex")

Neorn_matrix <- ReadCharacters("Neornithes_matrix.nex")
Neorn_tree <- read.tree("Neornithes_tree.nex")

Passe_matrix_a <- ReadCharacters("Passeriformes_matrix_a.nex")
Passe_a <- read.tree("Passeriformes_a_tree.phy")

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
tell_phydat <- morph.phydat(Tell_matrix)

#################

# Inflating transition rate using sim.morpho from package 'dispRity'

#Set the meanlog values to increase transition rate. 21 rates.
meanlog <- seq(-8, 2, by=0.5)

####

### Avialae: using tree from Avialae a ###

# Use acctran to populate tree with branch lengths based on character state changes in the matrix. Also randomly resolves polytomies.

av_tree <- acctran(Avial_a, av_phydat) 

# Save topology for Analyses 3 & 4

write.tree(av_tree, file="Avialae_tree_branchlengths_dichot.phy")

#Simulate 100 matrices per meanlog value with the empirical topology and same number of characters, but 0.9 binary chars and 0.1 3-state chars. This simulation method generates whole matrices using the sim.morpho() feature from the 'DispRity' package.

av_matrices <- c()
for(i in 1:length(meanlog)){
	av_matrices[[i]] <- replicate(100, sim.morpho(av_tree, characters=as.numeric(dim(Avial_matrix)[2]), states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog= meanlog[i], sdlog=0.1), invariant=FALSE))
	print(meanlog[i])	
}

#Convert simulated matrices from an array to a list
av_sims <- c()
for(j in 1:length(av_matrices)){

	sim <- list()
	for(i in 1:dim(av_matrices[[j]])[3]){
		sim[[i]] <- av_matrices[[j]][,,i]
		}
		av_sims[[j]] <- sim
}

# Convert matrices to phydat objects
av_phydats <- c()
for(i in 1:length(av_sims)){
	av_phydats[[i]] <- lapply(av_sims[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
ci <- c()
av_ci_all <- c()
for(i in 1:length(av_phydats)){
	for(j in 1:length(av_phydats[[i]])){
		ci[[j]] <- CI(av_tree, av_phydats[[i]][[j]]) #Calculate CI for each individual phydat
	}
	av_ci_all[[i]] <- ci #Make object of all CI values for all phydats
}

#Make numeric and calculate median etc

for(i in 1:length(av_ci_all)){
	av_ci_all[[i]] <- as.numeric(av_ci_all[[i]])
}

#Retention index
ri <- c()
av_ri_all <- c()
for(i in 1:length(av_phydats)){
	for(j in 1:length(av_phydats[[i]])){
		ri[[j]] <- RI(av_tree, av_phydats[[i]][[j]]) #Calculate RI for each individual phydat
	}
	av_ri_all[[i]] <- ri #Make object of all RI values for all phydats
}

for(i in 1:length(av_ri_all)){
	av_ri_all[[i]] <- as.numeric(av_ri_all[[i]])
}

# Relative homoplasy index
# n = 100
tic()
av_rhi <- c()
for(i in 1:length(av_sims)){
	av_rhi[[i]] <- lapply(av_phydats[[i]], RHI, av_tree, 100)
	print(i/length(av_sims)*100)
}
toc()

#Extract RHI values
av_rhi_val <- c()
av_rhi_vals <- c()

for(j in 1:length(av_rhi)){
	for(k in 1:length(av_rhi[[j]])){
		av_rhi_val[[k]] <- av_rhi[[j]][[k]][[2]]
	}
	av_rhi_vals[[j]] <- av_rhi_val
}

for(i in 1:length(av_rhi_vals)){
	av_rhi_vals[[i]] <- as.numeric(av_rhi_vals[[i]])
}

# L (tree length) - normalised H is identical to normalised L

av_L_val <- c()
av_L_vals <- c()
for(j in 1:length(av_rhi)){
	for(k in 1:length(av_rhi[[j]])){
		av_L_val[[k]] <- av_rhi[[j]][[k]][[4]]
	}
	av_L_vals[[j]] <- av_L_val
}

for(i in 1:length(av_L_vals)){
	av_L_vals[[i]] <- as.numeric(av_L_vals[[i]])
}

# Get median and quantiles

L_med <- as.numeric(lapply(av_L_vals, quantile, prob=0.5))
L_05 <- as.numeric(lapply(av_L_vals, quantile, prob=0.05))
L_95 <- as.numeric(lapply(av_L_vals, quantile, prob=0.95))

#Normalise values between 0 and 1
L_01 <- (L_med-min(L_med))/(max(L_med)-min(L_med))
L_01_05 <- (L_05-min(L_05))/(max(L_05)-min(L_05))
L_01_95 <- (L_95-min(L_95))/(max(L_95)-min(L_95))

rhi_med <- as.numeric(lapply(av_rhi_vals, quantile, prob=0.5))
rhi_05 <- as.numeric(lapply(av_rhi_vals, quantile, prob=0.05))
rhi_95 <- as.numeric(lapply(av_rhi_vals, quantile, prob=0.95))

ri_med <- as.numeric(lapply(av_ri_all, quantile, prob=0.5))
ri_05 <- as.numeric(lapply(av_ri_all, quantile, prob=0.05))
ri_95 <- as.numeric(lapply(av_ri_all, quantile, prob=0.95))

ci_med <- as.numeric(lapply(av_ci_all, quantile, prob=0.5))
ci_05 <- as.numeric(lapply(av_ci_all, quantile, prob=0.05))
ci_95 <- as.numeric(lapply(av_ci_all, quantile, prob=0.95))

# Plot

pdf(file="Analysis2_Avialae_2024-02-19.pdf")
plot(meanlog, ri_med, ylim=c(0,1.1), col="white")
lines(meanlog, 1-ri_med, col="blue")
lines(meanlog, 1-ri_05, col="light blue")
lines(meanlog, 1-ri_95, col="light blue")
lines(meanlog, 1-ci_med, col="red")
lines(meanlog, 1-ci_05, col="pink")
lines(meanlog, 1-ci_95, col="pink")
lines(meanlog, rhi_med, col="black")
lines(meanlog, rhi_05, col="gray")
lines(meanlog, rhi_95, col="gray")
lines(meanlog, L_01, col="dark green")
dev.off()

# Y axis range values:

min(1-ci_med)
max(1-ci_med)

min(1-ri_med)
max(1-ri_med)

min(rhi_med)
max(rhi_med)

# Calculate area under curve

trapz(meanlog, L_01)
trapz(meanlog, 1-ci_med)
trapz(meanlog, 1-ri_med)
trapz(meanlog, rhi_med)

#####

### Neornithes ###

# Use acctran to populate tree with branch lengths based on character state changes in the matrix. Also randomly resolves polytomies.

neo_tree <- acctran(Neorn_tree, neorn_phydat) 

# Save topology for Analyses 3 & 4

write.tree(neo_tree, file="Neornithes_tree_branchlengths_dichot.phy")

#Simulate 100 matrices per meanlog value with the empirical topology and same number of characters, but 0.9 binary chars and 0.1 3-state chars. This simulation method generates whole matrices using the sim.morpho() feature from the 'DispRity' package.

neo_matrices <- c()
for(i in 1:length(meanlog)){
	neo_matrices[[i]] <- replicate(100, sim.morpho(neo_tree, characters=as.numeric(dim(Neorn_matrix)[2]), states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog= meanlog[i], sdlog=0.1), invariant=FALSE))
	print(meanlog[i])	
}

#Convert simulated matrices from an array to a list
neo_sims <- c()
for(j in 1:length(neo_matrices)){

	sim <- list()
	for(i in 1:dim(neo_matrices[[j]])[3]){
		sim[[i]] <- neo_matrices[[j]][,,i]
		}
		neo_sims[[j]] <- sim
}

# Convert matrices to phydat objects
neo_phydats <- c()
for(i in 1:length(neo_sims)){
	neo_phydats[[i]] <- lapply(neo_sims[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
ci <- c()
neo_ci_all <- c()
for(i in 1:length(neo_phydats)){
	for(j in 1:length(neo_phydats[[i]])){
		ci[[j]] <- CI(neo_tree, neo_phydats[[i]][[j]]) #Calculate CI for each individual phydat
	}
	neo_ci_all[[i]] <- ci #Make object of all CI values for all phydats
}

#Make numeric and calculate median etc

for(i in 1:length(neo_ci_all)){
	neo_ci_all[[i]] <- as.numeric(neo_ci_all[[i]])
}

#Retention index
ri <- c()
neo_ri_all <- c()
for(i in 1:length(neo_phydats)){
	for(j in 1:length(neo_phydats[[i]])){
		ri[[j]] <- RI(neo_tree, neo_phydats[[i]][[j]]) #Calculate RI for each individual phydat
	}
	neo_ri_all[[i]] <- ri #Make object of all RI values for all phydats
}

for(i in 1:length(neo_ri_all)){
	neo_ri_all[[i]] <- as.numeric(neo_ri_all[[i]])
}

# Relative homoplasy index
# n = 100
tic()
neo_rhi <- c()
for(i in 1:length(neo_sims)){
	neo_rhi[[i]] <- lapply(neo_phydats[[i]], RHI, neo_tree, 100)
	print(i/length(neo_sims)*100)
}
toc()

#Extract RHI values
neo_rhi_val <- c()
neo_rhi_vals <- c()

for(j in 1:length(neo_rhi)){
	for(k in 1:length(neo_rhi[[j]])){
		neo_rhi_val[[k]] <- neo_rhi[[j]][[k]][[2]]
	}
	neo_rhi_vals[[j]] <- neo_rhi_val
}

for(i in 1:length(neo_rhi_vals)){
	neo_rhi_vals[[i]] <- as.numeric(neo_rhi_vals[[i]])
}

# L (tree length) - normalised H is identical to normalised L

neo_L_val <- c()
neo_L_vals <- c()
for(j in 1:length(neo_rhi)){
	for(k in 1:length(neo_rhi[[j]])){
		neo_L_val[[k]] <- neo_rhi[[j]][[k]][[4]]
	}
	neo_L_vals[[j]] <- neo_L_val
}

for(i in 1:length(neo_L_vals)){
	neo_L_vals[[i]] <- as.numeric(neo_L_vals[[i]])
}

# Get median and quantiles

L_med <- as.numeric(lapply(neo_L_vals, quantile, prob=0.5))
L_05 <- as.numeric(lapply(neo_L_vals, quantile, prob=0.05))
L_95 <- as.numeric(lapply(neo_L_vals, quantile, prob=0.95))

#Normalise values between 0 and 1
L_01 <- (L_med-min(L_med))/(max(L_med)-min(L_med))
L_01_05 <- (L_05-min(L_05))/(max(L_05)-min(L_05))
L_01_95 <- (L_95-min(L_95))/(max(L_95)-min(L_95))

rhi_med <- as.numeric(lapply(neo_rhi_vals, quantile, prob=0.5))
rhi_05 <- as.numeric(lapply(neo_rhi_vals, quantile, prob=0.05))
rhi_95 <- as.numeric(lapply(neo_rhi_vals, quantile, prob=0.95))

ri_med <- as.numeric(lapply(neo_ri_all, quantile, prob=0.5))
ri_05 <- as.numeric(lapply(neo_ri_all, quantile, prob=0.05))
ri_95 <- as.numeric(lapply(neo_ri_all, quantile, prob=0.95))

ci_med <- as.numeric(lapply(neo_ci_all, quantile, prob=0.5))
ci_05 <- as.numeric(lapply(neo_ci_all, quantile, prob=0.05))
ci_95 <- as.numeric(lapply(neo_ci_all, quantile, prob=0.95))

# Plot

pdf(file="Analysis2_Neornithes_2024-02-19.pdf")
plot(meanlog, L_01, ylim=c(0,1.1), col="white")
lines(meanlog, 1-ri_med, col="blue")
lines(meanlog, 1-ri_05, col="light blue")
lines(meanlog, 1-ri_95, col="light blue")
lines(meanlog, 1-ci_med, col="red")
lines(meanlog, 1-ci_05, col="pink")
lines(meanlog, 1-ci_95, col="pink")
lines(meanlog, rhi_med, col="black")
lines(meanlog, rhi_05, col="gray")
lines(meanlog, rhi_95, col="gray")
lines(meanlog, L_01, col="dark green")
dev.off()

# Y axis range values:

min(1-ci_med)
max(1-ci_med)

min(1-ri_med)
max(1-ri_med)

min(rhi_med)
max(rhi_med)

# Calculate area under curve

trapz(meanlog, L_01)
trapz(meanlog, 1-ci_med)
trapz(meanlog, 1-ri_med)
trapz(meanlog, rhi_med)

#####

### Passeriformes: using tree from Passeriformes a ###

# Use acctran to populate tree with branch lengths based on character state changes in the matrix. Also randomly resolves polytomies.

pas_tree <- acctran(Passe_a, passe_a_phydat) 

# Save topology for Analyses 3 & 4

write.tree(pas_tree, file="Passeriformes_tree_branchlengths_dichot.phy")

#Simulate 100 matrices per meanlog value with the empirical topology and same number of characters, but 0.9 binary chars and 0.1 3-state chars. This simulation method generates whole matrices using the sim.morpho() feature from the 'DispRity' package.

pas_matrices <- c()
for(i in 1:length(meanlog)){
	pas_matrices[[i]] <- replicate(100, sim.morpho(pas_tree, characters=as.numeric(dim(Passe_matrix_a)[2]), states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog= meanlog[i], sdlog=0.1), invariant=FALSE))
	print(meanlog[i])	
}

#Convert simulated matrices from an array to a list
pas_sims <- c()
for(j in 1:length(pas_matrices)){

	sim <- list()
	for(i in 1:dim(pas_matrices[[j]])[3]){
		sim[[i]] <- pas_matrices[[j]][,,i]
		}
		pas_sims[[j]] <- sim
}

# Convert matrices to phydat objects
pas_phydats <- c()
for(i in 1:length(pas_sims)){
	pas_phydats[[i]] <- lapply(pas_sims[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
ci <- c()
pas_ci_all <- c()
for(i in 1:length(pas_phydats)){
	for(j in 1:length(pas_phydats[[i]])){
		ci[[j]] <- CI(pas_tree, pas_phydats[[i]][[j]]) #Calculate CI for each individual phydat
	}
	pas_ci_all[[i]] <- ci #Make object of all CI values for all phydats
}

#Make numeric and calculate median etc

for(i in 1:length(pas_ci_all)){
	pas_ci_all[[i]] <- as.numeric(pas_ci_all[[i]])
}

#Retention index
ri <- c()
pas_ri_all <- c()
for(i in 1:length(pas_phydats)){
	for(j in 1:length(pas_phydats[[i]])){
		ri[[j]] <- RI(pas_tree, pas_phydats[[i]][[j]]) #Calculate RI for each individual phydat
	}
	pas_ri_all[[i]] <- ri #Make object of all RI values for all phydats
}

for(i in 1:length(pas_ri_all)){
	pas_ri_all[[i]] <- as.numeric(pas_ri_all[[i]])
}

# Relative homoplasy index
# n = 100
tic()
pas_rhi <- c()
for(i in 1:length(pas_sims)){
	pas_rhi[[i]] <- lapply(pas_phydats[[i]], RHI, pas_tree, 100)
	print(i/length(pas_sims)*100)
}
toc()

#Extract RHI values
pas_rhi_val <- c()
pas_rhi_vals <- c()

for(j in 1:length(pas_rhi)){
	for(k in 1:length(pas_rhi[[j]])){
		pas_rhi_val[[k]] <- pas_rhi[[j]][[k]][[2]]
	}
	pas_rhi_vals[[j]] <- pas_rhi_val
}

for(i in 1:length(pas_rhi_vals)){
	pas_rhi_vals[[i]] <- as.numeric(pas_rhi_vals[[i]])
}

# L (tree length) - normalised H is identical to normalised L

pas_L_val <- c()
pas_L_vals <- c()
for(j in 1:length(pas_rhi)){
	for(k in 1:length(pas_rhi[[j]])){
		pas_L_val[[k]] <- pas_rhi[[j]][[k]][[4]]
	}
	pas_L_vals[[j]] <- pas_L_val
}

for(i in 1:length(pas_L_vals)){
	pas_L_vals[[i]] <- as.numeric(pas_L_vals[[i]])
}

# Get median and quantiles

L_med <- as.numeric(lapply(pas_L_vals, quantile, prob=0.5))
L_05 <- as.numeric(lapply(pas_L_vals, quantile, prob=0.05))
L_95 <- as.numeric(lapply(pas_L_vals, quantile, prob=0.95))

#Normalise values between 0 and 1
L_01 <- (L_med-min(L_med))/(max(L_med)-min(L_med))
L_01_05 <- (L_05-min(L_05))/(max(L_05)-min(L_05))
L_01_95 <- (L_95-min(L_95))/(max(L_95)-min(L_95))

rhi_med <- as.numeric(lapply(pas_rhi_vals, quantile, prob=0.5))
rhi_05 <- as.numeric(lapply(pas_rhi_vals, quantile, prob=0.05))
rhi_95 <- as.numeric(lapply(pas_rhi_vals, quantile, prob=0.95))

ri_med <- as.numeric(lapply(pas_ri_all, quantile, prob=0.5))
ri_05 <- as.numeric(lapply(pas_ri_all, quantile, prob=0.05))
ri_95 <- as.numeric(lapply(pas_ri_all, quantile, prob=0.95))

ci_med <- as.numeric(lapply(pas_ci_all, quantile, prob=0.5))
ci_05 <- as.numeric(lapply(pas_ci_all, quantile, prob=0.05))
ci_95 <- as.numeric(lapply(pas_ci_all, quantile, prob=0.95))

# Plot

pdf(file="Analysis2_Passeriformes_2024-02-19.pdf")
plot(meanlog, L_01, ylim=c(0,1.1), col="white")
lines(meanlog, 1-ri_med, col="blue")
lines(meanlog, 1-ri_05, col="light blue")
lines(meanlog, 1-ri_95, col="light blue")
lines(meanlog, 1-ci_med, col="red")
lines(meanlog, 1-ci_05, col="pink")
lines(meanlog, 1-ci_95, col="pink")
lines(meanlog, rhi_med, col="black")
lines(meanlog, rhi_05, col="gray")
lines(meanlog, rhi_95, col="gray")
lines(meanlog, L_01, col="dark green")
dev.off()

# Y axis range values:

min(1-ci_med)
max(1-ci_med)

min(1-ri_med)
max(1-ri_med)

min(rhi_med)
max(rhi_med)

# Calculate area under curve

trapz(meanlog, L_01)
trapz(meanlog, 1-ci_med)
trapz(meanlog, 1-ri_med)
trapz(meanlog, rhi_med)

#####

### Telluraves ###

# Use acctran to populate tree with branch lengths based on character state changes in the matrix. Also randomly resolves polytomies.

tel_tree <- acctran(Tell_tree, tell_phydat) 

# Save topology for Analyses 3 & 4

write.tree(tel_tree, file="Telluraves_tree_branchlengths_dichot.phy")

#Simulate 100 matrices per meanlog value with the empirical topology and same number of characters, but 0.9 binary chars and 0.1 3-state chars. This simulation method generates whole matrices using the sim.morpho() feature from the 'DispRity' package.

tel_matrices <- c()
for(i in 1:length(meanlog)){
	tel_matrices[[i]] <- replicate(100, sim.morpho(tel_tree, characters=as.numeric(dim(Tell_matrix)[2]), states=c(0.9, 0.1), model="ER", rates=c(rlnorm, meanlog= meanlog[i], sdlog=0.1), invariant=FALSE))
	print(meanlog[i])	
}

#Convert simulated matrices from an array to a list
tel_sims <- c()
for(j in 1:length(tel_matrices)){

	sim <- list()
	for(i in 1:dim(tel_matrices[[j]])[3]){
		sim[[i]] <- tel_matrices[[j]][,,i]
		}
		tel_sims[[j]] <- sim
}

# Convert matrices to phydat objects
tel_phydats <- c()
for(i in 1:length(tel_sims)){
	tel_phydats[[i]] <- lapply(tel_sims[[i]], morph.phydat)
}

###

#Calculate metrics

#Consistency index
ci <- c()
tel_ci_all <- c()
for(i in 1:length(tel_phydats)){
	for(j in 1:length(tel_phydats[[i]])){
		ci[[j]] <- CI(tel_tree, tel_phydats[[i]][[j]]) #Calculate CI for each individual phydat
	}
	tel_ci_all[[i]] <- ci #Make object of all CI values for all phydats
}

#Make numeric and calculate median etc

for(i in 1:length(tel_ci_all)){
	tel_ci_all[[i]] <- as.numeric(tel_ci_all[[i]])
}

#Retention index
ri <- c()
tel_ri_all <- c()
for(i in 1:length(tel_phydats)){
	for(j in 1:length(tel_phydats[[i]])){
		ri[[j]] <- RI(tel_tree, tel_phydats[[i]][[j]]) #Calculate RI for each individual phydat
	}
	tel_ri_all[[i]] <- ri #Make object of all RI values for all phydats
}

for(i in 1:length(tel_ri_all)){
	tel_ri_all[[i]] <- as.numeric(tel_ri_all[[i]])
}

# Relative homoplasy index
# n = 100
tic()
tel_rhi <- c()
for(i in 1:length(tel_sims)){
	tel_rhi[[i]] <- lapply(tel_phydats[[i]], RHI, tel_tree, 100)
	print(i/length(tel_sims)*100)
}
toc()

#Extract RHI values
tel_rhi_val <- c()
tel_rhi_vals <- c()

for(j in 1:length(tel_rhi)){
	for(k in 1:length(tel_rhi[[j]])){
		tel_rhi_val[[k]] <- tel_rhi[[j]][[k]][[2]]
	}
	tel_rhi_vals[[j]] <- tel_rhi_val
}

for(i in 1:length(tel_rhi_vals)){
	tel_rhi_vals[[i]] <- as.numeric(tel_rhi_vals[[i]])
}

# L (tree length) - normalised H is identical to normalised L

tel_L_val <- c()
tel_L_vals <- c()
for(j in 1:length(tel_rhi)){
	for(k in 1:length(tel_rhi[[j]])){
		tel_L_val[[k]] <- tel_rhi[[j]][[k]][[4]]
	}
	tel_L_vals[[j]] <- tel_L_val
}

for(i in 1:length(tel_L_vals)){
	tel_L_vals[[i]] <- as.numeric(tel_L_vals[[i]])
}

# Get median and quantiles

L_med <- as.numeric(lapply(tel_L_vals, quantile, prob=0.5))
L_05 <- as.numeric(lapply(tel_L_vals, quantile, prob=0.05))
L_95 <- as.numeric(lapply(tel_L_vals, quantile, prob=0.95))

#Normalise values between 0 and 1
L_01 <- (L_med-min(L_med))/(max(L_med)-min(L_med))
L_01_05 <- (L_05-min(L_05))/(max(L_05)-min(L_05))
L_01_95 <- (L_95-min(L_95))/(max(L_95)-min(L_95))

rhi_med <- as.numeric(lapply(tel_rhi_vals, quantile, prob=0.5))
rhi_05 <- as.numeric(lapply(tel_rhi_vals, quantile, prob=0.05))
rhi_95 <- as.numeric(lapply(tel_rhi_vals, quantile, prob=0.95))

ri_med <- as.numeric(lapply(tel_ri_all, quantile, prob=0.5))
ri_05 <- as.numeric(lapply(tel_ri_all, quantile, prob=0.05))
ri_95 <- as.numeric(lapply(tel_ri_all, quantile, prob=0.95))

ci_med <- as.numeric(lapply(tel_ci_all, quantile, prob=0.5))
ci_05 <- as.numeric(lapply(tel_ci_all, quantile, prob=0.05))
ci_95 <- as.numeric(lapply(tel_ci_all, quantile, prob=0.95))

# Plot

pdf(file="Analysis2_Telluraves_2024-02-19.pdf")
plot(meanlog, L_01, ylim=c(0,1.1), col="white")
lines(meanlog, 1-ri_med, col="blue")
lines(meanlog, 1-ri_05, col="light blue")
lines(meanlog, 1-ri_95, col="light blue")
lines(meanlog, 1-ci_med, col="red")
lines(meanlog, 1-ci_05, col="pink")
lines(meanlog, 1-ci_95, col="pink")
lines(meanlog, rhi_med, col="black")
lines(meanlog, rhi_05, col="gray")
lines(meanlog, rhi_95, col="gray")
lines(meanlog, L_01, col="dark green")
dev.off()

# Y axis range values:

min(1-ci_med)
max(1-ci_med)

min(1-ri_med)
max(1-ri_med)

min(rhi_med)
max(rhi_med)

# Calculate area under curve

trapz(meanlog, L_01)
trapz(meanlog, 1-ci_med)
trapz(meanlog, 1-ri_med)
trapz(meanlog, rhi_med)






