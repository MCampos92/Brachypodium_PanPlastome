######################################################################
###################  DAPC Analysis Script ######################### 
################### Miguel Campos January 2025 #######################
######################################################################

library(adegenet)
library(dartR)

setwd("./DAPC")

# Load the plastomes data
Plastomes <- gl.read.fasta("PlastotypeS.fasta")

## PCA
pcoa <- gl.pcoa(Plastomes)
grp <- find.clusters(Plastomes, max.n.clust=10, method='ward', stat='AIC', glPca = pcoa)

head(grp$grp, 4)

# Proceed with DAPC
dapc1 <- dapc(Plastomes, grp$grp)
dapc1

# This step helps to determine the optimal number of PCs; if it's not optimal, 
# repeat the previous step using the selected PCs
temp <- optim.a.score(dapc1)

# Plot the DAPC
scatter(dapc1)

myCol <- c("blue", "purple", "darkgreen", "orange", "red", "pink", "brown")
scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4, cex=3, clab=0, leg=TRUE, txt.leg=paste("Cluster", 1:5))

# Plot with X representing the centers and the minimum spanning network between the groups
scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4, cex=4, clab=0, mstree=TRUE, scree.da=FALSE, posi.pca="bottomright", leg=TRUE, txt.leg=paste("Cluster", 1:5)) 
par(xpd=TRUE) 
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=2, lwd=5, col="black") 
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=2, lwd=2, col=myCol)

myInset <- function() { 
  temp <- dapc1$pca.eig 
  temp <- 100 * cumsum(temp) / sum(temp) 
  plot(temp, col=rep(c("black", "lightgrey"), c(dapc1$n.pca, 1000)), ylim=c(0, 100), xlab="PCAaxis", ylab="Cumulatedvariance(%)", cex=1, pch=20, type="h", lwd=2) 
} 
add.scatter(myInset(), posi="bottomright", inset=c(-0.03, -0.01), ratio=.28, bg=transp("white"))

# Structure-like plot
compoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:5), lab="", ncol=1, xlab="individuals", col=funky(5))

# Assignment of individuals
round(dapc1$posterior, 1)
write.csv(round(dapc1$posterior, 1), file="Itools.csv")

