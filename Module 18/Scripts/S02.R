#
# R script session S02 Multivariate Analysis for Genetic Data
# Jan Graffelman; July 2021. 
#

#
# R script session S02 Matrix decompositions 
#

#
# We use a subset of the gene expression data from the article:
#
#
# Zheng C, Xu, R. (2021) Molecular subtyping of Alzheimer's disease with consensus 
#             non-negative matrixfactorization. PLoSONE 16(5):e0250278.
#

rm(list=ls())

#
# The full data set
#
#X <- read.csv("http://www-eio.upc.es/~jan/SISG/ROSMAP_genexp_ad.csv",header=TRUE)
#X[1:5,1:5]
#rownames(X) <- X[,1]
#X <- X[,-1]
#X <- t(X)

#
# A subset of 100 genes
#

fn <- "http://www-eio.upc.es/~jan/SISG/SubsetRosmapGenExp.dat"
X <- fread(fn,data.table = FALSE)
rownames(X) <- X[,1]
X <- X[,-1]

X[1:5,1:5]

boxplot(X[,1:10])

boxplot(log(X[,1:10]))
boxplot(sqrt(X[,1:10]))

Xs <- scale(sqrt(X))

out <- svd(Xs)

#
# low-rank approximation to standardized expression data
#

U <- out$u
D <- diag(out$d)
V <- out$v

#
# Reconstitution
#

Reconst <- U%*%D%*%t(V)
max(abs(Xs-Reconst))

#
# Rank two approximation
#

U2 <- U[,1:2]
D2 <- D[1:2,1:2]
V2 <- V[,1:2]

Xhat <- U2%*%D2%*%t(V2)

for(i in 1:5) {
  plot(Xs[,i],Xhat[,i],xlab=colnames(Xs)[i],ylab="approximation")  
}

#
# Decay of squared singular values
#

d2 <- diag(D*D)
plot(d2)

#
# Goodness-of-fit
#

gof <- sum(d2[1:2])/sum(d2)
gof

gof <- sum(d2[1:3])/sum(d2)
gof

#
# Errors
#

E <- Xs-Xhat

sum(diag(t(E)%*%E))
sum(d2[3:length(d2)])

#
# Left singular vectors related to principal components (we will see later)
#

#
# Map of the individuals
#

Fp <- U%*%D
plot(Fp[,1],Fp[,2],asp=1)

#
# Spectral decomposition of R
#

R <- cor(Xs)
R[1:5,1:5]

Res <- eigen(R) 
Dlam <- diag(Res$values)

plot(diag(Dlam))

V <- Res$vectors
V

# verify the decomposition 

E <- R - V%*%Dlam%*%t(V)
max(abs(E))

# and the orthogonality of V

t(V[,1:3])%*%V[,1:3]
round(t(V[,1:3])%*%V[,1:3],digits=6)

#
# Goodness-of-fit
#

la <- diag(Dlam)
la
la2 <- la*la
sum(la2[1:2]/sum(la2))
plot(la2)

Rhat <- V[,1:2]%*%Dlam[1:2,1:2]%*%t(V[,1:2])

plot(R[lower.tri(R)],Rhat[lower.tri(Rhat)],
     xlab="Observed correlation",
     ylab="Approximated correlation")

abline(0,1,col="red",lwd=2)



