
# R script session S01 Multivariate Analysis for Genetic Data
# Jan Graffelman; July 2021. 
#
​
#
# We use allele intensites and called genotypes for two SNPs for calculation of basic
# matrices in MVA.
#
​
X <- read.table("http://www-eio.upc.es/~jan/SISG/Intensities.dat",header=TRUE)
head(X)
​
class(X)
X <- as.matrix(X)
​
Xi <- X[,1:4] # intensities only
​
n <- nrow(Xi)
p <- ncol(Xi)
n
p
​
m <- apply(Xi,2,mean)
m
​
colMeans(Xi)
​
s <- apply(Xi,2,sd)
s
​
Xc <- scale(Xi,scale=FALSE)
head(Xc)
​
Xc2 <- Xi - rep(1,n)%o%m
head(Xc2)
​
In <- diag(rep(1,n))
H <- In - (1/n)*rep(1,n)%o%rep(1,n)
Xc3 <- H%*%Xi
head(Xc3)
​
apply(Xc,2,mean)
​
apply(Xc,2,sd) 
s
​
Xs <- scale(Xi,scale=TRUE)
head(Xs)
​
apply(Xs,2,mean)
apply(Xs,2,sd)
​
Xs2 <- Xc%*%solve(diag(s))
head(Xs2)
​
Xs3 <- Xc/s
head(Xs3)
​
Xs4 <- sweep(Xc,2,s,FUN="/")
head(Xs4)
​
S <- cov(Xi)
S
​
(1/(n-1))*t(Xc)%*%Xc
​
R <- cor(Xi)
R
round(R,digits=2)
​
Ds <- diag(sqrt(diag(S)))
​
R2 <- solve(Ds)%*%S%*%solve(Ds)
R2
​
R3 <- (1/(n-1))*t(Xs)%*%Xs
R3
​
#
# Some distance calculations
#
​
?dist
De <- as.matrix(dist(Xi))
dim(De)
De[1:5,1:5]
​
De <- as.matrix(dist(Xc))
dim(De)
De[1:5,1:5]
​
Dwe <- as.matrix(dist(Xs))
dim(Dwe)
Dwe[1:5,1:5]
​
#
# Mahalanobis distance
#
​
Si <- solve(S)
Si
​
Si%*%S
​
round(Si%*%S,digits=8)
​
Dm <- matrix(0,nrow=n,ncol=n)
for(i in 1:n) {
  for(j in 1:n) {
    Dm[i,j] <- (Xi[i,]-Xi[j,])%*%Si%*%(Xi[i,]-Xi[j,])    
  }
}
​
Dm[1:5,1:5]
​
mahalanobis(Xi[1,],Xi[2,],cov=S)
​
?mahalanobis