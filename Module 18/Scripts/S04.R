#
# R script session S04 Principal component analysis.
#

rm(list=ls())

#install.packages(data.table)
#install.packages(calibrate)
library(data.table)
library(calibrate)


fn <- "http://www-eio.upc.es/~jan/SISG/CHD.raw"
X <- fread(fn,data.table=FALSE)

dim(X)
X[1:8,1:8]

X.fam <- X[,1:6]

X <- X[,7:ncol(X)]
rownames(X) <- X.fam[,1]

X <- as.matrix(X)

sum(is.na(X))

#
# We first analyse a small subset of the SNPs
#

X10 <- X[,1:10]

out.pca <- princomp(X10,cor=FALSE)

plot(out.pca)

biplot(out.pca)

S <- cov(X10)
V <- eigen(S)$vectors
Dl <- diag(eigen(S)$values)

Xc <- scale(X10,scale=FALSE)

Fp <- Xc%*%V
Gs <- V

bplot(Fp,Gs,cex.rowlab = 0.25)

bplot(Fp,5*Gs,cex.rowlab = 0.25,colch=NA,collab=colnames(X10),
      cex.collab=0.5,main="Form biplot")

#
# Let's do the covariance biplot
#

Fs <- Fp%*%sqrt(solve(Dl))
Gp <- Gs%*%sqrt(Dl)

bplot(Fs,5*Gp,cex.rowlab=0.25,colch=NA)

#
# Pretty similar....
#

diag(Dl)

#
# Now with more SNPs
#

dim(X)

out.pca <- princomp(X)

?prcomp
out.pca <- prcomp(X,scale=FALSE)

attributes(out.pca)

Fp <- out.pca$x

dim(Fp)

plot(Fp[,1],Fp[,2],asp=1,xlab="First principal axis",
     ylab="Second principal axis",main="CHD")
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")

la <- apply(Fp,2,var)
la
fr <- la/sum(la)
cu <- cumsum(fr)

Deco <- rbind(la,fr,cu)

Deco[,1:10]

Fp[abs(Fp[,1])>50,1:2]

plot(Fp[,1],Fp[,3],asp=1,xlab="First principal axis",
     ylab="Third principal axis")

Fp[abs(Fp[,3])>50,1:2]

#
# remove outliers
#

ii <- (abs(Fp[,1]) > 50) | (abs(Fp[,3]) >50)
sum(ii)

out.pca <- prcomp(X[!ii,],scale=FALSE)

attributes(out.pca)

Fp <- out.pca$x

dim(Fp)

plot(Fp[,1],Fp[,2],asp=1,xlab="First principal axis",
     ylab="Second principal axis",main="CHD")
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")

