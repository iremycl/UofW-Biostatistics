#
# R script session S06 Multidimensional scaling
#

library(calibrate)
rm(list=ls())

#
# First the classical geographic application
#

fn <- "http://www-eio.upc.es/~jan/SISG/USCityDistances.csv"
X <- read.csv(fn,header=FALSE)

head(X)
class(X)

city <- X[,1]
state <- X[,2]

X <- X[,c(-1,-2)]

X[1:5,1:5]
X <- as.matrix(X)

n <- nrow(X)
n

#
# Graffelman, J., (2019) Goodness-of-fit filtering in classical metric multidimensional 
# scaling with large datasets. 
# Journal of Applied Statistics. https://www.tandfonline.com/doi/full/10.1080/02664763.2019.1702929
#

# Set the cities to be the rownames of the matrix
rownames(X) <- city
colnames(X) <- city
out.mds <- PrinCoor(X)

# Check eigenvalues
out.mds$standard.decom

# Goodness-of-fit table with absolute eigenvalues
out.mds$absolute.decom
out.mds$absolute.decom[, 1] / sum(out.mds$absolute.decom[, 1])
cumsum(out.mds$absolute.decom[,2])

# Retrieve the principal coordinates
Fp <- out.mds$X
Fp <- -1*Fp

fr1 <- paste("(",toString(round(100*out.mds$absolute.decom[1,2],digits=2)),"%)",sep="")
fr2 <- paste("(",toString(round(100*out.mds$absolute.decom[2,2],digits=2)),"%)",sep="")

plot(Fp[,1],Fp[,2],xlim=c(-2000,2000),asp=1,cex=0.6,
     xlab=paste("First principal axis",fr1),ylab=paste("Second principal axis",fr2))
citycop <- city
city
citycop[3] <- ""
citycop[4] <- ""
citycop[5] <- ""
citycop[8] <- ""
citycop[11] <- ""
citycop[15] <- ""
citycop[13] <- ""
citycop[21] <- ""
citycop[26] <- ""
citycop[24] <- ""
citycop[28] <- ""
textxy(Fp[,1],Fp[,2],citycop,cex=0.6)
textxy(Fp[13,1],Fp[13,2],city[13],cex=0.6,offset=0.5)
origin(lty="dotted")
text(Fp[15,1],Fp[15,2],city[15],pos=1,cex=0.6)
text(Fp[24,1],Fp[24,2],city[24],pos=1,cex=0.6)
text(Fp[4,1],Fp[4,2],city[4],pos=3,cex=0.6)
text(Fp[11,1],Fp[11,2],city[11],pos=2,cex=0.6)
text(Fp[28,1],Fp[28,2],city[28],pos=1,cex=0.6,offset=0.2)
text(Fp[8,1],Fp[8,2],city[8],pos=1,cex=0.6,offset=0.2)
text(Fp[3,1],Fp[3,2],city[3],pos=3,cex=0.6)
text(Fp[21,1],Fp[21,2],city[21],pos=1,cex=0.6,offset=0.05)
text(Fp[5,1],Fp[5,2],city[5],pos=3,cex=0.6)
text(Fp[26,1],Fp[26,2],city[26],pos=1,cex=0.6)

#
# Goodness-of-fit
#

Dhat <- as.matrix(dist(Fp[,1:2]))
#De <- X

devec <- X[lower.tri(X)]
dhvec <- Dhat[lower.tri(Dhat)]

plot(devec,dhvec,xlab="Observed distance",ylab="Fitted distance")
abline(0,1,col="blue",lwd=2)

#
# Non-metric MDS for the same data
#

set.seed(123)
y <- matrix(runif(2*n),ncol=2)
out.nm <- isoMDS(X,y=y,maxit=100)
Fnm <- out.nm$points

Fnm[,2] <- -1*Fnm[,2]

plot(Fnm[,1],Fnm[,2],asp=1,xlim=c(-30,30),xlab="First dimension",ylab="Second dimension",
     main="US cities (non-metric MDS)")
textxy(Fnm[,1],Fnm[,2],city,cex=0.5)
origin(lty="dotted")    

#
# Now with some data from 1000G
#
rm(list=ls())
fn <- "http://www-eio.upc.es/~jan/SISG/DistCEUJPTYRI.rda"
load(url(fn))
ls()

dim(X.Gen)
X.Gen[1:5,1:5]
n <- nrow(X.Gen)
n
#
# Allele sharing distances (precomputed)
#

#
# Dals <- as.matrix(dist(X.Gen,method="manhattan"))
#
#?dist
#

dim(Dals)
Dals[1:5,1:5]
dim(Dals)

Dm <- as.dist(Dals)

Fp <- cmdscale(Dm,k=2)

table(Pop)
colvec <- rep(NA,n)
colvec[Pop=="CEU"] <- "red"
colvec[Pop=="JPT"] <- "blue"
colvec[Pop=="YRI"] <- "black"

plot(Fp[,1],Fp[,2],asp=1,xlab="First principal axis",ylab="Second principal axis",col=colvec)
origin(lty="dotted")
legend("topleft",c("CEU","JPT","YRI"),pch=1,col=c("red","blue","black"),cex=0.75)

#
# Diagnostic
#

Dhat <- as.matrix(dist(Fp[,1:2]))

dalsvec <- Dals[lower.tri(Dals)]
dhvec <- Dhat[lower.tri(Dhat)]

plot(dalsvec,dhvec,xlab="Observed distance",ylab="Estimated distance",xlim=c(0,1),asp=1)
abline(0,1,col="blue",lwd=2)
