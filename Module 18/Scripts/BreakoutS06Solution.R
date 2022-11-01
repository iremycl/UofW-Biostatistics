#
# Breakout S06 exercise; Day 02 session S06 Multidimensional scaling
#

rm(list=ls())

fn <- "http://www-eio.upc.es/~jan/SISG/RstYSTR.csv"

X <- read.csv(fn,header=TRUE)

X[1:5,1:5]

rownames(X) <- X[,1]
X <- X[,-1]

continent <- X[,1]
table(continent)

X <- X[,-1]

class(X)
X <- as.matrix(X)

X[is.na(X)] <- 0

D <- X + t(X)

out.mds <- cmdscale(D,eig=TRUE)
attributes(out.mds)
Y <- out.mds$points

colvec <- rep(NA,nrow(D))
table(continent)
colvec[continent=="SubSaharan"] <- "black"
colvec[continent=="Europe"]     <- "blue"
colvec[continent=="FarEast"]    <- "yellow"
colvec[continent=="America"]    <- "red"
colvec[continent=="MiddleEast"] <- "green"
colvec[continent=="Oceania"]     <- "orange"

plot(Y[,1],Y[,2],asp=1,col=colvec,pch=19,
     xlab="First principal axis",ylab="Second principal axis")
text(Y[,1],Y[,2],substr(rownames(Y),1,3),pos=1)

unique(continent)
legend("bottomleft",unique(continent),
       col=c("black","yellow","red","orange","blue","green"),pch=19)

plot(out.mds$eig)
out.mds$eig

out.mds <- cmdscale(D,eig=TRUE,k=2)
out.mds$GOF
