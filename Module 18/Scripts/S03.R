#
# R script session S03 Multivariate Analysis for Genetic Data; Biplots
# Jan Graffelman; July 2021. 
#

rm(list=ls())

#install.packages(data.table)
#install.packages(calibrate)
library(data.table)
library(calibrate)

fn <- "http://www-eio.upc.es/~jan/SISG/Data1000G_CEU_CHB_YRI.rda"

load(url(fn))

ls()

table(Pop)
table(Pop,useNA="always")

dim(Yn)

Xsup <- Yn[is.na(Pop),]
dim(Xsup)

X <- Yn[!is.na(Pop),]
dim(X)

Pop <- Pop[!is.na(Pop)]

X[1:10,1:10]

X <- X[,7:ncol(X)]
X <- as.matrix(X)

sum(is.na(X))

table(Pop,useNA="always")

colvec <- rep(NA,nrow(X))
colvec[Pop=="CEU"] <- "blue"
colvec[Pop=="CHB"] <- "yellow"
colvec[Pop=="YRI"] <- "black"

Xc <- scale(X,scale=FALSE)

out <- svd(Xc)

U <- out$u
D <- diag(out$d)
V <- out$v

Fp <- U%*%D

plot(Fp[,1],Fp[,2],asp=1,col=colvec,
     xlab="First principal axis",ylab="Second principal axis",
     main="1000G CEU-CHB-YRI")

legend("bottomright",c("CEU","CHB","YRI"),pch=1,col=c("blue","yellow","black"))

origin()

#
# Goodness of fit
#

d <- diag(D)
d2 <- d*d
fr <- d2/sum(d2)

goftable <- cbind(d2,fr,cumsum(fr))

head(goftable)

#
# Looking at 3 dimensions.
#

pairs(Fp[,1:3],col=colvec)

#
# Adding some supplementary individuals; first we verify our formulae
#

set.seed(123)
i <- sample(1:nrow(X),10)
i

Xadd <- X[i,]

m <- colMeans(X)

Xaddc <- Xadd - rep(1,10)%o%m

Fsup <- Xaddc%*%V%*%solve(t(V)%*%V)


plot(Fp[,1],Fp[,2],asp=1,col=colvec,
     xlab="First principal axis",ylab="Second principal axis",
     main="1000G CEU-CHB-YRI")

legend("bottomright",c("CEU","CHB","YRI"),pch=1,col=c("blue","yellow","black"))

origin()

points(Fsup[,1],Fsup[,2],asp=1,col="brown",pch=1,cex=2,lwd=2)

#
# Now we really map in the unknown ones.
#

dim(Xsup)
class(Xsup)

Xsup <- Xsup[,7:ncol(Xsup)]
Xsup <- as.matrix(Xsup)

Xsupc <- Xsup - rep(1,nrow(Xsup))%o%m

Fsup <- Xsupc%*%V%*%solve(t(V)%*%V)


plot(Fp[,1],Fp[,2],asp=1,col=colvec,
     xlab="First principal axis",ylab="Second principal axis",
     main="1000G CEU-CHB-YRI")

legend("bottomright",c("CEU","CHB","YRI"),pch=1,col=c("blue","yellow","black"))

origin()

points(Fsup[,1],Fsup[,2],asp=1,col="brown",pch=1,cex=2,lwd=2)



