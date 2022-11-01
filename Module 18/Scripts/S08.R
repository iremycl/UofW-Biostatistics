#
# R script session S08 on canonical correlation analysis 
#

rm(list=ls())

fn <- "http://www-eio.upc.es/~jan/SISG/Messina.csv"
X <- read.csv(fn)

head(X)

country <- substr(X$Country.Island,1,3)
country

n <- nrow(X)
n
colvec <- rep("black",n)
colvec[country=="Ita"] <- "green"
colvec[country=="Gre"] <- "blue"
colvec[country=="Tur"] <- "brown"
colvec[country=="Cze"] <- "red"
colvec[country=="Pal"] <- "black"

colnames(X)[6] <- "N"
colnames(X)

#
# Geographical map
#

plot(X$Long.E,X$Lat.N,asp=1,xlab="longitude",ylab="latitude",col=colvec)
countries <- c("Ita","Gre","Tur","Cze","Pal")
legend("topright",countries,col=c("green","blue","brown","red","black"),pch=1)

Xc <- scale(X[,4:5],scale=FALSE)
head(Xc)
colnames(Xc) <- c("Long (E)","Lat (N)")

Y <- X[,(8:ncol(X)-1)]
Y
sum(Y==0)

Y <- Y/rowSums(Y)
rowSums(Y)
head(Y)

Yn <- cbind(Y[,1]+Y[,2]+Y[,3],Y[,4:7],Y[,8]+Y[,7])
sum(Yn==0)

#install.packages("zCompositions")
library(zCompositions)

Ym <- cmultRepl(Yn)
head(Ym)
head(Yn)

colnames(Ym)[1] <- "A9_11_12"
colnames(Ym)[ncol(Ym)] <- "A17_18"


lX <- log(Ym)
Yclr <- scale(t(scale(t(lX),scale=FALSE)),scale=FALSE)

#
# We load a function for a covariance-based canonical analysis.
#

fn <- "http://www-eio.upc.es/~jan/SISG/canocov.R"
source(fn)

out <- canocov(Xc,Yclr)

out$ccor
out$fitRxy

library(calibrate)
bplot(out$Fs,25*out$Gp,rowarrow = TRUE,rowch=NA,colch=NA,rowlab = "",
      collab=colnames(Ym),cex.collab = 0.5)
text(out$Fs[,1],out$Fs[,2],colnames(Xc),cex=0.5,pos=1)

fa <- 1
Sa <- fa*out$U
points(Sa[,1],Sa[,2],pch=1,col=colvec)
legend("topright",countries,col=c("green","blue","brown","red","black"),
       pch=1,cex=0.75)

pct <- function(decom,di=1) {
  l1 <- paste("(",toString(round(100*decom[2,di],digits=2)),"%)",sep="")
  return(l1)
}


text(2.5,-0.2,pct(out$fitRxy),cex=0.5)
text(0.2,3,pct(out$fitRxy,2),cex=0.5,srt=90)

#
# Permutation test
# 

nsimul <- 10000

M <- matrix(NA,nrow=nsimul,ncol=2)
n
sample(n)
set.seed(123)
for(i in 1:nsimul) {
  index <- sample(n)
  Xc.scrambled <- Xc[index,]
  out.sim <- canocov(Xc.scrambled,Yclr)
  canvec <- diag(out.sim$ccor)
  M[i,] <- canvec
}

cc <- diag(out$ccor)
cc

head(M)
dim(M)


hist(M[,1],main="Permutation distribution",xlab="First canonical correlation")
abline(v=cc[1],col="red")


hist(M[,2],main="Permutation distribution",xlab="Second canonical correlation")
abline(v=cc[2],col="red")

# permutation p-values

pvals <- numeric(2)
for(i in 1:2) {
  pvals[i] <- sum(M[,i] >= cc[i])/nsimul
}
round(pvals,digits=10)




