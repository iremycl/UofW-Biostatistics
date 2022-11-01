#
# R script session S13; multivariate normal distribution
# Jan Graffelman; July 2021. 
#

rm(list=ls())

#
# We illustrate multivariate normality issues with MDS output of the NIST data
#

fn <- "http://www-eio.upc.es/~jan/SISG/NistDist.rda"
load(url(fn))

ls()
table(pop)

dim(Djaccard)
length(pop)

out.mds <- cmdscale(Djaccard,k=10)

n <- nrow(Djaccard)

colvec <- rep("black",n)
table(pop)
colvec[pop=="AA"] <- "black"
colvec[pop=="Asian"] <- "yellow"
colvec[pop=="Cauc"] <- "blue"


plot(out.mds[,1],out.mds[,2],asp=1,col=colvec,
     xlab="First principal axis",ylab="Second principal axis",main="MDS map")
origin(lty="dotted")
legend("topleft",c("African American","Asian","Caucasian"),col=c("black","yellow","blue"),
       pch=1,cex=0.75)

pop[pop=="AA"]   <- "African American"
pop[pop=="Cauc"] <- "Caucasian"

popf <- factor(pop)
length(popf)

#
# normality for a single variable 
#

qqnorm(out.mds[,1])
qqline(out.mds[,1],col="red",lwd=2)

qqnorm(out.mds[,1],col=colvec)
qqline(out.mds[,1],col="red",lwd=2)

table(popf)

qqnorm(out.mds[popf=="African American",1],main="African American")
qqline(out.mds[popf=="African American",1],col="red",lwd=2)


S <- cov(out.mds[,1:2])
m <- colMeans(out.mds[,1:2])

plot(out.mds[,1],out.mds[,2],asp=1,col=colvec,
     xlab="First principal axis",ylab="Second principal axis",main="MDS map")


library(ellipse)
Z1 <- ellipse(S,level=0.95,centre=m)
points(Z1,type="l",col="red",lwd=1)

#
# again stratify
#

Asi <- out.mds[popf=="Asian",1:2]
S <- cov(Asi)
mv <- colMeans(Asi)

Z1 <- ellipse(S,level=0.95,centre=mv)
points(Z1,type="l",col="yellow",lwd=1)

Z2 <- ellipse(S,level=0.50,centre=mv)
points(Z2,type="l",col="yellow",lwd=1)

Cau <- out.mds[popf=="Caucasian",1:2]
S <- cov(Cau)
mv <- colMeans(Cau)

Z1 <- ellipse(S,level=0.95,centre=mv)
points(Z1,type="l",col="blue",lwd=1)

Z2 <- ellipse(S,level=0.50,centre=mv)
points(Z2,type="l",col="blue",lwd=1)

Afr <- out.mds[popf=="African American",1:2]
S <- cov(Afr)
mv <- colMeans(Afr)

Z1 <- ellipse(S,level=0.95,centre=mv)
points(Z1,type="l",col="black",lwd=1)

Z2 <- ellipse(S,level=0.50,centre=mv)
points(Z2,type="l",col="black",lwd=1)

#
# Chi-square splot
#

Chisquareplot <- function (X, label = FALSE, main = paste("Chi-square plot", 
                                                          deparse(substitute(X))), ...) 
{
  n <- nrow(X)
  p <- ncol(X)
  m <- apply(X, 2, mean)
  S <- cov(X)
  Sinv <- solve(S)
  d2vec <- NULL
  for (i in 1:nrow(X)) {
    d2 <- (X[i, ] - m) %*% Sinv %*% (X[i, ] - m)
    d2vec <- c(d2vec, d2)
  }
  out <- sort(d2vec, index.return = TRUE)
  ind <- out$ix
  d2vecs <- out$x
  qtl <- (1:n - 0.5)
  theoquant <- qchisq(qtl/n, df = p)
  plot(theoquant, d2vecs, pch = 19, ylab = expression(d^2), 
       xlab = "Theoretical Chisquare quantiles", main=main, ...)
  if (label) 
    text(theoquant, d2vecs, ind, cex = 0.5, pos = 2)
  abline(a = 0, b = 1, col = "red")
  return(d2vec)
}

Chisquareplot(out.mds[popf=="Asian",1:10],main="Asian")

Chisquareplot(out.mds[popf=="Caucasian",1:10],main="Caucasian")

Chisquareplot(out.mds[popf=="African American",1:10],main="African American")

