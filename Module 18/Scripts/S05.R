#
# R script session S05 Log-ratio Principal component analysis.
#

rm(list=ls())

library(data.table)
library(calibrate)
#install.packages("robCompositions")
#library(robCompositions)

#
# Making a ternary diagram
#

rr <- c(AA=0.3,BB=0.1,AB=0.6)

#
# with package robCompositions
#


# Xm <- rbind(rr,rr) 
# Xm
# ternaryDiag(Xm,name=c(expression(f[AA]),expression(f[BB]),
#                      expression(f[AB])),main="",mcex=2,cex.main=2,
#            pch=19,col="red")

#
# Ternary diagram with package HardyWeinberg
#

library(HardyWeinberg)

rr <- rr[c(1,3,2)]

HWTernaryPlot(rr,n=100,vbounds=FALSE,region=0,hwcurve=FALSE)

#
# With HW parabole and acceptance region
#
              
HWTernaryPlot(rr,n=100,vbounds=FALSE,region=1,curvecols = c("blue","black"))

#
# A bunch of SNPs in a tenary diagram
#

data(Alzheimer)
head(Alzheimer)
class(Alzheimer)
?Alzheimer

Controls <- Alzheimer[Alzheimer$Group==0,]
Controls <- as.matrix(Controls)
head(Controls)
apply(Controls,1,sum)

max(apply(Controls,1,sum))
HWTernaryPlot(Controls[,1:3],n=715,vbounds=FALSE,
              region=1,curvecols = c("blue","black"),alpha = 0.05)
HWTernaryPlot(Controls[,1:3],n = 715,vbounds=FALSE,
              region=1,curvecols = c("blue","black"),alpha = 0.05/nrow(Controls))


#
# log-ratio PCA of Y chromosomal STR DYS448 from Purps et. al.
#

rm(list=ls())

fn <- "http://www-eio.upc.es/~jan/SISG/STRDYS448.rda"
load(url(fn))
ls()
dim(Z)
head(Z)

table(Z$Sample)

#
# The zero issue
#

zerofraction <- function(x) {
  n <- length(x)
  nozero <- sum(x==0)
  per <- nozero/n
  return(per)
}

colnames(Z)    

Tab <- table(Z$Sample,Z$DYS448)  
Tab
zf <- apply(Tab,2,zerofraction)
zf

dim(Tab)
sum(Tab==0)

pzero <- 100*sum(Tab==0)/(129*25)
pzero

zerolimit <- 0.05

ind <- zf < zerolimit

common <- Tab[,ind]
rares <- Tab[,!ind]
remain <- apply(rares,1,sum)

Tabn <- cbind(common,remain) # sum alleles that have many zeros.

head(Tabn)
dim(Tabn)

sum(Tabn==0)

library(zCompositions)

?cmultRepl
Tabc <- cmultRepl(Tabn,suppress.print=TRUE) # substitute zeros

sum(Tabc==0)


#
# clr transformation
#

lX <- log(Tabc)
Xclr <- t(scale(t(lX),scale=FALSE))
rownames(Xclr) <- rownames(Tab)
colnames(Xclr)[4] <- "Other"

head(Xclr)

#
# ordinary pca
#

out.pc <- princomp(Xclr,cor=FALSE)

#
# biplot coordinates
#

Fp <- out.pc$scores
Gs <- as.matrix(unclass(out.pc$loadings))
class(Gs)
Gs

#
# goodness-of-fit
#

res <- summary(out.pc)
attributes(res)
res$sdev
la <- res$sdev^2
la
fr <- la/sum(la)
cu <- cumsum(fr)
dec <- rbind(la,fr,cu)
dec
fr
fr1 <- round(100*fr[1],digits=2)
fr2 <- round(100*fr[2],digits=2)

#
# biplot
#

biplot(out.pc)

library(calibrate)

bplot(Fp,3*Gs,rowch=1,
      collab=colnames(Xclr),
      rowlab=rownames(Xclr),cex.rowlab = 0.5,
      colch = NA,main="DYS448")
text(3,-0.1,paste("PC1 (",fr1,"%)"),cex=0.5,col="blue")
text(-0.1,2.8,paste("PC2 (",fr2,"%)"),cex=0.5,col="blue",srt=90)






