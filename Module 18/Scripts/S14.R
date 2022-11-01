#
# R script session S14; multivariate inference
#

rm(list=ls())

#
# MVA inference
#

fn <- "http://www-eio.upc.es/~jan/SISG/hemophilia.dat"
X <- read.table(fn) 

#install.packages("ICSNP")
library(ICSNP)

colnames(X) <- c("Group","AHFact","AHFantigen")
head(X)

X$groupf <- factor(X$Group,labels=c("non-carrier","carrier"))

G1 <- X[X$groupf=="non-carrier",2:3]
G2 <- X[X$groupf=="carrier",2:3]

dim(G1)
dim(G2)

n1 <- nrow(G1)
n2 <- nrow(G2)

S1 <- cov(G1)
S2 <- cov(G2)
r1 <- c(n1,colMeans(G1))
r2 <- c(n2,colMeans(G2))

m1 <- colMeans(G1)
m2 <- colMeans(G2)

Ta <- rbind(r1,r2)
rownames(Ta) <- c("non-carrier","carrier")

Ta

colvec <- rep("black",nrow(X))
colvec[X$groupf=="carrier"] <- "red"
colvec[X$groupf=="non-carrier"] <- "blue"

plot(X$AHFact,X$AHFanti,col=colvec,xlab="AHF activity",ylab="AHF antigen")

nonc <- colMeans(G1)
nonc

ca <- colMeans(G2)

points(nonc[1],nonc[2],pch=17,col="blue")
points(ca[1],ca[2],pch=17,col="red")
legend("topleft",c("non-carrier","carrier"),pch=1,col=c("blue","red"))



plot(X$AHFact,X$AHFanti,col=colvec,xlab="AHF activity",ylab="AHF antigen")
points(nonc[1],nonc[2],pch=17,col="blue")
points(ca[1],ca[2],pch=17,col="red")
S1
S2
library(ellipse)
Z1 <- ellipse(S1,level=0.9,centre=m1)
points(Z1,type="l",col="blue",lwd=1)

Z2 <- ellipse(S2,level=0.9,centre=m2)
points(Z2,type="l",col="red",lwd=1)


legend("topleft",c("non-carrier","carrier"),pch=1,col=c("blue","red"))

S1 <- cov(G1)
S2 <- cov(G2)

#
# Equality of mean vectors
#


HotellingsT2(G1,G2,test="chi")
HotellingsT2(G1,G2,test="f")
HotellingsT2(G1,G2)

?HotellingsT2

boxplot(X$AHFact~X$Group)

#
# Homogeneity of covariances
#

#install.packages("biotools")
library(biotools)

head(X)
boxM(X[,2:3],grouping=X$groupf)

#
# MANOVA
#

rm(list=ls())


fn <- "http://www-eio.upc.es/~jan/SISG/NistDist.rda"

load(url(fn))

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


popf <- factor(pop)

colnames(out.mds) <- paste("PC",1:10,sep="")

Y <- data.frame(out.mds[,1:2],popf)

results1 <- manova(out.mds[,1:2]~popf, data = Y)
ss1 <- summary(results1, test="Wilks")
ss1

w1 <- ss1$stats[1,2]
w1
p1 <- ss1$stats[1,6]
p1

results2 <- manova(out.mds[,3:10]~popf, data = Y)
ss2 <- summary(results2, test="Wilks")
ss2

w2 <- ss2$stats[1,2]
p2 <- ss2$stats[1,6]



ZZZ <- cbind(c(w1,w2),c(p1,p2))
colnames(ZZZ) <- c("Wilks","pvalue")
rownames(ZZZ) <- c("Axes 1-2", "Axes 3-10")

ZZZ

Q <- data.frame(out.mds,popf)
colnames(Q)

anova(lm(PC1~popf,data=Q))
anova(lm(PC2~popf,data=Q))
anova(lm(PC3~popf,data=Q)) 
anova(lm(PC4~popf,data=Q))
anova(lm(PC5~popf,data=Q))
anova(lm(PC6~popf,data=Q)) 


