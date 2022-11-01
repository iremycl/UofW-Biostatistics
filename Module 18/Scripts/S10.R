#
# R script session S10 on cluster analysis; second part.
#

rm(list=ls())

library(dendextend)
library(cluster)
library(clusterSim)


fn <- "http://www-eio.upc.es/~jan/SISG/NistDist.rda"
load(url(fn))

ls()

table(pop)
#
# MDS of the NIST Data
#

dim(Djaccard)
length(pop)
table(pop)


out1 <- cmdscale(Djaccard,k=nrow(Djaccard)-1,eig=TRUE)

ev <- out1$eig 
ev
round(ev,digits=2) 

length(ev)
plot(out1$eig,type="b")
abline(h=0)

fit1 <- 100*ev[1]/sum(abs(ev))
fit2 <- 100*ev[2]/sum(abs(ev))
fit1
fit2

labx <- paste0("(",round(fit1,digits=1),"%)")
laby <- paste0("(",round(fit2,digits=1),"%)")

#
# MDS maps
#

n <- nrow(Djaccard)
colvec <- rep("grey",n)
table(pop)

colvec[pop=="Asian"] <- "yellow"
colvec[pop=="AA"]    <- "black"
colvec[pop=="Cauc"]  <- "blue"

col.anc  <- c("black","yellow","blue")
ancs     <- c("AA","Asian","Cauc")
ancsfull<- c("African American","Asian","Caucasian")

plot(out1$points[,1],out1$points[,2],asp=1,col=colvec,
     main="MDS map NIST STRs",
     xlab=paste0("First principal axis ",labx),
     ylab=paste0("Second principal axis ",laby))
legend("bottomleft",ancsfull,pch=1,col=col.anc,cex=0.75)  

#
# Model based clustering
#

library(mclust)
Fp <- out1$points
dim(Fp)
model1 <- Mclust(Fp[,1:2],G=3)

summary(model1,parameters=TRUE)
model1$parameters$pro
group.mb <- model1$classification
table(group.mb,pop)

plot(model1,what="classification",col=c("black","blue","yellow"),xlab=paste0("First principal axis ",labx),
     ylab=paste0("Second principal axis ",laby))
clu <- paste("Cluster",1:3)
legend("bottomleft",clu,pch=c(19,2,0),col=col.anc,cex=0.75)  


dis <- dist(Fp[,1:2])

#
# Silhouette scores
#

model3 <- Mclust(Fp[,1:2],G=3)
g.mb.3 <- model3$classification
sil3 = silhouette(g.mb.3, dis)
aaa <- summary(sil3)
attributes(aaa)
w3 <- summary(sil3)$avg.width
w3

AA <- as.matrix(sil3)
head(AA)

ms <- numeric(3)
for(i in 1:3) {
  ms[i] <- mean(AA[AA[,1]==i,3])
}
ms

model2 <- Mclust(Fp[,1:2],G=2)
g.mb.2 <- model2$classification
sil2 = silhouette(g.mb.2, dis)
w2 <- summary(sil2)$avg.width


model4 <- Mclust(Fp[,1:2],G=4)
g.mb.4 <- model4$classification
sil4 = silhouette(g.mb.4, dis)
w4 <- summary(sil4)$avg.width

model5 <- Mclust(Fp[,1:2],G=5)
g.mb.5 <- model5$classification
sil5 = silhouette(g.mb.5, dis)
w5 <- summary(sil5)$avg.width

c(s2=w2,s3=w3,s4=w4,s5=w5)

#
# F statistics
#

PseudoF <- numeric(5)
PseudoF.km <- numeric(5)
Silh <- numeric(5)
Silh.km <- numeric(5)
for(i in 2:5) {
  model3 <- Mclust(Fp[,1:2],G=i)
  group.mb <- model3$classification
  PseudoF[i] <- index.G1(Fp[,1:2],group.mb)
  Silh[i] <- index.S(dist(Fp[,1:2]),group.mb)
  
  outk <- kmeans(Fp[,1:2],i)
  group.km <- outk$cluster
  PseudoF.km[i] <- index.G1(Fp[,1:2],group.km)
  Silh.km[i] <- index.S(dist(Fp[,1:2]),group.km)
}

PseudoF <- PseudoF[-1]
Silh <- Silh[-1]
PseudoF.km <- PseudoF.km[-1]
Silh.km <- Silh.km[-1]


plot(2:5,PseudoF,type="b",xlab="Clusters",ylab="Pseudo F-statistic",
     main="F-statistics NIST STRs; model-based")

plot(2:5,PseudoF.km,type="b",xlab="Clusters",ylab="Pseudo F-statistic",
     main="F-statistics NIST STRs; k-means")

plot(2:5,Silh,type="b",xlab="Clusters",ylab="Silhouette coefficient")

#
# Silhouette plots
#

dis <- dist(Fp[,1:2])

model3 <- Mclust(Fp[,1:2],G=3)
group.mb3 <- model3$classification
sil = silhouette(group.mb3, dis)

plot(sil,col=c("green","red","blue"),
     max.strlen = 0,sub="",main="")

library(factoextra)
out <- fviz_silhouette(sil,main="")
plot(out)



