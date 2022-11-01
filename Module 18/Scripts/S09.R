#
# R script session S09 on cluster analysis 
#

rm(list=ls())

#install.packages("dendextend")
library(dendextend)

fn <- "http://www-eio.upc.es/~jan/SISG/DistCEUJPTYRI.rda"
load(url(fn))

ls()

table(Pop)
sum(table(Pop))
class(Dals)
dim(Dals)
Dals[1:5,1:5]

Dm <- as.dist(Dals)

#
# Hierarchical clustering
#

?hclust

hc.single <- hclust(Dm,method="single")

plot(hc.single,ylab="Distance",main="single linkage; allele-sharing distance",
     xlab="",hang=-1,las=1,cex.main=1,sub="")

hc.complete <- hclust(Dm,method="complete")

plot(hc.complete,ylab="Distance",main="complete linkage; allele-sharing distance",
     xlab="",hang=-1,las=1,cex.main=1,sub="")

clusters.complete <- cutree(hc.complete, k=3)

table(clusters.complete)

table(clusters.complete,Pop)

hc.average <- hclust(Dm,method="average")

plot(hc.average,ylab="Distance",main="average linkage; allele-sharing distance",
     xlab="",hang=-1,las=1,cex.main=1,sub="")


clusters.average <- cutree(hc.complete, k=3)

table(clusters.average)

table(clusters.average,Pop)

hc.ward <- hclust(Dm,method="ward.D2")
plot(hc.ward,ylab="Distance",main="Ward D; allele-sharing distance",
       xlab="",hang=-1,las=1,cex.main=1,sub="")

clusters.ward <- cutree(hc.ward, h=1)

table(clusters.ward,Pop)

#
# Embellishment of the dendrogram
#

dend <- as.dendrogram(hc.ward)

cluster.col <- c("grey","darkgoldenrod1","limegreen")

#install.packages("clusterSim")
library(clusterSim)

dend1 <- color_branches(dend, k = 3, 
                        col=cluster.col,groupLabels = c(1,2,3))

plot(dend1,leaflab="none",ylab="Distance",main="Ward; allele-sharing distance")

table(clusters.ward,Pop)


#
# But who is who?
#

?order.dendrogram
ii <- order.dendrogram(dend1)
ii

clusters.ward[ii]
Pop[ii]

nclusters <- clusters.ward
map <- unique(clusters.ward[ii])
map

for(i in 1:length(map)) {
  nclusters[clusters.ward==map[i]] <- i
}

table(clusters.ward,nclusters)
table(nclusters,Pop)



#
# Partitioning: K-means
#

dim(X.Gen)

nmis <- function (x) 
{
  y <- sum(is.na(x))
  return(y)
}

sum(is.na(X.Gen))
nmv <- apply(X.Gen,2,nmis)
table(nmv)

X.Gen2 <- X.Gen[,nmv==0]
sum(is.na(X.Gen2))


rm(X.Gen)
rm(Dals)

set.seed(123)

output3 <- kmeans(X.Gen2[,sample(1:ncol(X.Gen2),1e4)],3)

table(output3$cluster,Pop)

#output3
output3$totss
output3$withinss
output3$tot.withinss
output3$betweenss
output3$totss-sum(output3$withinss)
output3$size
output3$iter

Tab <- rbind(output3$size,output3$withinss)
Tab <- rbind(Tab,output3$withinss/(output3$size-1))
colnames(Tab) <- 1:ncol(Tab)
rownames(Tab) <- c("Size","SS","Var")

Tab

WSS <- output3$tot.withinss
BSS <- output3$betweenss
TSS <- output3$totss

c(WSS=WSS,BSS=BSS,TSS=TSS)

BSS/TSS


