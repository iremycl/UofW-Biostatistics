#
# R script session S07 Correspondence analysis 
#

rm(list=ls())


library(Ternary)


fn <- "http://www-eio.upc.es/~jan/SISG/Nist.rda"
load(url(fn))

ls()

head(Z)

dim(Z)

str1 <- paste(Z$F13B_1,Z$F13B_2,sep="/")

Pop <- Z$Pop
table(Pop,str1)

#
# We first do a three population analysis
#

str1 <- str1[Pop!="Hisp"]
Pop <- Pop[Pop!="Hisp"]

levels(str1)

#
# recast the factor !
#

Pop <- factor(Pop)
str1 <- factor(str1)

tt <- table(Pop,str1)
tt

tt <- as.matrix(tt)
tt <- t(tt)
tt

#
# Association?
#

chisq.test(tt)

fisher.test(tt,simulate.p.value = TRUE)

#
# Ternary diagram of row profiles
#

RP <- tt/rowSums(tt)


rowSums(RP)


rownames(RP)
RP
head(RP)
class(RP)
RP <- as.matrix(RP)

RP

TernaryPlot(atip = "Asian",
            btip = "Caucasian",
            ctip = "Afr-Ame"
)

TernaryPoints(RP[,c(2,3,1)],pch = 1)

AddToTernary(text, RP[,c(2,3,1)], rownames(RP), cex = 0.8, font = 2)

#
# the same with R package robCompositions
#

#ternaryDiag(RP[,c(1,3,2)],name=c("Afr-Ame","Caucasian",
#                                 "Asian"),
#            main="",mcex=1,cex.main=2,col="red",pch=1,text=rownames(RP),cex=2)

#dev.off()

library(ca)

out <- ca(tt)
summary(out)
plot(out) # not a biplot

plot(out,map="rowprincipal") # biplot of row profiles

plot(out,map="rowgreen") # contribution biplot

#
# Now the four populations (no longer fits in 2D)
#

str1 <- paste(Z$F13B_1,Z$F13B_2,sep="/")

Pop <- Z$Pop
tt <- table(Pop,str1)
tt <- t(tt)

chisq.test(tt)

fisher.test(tt,simulate.p.value = TRUE)

out <- ca(tt)
summary(out)
plot(out) # not a biplot

plot(out,map="rowprincipal") # biplot of row profiles

plot(out,map="rowgreen") # contribution biplot

Dc <- diag(out$colmass)
Dc

Gs <- sqrt(Dc)%*%out$colcoord

arrows(0,0,Gs[,1],Gs[,2],length=0.1,col="red")


