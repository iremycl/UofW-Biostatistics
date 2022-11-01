#
# R script session S12 on discriminant analysis (second part)

rm(list=ls())

#
# LDA and confusion matrix and classification rates
#

#
# We re-read the data and redo the LDA and QDA
#

library(calibrate)
library(HardyWeinberg)
data(Markers)
head(Markers)

#
# genotype coding
#

colvec <- rep("black",nrow(Markers))
colvec[Markers[,1]=="GG"] <- "blue"
colvec[Markers[,1]=="TT"] <- "red"
colvec[Markers[,1]=="GT"] <- "green"

#
# carrier coding
#

colvec <- rep("black",nrow(Markers))
colvec[Markers[,1]=="GG"] <- "red"
colvec[Markers[,1]=="TT"] <- "green"
colvec[Markers[,1]=="GT"] <- "green"


table(Markers$SNP1)

Markers$carrierT <- factor(Markers$SNP1=="TT" | Markers$SNP1=="GT",
                           labels=c("non T carrier","T carrier"))
table(Markers$carrierT,Markers$SNP1)

colnames(Markers)
Newdata <- Markers[is.na(Markers$SNP1),]
Newdata

colvec <- colvec[!is.na(Markers$SNP1)]
Train <- Markers[!is.na(Markers$SNP1),]

#
# Make linear classifier with the training data
#

out.lda <- lda(carrierT~iT+iG,data=Train)
out.lda
post.lda <- predict(out.lda)$posterior

out.qda <- qda(carrierT~iT+iG,data=Train)
out.qda

#
# classification/error rates
#

pre.lda <- predict(out.lda)
confusion.lda <- table(Train$carrierT,pre.lda$class)
confusion.lda

aper.lda <- (confusion.lda[1,2]+confusion.lda[2,1])/sum(confusion.lda)
aper.lda

#
# Applying jackknife/hold-one-out/cross validation
#

X <- Train
colnames(X)
n <- nrow(X)
nmisclas <- 0
for(i in 1:n) {
  ho <- X[i,]
  out <- lda(carrierT~iG+iT,data=X[-i,])
  yhat <- predict(out,newdata = ho)
  if(yhat$class!=ho$carrierT) nmisclas<-nmisclas+1
}

e_aer.lda <- nmisclas/n
e_aer.lda

out.qda <- qda(carrierT~iT+iG,data=Train)
out.qda

pre.qda <- predict(out.qda)

confusion.qda <- table(Train$carrierT,pre.qda$class)
confusion.qda

aper.qda <- (confusion.qda[1,2]+confusion.qda[2,1])/sum(confusion.qda)
aper.qda

n <- nrow(X)
nmisclas <- 0
for(i in 1:n) {
  ho <- X[i,]
  out <- qda(carrierT~iG+iT,data=X[-i,])
  yhat.qda <- predict(out,newdata = ho)
  if(yhat.qda$class!=ho$carrierT) nmisclas<-nmisclas+1
}

e_eae.qda <- nmisclas/n
e_eae.qda

#
# cross-validation with R's LDA function
#

out.cv <- lda(carrierT~iG+iT,data=X,CV=TRUE)
post.cv <- round(out.cv$posterior,digits=2)

#
# Compare
#

head(round(cbind(post.lda,post.cv),digits=3))

#
# Three group problem with LDA
#

rm(list=ls())
fn <- "http://www-eio.upc.es/~jan/SISG/NistDist.rda"
load(url(fn))

ls()

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

pop[pop=="AA"] <- "African American"
pop[pop=="Cauc"] <- "Caucasian"

popf <- factor(pop)
length(popf)
class(out.mds)
PC <- out.mds[,1:10]
head(PC)
colnames(PC) <- paste("PC",1:10,sep="")
out.lda <- lda(popf~PC)
attributes(out.lda)

out.lda$means

out.lda$prior

out.lda$counts

#
# Explained variability
#

sv <- out.lda$svd
la <- sv*sv
fr <- la/sum(la)
cf <- cumsum(fr)

Ta <- rbind(la,fr,cf)
rownames(Ta) <- c("Eigenvalue","Fraction","Cumulative")
Ta

lda.pre <- predict(out.lda)

LD <- lda.pre$x

fr1 <- paste("(",round(100*0.7389,digits=2),"%)",sep="") 
fr2 <- paste("(",round(100*0.2611,digits=2),"%)",sep="") 

plot(LD[,1],LD[,2],asp=1,col=colvec,
     xlab=paste("LD1",fr1),ylab=paste("LD2",fr2),
     main="LDA")
legend("topleft",c("African American","Asian","Caucasian"),col=c("black","yellow","blue"),
       pch=1,cex=0.75)
origin(lty="dotted")

C.lda <- table(popf,lda.pre$class)
C.lda
cr.lda <- sum(diag(C.lda))/sum(C.lda)
cr.lda

er.lda <- 1-cr.lda 
er.lda

100*cr.lda
100*er.lda




