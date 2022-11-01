#
# Breakout S12 exercise; Day 03 session S12 Discriminant analysis
#

rm(list=ls())

fn <- "http://www-eio.upc.es/~jan/SISG/BreakoutS12.dat"

X <- read.table(fn,header=TRUE)

n <- nrow(X)

head(X)

colvec <- rep(NA,n)
colvec[X$Genotype==0] <- "blue"
colvec[X$Genotype==1] <- "red"
colvec[X$Genotype==2] <- "green"
colvec[is.na(X$Genotype)] <- "black"

plot(X$iA,X$iB,xlab="Intensity A",ylab="Intensity B",col=colvec)

Xtrain <- X[!is.na(X$Genotype),]
Xtest  <- X[is.na(X$Genotype),]
dim(Xtrain)
colnames(Xtrain)
out.lda <- lda(Genotype~iA+iB,data=Xtrain)
out.lda
post.lda <- predict(out.lda)$posterior

pre.lda <- predict(out.lda)
confusion.lda <- table(Xtrain$Genotype,pre.lda$class)
confusion.lda

aper.lda <- (confusion.lda[1,2]+confusion.lda[2,1])/sum(confusion.lda)
aper.lda

class.rate <- 1- aper.lda
class.rate

lda.pre <- predict(out.lda)

LD <- lda.pre$x

s <- out.lda$svd
la <- s*s
fr <- la/sum(la)
rbind(la,fr)

fr1 <- paste("(",round(100*fr[1],digits=3),"%)",sep="") 
fr2 <- paste("(",round(100*fr[2],digits=3),"%)",sep="") 

plot(LD[,1],LD[,2],asp=1,col=colvec,
     xlab=paste("LD1",fr1),ylab=paste("LD2",fr2),
     main="LDA")

hat.lda <- predict(out.lda,newdata=Xtest)

LDAPred <- data.frame(hat.lda$class,hat.lda$x,hat.lda$posterior)

class(hat.lda$x)

cvec <- rep(NA,nrow(Xtest))
cvec[hat.lda$class==0] <- "blue"
cvec[hat.lda$class==1] <- "red"
cvec[hat.lda$class==2] <- "green"
table(cvec)
points(hat.lda$x[,1],hat.lda$x[,2],pch=2,col=cvec,cex=0.5)
head(LDAPred)

class(hat.lda$posterior)

maxpos <- apply(hat.lda$posterior,1,max)
imaxpos <- maxpos < 0.95

points(hat.lda$x[imaxpos,1],hat.lda$x[imaxpos,2],pch=2,
       col=cvec[imaxpos],cex=2)



