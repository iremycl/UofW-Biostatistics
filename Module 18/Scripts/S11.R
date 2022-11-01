#
# R script session S11 on discriminant analysis.

rm(list=ls())

#
# LDA
#

library(MASS)
library(HardyWeinberg)
data(Markers)
head(Markers)

colvec <- rep("black",nrow(Markers))
colvec[Markers[,1]=="GG"] <- "blue"
colvec[Markers[,1]=="TT"] <- "red"
colvec[Markers[,1]=="GT"] <- "green"

plot(Markers[,2],Markers[,3],col=colvec,xlab="Intensity G",ylab="Intensity T")
legend("topright",c("GG","GT","TT","NA"),pch=1,col=c("blue","green","red","black"))

colvec <- rep("black",nrow(Markers))
colvec[Markers[,1]=="GG"] <- "red"
colvec[Markers[,1]=="TT"] <- "green"
colvec[Markers[,1]=="GT"] <- "green"

plot(Markers[,2],Markers[,3],col=colvec,xlab="Intensity G",ylab="Intensity T")
legend("topright",c("T carrier","non-carrier","NA"),pch=1,col=c("green","red","black"))


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

out <- lda(carrierT~iT+iG,data=Train)
out

pp <- out$prior
out$means
out$scaling

stats <- cbind(pp,out$means)
stats

plot(out)

results <- predict(out)

table(results$class)

m1 <- mean(results$x[results$class=="T carrier"])
m1

m2 <- mean(results$x[results$class=="non T carrier"])
m2


plot(results$x,jitter(as.numeric(results$class),),xlab="Linear discriminant",ylim=c(0,3),
     ylab="Group",col=colvec)
points(m1,2,pch=19,col="green")
points(m2,1,pch=19,col="red")
abline(h=2,lty="dotted")
abline(h=1,lty="dotted")
legend("topleft",c("Non-carrier","Carrier"),pch=1,col=c("red","green"))

#
# Prediction
#

hat.lda <- predict(out,newdata=Newdata)

LDAPred <- data.frame(hat.lda$class,hat.lda$x,hat.lda$posterior)

head(LDAPred)

Sc <- cov(Train[Train$carrierT=="non T carrier",2:3])
Sc
cor(Train[Train$carrierT=="non T carrier",2:3])

Snc <- cov(Train[Train$carrierT=="T carrier",2:3])
Snc
cor(Train[Train$carrierT=="T carrier",2:3])

#
# Quadratic discriminant analysis
#

out.qda <- qda(carrierT~iT+iG,data=Train)
out.qda

hat.qda <- predict(out.qda)
table(hat.qda$class)

newpr.qda <- predict(out.qda,newdata=Newdata)
newpr.qda

QDAPred <- data.frame(newpr.qda$class,newpr.qda$posterior)
head(QDAPred)

ZZZ <- data.frame(LDAPred,QDAPred)
head(ZZZ)

table(hat.qda$class)

qprd <- hat.qda$class

symvec <- rep(1,nrow(Train))
symvec[qprd=="T carrier"] <- 2
symvec[qprd=="non T carrier"] <- 3

symvecnew <- rep(1,nrow(Train))
symvecnew[newpr.qda$class=="T carrier"] <- 2
symvecnew[newpr.qda$class=="non T carrier"] <- 3

plot(Train[,2],Train[,3],col=colvec,xlab="Intensity G",ylab="Intensity T",pch=symvec)
legend("topright",c("T carrier","non-carrier","NA"),pch=c(2,3),col=c("green","red","black"))

points(Newdata[,2],Newdata[,3],col="black",pch=symvecnew)




