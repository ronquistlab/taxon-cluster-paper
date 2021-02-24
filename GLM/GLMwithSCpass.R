MatrixSCpass <- read.csv("MatrixSCpass.csv", header=T)

#Visualise collinearity
library('corrplot')
MP<-cor(MatrixSCpass[,2:8])
corrplot.mixed(MP)
#At first look, Haplo/spec, stability/max_p, and zones/sites are strongly correlated so collinearity must be examined when creating the model

#Create the first model with all variables included
modP <- glm(verified ~ haplo+spec+stability+max_p+zones+sites, data=MatrixSCpass, family=quasibinomial, maxit=100)
par(mfrow=c(2,2))
plot(modP)
summary(modP)

#Use Farrar-Glauber test 1 to detect collinearity
library(mctest)
omcdiag(modP)

#All tests detect collinearity, Chi-square very high, move to F-G step 2 to locate collinearity
imcdiag(modP)

#VIF of variable haplo is high, remove from model
modP2<-update(modP,~. -haplo)

#Re-run F-G test 1 on new model
omcdiag(modP2)
#Collinearity still detected, Chi-square still high

#Re-run F-G test 2 to locate collinearity
imcdiag(modP2)
#VIF values are all under 10, but collinearity is still detected for ALL VARIABLES

#max_p has highest VIF value, remove from model
modP3<-update(modP2,~. - max_p)

#Re-run F-G test 1 on new model
omcdiag(modP3)
#Collinearity still detected 

#Re-run F-G test 2 to locate collinearity
imcdiag(modP3)
#Collinearity still detected for zones and sites

#Sites has higher VIF value, remove from model
modP4<-update(modP3,~. -sites)

#Re-run F-G test 1 on new model
omcdiag(modP4)
#Collinearity still detected, run F-G part 2 to locate
imcdiag(modP4)
#VIF values are all low (<1.5), no collinearity detected for spec, stability, or zones. 

summary(modP4)
#shows that only stability is statistically significant