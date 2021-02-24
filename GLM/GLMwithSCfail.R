MatrixSCfail <- read.csv("MatrixSCfail.csv", header=T)

#Collinearity suspected, visualise:
library('corrplot')
M<-cor(MatrixSCfail[,2:8])
corrplot.mixed(M)
#At first look, Haplo/spec, stability/max_p, and zones/sites are strongly correlated so collinearity must be examined when creating the model

#Create the first model with all variables included
modF <- glm(verified ~ haplo+spec+stability+max_p+zones+sites, data=MatrixSCfail, family=quasibinomial, maxit=100)
par(mfrow=c(2,2))
plot(modF)
summary(modF)

#Use Farrar-Glauber test 1 to detect collinearity
library(mctest)
omcdiag(modF)

#All tests detect collinearity, Chi-square very high, move to F-G step 2 to locate collinearity
imcdiag(modF)

#VIF of variable haplo is high, remove from model
modF2<-update(modF,~. -haplo)

#Re-run F-G test 1 on new model
omcdiag(modF2)
#Collinearity still detected, Chi-square still high

#Re-run F-G test 2 to locate collinearity
imcdiag(modF2)
#VIF values are all under 10, but collinearity is detected for stability, max_p, and sites 

#max_p has highest VIF value, remove from model
modF3<-update(modF2,~. - max_p)

#Re-run F-G test 1 on new model
omcdiag(modF3)
#Collinearity still detected by Chi-square and value is still high

#Re-run F-G test 2 to locate collinearity
imcdiag(modF3)
#Collinearity still detected for zones and sites

#Sites has higher VIF value, remove from model
modF4<-update(modF3,~. -sites)

#Re-run F-G test 1 on new model
omcdiag(modF4)
#Collinearity still detected, tun F-G part 2 to locate
imcdiag(modF4)
#VIF values are all low (<1.5), no collinearity detected for spec, stability, or zones

summary(modF4)
#shows that only stability is statistically significant