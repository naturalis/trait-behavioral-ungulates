# SocialHierarchy
responsev <- subset(binData, select=c(SocialHierarchy))

glmData <- cbind(responsev, Predictors, ungulatesTraits$CanonicalName)
colnames(glmData)[colnames(glmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
# Only the rows without any missing values are selected
glmData <- glmData[complete.cases(glmData),] 
# Dropping species from the tree
dropTips <- setdiff(tree$tip.label, glmData$CanonicalName) 
# Phylogenetic tree for THIS model
modelTree <- drop.tip(tree, dropTips) 
glmData$CanonicalName <- NULL

formula <-  SocialHierarchy ~ bio18 + bio19

## MODEL ANALYSIS

# Using the phyloglm
phyloGLM <- phyloglm(formula, data = glmData, phy = modelTree,
                     method = "logistic_IG10", btol = 36, log.alpha.bound = 4)
# print(summary(phyloGLM))

glmData$residuals <- residuals(phyloGLM)
glmData$fitted <- fitted.values(phyloGLM)
plot(phyloGLM)
lines(lowess(fitted,residuals),col="black",lwd=2)

# bio18
plot(glmData$bio18,residuals(phyloGLM),col=c("blue","red")) + text(glmData$bio18,residuals(phyloGLM), cex=0.6, pos=1, col="red")
lines(lowess(glmData$bio18,residuals(phyloGLM)),col="black",lwd=2)
lines(lowess(glmData$bio18[glmData$SocialHierarchy==0],residuals(phyloGLM)[glmData$SocialHierarchy==0]),col="blue")
lines(lowess(glmData$bio18[glmData$SocialHierarchy==1],residuals(phyloGLM)[glmData$SocialHierarchy==1]),col="red")
abline(h=0,lty=2,col="grey")
