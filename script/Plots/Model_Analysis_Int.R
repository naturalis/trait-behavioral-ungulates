### Model Analysis for Integer Traits

## X6.1_DietBreadth
responsev <- subset(intData, select=c(X6.1_DietBreadth))

lmData <- cbind(responsev, ungulatesTraits$CanonicalName, Predictors)
colnames(lmData)[colnames(lmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
# Only the rows without any missing values are selected
lmData <- lmData[complete.cases(lmData),] 
# Dropping species from the tree
dropTips <- setdiff(tree$tip.label, lmData$CanonicalName) 
# Phylogenetic tree for THIS model
modelTree <- drop.tip(tree, dropTips) 

formula <- X6.1_DietBreadth ~ Aspect + BulkDensity + ClayPercentage + PETDriestQuarter + 
  PETseasonality + PETWettestQuarter + OrganicCarbon + PhCaCL + Slope + bio2 + bio14 +
  bio15 + bio18 + bio19

phyloGLM <- phyloglm(formula, data = lmData, phy=modelTree, 
                      method= "poisson_GEE", start.beta=NULL, boot = 0, full.matrix = TRUE)
print(summary(phyloGLM))

# print(head(phyloGLM$residuals))
# print(head(phyloGLM$fitted.values))

wn = data.frame(CanonicalName = names(phyloGLM$residuals), res = phyloGLM$residuals,
                fitted = exp(phyloGLM$linear.predictors))
lmdatamrg <- merge(x = lmData, y = wn, by = "CanonicalName", all = T)
lmdatamrg$res <- lmdatamrg$X6.1_DietBreadth - lmdatamrg$fitted

# Residual Plot
binnedplot(lmdatamrg$fitted, 
           lmdatamrg$redisuals, 
           nclass = 16, 
           xlab = "Expected Values", 
           ylab = "Average Residual", 
           cex.pts = 0.7, 
           col.pts = 1, 
           col.int = "gray")

rm(responsev, dropTips, formula, lmData, modelTree, phyloGLM, wn, lmdatamrg)


## NaturalPredators
responsev <- subset(ungulatesTraits, select=c(NaturalPredators))
responsev$NaturalPredators <- as.integer(responsev$NaturalPredators)

lmData <- cbind(responsev, Predictors, ungulatesTraits$CanonicalName)
colnames(lmData)[colnames(lmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
# Only the rows without any missing values are selected
lmData <- lmData[complete.cases(lmData),] 
# Dropping species from the tree
dropTips <- setdiff(tree$tip.label, lmData$CanonicalName) 
# Phylogenetic tree for THIS model
modelTree <- drop.tip(tree, dropTips) 

## MODEL ANALYSIS
formula <- NaturalPredators ~ Aspect + BulkDensity + ClayPercentage + OrganicCarbon + bio2 + bio15

phyloLM <- phylolm(formula, data = lmData, phy = modelTree, model = "BM",
                   lower.bound = NULL, upper.bound = NULL, starting.value = NULL, 
                   measurement_error = FALSE, boot=0, full.matrix = TRUE)
print(summary(phyloLM))
# print(head(phyloLM$residuals))
# print(head(phyloLM$fitted.values))
# plot(phyloLM)

wn = data.frame(CanonicalName = names(phyloLM$residuals), redisuals = phyloLM$residuals,
                fitted = phyloLM$fitted.values)
lmdatamrg <- merge(x = lmData, y = wn, by = "CanonicalName", all = T)

# Residual Plot
binnedplot(lmdatamrg$fitted, 
           lmdatamrg$redisuals, 
           nclass = 16, 
           xlab = "Expected Values", 
           ylab = "Average Residual", 
           cex.pts = 0.7, 
           col.pts = 1, 
           col.int = "gray")

rm(responsev, dropTips, formula, lmData, modelTree, phyloLM, wn, lmdatamrg)
