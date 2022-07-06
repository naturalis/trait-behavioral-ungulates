# Sociality
install.packages("rr2")
library(rr2)

# Prepare dataset
responsev <- subset(binData, select=c(Sociality))
glmData <- cbind(responsev, Predictors, ungulatesTraits$CanonicalName)
colnames(glmData)[colnames(glmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
# Only the rows without any missing values are selected
glmData <- glmData[complete.cases(glmData),] 
# Dropping species from the tree
dropTips <- setdiff(tree$tip.label, glmData$CanonicalName) 
# Phylogenetic tree for THIS model
modelTree <- drop.tip(tree, dropTips) 
glmData$CanonicalName <- NULL

# Full model
formula <-  Sociality ~ bio14 + bio15

## MODEL ANALYSIS
  
# Using the phyloglm
phyloGLM <- phyloglm(formula, data = glmData, phy = modelTree,
                      method = "logistic_IG10", btol = 36, log.alpha.bound = 4)
print(summary(phyloGLM))


## Partial R2
# Compare full model and reduced models

R2.lik(phyloGLM) # 0.3030795

# Reduced model: without ecological variables
r_formula <- Sociality ~ 1
r_phyloGLM <- phyloglm(r_formula, data = glmData, phy = modelTree,
                       method = "logistic_IG10", btol = 36, log.alpha.bound = 4)
R2.lik(phyloGLM, r_phyloGLM) # 0.2252159

# Reduced model: without phylogeny
r_glm_formula <- Sociality ~ bio14 + bio15
r_glm <- glm(r_glm_formula, family=binomial(link='logit'), data = glmData)
R2.lik(phyloGLM, r_glm) # 0.0420794

