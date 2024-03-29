### Model Analysis for Numeric Traits
# LitterSize: doesn't need clr-transformation
# AVGFoodConsumption: needs clr-transformation but contains 0 value
# PopulationDensity: needs clr-transformation
# thus, these three traits will be treated respectively
## The general process for general data is at the end of the file 

# X15.1_LitterSize
responsev <- subset(numData, select=c(X15.1_LitterSize))

lmData <- cbind(responsev, ungulatesTraits$CanonicalName, Predictors)
colnames(lmData)[colnames(lmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
# Only the rows without any missing values are selected
lmData <- lmData[complete.cases(lmData),] 
# Dropping species from the tree
dropTips <- setdiff(tree$tip.label, lmData$CanonicalName) 
# Phylogenetic tree for THIS model
modelTree <- drop.tip(tree, dropTips) 

f <- X15.1_LitterSize ~ BulkDensity + ClayPercentage + PhCaCL + bio2 + bio14 + bio15

## MODEL ANALYSIS
# Using the phyloglm
phyloLM <- phylolm(formula = f, data = lmData, phy = modelTree, model = "BM",
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
plot(lmdatamrg$fitted, lmdatamrg$redisuals, xlab = "Fitted value", ylab = "Residuals")

rm(dropTips, f, lmData, modelTree,  phyloLM, wn, lmdatamrg)


# X21.1_PopulationDensity_n_km2
responsev <- subset(numData_clr, select=c(X21.1_PopulationDensity_n_km2))

lmData <- cbind(responsev, ungulatesTraits$CanonicalName, Predictors)
colnames(lmData)[colnames(lmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
# Only the rows without any missing values are selected
lmData <- lmData[complete.cases(lmData),] 
# Dropping species from the tree
dropTips <- setdiff(tree$tip.label, lmData$CanonicalName) 
# Phylogenetic tree for THIS model
modelTree <- drop.tip(tree, dropTips) 

f <- X21.1_PopulationDensity_n_km2 ~ PETDriestQuarter + PETseasonality + PETWettestQuarter + bio14 + bio15 + bio18

## MODEL ANALYSIS
# Using the phyloglm
phyloLM <- phylolm(formula = f, data = lmData, phy = modelTree, model = "BM",
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
plot(lmdatamrg$fitted, lmdatamrg$redisuals, xlab = "Fitted value", ylab = "Residuals")

rm(dropTips, f, lmData, modelTree,  phyloLM, wn, lmdatamrg)


# AVGFoodConsumption
responsev <- subset(numData, select=c(AVGFoodConsumption))

lmData <- cbind(responsev, ungulatesTraits$CanonicalName, Predictors)
colnames(lmData)[colnames(lmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
# Only the rows without any missing values are selected
lmData <- lmData[complete.cases(lmData),] 
lmData <- lmData[!(row.names(lmData) %in% c("Cephalophus_natalensis","Tragulus_napu")), ]
# Dropping species from the tree
dropTips <- setdiff(tree$tip.label, lmData$CanonicalName) 
# Phylogenetic tree for THIS model
modelTree <- drop.tip(tree, dropTips) 

# Centred log-ratio transformation
# Geometrical Mean 3.167097
geoMean <- exp(mean(log(lmData$AVGFoodConsumption)))
lmData$AVGFoodConsumption <- log(lmData$AVGFoodConsumption/geoMean)

f <- AVGFoodConsumption ~ Aspect + PETDriestQuarter + PETseasonality + PETWettestQuarter + PhCaCL + bio2 + bio14 + bio18 + bio19

## MODEL ANALYSIS
# Using the phyloglm
phyloLM <- phylolm(formula = f, data = lmData, phy = modelTree, model = "BM",
                   lower.bound = NULL, upper.bound = NULL, starting.value = NULL, 
                   measurement_error = FALSE, boot=0, full.matrix = TRUE)
print(summary(phyloLM))
# plot(phyloLM)
# print(head(phyloLM$residuals))
# print(head(phyloLM$fitted.values))

wn = data.frame(CanonicalName = names(phyloLM$residuals), redisuals = phyloLM$residuals,
                fitted = phyloLM$fitted.values)
lmdatamrg <- merge(x = lmData, y = wn, by = "CanonicalName", all = T)

# Residual Plot
plot(lmdatamrg$fitted, lmdatamrg$redisuals, xlab = "Fitted value", ylab = "Residuals")

rm(dropTips, f, lmData, modelTree,  phyloLM, wn, lmdatamrg)



## General process
# Used to iterate over response variable
i <- 1
# Output the summary information
sink("model_analysis_num.txt", append = TRUE)

mAnalysis_num <- function(y, c1, c2, c3){ 
  # 'y' is for every formula,
  # 'c1' is for dataset that contains every dependent variable,
  # "c2" is for dataset that contains the independent variables, 
  # "c3" is for the phylogeny
  print(i)
  
  ## Data Preparation
  
  # Change the type of y to formula
  f <- as.formula(y)
  # Construct a data set "lmData" that is specific to THIS model (i.e. for this response trait)
  lmData <- cbind(c1[,i], ungulatesTraits$CanonicalName, c2)
  colnames(lmData)[1] <- 'x'
  colnames(lmData)[colnames(lmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
  # Only the rows without any missing values are selected
  lmData <- lmData[complete.cases(lmData),] 
  # Dropping species from the tree
  dropTips <- setdiff(c3$tip.label, lmData$CanonicalName) 
  # Phylogenetic tree for THIS model
  modelTree <- drop.tip(c3, dropTips) 
  lmData$CanonicalName <- NULL
  print(dim(lmData))
  
  ## MODEL ANALYSIS
  
  # Using the phyloglm
  phyloLM <- phylolm(formula = f, data = lmData, phy = modelTree, model = "BM",
                     lower.bound = NULL, upper.bound = NULL, starting.value = NULL, 
                     measurement_error = FALSE, boot=0, full.matrix = TRUE)
  print(summary(phyloLM))
  
  i <<- i+1
  
  rm(dropTips, f, lmData, modelTree,  phyloLM)
  
}

apply(numFormulas, 1, mAnalysis_num, c1=sig_numData, c2=Predictors, c3=tree)
print(i)
print('Finish')
sink()
rm(i)