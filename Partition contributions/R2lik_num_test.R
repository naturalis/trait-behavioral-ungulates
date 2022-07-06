### Model Analysis for Numeric Traits
# LitterSize: doesn't need clr-transformation
# AVGFoodConsumption: needs clr-transformation but contains value '0'
# PopulationDensity: needs clr-transformation
# thus, these three traits will be treated respectively

# R2.lik Function (using StarTree)
var_par_num_stree <- function(y, c1, c2, c3){
  ## Partial R2
  # y - formula for full modedl
  # c1 - dataset
  # c2 - full tree
  # c3 - full model
  print('VARIATION PARTITION')
  
  # formula
  # Change the type of y to formula
  f <- as.formula(y)
  # formula for reduced model
  r_formula <- as.formula(paste(colnames(c1)[1], "~ 1"))
  
  # phylogeny for reduced model
  # Build a star phylogeny
  b_length <- rep(10, times = length(c2$tip.label))
  r_tree <- starTree(c1$CanonicalName, branch.lengths = b_length)
  
  # DV ~ PHY
  rmodel_WITHOUTeco <- phylolm(r_formula, data = c1, phy = c2, model = "BM",
                               lower.bound = NULL, upper.bound = NULL, starting.value = NULL, 
                               measurement_error = FALSE, boot=0, full.matrix = TRUE)
  # print('Reduced Model Summary Info: ')
  # print(summary(rmodel_WITHOUTeco))
  
  print('***Variation that ecological variables explain:')
  print(R2.lik(c3, rmodel_WITHOUTeco))
  
  # DV ~ E
  rmodel_WITHOUTphy <- phylolm(formula = f, data = c1, phy = r_tree, model = "BM",
                               lower.bound = NULL, upper.bound = NULL, starting.value = NULL, 
                               measurement_error = FALSE, boot=0, full.matrix = TRUE)
  # print('Reduced Model Summary Info: ')
  # print(summary(rmodel_WITHOUTphy))
  
  print('***Variation that phylogenetic conservatism explain:')
  print(R2.lik(c3, rmodel_WITHOUTphy))
  
  rm(f, r_formula, b_length, r_tree, rmodel_WITHOUTeco, rmodel_WITHOUTphy)
  
}

# R2.lik Function (using normal lm)
var_par_num_lm <- function(y, c1, c2, c3){
  ## Partial R2
  # y - formula for full modedl
  # c1 - dataset
  # c2 - full tree
  # c3 - full model
  print('VARIATION PARTITION')
  
  # formula
  # Change the type of y to formula
  f <- as.formula(y)
  # formula for reduced model
  r_formula <- as.formula(paste(colnames(c1)[1], "~ 1"))
  
  # DV ~ PHY
  rmodel_WITHOUTeco <- phylolm(r_formula, data = c1, phy = c2, model = "BM",
                               lower.bound = NULL, upper.bound = NULL, starting.value = NULL, 
                               measurement_error = FALSE, boot=0, full.matrix = TRUE)
  # print('Reduced Model Summary Info: ')
  # print(summary(rmodel_WITHOUTeco))
  
  print('***Variation that ecological variables explain:')
  print(R2.lik(c3, rmodel_WITHOUTeco))
  
  # DV ~ E
  rmodel_WITHOUTphy <- lm(formula = f, data = c1)
  # print('Reduced Model Summary Info: ')
  # print(summary(rmodel_WITHOUTphy))
  
  print('***Variation that phylogenetic conservatism explain:')
  print(R2.lik(c3, rmodel_WITHOUTphy))
  
  rm(f, r_formula, rmodel_WITHOUTeco, rmodel_WITHOUTphy)
}

#### X15.1_LitterSize
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

# Phylogenetic Signal
PhySig <- phylolm(formula = f, data = lmData, phy = modelTree, model = "lambda",
                  lower.bound = NULL, upper.bound = NULL, starting.value = NULL, 
                  measurement_error = FALSE, boot=0, full.matrix = TRUE)
PhySig$optpar

# R2.lik
var_par_num_stree(f, lmData, modelTree, phyloLM)
var_par_num_lm(f, lmData, modelTree, phyloLM)

rm(responsev, lmData, dropTips, modelTree, phyloLM, f, PhySig)

####################################################################################
#### X21.1_PopulationDensity_n_per_km2
responsev <- subset(numData, select=c(X21.1_PopulationDensity_n_per_km2))

lmData <- cbind(responsev, ungulatesTraits$CanonicalName, Predictors)
colnames(lmData)[colnames(lmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
# Only the rows without any missing values are selected
lmData <- lmData[complete.cases(lmData),] 
# Dropping species from the tree
dropTips <- setdiff(tree$tip.label, lmData$CanonicalName) 
# Phylogenetic tree for THIS model
modelTree <- drop.tip(tree, dropTips) 

# Centred log-ratio transformation
geoMean <- exp(mean(log(lmData$X21.1_PopulationDensity_n_per_km2)))
lmData$X21.1_PopulationDensity_n_per_km2 <- log(lmData$X21.1_PopulationDensity_n_per_km2/geoMean)

f <- X21.1_PopulationDensity_n_per_km2 ~ PETDriestQuarter + PETseasonality + PETWettestQuarter + bio14 + bio15 + bio18

## MODEL ANALYSIS
# Using the phyloglm
phyloLM <- phylolm(formula = f, data = lmData, phy = modelTree, model = "BM",
                   lower.bound = NULL, upper.bound = NULL, starting.value = NULL, 
                   measurement_error = FALSE, boot=0, full.matrix = TRUE)
print(summary(phyloLM))

# Phylogenetic Signal
PhySig <- phylolm(formula = f, data = lmData, phy = modelTree, model = "lambda",
                  lower.bound = NULL, upper.bound = NULL, starting.value = NULL, 
                  measurement_error = FALSE, boot=0, full.matrix = TRUE)
PhySig$optpar

# Variation Partition
var_par_num_stree(f, lmData, modelTree, phyloLM)
var_par_num_lm(f, lmData, modelTree, phyloLM)

rm(responsev, lmData, dropTips, modelTree, geoMean, phyloLM, f, PhySig)

##########################################################################
#### AVGFoodConsumption
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

# Phylogenetic Signal
PhySig <- phylolm(formula = f, data = lmData, phy = modelTree, model = "lambda",
                  lower.bound = NULL, upper.bound = NULL, starting.value = NULL, 
                  measurement_error = FALSE, boot=0, full.matrix = TRUE)
PhySig$optpar

# Variation Partition
var_par_num_stree(f, lmData, modelTree, phyloLM)
var_par_num_lm(f, lmData, modelTree, phyloLM)

rm(responsev, lmData, dropTips, modelTree, geoMean, phyloLM, f, PhySig)
