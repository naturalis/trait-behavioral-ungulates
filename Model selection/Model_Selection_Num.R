### Model Selection for Numeric Traits


## 1. Data that need transformation

# The gap between the maximum and minimum values of some traits is too large (at least 10^2). 
# In order to make it closer to a normal distribution to get a better regression result, we did a
# centred log-ratio tranformation.
numData_clr <- subset(numData, select = c(X5.1_AdultBodyMass_g, X21.1_PopulationDensity_n_km2, 
                                          X22.1_HomeRange_km2, CarryWeight, AVGTravelDistance,
                                          X10.2_SocialGrpSize))
# Since centred log-ratio transformation requires there is no NA and 0, 
# AVGFoodConsumption will be calculated seperately at the end.

# Output the summary information
sink("model_selection_num_clr.txt", append = TRUE)

mSelection_num_clr <- function(x, c1, c2){ 
  # "x" is for dataset that contains every dependent variable,
  # "c1" is for dataset that contains the independent variables, 
  # "c2" is for the phylogeny
  
  ## Data Preparation
  
  # The result of 'apply' is type list. Change x to dataframe with type numeric
  dim(x)<-c(148,1)
  x<- as.numeric(x)
  # Construct a data set "lmData" that is specific to THIS model (i.e. for this response trait)
  lmData <- cbind(x, c1, ungulatesTraits$CanonicalName)
  colnames(lmData)[colnames(lmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
  # Only the rows without any missing values are selected
  lmData <- lmData[complete.cases(lmData),] 
  # Dropping species from the tree
  dropTips <- setdiff(c2$tip.label, lmData$CanonicalName) 
  # Phylogenetic tree for THIS model
  modelTree <- drop.tip(c2, dropTips) 
  lmData$CanonicalName <- NULL
  
  # Centred log-ratio transformation
  geoMean <- exp(mean(log(lmData$x)))
  lmData$x <- log(lmData$x/geoMean)
  
  formula <- x ~ Aspect + BulkDensity + ClayPercentage + PETDriestQuarter + PETseasonality + PETWettestQuarter + OrganicCarbon + PhCaCL + Slope + bio2 + bio14 +bio15 + bio18 + bio19
  
  ## MODEL SELECTION
  
  #Using the phyloglmstep
  modelstep <- phylostep(formula, starting.formula = NULL, data = lmData,
                         phy = modelTree, model = "BM", direction = "forward", trace = 2, 
                         lower.bound = NULL, upper.bound = NULL, starting.value = NULL, k=2)
  print(summary(modelstep))

  rm(dropTips, formula, lmData, modelTree, geoMean, modelstep)
  
}

apply(numData_clr, 2, mSelection_num_clr, c1=Predictors, c2=tree)
sink()


## 2. Data that don't need a transformation

numData_orign <- subset(numData, select = -c(X5.1_AdultBodyMass_g, X21.1_PopulationDensity_n_km2, 
                                          X22.1_HomeRange_km2, AVGFoodConsumption, 
                                          CarryWeight, AVGTravelDistance, X10.2_SocialGrpSize))
# Output the summary information
sink("model_selection_num_origin.txt", append = TRUE)

mSelection_num <- function(x, c1, c2){ 
  # "x" is for dataset that contains every dependent variable,
  # "c1" is for dataset that contains the independent variables, 
  # "c2" is for the phylogeny
  
  ## Data Preparation
  
  # The result of 'apply' is type list. Change x to dataframe with type numeric
  dim(x)<-c(148,1)
  x<- as.numeric(x)
  # Construct a data set "lmData" that is specific to THIS model (i.e. for this response trait)
  lmData <- cbind(x, c1, ungulatesTraits$CanonicalName)
  colnames(lmData)[colnames(lmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
  # Only the rows without any missing values are selected
  lmData <- lmData[complete.cases(lmData),] 
  # Dropping species from the tree
  dropTips <- setdiff(c2$tip.label, lmData$CanonicalName) 
  # Phylogenetic tree for THIS model
  modelTree <- drop.tip(c2, dropTips) 
  lmData$CanonicalName <- NULL
  
  formula <- x ~ Aspect + BulkDensity + ClayPercentage + PETDriestQuarter + PETseasonality + PETWettestQuarter + OrganicCarbon + PhCaCL + Slope + bio2 + bio14 +bio15 + bio18 + bio19
  
  ## MODEL SELECTION
  
  #Using the phyloglmstep
  modelstep <- phylostep(formula, starting.formula = NULL, data = lmData,
                         phy = modelTree, model = "BM", direction = "forward", trace = 2, 
                         lower.bound = NULL, upper.bound = NULL, starting.value = NULL, k=2)
  print(summary(modelstep))
  
  rm(dropTips, formula, lmData, modelTree, modelstep)
  
}

apply(numData_orign, 2, mSelection_num, c1=Predictors, c2=tree)
sink()


## AVGFoodConsumption

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

formula <- AVGFoodConsumption ~ Aspect + BulkDensity + ClayPercentage + PETDriestQuarter + PETseasonality + PETWettestQuarter + OrganicCarbon + PhCaCL + Slope + bio2 + bio14 +bio15 + bio18 + bio19

## MODEL SELECTION

#Using the phyloglmstep
modelstep <- phylostep(formula, starting.formula = NULL, data = lmData,
                       phy = modelTree, model = "BM", direction = "forward", trace = 2, 
                       lower.bound = NULL, upper.bound = NULL, starting.value = NULL, k=2)
print(summary(modelstep))

rm(responsev, dropTips, geoMean, formula, lmData, modelTree, modelstep)
