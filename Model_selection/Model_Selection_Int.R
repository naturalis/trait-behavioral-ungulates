### Model Selection for Integer Traits

## 1. Count data
# For count data, they follow the Poisson distribution. 
# Since there is no poisson method provided in phyloglmstep, we will use all the predictors
# and do the model analysis directly for these traits.

intData_count <- subset(intData, select=c(X6.1_DietBreadth, X12.1_HabitatBreadth,
                                            X24.1_TeatNumber))

# Output the summary information
sink("model_selection_int_count.txt", append = TRUE)

mAnalysis_int_count <- function(x, c1, c2){ 
  # "x" is for dataset that contains every dependent variable,
  # "c1" is for dataset that contains the independent variables, 
  # "c2" is for the phylogeny
  
  ## Data Preparation
  
  # The result of 'apply' is type list. Change x to dataframe with type numeric
  dim(x)<-c(148,1)
  x<- as.integer(x)
  # Construct a data set "glmData" that is specific to THIS model (i.e. for this response trait)
  glmData <- cbind(x, c1, ungulatesTraits$CanonicalName)
  colnames(glmData)[colnames(glmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
  # Only the rows without any missing values are selected
  glmData <- glmData[complete.cases(glmData),] 
  # Dropping species from the tree
  dropTips <- setdiff(c2$tip.label, glmData$CanonicalName) 
  # Phylogenetic tree for THIS model
  modelTree <- drop.tip(c2, dropTips) 
  glmData$CanonicalName <- NULL
  
  formula <- x ~ Aspect + BulkDensity + ClayPercentage + PETDriestQuarter + 
    PETseasonality + PETWettestQuarter + OrganicCarbon + PhCaCL + Slope + bio2 + bio14 +
    bio15 + bio18 + bio19
  
  ## MODEL ANALYSIS
  
  #Using the phyloglmstep
  phyloGLM <- phyloglm(formula, data = glmData, phy=modelTree, method= "poisson_GEE", 
                       start.beta=NULL, full.matrix = TRUE)
  print(summary(phyloGLM))
  # print("One trait is finished.")
  
  rm(dropTips, formula, glmData, modelTree, modelstep)
  
}

apply(intData_count, 2, mAnalysis_int_count, c1=Predictors, c2=tree)
sink()



## 2. Ordinal data

intData_ordn <- subset(intData, select= -c(X6.1_DietBreadth, X12.1_HabitatBreadth,
                                          X24.1_TeatNumber))

# Output the summary information
sink("model_selection_int_ordinal.txt", append = TRUE)

mSelection_int <- function(x, c1, c2){ 
  # "x" is for dataset that contains every dependent variable,
  # "c1" is for dataset that contains the independent variables, 
  # "c2" is for the phylogeny
  
  ## Data Preparation
  
  # The result of 'apply' is type list. Change x to dataframe with type numeric
  dim(x)<-c(148,1)
  x<- as.integer(x)
  # Construct a data set "glmData" that is specific to THIS model (i.e. for this response trait)
  glmData <- cbind(x, c1, ungulatesTraits$CanonicalName)
  colnames(glmData)[colnames(glmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
  # Only the rows without any missing values are selected
  glmData <- glmData[complete.cases(glmData),] 
  # Dropping species from the tree
  dropTips <- setdiff(c2$tip.label, glmData$CanonicalName) 
  # Phylogenetic tree for THIS model
  modelTree <- drop.tip(c2, dropTips) 
  glmData$CanonicalName <- NULL
  
  formula <- x ~ Aspect + BulkDensity + ClayPercentage + PETDriestQuarter + 
    PETseasonality + PETWettestQuarter + OrganicCarbon + PhCaCL + Slope + bio2 + bio14 +
    bio15 + bio18 + bio19
  
  ## MODEL SELECTION
  
  #Using the phyloglmstep
  modelstep <- phylostep(formula, starting.formula = NULL, data = glmData,
                         phy = modelTree, model = "BM", direction = "forward", trace = 2, 
                         lower.bound = NULL, upper.bound = NULL, starting.value = NULL, k=2)
  print(summary(modelstep))
  # print("One trait is finished.")
  
  rm(dropTips, formula, glmData, modelTree, modelstep)
  
}

apply(intData_ordn, 2, mSelection_int, c1=Predictors, c2=tree)
sink()