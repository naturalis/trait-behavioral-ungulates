### Model Selection for Numeric Traits

mSelection_num <- function(x, c1, c2){ 
  # "x" is for dataset that contains every dependent variable,
  # "c1" is for dataset that contains the independent variables, 
  # "c2" is for the phylogeny
  # Output the summary information
  sink("model_selection_num_0510.txt", append = TRUE)
  
  ## Data Preparation
  
  # The result of 'apply' is type list. Change x to dataframe with type numeric
  dim(x)<-c(148,1)
  x<- as.numeric(x)
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

  rm(dropTips, formula, glmData, modelTree, modelstep)
  
}

apply(numData, 2, mSelection_num, c1=Predictors, c2=tree)
sink()