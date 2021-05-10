### Model Selection for Factor Traits that Contain 2 Levels

mSelection_binary <- function(x, c1, c2){ 
  # "x" is for dataset that contains every dependent variable,
  # "c1" is for dataset that contains the independent variables, 
  # "c2" is for the phylogeny
  
  # Output the summary information
  sink("model_selection_factor_binary.txt", append = TRUE)
  
  ## Data Preparation
  
  # The result of 'apply' is type list. Change x to dataframe with type factor
  dim(x)<-c(148,1)
  x<- as.factor(x)
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
  
  # Using the phyloglmstep
  modelstep <- phyloglmstep(formula, starting.formula = NULL, data = glmData, phy=modelTree, 
                              method= "logistic_IG10", direction = "forward", trace = 2, 
                              btol = 20, log.alpha.bound = 4, start.beta=NULL, 
                              start.alpha=NULL, boot = 0, full.matrix = TRUE, k=2)
  summary(modelstep)
  
  # Remove the intermediate variable in preparation for the next model
  rm(glmData, dropTips, modelTree, formula, modelstep)
  
}

apply(binData, 2, mSelection_binary, c1=Predictors, c2=tree)
sink()
