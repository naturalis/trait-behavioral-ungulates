### Model Analysis for Factor Traits Contain two levels

# Used to iterate over response variable
i <- 1
# Output the summary information
sink("model_analysis_factor_binary_IG10_0527.txt", append = TRUE)
# sink("model_analysis_factor_binary_MPLE_0525.txt", append = TRUE)

mAnalysis_binary <- function(y, c1, c2, c3){ 
  # 'y' is for every formula,
  # 'c1' is for dataset that contains every dependent variable,
  # "c2" is for dataset that contains the independent variables, 
  # "c3" is for the phylogeny
  print(i)
  
  ## Data Preparation
  
  # Change the type of y to formula
  f <- as.formula(y)
  # Construct a data set "glmData" that is specific to THIS model (i.e. for this response trait)
  glmData <- cbind(c1[,i], ungulatesTraits$CanonicalName, c2)
  colnames(glmData)[1] <- 'x'
  colnames(glmData)[colnames(glmData)=="ungulatesTraits$CanonicalName"]<-c("CanonicalName")
  # Only the rows without any missing values are selected
  glmData <- glmData[complete.cases(glmData),] 
  # Dropping species from the tree
  dropTips <- setdiff(c3$tip.label, glmData$CanonicalName) 
  # Phylogenetic tree for THIS model
  modelTree <- drop.tip(c3, dropTips) 
  glmData$CanonicalName <- NULL
  print(dim(glmData))
  
  ## MODEL ANALYSIS
  
  # Using the phyloglm
  phyloGLM <- phyloglm(formula = f, data = glmData, phy = modelTree, 
                         method = "logistic_IG10", btol = 36, log.alpha.bound = 4)
  # phyloGLM <- phyloglm(formula = f, data = glmData, phy = modelTree, 
                         method = "logistic_MPLE", btol = 36, log.alpha.bound = 4)
  print(summary(phyloGLM))
    
  i <<- i+1
  
  rm(dropTips, f, glmData, modelTree,  phyloGLM)
  
}

apply(binFormulas, 1, mAnalysis_binary, c1=sig_binData, c2=Predictors, c3=tree)
print(i)
print('Finish')
sink()
rm(i)
