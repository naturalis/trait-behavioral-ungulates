### Model Analysis for Numeric Traits

# Used to iterate over response variable
i <- 1
# Output the summary information
sink("model_analysis_factor_num.txt", append = TRUE)

mAnalysis_binary <- function(y, c1, c2, c3){ 
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
                     measurement_error = FALSE, boot=0,full.matrix = TRUE)
  print(summary(phyloLM))
  
  i <<- i+1
  
  rm(dropTips, f, lmData, modelTree,  phyloLM)
  
}

apply(numFormulas, 1, mAnalysis_num, c1=numData, c2=Predictors, c3=tree)
print(i)
print('Finish')
sink()
rm(i)