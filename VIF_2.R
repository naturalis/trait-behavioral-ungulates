## BEGIN VALUES
# Do the VIF analysis for "ecology"
Predictors <- ecology

# The first correlation matrix is a visualization of the dataset without the removal of any 
# traits. There is a lot of collinearity present. The corrplot package is used to 
# visualize the collinearity between the traits, using the correlation matrix.
matrix <- cor(Predictors, use = "pairwise.complete.obs")
corrplot::corrplot(matrix, type="lower", order = "hclust", tl.pos = "l", tl.col = "black", tl.cex = 0.6)


## CHECK CORRELATION MATRIX
# With the findCorrelation() function, the columns with the highest collinearity are found.
# By visualizing the matrix, the outcome of the findCorrelation function can be easily checked.
# The VIF function can be run to easily check if the results are still infinite.
# A cutoff value of .83 is chosen, using trial and error.

# matrix <- cor(Predictors, use = "pairwise.complete.obs")
delColumn <- findCorrelation(matrix, cutoff = .83, verbose = FALSE, exact=FALSE)
Predictors <- Predictors[,-c(delColumn)]
vif(Predictors)

# Check CORRELATION MATRIX
# A visualization of the dataset after the cutoff value is implemented.
matrix <- cor(Predictors, use = "pairwise.complete.obs")
corrplot::corrplot(matrix, type="lower", order = "hclust", tl.pos = "l", tl.col = "black", tl.cex = 0.8)


## VIF ANALYSIS
# Examine the VIF of the variables in the filtered Predictor in turn. 
# If the value of VIF is greater than 10, then remove the variable with the highest VIF. 
# Repeat this operation until all variables have VIF values less than 10.
# Check VIF -> PETWarmestQuarter has highest VIF value above 10
# Remove PETWarmestQuarter
vif(Predictors)
Predictors$PETWarmestQuarter <- NULL

# Check VIF -> minTempWarmest has highest VIF value above 10
# Remove minTempWarmest
vif(Predictors)
Predictors$minTempWarmest <- NULL

# Check VIF -> bio3 has highest VIF value above 10
# Remove bio3
vif(Predictors)
Predictors$bio3 <- NULL

## FINAL PREDICTOR SET
# Final set of predictors, containing: Aspect, BulkDensity, ClayPercentage, 
# PETDriestQuarter, PETseasonality, PETWettestQuarter, OrganicCarbon, PhCaCL, 
# Slope, bio2, bio14, bio15, bio18, bio19
Predictors

## FINAL VALUES
# The final correlation matrix is the final product after removing all highly correlated
# traits.
matrix <- cor(Predictors, use = "pairwise.complete.obs")
corrplot::corrplot(matrix, type="lower", order = "hclust", tl.pos = "l", tl.col = "black", tl.cex = 0.9)

## REMOVE VARIABLES
rm(matrix, delColumn)
