# Check that lm/lmer model is valid
# Throw warning if
#	1) Intercept is ommited
#	2) Any coefficient is NA
#	3) a categorical variable is modeled as a fixed effect
setGeneric("checkModelStatus", signature="fit",
           function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999 )
             standardGeneric("checkModelStatus")
)

setMethod("checkModelStatus", "lm",
          function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999 )
          {
            # if no intercept is specified, give warning
            if( showWarnings && length(which(names(coef(fit)) == "(Intercept)")) == 0 ){
              warning("No Intercept term was specified in the formula:\nThe results will not behave as expected and may be very wrong!!")
            }
            
            # if any coefficient is NA
            if( showWarnings && any(is.na(coef(fit))) ){
              stop("The variables specified in this model are redundant,\nso the design matrix is not full rank")
            }
            
            # check colinearity
            score = colinearityScore(fit)
            if( score > colinearityCutoff ){
              stop(paste("Colinear score =", format(score, digits=4), ">", colinearityCutoff,"\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
            }
          }
)

setMethod("checkModelStatus", "lmerMod",
          function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999 )
          {
            # if no intercept is specified, give warning
            if( !dream && showWarnings && length(which(colnames(fit@pp$X) == "(Intercept)")) == 0 ){
              warning("No Intercept term was specified in the formula:\nThe results will not behave as expected and may be very wrong!!")
            }
            
            # if any coefficient is NA
            if( ( showWarnings | dream) && any(is.na(coef(fit))) ){
              stop("The variables specified in this model are redundant,\nso the design matrix is not full rank")
            }
            
            # check colinearity
            ###################
            score = colinearityScore(fit)
            
            if( score > colinearityCutoff ){
              stop(paste("Colinear score =", format(score, digits=4), ">", colinearityCutoff,"\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
            }
            
            # check that factors are random and continuous variables are fixed
            ###################################################################
            
            # remove backticks with gsub manually
            # solve issue that backticks are conserved is some but not all parts of lmer()
            
            # random variables
            randVar = as.character(attr(attr(fit@frame, "terms"), "predvars.random"))[-c(1,2)]
            randVar = gsub("`", "", randVar)
            
            # fixed variables
            fixedVar = as.character(attr(attr(fit@frame, "terms"), "predvars.fixed"))[-c(1,2)]
            fixedVar = gsub("`", "", fixedVar)
            
            varType = attr(attr(fit@frame, "terms"), "dataClasses")[-1]
            
            # variables fit by regression
            testVar = attr(attr(fit@frame, "terms"), "term.labels")
            testVar = gsub("`", "", testVar)
            
            # keep only tested variables
            varType = varType[testVar]
            
            for( i in 1:length(varType) ){
              
              # if factor is not random
              if( (showWarnings && ! dream) && varType[i] %in% c("factor", "character") && (! names(varType)[i] %in% randVar) ){
                stop(paste("Categorical variables modeled as fixed effect:", paste(names(varType)[i], collapse=', '), "\nThe results will not behave as expected and may be very wrong!!"))		
              }
              
              # If numeric/double is not fixed
              if( (showWarnings && ! dream) && varType[i] %in% c("numeric", "double") && (!names(varType)[i] %in% fixedVar) ){
                stop(paste("Continuous variable cannot be modeled as a random effect:", names(varType)[i]))		
              }
            }
            
            if( (! dream) && showWarnings && isVaryingCoefficientModel( fit ) ){
              stop(paste("Random slope models {i.e. ~ (var1 | var2) } are no longer supported\nfor estimating variance fractions.\nThey produced results that were not interpretable."))			
            }
            
            # show convergance message
            if( showWarnings && !is.null(fit@optinfo$conv$lme4$messages) ){
              stop(fit@optinfo$conv$lme4$messages)
            }
          }
)