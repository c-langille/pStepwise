#' Easily access p-values from linear models
#'
#' This function gives the user easy access to p-values for specific predictors from a linear model.  It is mainly used
#' to easily pass p-values to other functions.
#'
#' @param pred the predictor from the linear model whose p-value is desired
#' @param fit a linear model of type "lm" containing the desired predictor
#' @return the p-value for the desired predictor from the linear model
#' @author Cory Langille <lang1729@gmail.com>
#' @examples
#' #Using the leafshape dataset from the DAAG package
#' data(leafshape)
#' attach(leafshape)
#' currentModel <- lm(bladelen~., data = leafshape)
#' extractp("bladewid", currentModel)
#' #    bladewid
#' #    2.579703e-43
#' @seealso \code{\link{lm}}, \code{\link{pStepwise}}
#' @import stats
#' @export
#'
extractp <- function(pred, fit) {
  predClasses <- attr(fit$terms, "dataClasses")        #finds the predictors that are categorical variables
  fact <- predClasses[names(predClasses) == pred]
  if(fact == "factor") {
    pvalues <- anova(fit)[,5]
    predList <- rownames(anova(fit))
    index <- which(predList == pred)
    return(pvalues[index])
  }else {
    pvalues <- summary(fit)$coefficients[,4]
    return(pvalues[names(pvalues)==pred])
  }
}



#' Create formulas from strings to be used in linear models
#'
#' This function takes a string and a linear model and either adds the predictor to the model, or removes it.
#'
#' @param pred the predictor to be added or removed from the current model
#' @param fitCurrent the current model to be updated
#' @param add by default adds the predictor to to the model.  add=F removes the predictor from the model
#' @seealso \code{\link{pStepwise}}
#' @return the updated model of type "lm"
#' @examples
#' #Using the leafshape data set from the DAAG package
#'
#' #Adding a predictor to a linear model
#' data(leafshape)
#' attach(leafshape)
#' currentModel <- lm(bladelen~bladewid + latitude, data = leafshape)
#' newModel <- fMaker("petiole", currentModel)
#' newModel$call
#' #lm(formula = bladelen ~ bladewid + latitude + petiole, data = leafshape)
#'
#' #Removing a predictor from a linear model
#' currentModel <- lm(bladelen~bladewid + latitude, data = leafshape)
#' newModel <- fMaker("latitude", currentModel, add = FALSE)
#' newModel$call
#' #lm(formula = bladelen ~ bladewid, data = leafshape)
#' @export
#'
fMaker <- function(pred, fitCurrent, add=T) {
  addNew <- as.formula(paste(".~.+", pred))
  remNew <- as.formula(paste(".~.-", pred))
  ifelse(add, return(update(fitCurrent, addNew)), return(update(fitCurrent, remNew)))
}



#' Adds a single predictor to a linear model based on its p-value
#'
#' This function will try and add a new predictor to a current model.  A predictor will be added if it has minimum p-value among all predictors and its p-value is below a certain threshold
#'
#' @param fitCurrent the current model of type "lm"
#' @param fullmodel a linear model containing all possible predictors.  Typically of the form lm(y~., data=data)
#' @param aEnter the threshold for adding the predictor, set to 0.1 be default
#' @param forcedOut vector of predictors that will be forced out of the final model
#' @seealso \code{\link{pStepwise}}
#' @return an updated linear model of type "lm"
#' @export
#'
stepfwd <- function(fitCurrent, fullmodel, aEnter = 0.1, forcedOut = NULL) {
  predsInModel <- rownames(anova(fitCurrent))                                                   #list of predictors in current model
  predsFull <- rownames(anova(fullmodel))                                                       #list of predictors in full model
  predsNotInModel <- setdiff(predsFull, predsInModel)                                           #list of predictors not in current model, not including forced out predictors
  pvalues <- unlist(sapply(predsNotInModel, function(x) as.numeric(extractp(x, fMaker(x, fitCurrent)))))#takes each predictor not in the current model, creates a new lm which includes it, and stores its respective p-value. Really ugly, but the only way I could make it work.
  toAdd <- pvalues[which(pvalues==min(pvalues))]                                                #possible new predictor
  if(length(toAdd)==0) return(fitCurrent)                                                       #returns original model if no new predictors can be added
  if(toAdd <= aEnter) {
    cat("Add predictor", names(toAdd), "\n")
    print(summary(fMaker(names(toAdd), fitCurrent)))
    return(fMaker(names(toAdd), fitCurrent))                                                    #updates and returns new model with additional predictor
  }
  return(fitCurrent)

}



#' Removes a single predictor from a linear model based on its p-value
#'
#' This function will try and remove a single predictor from a current linear model.  A predictor will be removed if it has maximal p-value and its p-value is greater than a certain threshold.
#'
#' @param fitCurrent the current model of type "lm"
#' @param fullmodel a linear model containing all possible predictors.  Typically of the form lm(y~., data=data)
#' @param forcedIn vector of predictors that will be forced into the final model
#' @param aRemove the threshold for removing the predictor, set to 0.1 by default
#' @return an updated linear model of type "lm"
#' @seealso \code{\link{pStepwise}}
#' @export
#'
stepbwd <- function(fitCurrent, fullmodel, aRemove = 0.1, forcedIn = NULL) {
  predsIncluded <- rownames(anova(fitCurrent))                                               #predictors in current model
  predsIncluded <- predsIncluded[(predsIncluded != "Residuals")]                             #removes "residuals" from predictors
  predsIncluded <- setdiff(predsIncluded, intersect(predsIncluded, forcedIn))                #makes sure no forced in predictors get removed
  pvalues <- unlist(sapply(predsIncluded, function(x) as.numeric(extractp(x, fitCurrent))))  #checks the p-value for each predictor in current model
  if(length(pvalues)==0) return(fitCurrent)                                                  #returns current model if there are no more possible predictors to remove
  toRemove <- pvalues[which(pvalues == max(pvalues))]                                        #selects the predictor with maximal p-value
  if(length(toRemove)==0) return(fitCurrent)
  if(toRemove > aRemove){
    cat("Remove predictor", names(toRemove), "\n")
    return(fMaker(names(toRemove), fitCurrent, add=F))                                       #returns an updated model if the p-value is above the threshold
  }
    return(fitCurrent)                                                                       #else, returns original model
}

#' Selects the best predictors for a linear model based on p-values
#'
#' This function will attempt to create the "best' linear model by finding the most significant predictors. A predictor will be included/excluded in the final model if when it is added/removed its p-value is below/above a certain threshold.
#' @usage pStepwise(response, fullmodel, aEnter = 0.1,
#'                  aRemove = 0.1, forcedIn = NULL, forcedOut = NULL, method = "both")
#' @param response the response variable of interest in the model
#' @param fullmodel a linear model containing all possible predictors, typically of the form lm(y ~ ., data = data)
#' @param aEnter the threshold for adding new predictors, set to 0.1 by default
#' @param aRemove the threshold for removing predictors from the current model, set to 0.1 by default
#' @param forcedIn a vector of predictors that will be forced into the final model regardless of their p-values
#' @param forcedOut a vector of predictors that will not be included in the final model regardless of their p-values
#' @param method "forward" will only add predictors, "backward" will only remove predictors and "both" will perfrom stepwise.  "both" by default
#' @return a linear model of type lm containing the "best" predictors
#' @author Cory Langille <lang1729@gmail.com>
#' @seealso \code{\link{extractp}}, \code{\link{stepfwd}}, \code{\link{stepbwd}}, \code{\link{fMaker}}
#' @examples
#' #Using the leafshape dataset from the DAAG package
#' data(leafshape)
#' attach(leafshape)
#' response <- "bladelen"
#' fullmodel <- lm(bladelen ~ . , data = leafshape)
#' forcedOut <- c("loglen", "logwid", "logpet" )
#' pStepwise(response, fullmodel, forcedOut = forcedOut)
#' @export
#'
pStepwise <- function(response, fullmodel, aEnter = 0.1, aRemove = 0.1,
                      forcedIn = NULL, forcedOut = NULL, method = "both") {
  fitBwd <- lm(as.formula(paste(response, "~1")))                                  #creates an empty model to begin with (poor name choice but it makes the while loop below easier)
  for(pred in forcedIn) fitBwd <- fMaker(pred, fitBwd)                             #adds in forced predictors to initial model
  for(pred in forcedOut) fullmodel <- fMaker(pred, fullmodel, add=F)               #removes forced out predictors from fullmodel
  if(method=="forward") aRemove <- 1                                               #makes it impossible to remove predictors if using "forward" method
  if(method=="backward") {                                                         #section for backward selection
    fitFwd <- fullmodel                                                                         #starts with a full model as the inital model
    cat("Initial model: ")
    print(summary(fullmodel))
    while(TRUE) {
      fitBwd <- stepbwd(fitFwd, fullmodel, aRemove = aRemove, forcedIn = forcedIn)              #continously tries removing predictors until no longer possible
      if(identical(fitFwd, fitBwd)==T) {
        cat("Predictors forced in: ", forcedIn, "\n")
        cat("Predictors forced out: ", forcedOut, "\n")
        cat("--------Final model--------", "\n")
        print(summary(fitFwd))
        return(fitFwd)
      } else {
        fitFwd <- fitBwd
      }
    }
  }
  cat("Initial model: ")
  print(summary(fitBwd))
  while(TRUE){
    fitFwd <- stepfwd(fitBwd, fullmodel, aEnter = aEnter, forcedOut = forcedOut)   #function that tries to add a predictor to current model
    if(identical(fitFwd, fitBwd) == T) {                                           #if new model is the same as old model, then no new predictors were added, it will stop
      cat("Predictors forced in: ", forcedIn, "\n")
      cat("Predictors forced out: ", forcedOut, "\n")
      cat("--------Final model--------", "\n")
      print(summary(fitFwd))
      return(fitFwd)
    }else {                                                                        #function that tries to remove a predictor from current model
      fitBwd <- stepbwd(fitFwd, fullmodel, forcedIn = forcedIn, aRemove = aRemove)
    }
  }
}
