#' wrap_cvAUC
#' 
#' This function is a helper wrapper function for the \code{cvAUC} function
#' included in the \code{cvAUC} package by Erin LeDell. The function allows
#' the user to use the data splitting options available in the \code{SuperLearner}
#' package and provides a specific structure for different learners to be used 
#' to generate predictions. The biggest addition with this function is that influence
#' functions are returned, which can be used to develop hypothesis tests comparing
#' the CV-AUC between two different learners. The function \code{diff_cvAUC} 
#' performs these tests. 
#' 
#' @param Y A \code{numeric} vector of class labels
#' @param X A \code{data.frame} of variables that \code{learner} will use
#' to predict. It is assumed that the format at code{X} will place nicely with
#' the function specified by \code{learner}.
#' @param learner A \code{character} name of a function that generates predictions.
#' The function should take as input \code{Y}, \code{X}, and \code{newX}, use \code{X}
#' to predict \code{Y} and return predictions for \code{newX}. See examples below. 
#' @param confidence A \code{numeric} between 0 and 1 specifying the nominal coverage
#' probability for the confidence interval. Default is 0.95. 
#' @param seed A \code{numeric} specifying what seed to set prior to data splitting. 
#' If \code{diff_cvAUC} is to be used afterwards to compare CV-AUCs for different fits, 
#' be sure to specify the same \code{seed} so that the sample splits are the same.
#' @param id A \code{numeric} vector of observation identifiers. Only used for splitting
#' data and should probably be ignored for now as the CV-AUC calculations do not account for
#' dependent data in any other way. 
#' @param cvControl A \code{list} of a specific form. See \code{?SuperLearner.CV.control} for
#' more information. 
#' @param returnFits A \code{boolean} indicating whether or not to return the model fit objects
#' for each fold.
#' @param parallel A \code{boolean} indicating whether to perform the model fitting across folds in 
#' parallel. If \code{TRUE} then \code{foreach} is used to parallelize the fitting. 
#' @param ... Not currently used
#' 
#' @return An object of class \code{wrap_cvAUC} with the following entries: 
#' \item{cvAUC}{The estimated cross-validated AUC.}
#' \item{se}{The standard error for the estimated CV-AUC.}
#' \item{ci}{A \code{100*confidence} percent confidence interval.}
#' \item{confidence}{The level of confidence for the interval.}
#' \item{ic}{The estimated influence function evaluated on the observations.}
#' \item{folds}{The row indices for each validation sample.}
#' \item{fitLibrary}{The fit objects from \code{learner}.}
#' \item{learner}{The learner that was used to generate predictions.}
#' \item{p}{The one-sided p-value testing the null hypothesis that CV-AUC = 0.5 
#' against the alternative that CV-AUC > 0.5. }
#' 
#' @export
#' 
#' @importFrom SuperLearner CVFolds SuperLearner.CV.control
#' 
#' @examples
#' n <- 1000
#' X <- data.frame(x1=rnorm(n),x2=rnorm(n))
#' Y <- rbinom(n,1,plogis(X$x1 + X$x2))
#' myglm1 <- function(Y,X,newX){
#'    fm <- glm(Y~.,data=X,family=binomial())
#'    pred <- predict(fm,newdata=newX,type="response")
#'    return(list(fit = fm, pred = pred))
#' }
#' myglm2 <- function(Y,X,newX){
#'   fm <- glm(Y~x1,data=X,family=binomial())
#'   pred <- predict(fm,newdata=newX,type="response")
#'   return(list(fit = fm, pred = pred))
#' }
#' out1 <- wrap_cvAUC(Y = Y, X=X, learner = "myglm1")
#' out2 <- wrap_cvAUC(Y = Y, X=X, learner = "myglm2")

wrap_cvAUC <- function(
    Y, X, learner, confidence = 0.95, seed=1234, id=NULL, cvControl=list(V = 10L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL), 
    returnFits=FALSE, parallel = FALSE, ...
){
  N <- length(Y)
  cvControl <- do.call("SuperLearner.CV.control",cvControl)
  set.seed(seed)
  folds <- SuperLearner::CVFolds(N = N, id = id, Y = Y, cvControl = cvControl)

  predFitList <- .getPredictions(learner=learner, Y=Y, X=X, V=length(folds), folds=folds, returnFits=returnFits, parallel=parallel)

  pred <- predFitList$pred
  fitList <- predFitList$fit

  ciOut <- ci.cvAUC_withIC(prediction = pred, labels = Y, folds = folds, confidence = confidence)
  p <- pnorm(-abs((ciOut$cvAUC-0.5)/ciOut$se))
  out <- c(ciOut, folds = list(folds), fitLibrary = fitList, learner = learner, confidence = confidence, p = p)
  class(out) <- "wrap_cvAUC"
  return(out)
}

#' diff_cvAUC
#' 
#' This function estimates the difference between two cross-validated AUCs
#' fit using the \code{wrap_cvAUC} function. It also provides an influence 
#' function-based confidence interval estimate and hypothesis test of the null
#' hypothesis that the two CV-AUCs are equal. 
#' 
#' @param fit1 An object of class \code{wrap_AUC}
#' @param fit2 An object of class \code{wrap_AUC}
#' @param confidence A \code{numeric} between 0 and 1 indicating the nominal coverage
#' probability for the confidence interval. 
#' 
#' @return An object of class \code{diff_cvAUC} with the following entries:
#' \item{diff}{The difference in CV-AUC between the two fits.}
#' \item{ci}{The confidence interval for the difference between the two fits. }
#' \item{p}{The two-sided p-value for the test that the two CV-AUCs are equal. }
#' \item{folds}{The number of folds used by \cote{fit1} and \code{fit2}.}
#' \item{learner1}{The name of the learner used for \code{fit1}.}
#' \item{learner2}{The name of the learner used for \code{fit2}.}
#' \item{confidence}{The confidence interval level.}
#' 
#' @export

diff_cvAUC <- function(fit1, fit2, confidence = 0.95){
  # check if folds are identical
  sameFolds <- all(unlist(fit1$folds) == unlist(fit2$folds))
  if(!sameFolds) stop("Different folds used for fit1 and fit2. Be sure to set same seed.")

  diff <- fit1$cvAUC - fit2$cvAUC
  a <- matrix(c(1,-1))
  icMat <- rbind(fit1$ic, fit2$ic)
  seDiff <- sqrt( tcrossprod(crossprod(a,icMat))) / dim(icMat)[2]
  z <- qnorm(1- (1-confidence)/2)
  ci <- c(diff - z * seDiff, diff + z * seDiff)
  p <- 2*pnorm(-abs(diff/seDiff))
  out <- list(diff = diff, ci = ci, p = p, folds = length(fit1$folds), learner1 = fit1$learner, learner2 = fit2$learner, confidence = confidence)
  class(out) <- "diff_cvAUC"
  return(out)
}

#' print.wrap_cvAUC
#' 
#' A print method for objects of class \code{wrap_cvAUC}
#' 
#' @param object An object of class \code{wrap_cvAUC}.
#' @param digits A \code{numeric} providing the number of digits to round to.
#' 
#' @export 
print.wrap_cvAUC <- function(object, digits = 2){
  cat(paste0(length(object$folds),"-fold CV-AUC for ",
             object$learner, "      : ", 
             sprintf(paste0("%1.",digits,"f"),round(object$cvAUC,digits))))
  cat("\n")
  cat(paste0(object$confidence*100,"% confidence interval        : ", 
               paste0(sprintf(paste0("%1.",digits,"f"),round(object$ci, digits)),collapse=" , ")))
  cat("\n")
  cat(paste0("One-sided p-value CV-AUC > 0.5 : ", ifelse(object$p < 10^(-digits), paste0("< ",10^(-digits)),
                                                         paste0("= ", sprintf(paste0("%1.",digits,"f"),round(object$p, digits))))))
  cat("\n")
}


#' print.diff_cvAUC
#' 
#' A print method for objects of class \code{diff_cvAUC}
#' 
#' @param object An object of class \code{diff_cvAUC}.
#' @param digits A \code{numeric} providing the number of digits to round to.
#' 
#' @export 
print.diff_cvAUC <- function(object, digits = 2){
  cat(paste0("Diff ",object$folds,"-fold CV-AUC (", 
             object$learner1, " - ", object$learner2, ") : ", sprintf(paste0("%1.",digits,"f"),round(object$diff,digits))))
  cat("\n")
  cat(paste0(object$confidence*100,"% confidence interval               : ", 
               paste0(sprintf(paste0("%1.",digits,"f"),round(object$ci, digits)),collapse=" , ")))
  cat("\n")
  cat(paste0("Two-sided p-value diff not 0          : ", ifelse(object$p < 10^(-digits), paste0("< ",10^(-digits)),
                                                         paste0("= ", sprintf(paste0("%1.",digits,"f"),round(object$p, digits))))))
  cat("\n")
}

#' ci.cvAUC_withIC
#' 
#' This function is taken nearly verbatim from the \code{cvAUC} package. The only difference
#' is that it additionally returns estimated influence functions, which allows
#' us to do later comparisons of CV-AUC between different learners. 
#'  

ci.cvAUC_withIC <- function(predictions, labels, label.ordering = NULL, folds = NULL, confidence = 0.95) {
  
  # Pre-process the input
  clean <- cvAUC:::.process_input(predictions = predictions, labels = labels, 
                          label.ordering = label.ordering, folds = folds,
                          ids = NULL, confidence = confidence)
  
  predictions <- clean$predictions  # Length-V list of predicted values
  labels <- clean$labels  # Length-V list of true labels
  pos <- levels(labels[[1]])[[2]]  # Positive class label
  neg <- levels(labels[[1]])[[1]]  # Negative class label
  n_obs <- length(unlist(labels))  # Number of observations
  
  # Inverse probability weights across entire data set
  w1 <- 1/(sum(unlist(labels) == pos)/n_obs)  # Inverse weights for positive class
  w0 <- 1/(sum(unlist(labels) == neg)/n_obs)  # Inverse weights for negative class

  # This is required to cleanly get past R CMD CHECK
  # https://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
  # pred <- label <- NULL
  # fracNegLabelsWithSmallerPreds <- fracPosLabelsWithLargerPreds <- icVal <- NULL  

  .IC <- function(fold_preds, fold_labels, pos, neg, w1, w0) {
      n_rows <- length(fold_labels)
      n_pos <- sum(fold_labels == pos)
      n_neg <- n_rows - n_pos
      auc <- cvAUC:::AUC(fold_preds, fold_labels)
      DT <- data.table(pred = fold_preds, label = fold_labels)
      DT <- DT[order(pred, -xtfrm(label))]
      DT[, `:=`(fracNegLabelsWithSmallerPreds, cumsum(label == 
                                                        neg)/n_neg)]
      DT <- DT[order(-pred, label)]
      DT[, `:=`(fracPosLabelsWithLargerPreds, cumsum(label == 
                                                       pos)/n_pos)]
      DT[, `:=`(icVal, ifelse(label == pos, w1 * (fracNegLabelsWithSmallerPreds - 
                                                    auc), w0 * (fracPosLabelsWithLargerPreds - auc)))]
      return(DT$icVal)
  }

  icOut <- mapply(FUN = .IC, SIMPLIFY = FALSE, fold_preds = predictions, 
    fold_labels = labels, MoreArgs = list(pos = pos, neg = neg, w1 = w1, w0 = w0))
  ic <- rep(NA, n_obs)
  ic[unlist(folds)] <- unlist(icOut)
  # Estimated variance
  sighat2 <- mean(unlist(lapply(icOut, function(x){mean(x^2)})))
  se <- sqrt(sighat2/n_obs)  
  cvauc <- cvAUC::cvAUC(predictions, labels)$cvAUC
  z <- qnorm(confidence + (1 - confidence)/2)
  ci_cvauc <- c(cvauc - (z * se), cvauc + (z * se))
  ci_cvauc[1] <- ifelse(ci_cvauc[1] < 0, 0, ci_cvauc[1])  #Truncate CI at [0,1]
  ci_cvauc[2] <- ifelse(ci_cvauc[2] > 1, 1, ci_cvauc[2]) 
  
  return(list(cvAUC = cvauc, se = se, ci = ci_cvauc, confidence = confidence, ic = ic))
}

#' .getPredictions
#' 
#' This function calls the specified \code{learner} over each fold and
#' returns a list of cross-validated predictions and model fits. 
#' 

.getPredictions <- function(learner, Y, X, V, folds, returnFits, parallel){

  .doFit <- function(x, tmpX, Y, folds, learner){
    out <- do.call(learner, args=list(Y=Y[-folds[[x]]], X=tmpX[-folds[[x]],,drop=FALSE], newX=tmpX[folds[[x]],,drop=FALSE]))
    return(out)
  }

  if(parallel){
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    predFitList <- foreach(v = 1:length(folds), .export=learner) %dopar% 
      .doFit(v, tmpX = X, Y = Y, folds = folds, learner = learner)
    stopCluster(cl)
  }else{
    predFitList <- lapply(split(seq(V),factor(seq(V))),FUN=.doFit, tmpX = X, Y=Y, folds = folds, learner = learner)
  }
  
  # separate predictions from model fits
  predList <- lapply(predFitList, function(x){x$pred})
  if(returnFits){
    fitList <- lapply(predFitList, function(x){x$fit})
  }else{
    fitList <- NULL
  }
  # re-order predictions
  tmp <- unlist(predList)
  pred <- rep(NA, length(tmp))
  pred[unlist(folds)] <- tmp
  
  # return results
  return(list(pred=pred,fit=fitList))
}

