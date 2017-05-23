#' plot_cvROC
#'
#' Plot cross-validated ROC curves and the averaged ROC curve. 
#' 
#' @param Y A \code{numeric} vector of class labels
#' @param X A \code{data.frame} of variables that was used in the call to \code{wrap_cvAUC}
#' @param wrap_cvAUC_fit An object of class \code{wrap_cvAUC}. It is assumed that \code{Y}
#' is in the same order as when \code{wrap_cvAUC} was called. 
#' @param plotAverage A \code{boolean} indicating whether to plot the average ROC curve
#' over cross-validation folds.
#' @param ... Other options passed to \code{plot}.
#' @importFrom stats predict
#' @importFrom ROCR prediction performance

plot_cvROC <- function(Y, X, wrap_cvAUC_fit, plotAverage, ...){
	# Work flow
	# ---------
	# In each fold: 
	# 	1. Make prediction object using ROCR::prediction
	# 	2. Get
	pred_valid <- mapply(ind = wrap_cvAUC_fit$folds, 
	                     fit = wrap_cvAUC_fit$fitLibrary, 
	                     FUN = function(ind, fit, X){ 
	                     	stats::predict(fit, newdata = X[ind,])
	                     }, MoreArgs = list(X=X))

	labels_valid <- lapply(wrap_cvAUC_fit$folds, function(i,Y){ Y[i] }, Y=Y)

	pred <- ROCR::prediction(pred_valid, labels_valid)
	perf <- performance( pred, "tpr", "fpr")
	plot(perf, avg="threshold", ...)
}