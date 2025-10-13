#' Print method for class "varshrinkest"
#'
#' Print method for an object of class "varshrinkest"
#' @param x An object of class "varshrinkest"
#' @param digits,... Other arguments for print() method
#' @importFrom stats coef
#' @examples
#' data(Canada, package = "vars")
#' y <- diff(Canada)
#' estim <- VARshrink(y, p = 2, type = "const", method = "ridge")
#' print(estim)
#' @export
print.varshrinkest <- function(x, digits = max(3, getOption("digits") - 3),
                               ...) {
  K <- length(x$varresult)
  tsnames <- names(x$varresult)
  text1 <- "VAR Shrinkage Estimation Results:"
  cat(paste("\n", text1, "\n", sep = ""))
  row <- paste(rep("=", nchar(text1)), collapse = "")
  cat(row, "\n")
  cat("\n")
	## print coefficients ##
  for (i in 1:K) {
      result <- coef(x$varresult[[i]])
      text1 <- paste("Estimated coefficients for equation ",
                     tsnames[i], ":", sep = "")
      cat(text1, "\n")
      row <- paste(rep("=", nchar(text1)), collapse = "")
      cat(row, "\n")
      text2 <- paste("Call:\n", tsnames[i], " = ",
                     paste(names(result), collapse = " + "), sep = "")
      cat(text2, "\n\n")
      print(result, ...)
      cat("\n\n")
  }
	## print other parameters ##
	if (!is.null(x$Sigma)) {
		cat("Sigma for noise:\n")
		print(x$Sigma, ...)
	}
	if (!is.null(x$dof)) {
		text1 <- paste("dof for noise: ", x$dof,
		               " (estimated: ",
		               ifelse(is.null(e <- x$dof.estimated), FALSE, e), ")",
		               sep = "")
		cat(text1, "\n")
	}
	if (!is.null(x$lambda)) {
		text1 <- paste("lambda: ", paste(x$lambda, collapse = " "),
		               " (estimated: ",
		               ifelse(is.null(e <- x$lambda.estimated), FALSE, e),
		               ")",
		               sep = "")
		cat(text1, "\n")
	}
	if (!is.null(x$lambda_var)) {
		text1 <- paste("lambda_var: ", paste(x$lambda_var, collapse = " "),
		               " (estimated: ",
		               ifelse(is.null(e <- x$lambda_var.estimated), FALSE, e),
		               ")",
		               sep = "")
		cat(text1, "\n")
	}
	if (!is.null(x$GCV)) {
		text1 <- paste("GCV: ", paste(round(x$GCV, digits), collapse = " "))
		cat(text1, "\n")
	}
	## return value
	invisible(x)
}
