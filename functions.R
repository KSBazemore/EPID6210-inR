#Functions for longitudinal and clustered data analysis
#Author: Katrina Bazemore, MPH
#Have linked to sources of functions where relevant

#Loading libraries
library(gee)
library(geepack)
library(geeM)

#Function to get p-values in coefficients table from gee package models (written by KB)
#Example P_table(gee_ar1, 4)

P_table <- function(fitted.model, round.to) {
  coeffs <- round(summary(fitted.model)$coefficients, round.to)
  model_based_p <- round((2*pnorm(abs(coef(summary(fitted.model))[,3]), lower.tail = F)), round.to)
  sandwich_based_p <- round((2*pnorm(abs(coef(summary(fitted.model))[,5]), lower.tail = F)), round.to)
  table <- cbind(coeffs, model_based_p, sandwich_based_p)
  return(table)
}



#Function for working correlation matrix from: https://stackoverflow.com/questions/78003795/extracting-working-correlation-matrix-from-a-geeglm-object-produced-using-geepac

#Use: specify model, correclation structure, and package
#Accepted correlation structures: "independent", "exchangeable", "ar1", "unstructured"
#Accepted packages: "geepack", "geem", "gee"

#Example: createCorrelationMatrixR(model.1, model.corstr = "ar1", fitting.package = "geepack")

createCorrelationMatrixR <- function(fitted.model = NULL, 
                                     model.corstr = NULL, 
                                     fitting.package = NULL) {
  # check inputs 
  if (is.null(fitted.model) || is.null(model.corstr) || is.null(fitting.package)) {
    stop("Some arguments to createCorrelationMatrixR() are missing.")
  }
  fitting.package <- tolower(fitting.package) 
  if (!fitting.package %in% c("geepack", "geem", "gee")) {
    stop("Unrecognized fitting.package in createCorrelationMatrixR().")
  }
  model.corstr <- tolower(model.corstr)
  if (!model.corstr %in% c("independent", "exchangeable", "ar1", "unstructured")) {
    stop("Unsupported model.corstr in createCorrelationMatrixR().")
  }
  # create working correlation matrix
  model_sumy <- summary(fitted.model)
  if (fitting.package == "geepack") {
    est_alpha <- model_sumy$corr$Estimate
    if (model.corstr == "independent") {
      working_corrmat <- diag(1, 
                              nrow = max(table(fitted.model$id)), 
                              ncol = max(table(fitted.model$id)))
    } else if (model.corstr == "exchangeable") {
      working_corrmat <- matrix(est_alpha, 
                                nrow = max(table(fitted.model$id)), 
                                ncol = max(table(fitted.model$id)))
      diag(working_corrmat) <- 1
    } else if (model.corstr == "ar1") {
      working_corrmat <- matrix(NA_real_, 
                                nrow = max(table(fitted.model$id)), 
                                ncol = max(table(fitted.model$id)))
      for (i in seq(ncol(working_corrmat))) {
        for (j in seq(nrow(working_corrmat))) {
          working_corrmat[i, j] <- est_alpha^abs(i - j)
        }
      }
    } else if (model.corstr == "unstructured") {
      names(est_alpha) <- rownames(model_sumy$corr)
      working_corrmat <- matrix(NA_real_, 
                                nrow = max(table(fitted.model$id)), 
                                ncol = max(table(fitted.model$id)))
      for (i in seq(est_alpha)) {
        alpha_dims <- gsub("alpha.", "", names(est_alpha)[i])
        alpha_dims <- stringr::str_split(alpha_dims, 
                                         pattern = ":", 
                                         simplify = TRUE)
        alpha_dims <- as.numeric(alpha_dims)
        alpha_dims_rev <- rev(alpha_dims)
        working_corrmat[alpha_dims[1], alpha_dims[2]] <- as.numeric(est_alpha[i])
        working_corrmat[alpha_dims_rev[1], alpha_dims_rev[2]] <- as.numeric(est_alpha[i])
      }
      diag(working_corrmat) <- 1
    }
  } else if (fitting.package == "geem") {
    working_corrmat <- as.matrix(model_sumy$biggest.R.alpha)
  } else if (fitting.package == "gee") {
    working_corrmat <- model_sumy$working.correlation
  }
  rownames(working_corrmat) <- colnames(working_corrmat) <- paste0("t", seq(ncol(working_corrmat)))
  return(working_corrmat)
}


#Function for calculating QIC and CIC on GEE package models from: https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture23.htm

QIC.binom.gee <- function(model.R,model.independence)
{
  #calculates binomial QAIC of Pan (2001)
  #obtain trace term of QAIC
  AIinverse <- solve(model.independence$naive.variance)
  V.msR <- model.R$robust.variance
  trace.term <- sum(diag(AIinverse%*%V.msR))
  #estimated mean and observed values
  mu.R <- model.R$fitted.values
  y <- model.R$y
  #scale for binary data
  scale <- 1
  #quasilikelihood for binomial model
  quasi.R <- sum(y*log(mu.R/(1-mu.R))+log(1-mu.R))/scale
  QIC <- (-2)*quasi.R + 2*trace.term
  output <- c(QIC,trace.term)
  names(output) <- c('QIC','CIC')
  output
}


#Function for calculating QIC and CIC on geepack package models from: https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture23.htm
QIC.binom.geeglm <- function(model.geeglm, model.independence)
{
  #calculates binomial QAIC of Pan (2001)
  #obtain trace term of QAIC
  AIinverse <- solve(model.independence$geese$vbeta.naiv)
  V.msR <- model.geeglm$geese$vbeta
  trace.term <- sum(diag(AIinverse%*%V.msR))
  #estimated mean and observed values
  mu.R <- model.geeglm$fitted.values
  y <- model.geeglm$y
  #scale for binary data
  scale <- 1
  #quasilikelihood for binomial model
  quasi.R <- sum(y*log(mu.R/(1-mu.R))+log(1-mu.R))/scale
  QIC <- (-2)*quasi.R + 2*trace.term
  output <- c(QIC,trace.term)
  names(output) <- c('QIC','CIC')
  output
}