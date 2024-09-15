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

#Function for viewing patterns of observations in a frequency table (written by KB)
view_patterns <- function(data, id.column, time.column, up.through) {
  #Filter data to only include observations up through the specified month
  data_filter <- data %>% filter({{time.column}} <= up.through)  
  
  #Making a new dataframe to store patterns for each ID
  pattern_set <- data.frame(
    id <- rep(NA, max(data_filter[[deparse(substitute(id.column))]])),
    pattern <- rep(NA, max(data_filter[[deparse(substitute(id.column))]]))
  )
  
  #Getting a list of ids
  ids <- unique(data_filter[[deparse(substitute(id.column))]])
  
  #Loop over the ids to get the observation pattern for each
  for (i in ids) {
    pattern_set$id[i] <- i
    id_set <- dplyr::filter(data_filter, {{id.column}} == i)
    pattern_set$pattern[i] <- paste(id_set[[deparse(substitute(time.column))]], collapse = ", ")
  }
  
  #Cleaning up the pattern variable
  pattern_set$pattern <- gsub(" ", "", pattern_set$pattern)
  
  #Selecting only the relevent columns from pattern_set
  pattern_set <- pattern_set %>% dplyr::select(id, pattern) %>% filter(!is.na(pattern) & pattern != "")
  
  #Tabulating the frequencies of each pattern
  pattern_counts <- pattern_set %>% group_by(pattern) %>% count() %>% mutate(percent = (n/length(ids))*100) %>% arrange(desc(n))
  
  #Returning the pattern counts table
  return(pattern_counts)
}


#Updated version of function for correlation matrices that takes "markov" as an input correlation structure - adapted original version of function with help from Claude
createCorrMatrixR <- function(fitted.model = NULL, 
                              model.corstr = NULL, 
                              fitting.package = NULL,
                              time_var = NULL,
                              use_time_labels = TRUE) {
  # check inputs 
  if (is.null(fitted.model) || is.null(model.corstr) || is.null(fitting.package)) {
    stop("Some arguments to createCorrelationMatrixR() are missing.")
  }
  fitting.package <- tolower(fitting.package) 
  if (!fitting.package %in% c("geepack", "geem", "gee", "qlspack")) {
    stop("Unrecognized fitting.package in createCorrelationMatrixR().")
  }
  model.corstr <- tolower(model.corstr)
  if (!model.corstr %in% c("independent", "exchangeable", "ar1", "unstructured", "markov")) {
    stop("Unsupported model.corstr in createCorrelationMatrixR().")
  }
  
  # create working correlation matrix
  model_sumy <- summary(fitted.model)
  
  if (fitting.package %in% c("geepack", "qlspack")) {
    est_alpha <- model_sumy$corr$Estimate
    max_cluster_size <- max(table(fitted.model$id))
    
    if (model.corstr == "independent") {
      working_corrmat <- diag(1, nrow = max_cluster_size, ncol = max_cluster_size)
    } else if (model.corstr == "exchangeable") {
      working_corrmat <- matrix(est_alpha, nrow = max_cluster_size, ncol = max_cluster_size)
      diag(working_corrmat) <- 1
    } else if (model.corstr == "ar1") {
      working_corrmat <- matrix(NA_real_, nrow = max_cluster_size, ncol = max_cluster_size)
      for (i in seq(ncol(working_corrmat))) {
        for (j in seq(nrow(working_corrmat))) {
          working_corrmat[i, j] <- est_alpha^abs(i - j)
        }
      }
    } else if (model.corstr == "unstructured") {
      names(est_alpha) <- rownames(model_sumy$corr)
      working_corrmat <- matrix(NA_real_, nrow = max_cluster_size, ncol = max_cluster_size)
      for (i in seq(est_alpha)) {
        alpha_dims <- gsub("alpha.", "", names(est_alpha)[i])
        alpha_dims <- stringr::str_split(alpha_dims, pattern = ":", simplify = TRUE)
        alpha_dims <- as.numeric(alpha_dims)
        alpha_dims_rev <- rev(alpha_dims)
        working_corrmat[alpha_dims[1], alpha_dims[2]] <- as.numeric(est_alpha[i])
        working_corrmat[alpha_dims_rev[1], alpha_dims_rev[2]] <- as.numeric(est_alpha[i])
      }
      diag(working_corrmat) <- 1
    } else if (model.corstr == "markov") {
      # Extract time variable
      if (is.null(time_var)) {
        stop("Time variable name must be provided for Markov correlation structure.")
      }
      
      if (!time_var %in% names(fitted.model$data)) {
        stop(paste("Time variable", time_var, "not found in the model data."))
      }
      
      time_points <- unique(fitted.model$data[[time_var]])
      time_points <- sort(time_points)
      
      # Create time difference matrix
      time_diff_mat <- abs(outer(time_points, time_points, "-"))
      
      # Calculate correlation matrix
      working_corrmat <- est_alpha^time_diff_mat
      
      if (use_time_labels) {
        # Use actual time values (e.g., months) as labels
        labels <- as.character(time_points)
      } else {
        # Use t1, t2, t3, etc. as labels
        labels <- paste0("t", seq_along(time_points))
      }
      
      rownames(working_corrmat) <- colnames(working_corrmat) <- labels
    }
  } else if (fitting.package == "geem") {
    working_corrmat <- as.matrix(model_sumy$biggest.R.alpha)
  } else if (fitting.package == "gee") {
    working_corrmat <- model_sumy$working.correlation
  }
  
  return(working_corrmat)
}


#Function for subsetting an estimated correlation matrix based on timepoints of interest (takes output for createCorrMatrix) - written with help from Claude
  #example subsetCorrMatrix(full_matrix, c("1", "6", "12"))

subsetCorrMatrix <- function(corr_matrix, time_points) {
  # Check if all specified time points are in the matrix
  if (!all(time_points %in% rownames(corr_matrix))) {
    missing_points <- setdiff(time_points, rownames(corr_matrix))
    stop(paste("The following time points are not in the correlation matrix:", 
               paste(missing_points, collapse = ", ")))
  }
  
  # Subset the correlation matrix
  subset_matrix <- corr_matrix[time_points, time_points]
  
  # If only one time point is selected, ensure it's still a matrix
  if (length(time_points) == 1) {
    subset_matrix <- matrix(subset_matrix, nrow = 1, ncol = 1)
    rownames(subset_matrix) <- colnames(subset_matrix) <- time_points
  }
  
  return(subset_matrix)
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