# nFGS (a covariate selection algorithm in non-Gaussianity models in r)

## 1. get best single covariate, augment that to best two covariates
## 2. get best two covariates, augment that to best triple covariates and so on
## [1 check p value (significance); 2 check normality of rx; 3 check independence of rx and ry]

library(nortest)
library(dHSIC)

## x, y, W need to have same number of observations
nGFS <-
  function(x, # the focal predictor
           y, # the outcome
           W, # a set of potential covariates
           fixed = NULL, # fixed effects
           is.factor = TRUE, # if fixed effects indicator is a factor
           alpha_Gauss = 0.05, # normality check
           alpha_indep = 0.05, # independence check
           alpha_p = 0.05, # significance check
           independence = "hsic", # independence test, default "hsic", alternative "corr"
           predefine = NULL) # a set of predefined confounders
  {
    # Data preparation
    x <- as.matrix(x)
    y <- as.matrix(y)
    W <- as.matrix(W)
    nobs <- NROW(x)
    
    ## check dimensions of the variables
    if (NCOL(y) != 1)
      stop("forward_selection works for univariate linear model")
    if (NCOL(x) != 1)
      stop("forward_selection works for one interested predictor")
    
    fix_name <- NULL
    ## additional terms: fixed effects and predefined confouders if possible
    if (is.null(fixed) & is.null(predefine)) {
      add_term <- NULL
    } else {
      if (is.factor & is.null(fixed) == FALSE) {
        add_term <- cbind(predefine, factor(fixed))
        fix_name <- NULL
        for (factor_n in 1:(length(unique(fixed)) - 1)) {
          factor1 <- paste('fixed_factor', factor_n, sep = '')
          fix_name <- cbind(fix_name, factor1)
        }
      } else {
        if (is.null(fixed)) {
          add_term <- predefine
        } else {
          add_term <- cbind(predefine, fixed)
          fix_name <- "fixed"
        }
      }
      if (is.list(add_term)) {
        add_term <- matrix(unlist(add_term), ncol = length(add_term))
      }
    }
    
    ## work with variable names
    if (is.null(colnames(x))) {
      zname <- c('intercept', 'x')
    } else {
      xname <- colnames(x) # with x name
      zname <- c('intercept', colnames(x))
    }
    
    pre_name <- NULL
    for (pre_n in 1:NCOL(predefine)) {
      pre1 <- paste('predifined', pre_n, sep = '')
      pre_name <- c(pre_name, pre1)
    }
    
    if (is.null(predefine)) {
      add_name <- c(fix_name)
    } else {
      if (is.null(fixed)) {
        add_name <- c(pre_name)
      } else {
        add_name <- c(pre_name, fix_name)
      }
    }
    
    # linear regression: unconditional model
    if (is.null(add_term)) {
      m0 <- lm(y ~ x)
      uncond <- summary(m0)$coefficients
      rownames(uncond) <- zname
    } else {
      m0 <- lm(y ~ x + add_term)
      uncond <- summary(m0)$coefficients
      rownames(uncond) <- c(zname, add_name)
    }
    
    ## print unconditional model
    cat("Unconditional model (no covariates)", "\n")
    print(uncond)
    cat("\n")
    cat("R-squared and Adjusted R-squared:", "\n")
    print(summary(m0)$r.squared)
    print(summary(m0)$adj.r.squared)
    cat("\n")
    
    # find best set of z from W
    best_z <- matrix(NA, nrow = nobs, ncol = 0) # empty z
    ind <- list() # empty index
    
    # check prefect linearity with x
    ind_linear <- list()
    for (col in 1:NCOL(W)) {
      if (cor.test(x, W[, col])$esti > 0.99) {
        ind_linear <- cbind(ind_linear, col)
      }
    }
    
    for (newCol in 1:NCOL(W)) {
      pval <- double(NCOL(W)) # save p-value of independence test
      for (i in 1:NCOL(W)) {
        if (i %in% ind == TRUE | i %in% ind_linear == TRUE) {
          next
        }
        
        if (is.null(add_term)) {
          # without additional terms
          new_W <- cbind(best_z, W[, i])
          new_x <- cbind(x, new_W)
        } else {
          # with additional terms
          new_W <- cbind(add_term, best_z, W[, i])
          new_x <- cbind(x, add_term, new_W)
        }
        
        rx <- lm.fit(new_W, x)$residuals
        ry <- lm.fit(new_x, y)$residuals
        
        # drop selected z if non-significance
        m01 <- lm(y ~ new_x)
        p <- summary(m01)$coefficients[, 4]
        if (tail(p, n = 1) > alpha_p) {
          next
        }
        
        p_Gauss <-
          lillie.test(rx)$p #the Lilliefors (Kolmogorov-Smirnov) normality test
        if (p_Gauss > alpha_Gauss) {
          next
        }
        
        # Independence test: two ways (HSIC or Nonlinear corr)
        ## HSIC
        if (independence == "hsic") {
          hsic.p <-
            dhsic.test(ry, rx, method = "gamma", kernel = "gaussian")$p.value
          pval[i] <- hsic.p
        }
        
        ## Nonlinear corr
        if (independence == "corr") {
          pval_corr <- c(
            cor.test(tanh(rx), ry)$p.value,
            cor.test(rx, tanh(ry))$p.value,
            cor.test(rx ^ 2, ry)$p.value,
            cor.test(rx, ry ^ 2)$p.value,
            cor.test(tanh(rx), tanh(ry))$p.value,
            cor.test(tanh(rx), ry ^ 2)$p.value,
            cor.test(rx ^ 2, tanh(ry))$p.value,
            cor.test(rx ^ 2, ry ^ 2)$p.value
          )
          pval[i] <- min(pval_corr)
        }
      }
      
      # find best pval (max) from independence test
      if (all(pval == 0)) {
        break
      }
      
      if (max(pval) > alpha_indep) {
        index <- which.max(pval)
        #save index
        ind <- cbind(ind, index)
        best_z = cbind(best_z, W[, index])
      }
      else {
        break
      }
    }
    
    if (NCOL(best_z) > 0) {
      if (is.null(add_term)) {
        m1 <- lm(y ~ x + best_z)
      } else {
        m1 <- lm(y ~ x + best_z + add_term)
      }
      
      cat("Best z index order: ", "\n")
      print (as.numeric(ind)) # print index from the best single z, the best two...
      cat("\n")
      
      output <- summary(m1)$coefficients
      
      ## variable names for z
      for (index in ind) {
        if (is.null(colnames(W))) {
          z <- paste('z', index, sep = '')
        } else {
          z <- colnames(W)[index]
        }
        zname <- c(zname, z)
      }
      
      if (is.null(add_term)) {
        rownames(output) <- zname
      } else {
        rownames(output) <- c(zname, add_name)
      }
      
      cat("Conditional model (with selected z)", "\n")
      print(output)
      cat("\n")
      cat("R-squared and Adjusted R-squared:", "\n")
      print(summary(m1)$r.squared)
      print(summary(m1)$adj.r.squared)
      
      # return index and regression coefficients table
      invisible(list(ind, output))
      
    } else {
      cat("No z (covariates) found")
      cat("\n")
    }
  }