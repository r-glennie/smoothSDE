#' R6 class for stochastic differential equation with state-switching 
#' 
#' Contains the model formulas and data.
#' 
#' @importFrom R6 R6Class
#' @importFrom mgcv gam rmvn
#' @importFrom ggplot2 ggplot aes theme_light geom_line theme scale_colour_manual
#' facet_wrap label_bquote xlab ylab ggtitle element_blank element_text
#' @importFrom TMB MakeADFun sdreport
#' @importFrom MASS ginv
#' 
#' @useDynLib smoothSDE, .registration = TRUE
#' 
#' @export
#' 
SDE_HMM <- R6Class("SDE_HMM", inherit = SDE, 
    public = list(
        #################
        ## Constructor ##
        #################
        #' @description Create a SDE_HMM object
        #' 
        #' @param formulas List of formulas for model parameters
        #' @param data Data frame with covariates, response variable,
        #' time, coarse_time, and ID
        #' @param type Type of SDE ("BM", "OU"...)
        #' @param response Vector of names of response variables
        #' @param number of coarse scale states 
        #' @param par0 Vector of initial values for SDE parameters
        #' @param fixpar Vector of names of fixed SDE parameters
        #' @param other_data Named list of data objects to pass to
        #' likelihood
        #' 
        #' @return A new SDE object
        initialize = function(formulas, data, type, response, nstates = 2, par0 = NULL, 
                              tpm0 = NULL, 
                              fixpar = NULL, other_data = NULL) {
            private$formulas_ <- formulas
            private$type_ <- type
            private$response_ <- response
            private$fixpar_ <- fixpar
            
            if(any(!response %in% colnames(data)))
                stop("'response' not found in 'data'")
            
            # Link functions for SDE parameters
            n_dim <- length(response)
            link <- switch (type,
                            "BM_HMM" = list(mu = identity, sigma = log),
                            "BM-t_HMM" = list(mu = identity, sigma = log),
                            "OU_HMM" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                             tau = log, kappa = log)),
                            "CTCRW_HMM" = list(beta = log, sigma = log))
            
            # Inverse link functions for SDE parameters
            invlink <- switch (type,
                               "BM_HMM" = list(mu = identity, sigma = exp),
                               "BM-t_HMM" = list(mu = identity, sigma = exp),
                               "OU_HMM" = as.list(c(mu = lapply(1:n_dim, function(i) identity), 
                                                tau = exp, kappa = exp)),
                               "CTCRW_HMM" = list(beta = exp, sigma = exp))
            
            private$link_ <- link
            private$invlink_ <- invlink
            
            # Check that "formulas" is of the right length and with right names
            if(length(formulas) != length(invlink)) {
                err <- paste0("'formulas' should be a list of length ", 
                              length(invlink), " for the model ", type,
                              ", with components ", 
                              paste(names(invlink), collapse = ", "))
                stop(err)
            } else if(any(names(formulas) != names(invlink))) {
                err <- paste0("'formulas' should be a list with components ", 
                              paste(names(invlink), collapse = ", "))
                stop(err)
            }
            if(any(sapply(self$formulas()[fixpar], function(f) f != ~1))) {
                stop("formulas should be ~1 for fixed parameters")
            }
            
            # Check that data has an "ID" column, and that it's a factor
            if(!any(colnames(data) == "ID")) {
                warning(paste("No ID column found in 'data',",
                              "assuming same ID for all observations"))
                data$ID <- factor(1)
            } else {
                data$ID <- factor(data$ID)
            }
            
            # Check that data has a "time" column
            if(!any(colnames(data) == "time")) {
                stop("'data' should have a time column")
            }
            private$data_ <- data
            
            # Save terms of model formulas and model matrices
            mats <- self$make_mat()
            # Convert design matrices to multi-state design matrices 
            X_fe_list <- lapply(1:nstates, function(i) mats$X_fe)
            X_re_list <- lapply(1:nstates, function(i) mats$X_re)
            S_list <- lapply(1:nstates, function(i) mats$S)
            if (!is.null(mats$X_fe)) {
                X_fe <- bdiag(X_fe_list)
            } else {
                X_fe <- NULL
            }
            if(!is.null(mats$X_re)) {
                X_re <- bdiag(X_re_list)
            } else {
                X_re <- NULL
            }
            if (!is.null(mats$S)) {
                S <- bdiag(S_list)
            } else {
                S <- NULL
            }
            private$mats_ <- list(X_fe = X_fe, X_re = X_re, S = S)
            ncol_fe <- rep(mats$ncol_fe, nstates) 
            ncol_re <- rep(mats$ncol_re, nstates) 
            private$terms_ <- list(ncol_fe = ncol_fe,
                                   ncol_re = ncol_re,
                                   names_fe = paste0(rep(colnames(mats$X_fe), nstates), rep(1:nstates, each = ncol(mats$X_fe))),
                                   names_re_all = paste0(rep(colnames(mats$X_re), nstates), rep(1:nstates, each = ncol(mats$X_re))),
                                   names_re = rep(names(mats$ncol_re), nstates))
            # Initial parameters (zero if par0 not provided)
            self$update_coeff_fe(rep(0, sum(ncol_fe)))
            self$update_coeff_re(rep(0, sum(ncol_re)))
            self$update_lambda(rep(1, length(ncol_re)))
            
            # Set initial fixed effect coefficients if provided (par0)
            if(!is.null(par0)) {
                # Number of SDE parameters
                n_par <- length(self$formulas())  
                
                if(length(par0) != n_par) {
                    stop("'par0' should be of length ", n_par,
                         " with one entry for each SDE parameter (",
                         paste0(names(self$formulas()), collapse = ", "), ")")
                }
                
                # First column of X_fe for each SDE parameter
                i0 <- c(1, cumsum(mats$ncol_fe)[-n_par] + 1)
                
                # Apply link to get parameters on working scale
                for (b in 1:nstates) {
                    private$coeff_fe_[i0 + (b - 1) * sum(mats$ncol_fe)] <- sapply(1:n_par, function(i) {
                        self$link()[[i]](par0[i])
                    })
                }
            }
            if (!is.null(tpm0)) {
                private$tpm_ <- tpm0
            } else {
                private$tpm_ <- matrix(1/nstates, nr = nstates, nc = nstates)
            }
            private$other_data_ <- other_data
            private$nstates_ <- nstates
        }, 
        
        #' @description TMB setup
        #'  
        #' @details This creates an attribute \code{tmb_obj}, which can be used to 
        #' evaluate the negative log-likelihood function.
        #' 
        #' @param silent Logical. If TRUE, all tracing outputs are hidden (default).
        #' @param map List passed to MakeADFun to fix parameters. (See TMB documentation.)
        setup = function(silent = TRUE, map = NULL) {
            # Number of time steps
            n <- nrow(self$data())
            
            # Create model matrices
            X_fe <- self$mats()$X_fe
            X_re <- self$mats()$X_re
            S <- self$mats()$S
            ncol_fe <- self$terms()$ncol_fe
            ncol_re <- self$terms()$ncol_re
            
            # Get tpm parameters
            log_tpm <- log(self$tpm() / diag(self$tpm()))
            log_tpm <- as.vector(log_tpm[!diag(self$nstates())])
            
            # Format initial parameters for TMB
            # (First fixed effects, then random effects)
            tmb_par <- list(coeff_fe = self$coeff_fe(),
                            log_lambda = 0,
                            coeff_re = 0, 
                            log_tpm = log_tpm)
            
            # Setup random effects
            random <- NULL
            if(is.null(S)) {
                # If there are no random effects, 
                # coeff_re and log_lambda are not estimated
                map <- c(map, list(coeff_re = factor(NA),
                                   log_lambda = factor(NA)))
                S <- as(matrix(0, 1, 1), "sparseMatrix")
                ncol_re <- 0
                X_re <- as(rep(0, nrow(X_fe)), "sparseMatrix")
            } else {
                # If there are random effects, 
                # set initial values for coeff_re and log_lambda
                random <- c(random, "coeff_re")
                tmb_par$coeff_re <- self$coeff_re()
                tmb_par$log_lambda <- log(self$lambda())
            }
            
            # Setup fixed parameters
            if(!is.null(self$fixpar())) {
                # Indices of fixed coefficients in coeff_fe
                ind_fixcoeff <- self$ind_fixcoeff()
                
                # Define vector with a different integer for each coefficient
                # to be estimated, and NA for each fixed coefficient
                coeff_fe_map <- 1:ncol(X_fe)
                coeff_fe_map[ind_fixcoeff] <- NA
                
                # Update map (to be passed to TMB)
                map <- c(map, list(coeff_fe = factor(coeff_fe_map)))
            }
            
            # TMB data object
            tmb_dat <- list(type = self$type(),
                            ID = self$data()$ID,
                            times = self$data()$time,
                            coarse_time = self$data()$coarse_time, 
                            n_coarse = length(unique(self$data()$coarse_time)),
                            n_states = self$nstates(), 
                            obs = as.matrix(self$obs()),
                            X_fe = X_fe,
                            X_re = X_re,
                            S = S,
                            ncol_re = ncol_re)
            
            # Model-specific data objects
            if(self$type() == "BM-t_HMM") {
                # Pass degrees of freedom for BM-t model
                tmb_dat$other_data <- self$other_data()$df
            } else if(self$type() == "BM_HMM" | self$type() == "OU_HMM") {
                # No extra data needed for BM and OU models
                tmb_dat$other_data <- 0                
            } else if(self$type() == "CTCRW_HMM") {
                # Number of dimensions
                n_dim <- ncol(self$obs())
                # Define initial state and covariance for Kalman filter
                # First index for each ID
                i0 <- c(1, which(self$data()$ID[-n] != self$data()$ID[-1]) + 1)
                # Initial state = (x1, 0, y1, 0, ...)
                a0 <- matrix(0, length(i0), 2*n_dim)
                for(i in 1:n_dim) {
                    a0[, 2*(i-1)+1] <- self$obs()[i0, i]
                }
                tmb_dat$a0 <- a0
                # Initial state covariance
                if(is.null(self$other_data()$P0)) {
                    # Default if P0 not provided by user
                    tmb_dat$P0 <- diag(rep(c(1, 10), n_dim))                    
                } else {
                    tmb_dat$P0 <- self$other_data()$P0
                }
            } 
            
            # Create TMB object
            tmb_obj <- MakeADFun(data = tmb_dat, parameters = tmb_par, 
                                 dll = "smoothSDE", silent = silent,
                                 map = map, random = random)
            
            # Negative log-likelihood function
            private$tmb_obj_ <- tmb_obj
        }, 
        
        nstates = function() {return(private$nstates_)}, 
        
        #' @description Get SDE parameters
        #' 
        #' @param t Time points for which the parameters should be returned.
        #' If "all", returns parameters for all time steps. Default: 1.
        #' @param X_fe Optional design matrix for fixed effects, as returned
        #' by \code{make_mat}. By default, uses design matrix from data.
        #' @param X_re Optional design matrix for random effects, as returned
        #' by \code{make_mat}. By default, uses design matrix from data.
        #' @param coeff_fe Optional vector of fixed effect parameters
        #' @param coeff_re Optional vector of random effect parameters
        #' @param resp Logical (default: TRUE). Should the output be on 
        #' the response scale? If FALSE, the output is on the linear 
        #' predictor scale.
        #' @param T total number of time points 
        #' 
        #' @return Matrix with one row for each time point in t, and one
        #' column for each SDE parameter
        par = function(t = 1, state = 1, X_fe = NULL, X_re = NULL, 
                       coeff_fe = NULL, coeff_re = NULL, 
                       resp = TRUE, n = NULL) {
            # Use design matrices from data if not provided
            if(is.null(X_fe)) {
                X_fe <- self$mats()$X_fe
            }
            
            if(is.null(X_re)) {
                # Apply decay if necessary
                if(is.null(self$other_data()$t_decay)) {
                    X_re <- self$mats()$X_re
                } else {
                    X_re <- self$X_re_decay()
                }
            }
            
            # Use estimated coeff if not provided
            if(is.null(coeff_fe))
                coeff_fe <- self$coeff_fe()
            if(is.null(coeff_re))
                coeff_re <- self$coeff_re()
            
            # if null, use observations 
            if (is.null(n)) {
                n <- nrow(self$data())
            }
            
            # Get linear predictor and put into matrix where each row
            # corresponds to a time step and each column to a parameter
            lp <- X_fe %*% coeff_fe + X_re %*% coeff_re
            lp_mat <- matrix(lp, ncol = length(self$formulas()))
            
            # Apply inverse link to get parameters on natural scale
            if(resp) {
                par_mat <- matrix(NA, nrow = nrow(lp_mat), ncol = ncol(lp_mat))
                for(i in 1:ncol(lp_mat)) {
                    par_mat[,i] <- self$invlink()[[i]](lp_mat[,i])
                }                
            } else {
                par_mat <- lp_mat
            }
            colnames(par_mat) <- names(self$invlink())
            
            # Keep rows of par_mat given in 't'
            if(length(t) == 1) {
                if(t == "all")
                    t <- 1:n
            }
            if(any(t < 1 | t > nrow(par_mat))) {
                stop("'t' should be between 1 and", nrow(par_mat))
            }
            par_mat <- par_mat[t + (state - 1) * n,, drop = FALSE]
            
            return(par_mat)
        }, 
        
        tpm = function() {return(private$tpm_)}, 
        
        update_tpm = function(new_tpm) {private$tpm_ <- new_tpm}
        
    ), 
                    
    private = list(
        nstates_ = NULL, 
        tpm_ = NULL
                        
    )
)
