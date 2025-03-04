# aplicar n.min a dat
# Credit: Joaquim Radua, 2024



################################################################################
#                                                                              #
# combat_fit and combat_apply functions                                        #
# =====================================                                        #
#                                                                              #
# Authors:                                                                     #
#                                                                              #
#  1) The original ComBat function was in the sva package that can be found at #
#     https://bioconductor.org/packages/release/bioc/html/sva.html             #
#                                                                              #
#  2) First modification by Jean-Philippe Fortin for the harmonization of MRI  #
#     data. Please cite https://10.1016/j.neuroimage.2017.08.047               #
#                                                                              #
#  3) Second modification by Joaquim Radua to separate functions for fitting   #
#     and applying the harmonization, allow missings and constant rows and mi- #
#     nor changes in the arguments of the functions to facilitate their use.   #
#     Please cite "Increased power by harmonizing structural MRI site diffe-   #
#     rences with the ComBat batch adjustment method in ENIGMA".               #
#                                                                              #
#  4) Third modification by Joaquim Radua to include xxx
#                                                                              #
# The original and present code is under the Artistic License 2.0. If using    #
# this code, make sure you agree and accept this license.                      #
#                                                                              #
# Additions by Joaquim Radua: use lme instead of ComBat                        #
#                                                                              #
################################################################################


# Calculate some characteristics on the batches
.combat_tmp1 <- function (dat, batch, levels_batch, mod) {
  batchmod <- model.matrix(~ -1 + batch)
  # A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels_batch[i])
  }
  # List of samples in each batch
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  # Combine batch variable and covariates
  design <- cbind(batchmod, mod)
  # Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function (x) all(x == 1))
  design <- as.matrix(design[, !check])
  batch.design <- design[,1:n.batch]
  return(list(
    dat = dat,
    batchmod = batchmod,
    n.batch = n.batch,
    batches = batches,
    n.batches = n.batches,
    n.array = n.array,
    design = design,
    batch.design = batch.design
  ))
}
# Estimate B.hat, grand.mean and var.pooled
.combat_tmp2 <- function (tmp1, verbose) {
  # Number of covariates or covariate levels
  if (verbose) {
    cat(
      "[combat.enigma] Adjusting for",
      ncol(tmp1$design) - ncol(tmp1$batchmod),
      "covariate(s) or covariate level(s)\n"
    )
  }
  # Check if the design is confounded
  if (qr(tmp1$design)$rank < ncol(tmp1$design)) {
    if (ncol(tmp1$design) == (tmp1$n.batch + 1)) {
      stop("[combat.enigma] The covariate is confounded with batch. Remove the covariate and rerun ComBat.")
    }
    if (ncol(tmp1$design) > (tmp1$n.batch + 1)) {
      if ((qr(tmp1$design[,-c(1:tmp1$n.batch)])$rank < ncol(tmp1$design[,-c(1:tmp1$n.batch)]))) {
        stop("The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.")
      } else {
        stop("At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.")
      }
    }
  }
  # Standardize data across features
  B.hat <- solve(t(tmp1$design) %*% tmp1$design) %*% t(tmp1$design) %*% t(as.matrix(tmp1$dat))
  # Standarization Model
  grand.mean <- t(tmp1$n.batches / tmp1$n.array) %*% B.hat[1:tmp1$n.batch,]
  var.pooled <- ((tmp1$dat - t(tmp1$design %*% B.hat))^2) %*% rep(1 / tmp1$n.array, tmp1$n.array)
  return(list(
    B.hat = B.hat,
    grand.mean = grand.mean,
    var.pooled = var.pooled
  ))
}
# Standardize data
.combat_tmp3 <- function (dat, tmp1, tmp2, verbose) {
  if (verbose) {
    cat("[combat.enigma] Standardizing data across features\n")
  }
  stand.mean <- t(tmp2$grand.mean) %*% t(rep(1, tmp1$n.array))
  if (!is.null(tmp1$design)) {
    tmp <- tmp1$design;tmp[,c(1:tmp1$n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% tmp2$B.hat)
  } 
  s.data <- (dat-stand.mean) / (sqrt(tmp2$var.pooled) %*% t(rep(1, tmp1$n.array)))
  return(list(
    stand.mean = stand.mean,
    s.data = s.data
  ))
}
# Fit L/S model 
.combat_tmp4 <- function (tmp1, tmp2, tmp3, eb, verbose) {
  # Get regression batch effect parameters
  if (eb) {
    if (verbose) {
      cat("[combat.enigma] Fitting L/S model and finding priors\n")
    }
  } else {
    if (verbose) {
      cat("[combat.enigma] Fitting L/S model\n")
    }
  }
  gamma.hat <- solve(t(tmp1$batch.design) %*% tmp1$batch.design) %*% t(tmp1$batch.design) %*% t(as.matrix(tmp3$s.data))
  delta.hat <- NULL
  for (i in tmp1$batches) {
    delta.hat <- rbind(delta.hat, apply(tmp3$s.data[,i], 1, var, na.rm = T))
  }
  # Empirical Bayes correction:
  gamma.star <- delta.star <- NULL
  gamma.bar <- t2 <- a.prior <- b.prior <- NULL
  if (eb) {
    # Find Priors
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, .aprior)
    b.prior <- apply(delta.hat, 1, .bprior)
    # Find EB batch adjustments
    if (verbose) {
      cat("[combat.enigma] Finding parametric adjustments\n")
    }
    for (i in 1:tmp1$n.batch) {
      temp <- .it.sol(tmp3$s.data[,tmp1$batches[[i]]], gamma.hat[i,], delta.hat[i,], gamma.bar[i], t2[i], a.prior[i], b.prior[i])
      gamma.star <- rbind(gamma.star, temp[1,])
      delta.star <- rbind(delta.star, temp[2,])
    }
  } 
  return(list(
    gamma.hat = gamma.hat,
    delta.hat = delta.hat, 
    gamma.star = gamma.star,
    delta.star = delta.star, 
    gamma.bar = gamma.bar,
    t2 = t2,
    a.prior = a.prior,
    b.prior = b.prior
  ))
}
# Adjust the data
.combat_tmp5 <- function (tmp1, tmp2, tmp3, tmp4, eb, verbose) {
  # Normalize the data
  if (verbose) {
    cat("[combat.enigma] Adjusting the data\n")
  }
  bayesdata <- tmp3$s.data
  j <- 1
  for (i in tmp1$batches) {
    if (eb) {
      bayesdata[,i] <- (bayesdata[,i] - t(tmp1$batch.design[i,] %*% tmp4$gamma.star)) / (sqrt(tmp4$delta.star[j,]) %*% t(rep(1, tmp1$n.batches[j])))
    } else {
      bayesdata[,i] <- (bayesdata[,i] - t(tmp1$batch.design[i,] %*% tmp4$gamma.hat )) / (sqrt(tmp4$delta.hat[j,] ) %*% t(rep(1, tmp1$n.batches[j])))
    }
    j <- j + 1
  }
  return((bayesdata * (sqrt(tmp2$var.pooled) %*% t(rep(1, tmp1$n.array)))) + tmp3$stand.mean)
}
# Fit lme model
.lme_tmp1 <- function (dat, batch, levels_batch, mod) {
  random = NULL
  for (i in 1:nrow(dat)) {
    dat_i = dat[i,]
    if (!is.null(mod)) {
      m = tryCatch({
        lme(dat_i ~ mod, random = ~ 1 | batch, na.action = na.omit)
      }, error = function (e) {
          warning("lme with cov failed, using lme without cov")
          m = lme(dat_i ~ 1, random = ~ 1 | batch, na.action = na.omit)
      })
    } else {
      m = lme(dat_i ~ 1, random = ~ 1 | batch, na.action = na.omit)
    }
    random_i = m$coefficients$random$batch
    random_i = t(t(random_i[match(levels_batch, rownames(random_i)),]))
    colnames(random_i) = rownames(dat)[i]
    random = cbind(random, random_i)
  }
  return(list(
    random = random
  ))
}
# Adjust the data
.lme_tmp2 <- function (dat, batch, levels_batch, tmp1, verbose) {
  if (verbose) {
    cat("[combat.enigma] Adjusting the data\n")
  }
  return(dat - t(tmp1$random)[,match(batch, levels_batch)])
}
# Impute missings in cov
.impute_cov = function(cov, site, cov_means) {
  for (batch_i in sort(unique(site))) {
    i <- which(site == batch_i)
    cov_means_i = cov_means[which(rownames(cov_means) == batch_i),]
    for (j in 1:ncol(cov)) {
      is_na <- which(is.na(cov[i,j]))
      if (length(is_na) > 0) {
        if (!is.na(cov_means_i[j])) {
          cov[i[is_na],j] <- cov_means_i[j]
        } else {
          cov[i[is_na],j] <- median(cov_means[,j], na.rm = TRUE)
        }
      }
    }
  }
  cov
}
combat_fit <- function (dat, site, cov = NULL, method = "combat",
                        n.min = 10,
                        impute_missing_cov = FALSE,
                        prescaling = FALSE,
                        impute_empty_sites = FALSE,
                        eb = TRUE,
                        verbose = TRUE) {
  
  # CHECK ARGUMENTS ============================================================
  # Check site and verbose before dat because they are used when checking dat
  if (!is.factor(site)) {
    stop("site must be a factor")
  }
  if (!(verbose %in% c(TRUE, FALSE))) {
    stop("verbose must be TRUE or FALSE")
  }
  if (is.data.frame(dat)) {
    dat <- as.matrix(dat)
  } else if (!is.matrix(dat)) {
    stop("dat must be a matrix")
  }
  if (ncol(dat) == length(site)) {
    transpose <- FALSE
    if (verbose) {
      cat("[combat.enigma] Subjects are COLUMNS\n")
    }
  } else if (nrow(dat) == length(site)) {
    transpose <- TRUE
    dat <- t(dat)
    if (verbose) {
      cat("[combat.enigma] Subjects are ROWS\n")
    }
  } else {
    stop("dat must have the same number of columns or rows than the length of site")
  }
  # Check impute_missing_cov before cov because it is used when checking cov
  if (!(impute_missing_cov %in% c(TRUE, FALSE))) {
    stop("impute_missing_cov must be TRUE or FALSE")
  }
  if (is.data.frame((cov))) {
    cov <- as.matrix(cov)
  } else if (!(is.matrix(cov) || is.null(cov))) {
    stop("cov must be a matrix or NULL")
  }
  if (any(is.na(cov)) & !impute_missing_cov) {
    stop("missing values in cov, use impute_missing_cov")
  }
  if (!(method %in% c("combat", "lme"))) {
    stop("method must be combat or lme")
  }
  if (!(prescaling %in% c(TRUE, FALSE))) {
    stop("prescaling must be TRUE or FALSE")
  }
  if (method == "combat") {
    if (!(impute_empty_sites %in% c(TRUE, FALSE))) {
      stop("impute_empty_sites must be TRUE or FALSE")
    }
    if (!(eb %in% c(TRUE, FALSE))) {
      stop("eb must be TRUE or FALSE")
    }
  }
  
  # COV IMPUTATION AND PRESCALING ==============================================
  # Pre-scaling find scaling factor for all (ROI or voxel)s
  if (impute_missing_cov) {
    cov_means = apply(cov, 2, by, site, mean, na.rm = TRUE)
    if (any(is.na(cov))) {
      if (verbose) {
        cat("[combat.enigma] Imputing missing covariates\n")
      }
      cov = .impute_cov(cov, site, cov_means)
    }
  }
  if (prescaling) {
    if (verbose) {
      cat("[combat.enigma] Prescaling\n")
    }
    # Calculate SD for each site
    dat_sd = NULL
    for (batch_i in sort(unique(site))) {
      i <- which(site == batch_i)
      cov_i = cov[i,]
      dat_sd_i = matrix(apply(dat[,i], 1, function (x) {
        if (sum(!is.na(x)) < n.min) {
          return(NA)
        }
        sd(lm(x ~ cov_i)$residuals)
      }), nrow = 1)
      rownames(dat_sd_i) = batch_i
      dat_sd = rbind(dat_sd, dat_sd_i)
    }
    # Get (ROI or voxel)s with information from all sites
    common_dat_sd = dat_sd[, which(apply(dat_sd, 2, function (x) {all(!is.na(x))}))]
    # Find prescaling factors
    prescaling_dat = apply(common_dat_sd, 2, function (x) {
      exp(median(log(x)) - log(x))
    })
    prescaling_factors = apply(prescaling_dat, 1, median)
    # Apply prescaling factors
    for (batch_i in sort(unique(site))) {
      i <- which(site == batch_i)
      dat[,i] = prescaling_factors[which(names(prescaling_factors) == batch_i)] * dat[,i]
    }
  }
  
  # DAT IMPUTATION =============================================================
  # Imputation of missing values, separately for each site and (ROI or voxel).
  # (only needed for ComBat)
  if (method == "combat") {
    if (any(is.na(dat))) {
      if (verbose) {
        cat("[combat.enigma] Imputing missing data (only for fit)\n")
      }
      for (batch_i in sort(unique(site))) {
        i <- which(site == batch_i)
        dat_i <- dat[,i]
        if (any(is.na(dat_i))) {
          for (j in 1:nrow(dat)) {
            dat_ji <- dat_i[j,]
            is_na <- which(is.na(dat_ji))
            # Some missing, impute from other subjects of the site
            if (length(is_na) > 0 && length(is_na) < length(i)) {
              if (!is.null(cov)) {
                if (length(is_na) == 1) {
                  mod_i_is_na <- matrix(cov[i[is_na],], nrow = 1)
                } else {
                  mod_i_is_na <- cov[i[is_na],]
                }
                beta <- matrix(coef(lm(dat_ji ~ cov[i,])))
                beta[which(is.na(beta))] <- 0
                dat[j,i[is_na]] <- cbind(1, mod_i_is_na) %*% beta
              } else {
                dat[j,i[is_na]] <- mean(dat_ji, na.rm = TRUE)
              }
            }
          }
        }
      }
      # If there are still missing data, they are because a site has all missing
      # data for a (ROI or voxel)
      for (batch_i in sort(unique(site))) {
        i <- which(site == batch_i)
        if (any(is.na(dat[,i]))) {
          if (!impute_empty_sites) {
            stop("empty site, use impute_empty_sites or lme method")
          }
          for (j in 1:nrow(dat)) {
            dat_j <- dat[j,]
            if (is.na(dat_j[i[1]])) {
              if (!is.null(cov)) {
                beta <- matrix(coef(lm(dat_j ~ cov)))
                beta[which(is.na(beta))] <- 0
                dat[j,i] <- cbind(1, cov[i,]) %*% beta
              } else {
                dat[j,i] <- mean(dat_j, na.rm = TRUE)
              }
            }
          }
        }
      }
    }
  }
  
  # COMBAT AND OUTPUT ==========================================================
  not_constant <- which(apply(dat, 1, function (x) {var(x, na.rm = TRUE) > 0}))
  dat <- dat[not_constant,]
  levels_batch <- levels(site)
  if (method == "combat") {
    if (verbose) {
      if (eb) {
        cat("[combat.enigma] Fitting ComBat with empirical Bayes\n")
      } else {
        cat("[combat.enigma] Fitting ComBat without empirical Bayes\n")
      }
    }
    tmp1 <- .combat_tmp1(dat, site, levels_batch, cov)
    tmp2 <- .combat_tmp2(tmp1, verbose)
    tmp3 <- .combat_tmp3(dat, tmp1, tmp2, verbose)
    tmp4 <- .combat_tmp4(tmp1, tmp2, tmp3, eb, verbose)
    combat_parameters = list(
      levels_batch = levels_batch,
      transpose = transpose,
      not_constant = not_constant,
      method = "combat",
      impute_missing_cov = impute_missing_cov,
      prescaling = prescaling,
      impute_empty_sites = impute_empty_sites,
      eb = eb,
      tmp2 = tmp2,
      tmp4 = tmp4
    )
  }
  if (method == "lme") {
    if (verbose) {
      cat("[combat.enigma] Fitting linear mixed-effects models\n")
    }
    tmp1 <- .lme_tmp1(dat, site, levels_batch, cov)
    combat_parameters = list(
      levels_batch = levels_batch,
      transpose = transpose,
      not_constant = not_constant,
      method = "lme",
      impute_missing_cov = impute_missing_cov,
      prescaling = prescaling,
      tmp1 = tmp1
    )
  }
  if (impute_missing_cov) {
    combat_parameters$cov_means = cov_means
  }
  if (prescaling) {
    combat_parameters$prescaling_factors = prescaling_factors
  }
  class(combat_parameters) = "combat.enigma"
  return(combat_parameters)
}

combat_apply <- function (combat_parameters, dat, site, cov = NULL,
                          verbose = TRUE) {
  
  # CHECK ARGUMENTS ============================================================
  if (class(combat_parameters) != "combat.enigma") {
    stop("combat_parameters must be a list of class combat.enigma")
  }
  # Check site and verbose before dat because they are used when checking dat
  if (!is.factor(site) || nlevels(site) != length(combat_parameters$levels_batch) || any(levels(site) != combat_parameters$levels_batch)) {
    stop("site must be a factor with the same levels than when fitting combat")
  }
  if (!(verbose %in% c(TRUE, FALSE))) {
    stop("verbose must be TRUE or FALSE")
  }
  if (is.data.frame(dat)) {
    dat <- as.matrix(dat)
  } else if (!is.matrix(dat)) {
    stop("dat must be a matrix")
  }
  if (combat_parameters$transpose) {
    if (nrow(dat) != length(site)) {
      stop("dat must have the same number of rows than the length of site")
    }
    dat <- t(dat)
  } else {
    if (ncol(dat) != length(site)) {
      stop("dat must have the same number of columns than the length of site")
    }
  }
  if (is.data.frame((cov))) {
    cov <- as.matrix(cov)
  } else if (!(is.matrix(cov) || is.null(cov))) {
    stop("cov must be a matrix or NULL")
  }
  if (any(is.na(cov)) && !combat_parameters$impute_missing_cov) {
    stop("missing values in cov, use impute_missing_cov in combat_fit")
  }
  
  # COV IMPUTATION AND PRESCALING ==============================================
  if (combat_parameters$impute_missing_cov & any(is.na(cov))) {
    if (verbose) {
      cat("[combat.enigma] Imputing missing covariates\n")
    }
    cov = .impute_cov(cov, site, combat_parameters$cov_means)
  }
  if (combat_parameters$prescaling) {
    if (verbose) {
      cat("[combat.enigma] Prescaling\n")
    }
    for (batch_i in sort(unique(site))) {
      i <- which(site == batch_i)
      dat[,i] = combat_parameters$prescaling_factors[which(names(combat_parameters$prescaling_factors) == batch_i)] * dat[,i]
    }
  }
  
  # COMBAT AND OUTPUT ==========================================================
  dat.combat <- dat
  dat <- dat[combat_parameters$not_constant,]
  if (combat_parameters$method == "combat") {
    tmp1 <- .combat_tmp1(dat, site, combat_parameters$levels_batch, cov)
    tmp3 <- .combat_tmp3(dat, tmp1, combat_parameters$tmp2, verbose)
    tmp5 <- .combat_tmp5(tmp1, combat_parameters$tmp2, tmp3, combat_parameters$tmp4, combat_parameters$eb, verbose)
    dat.combat[combat_parameters$not_constant,] <- tmp5
  }
  if (combat_parameters$method == "lme") {
    tmp2 <- .lme_tmp2(dat, site, combat_parameters$levels_batch, combat_parameters$tmp1, verbose)
    dat.combat[combat_parameters$not_constant,] <- tmp2
  }
  if (combat_parameters$transpose) {
    dat.combat <- t(dat.combat)
  }
  if (combat_parameters$method == "combat") {
    return(list(
      dat.combat = dat.combat, 
      gamma.hat = combat_parameters$tmp4$gamma.hat,
      delta.hat = combat_parameters$tmp4$delta.hat, 
      gamma.star = combat_parameters$tmp4$gamma.star,
      delta.star = combat_parameters$tmp4$delta.star, 
      gamma.bar = combat_parameters$tmp4$gamma.bar,
      t2 = combat_parameters$tmp4$t2,
      a.prior = combat_parameters$tmp4$a.prior,
      b.prior = combat_parameters$tmp4$b.prior,
      stand.mean = tmp3$stand.mean,
      stand.sd = sqrt(combat_parameters$tmp2$var.pooled)[,1]
    ))
  }
  if (combat_parameters$method == "lme") {
    return(list(
      dat.combat = dat.combat,
      random = combat_parameters$tmp1$random
    ))
  }
}


################################################################################
#                                                                              #
# This is a copy of the original code from the standard version of the sva     #
# package that can be found at                                                 #
# https://bioconductor.org/packages/release/bioc/html/sva.html                 #
# The original and present code is under the Artistic License 2.0.             #
# If using this code, make sure you agree and accept this license.             #
#                                                                              #
################################################################################


# Following four find empirical hyper-prior values
.aprior <- function (gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (2 * s2 + m^2) / s2
}
.bprior <- function (gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (m * s2 + m^3) / s2
}
.postmean <- function (g.hat, g.bar, n, d.star, t2) {
  (t2 * n * g.hat + d.star * g.bar) / (t2 * n + d.star)
}
.postvar <- function (sum2, n, a, b) {
  (0.5 * sum2 + b) / (n / 2 + a - 1)
}
# Pass in entire data set, the design matrix for the entire data, the batch
# means, the batch variances, priors (m, t2, a, b), columns of the data 
# matrix for the batch. Uses the EM to find the parametric batch adjustments
.it.sol <- function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 0.0001) {
  n <- apply(!is.na(sdat), 1, sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while (change > conv) {
    g.new <- .postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- apply((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, 1, sum, na.rm = T)
    d.new <- .postvar(sum2, n, a, b)
    change <- max(abs(g.new - g.old) / g.old, abs(d.new - d.old) / d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count + 1
  }
  # cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star", "d.star")
  adjust
}
