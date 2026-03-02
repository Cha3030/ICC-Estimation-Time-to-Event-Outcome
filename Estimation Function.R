library(survival)
library(tidyverse)
library(lme4)
library(coxme)
library(parfm)


# =============================================================================
# ICC ESTIMATION
# =============================================================================
# The following 10 methods for estimating ICC are presented. 
# The dataset must consist of clustered time-to-event data. In our examples, 
# clusters are represented by hospital, the event indicator is status 
# (1 = event, 0 = censored), event times are stored in eventtime, and 
# interval defines discrete time steps for methods that require discretization.

icc_estimation <- function(data) {
  
  # --- Method 1: Weibull model with log-normal frailty (Williams, 1995) 
  # This parametric approach utilizes a Weibull survival model with log-normal 
  # frailty. The ICC is computed analytically as the ratio of the cluster-level 
  # variance to the sum of the cluster and residual variances, defining the 
  # residual variance as pi^2 / (6 * gamma^2).
  icc_weibull <- tryCatch({
    fit_parfm <- parfm(
      Surv(eventtime, status) ~ 1,
      cluster = "hospital",
      data = data,
      dist = "weibull",
      frailty = "lognormal",
      maxit = 100,  
      correct = 0   
    )
    a <- as_tibble(fit_parfm)
    gamma <- as.numeric(a[2,1])
    var_cluster2 <- as.numeric(a[1,1])
    var_resid <- pi^2 / (6 * gamma^2)
    var_cluster2 / (var_cluster2 + var_resid)
  }, error = function(e) NA)
  
  # --- Method 2: Proportion of event (Xie & Waksman, 2003)
  # ICC estimated with a closed-form formula directly from clustered 
  # survival data using the overall event rate (p). The numerator sums the products 
  # of deviations of pairs within clusters, while the denominator normalizes 
  # by p(1-p) and the total number of pairs. 
  p_hat <- mean(data$status, na.rm = TRUE)
  data$dstatus <- data$status - p_hat
  
  agg <- by(data$dstatus, data$hospital, function(x) {
    s <- sum(x)
    ssq <- sum(x^2)
    list(K = length(x), pair_sum = (s^2 - ssq))
  })
  
  agg_list <- lapply(agg, unclass)
  Ks <- vapply(agg_list, `[[`, numeric(1), "K")
  pair_sums <- vapply(agg_list, `[[`, numeric(1), "pair_sum")
  
  numerator <- sum(pair_sums)
  denom_pairs <- sum(Ks * (Ks - 1))
  denominator <- p_hat * (1 - p_hat) * denom_pairs
  icc_proportion <- if (denominator != 0) numerator / denominator else NA
  
  # --- Method 3: Martingale residual (Gagnon, 2004)
  # This nonparametric method estimates the ICC using martingale residuals 
  # from Cox models. It quantifies the proportion of variation attributable 
  # to clustering by decomposing between- and within-cluster variability.
  cox_fit <- coxph(Surv(eventtime, status) ~ 1, data = data)
  data$mart_resid <- residuals(cox_fit, type = "martingale")
  
  cluster_means_m <- data %>%
    group_by(hospital) %>%
    summarise(mean_resid = mean(mart_resid, na.rm = TRUE),
              n = n(), .groups = "drop")
  
  grand_mean_m <- mean(data$mart_resid, na.rm = TRUE)
  MSB_m <- sum(cluster_means_m$n * (cluster_means_m$mean_resid - grand_mean_m)^2) /
    (nrow(cluster_means_m) - 1)
  
  data_m <- left_join(data, cluster_means_m, by = "hospital")
  MSW_m <- sum((data_m$mart_resid - data_m$mean_resid)^2) /
    (nrow(data_m) - nrow(cluster_means_m))
  
  icc_martingale <- (MSB_m - MSW_m) / (MSB_m + (mean(cluster_means_m$n) - 1) * MSW_m)
  
  # --- Method 4: Censoring indicators (Kalia, Klar & Donner, 2016)
  # ICC estimated with closed-form formula from the between and within-cluster 
  # variance of censoring indicators. A binary variable (1 event observed, 0 censored) 
  # is created, and the ICC is computed using ANOVA-type variance decomposition.
  # It assumes approximately equal cluster sizes and ignores exact event times.
  grand_mean_c <- mean(data$status, na.rm = TRUE)
  cluster_n_c <- tapply(data$status, data$hospital, length)
  cluster_mean_c <- tapply(data$status, data$hospital, mean)
  
  r_c <- length(cluster_n_c)
  MSB_c <- sum(cluster_n_c * (cluster_mean_c - grand_mean_c)^2) / (r_c - 1)
  
  cluster_mean_ind_c <- cluster_mean_c[data$hospital]
  MSW_c <- sum((data$status - cluster_mean_ind_c)^2) / (nrow(data) - r_c)
  
  m_bar_c <- mean(cluster_n_c)
  icc_censoring <- (MSB_c - MSW_c) / (MSB_c + (m_bar_c - 1) * MSW_c)
  
  # --- Method 5: Observed event times (Kalia, Klar & Donner, 2016)
  # ICC estimated with closed-form formula directly from observed 
  # survival times using regression-based decomposition of variance, 
  # considering only uncensored events (status = 1).
  obs_data <- data[data$status == 1, ]
  obs_data$hospital <- as.factor(obs_data$hospital)
  obs_data$hospital <- droplevels(obs_data$hospital)
  
  if (nrow(obs_data) > 0 && length(unique(obs_data$hospital)) > 1) {
    cluster_means_e <- tapply(obs_data$eventtime, obs_data$hospital, mean)
    overall_mean_e <- mean(obs_data$eventtime)
    mi_e <- table(obs_data$hospital)
    r_e <- length(mi_e)
    N_star_e <- sum(mi_e)
    m_o_e <- (1 / (r_e - 1)) * (N_star_e - sum(mi_e^2) / N_star_e)
    
    MSB_e <- (1 / (r_e - 1)) * sum(mi_e * (cluster_means_e - overall_mean_e)^2)
    MSW_e <- (1 / (N_star_e - r_e)) * sum(tapply(obs_data$eventtime, obs_data$hospital,
                                                 function(x) sum((x - mean(x))^2)))
    icc_event <- (MSB_e - MSW_e) / (MSB_e + (m_o_e - 1) * MSW_e)
  } else { icc_event <- NA }
  
  # --- Method 6: Exponential-Normal Approximation (Oliveira et al., 2016)
  # Special case of the Weibull–Normal model. Fits a  survival model with a single
  # log-normal frailty for clustering. The ICC is computed using a theoretical
  # formula based on frailty variance to determine the cluster effect proportion.
  icc_exp_norm_fit <- tryCatch({
    coxme(Surv(eventtime, status) ~ 1 + (1|hospital), data = data)
  }, error = function(e) NULL)
  
  icc_exp_norm <- if (!is.null(icc_exp_norm_fit)) {
    D <- VarCorr(icc_exp_norm_fit)$hospital[[1]]
    1 - (1 / (exp(D) * (2 * exp(D) - 1)))
  } else { NA }
  
  # --- Method 7: CoxME with Gaussian frailty (McCune, 2025)
  # A Cox proportional hazards model with a Gaussian random intercept 
  # for hospital. ICC is estimated as the fraction of cluster-level random 
  # effect variance, with residual variance approximated as pi^2/6. 
  mod_coxme <- tryCatch({
    coxme(Surv(eventtime, status) ~ 1 + (1|hospital), data = data)
  }, error = function(e) NULL)
  
  icc_coxme <- if (!is.null(mod_coxme)) {
    var_coxme <- VarCorr(mod_coxme)$hospital[[1]]
    var_coxme / (var_coxme + pi^2 / 6)
  } else { NA }
  
  # --- Method 8: Cox with Gamma frailty (McCune, 2025)
  # A Cox proportional hazards model with a gamma-distributed frailty. 
  # ICC is estimated from the cluster-level random effect variance, 
  # with the residual variance estimated numerically via Laplace approximation.
  cox_frailty_gamma <- tryCatch({
    coxph(Surv(eventtime, status) ~ 1 + frailty(hospital, distribution = "gamma"), data = data)
  }, error = function(e) NULL)
  
  var_gamma <- if (!is.null(cox_frailty_gamma)) {
    as.numeric(cox_frailty_gamma$history$`frailty(hospital, distribution = "gamma")`[1])
  } else { NA }
  
  if (!is.na(var_gamma)) {
    g <- function(w, k, s, sigma2) { -k * w + exp(w) * s + w^2 / (2 * sigma2) }
    g2 <- function(w, k, s, sigma2) { exp(w) * s + 1 / sigma2 }
    Lapl <- Vectorize(function(s, k, sigma2) {
      wTilde <- optimize(f = g, c(-1e10, 1e10), maximum = FALSE, k=k, s=s, sigma2=sigma2)$minimum
      (-1)^k * exp(-g(wTilde, k, s, sigma2)) / sqrt(sigma2 * g2(wTilde, k, s, sigma2))
    }, "s")
    fr.lognormal <- function(k, s, sigma2) {
      intTau <- Vectorize(function(x, i.sigma2 = sigma2) {
        x * Lapl(x, 0, i.sigma2) * Lapl(x, 2, i.sigma2)
      }, "x")
      4 * integrate(intTau, 0, Inf, i.sigma2 = sigma2)$value - 1
    }
    icc_gamma_np <- fr.lognormal(k = 1, s = 1, sigma2 = var_gamma)
  } else { icc_gamma_np <- NA }
  
  # --- Method 9: Binomial GLMM on discretized time (McCune, 2025; Lam & Ip, 2003)
  # Survival times are discretized into an expanded dataset. A GLMM 
  # with a cloglog link is fitted with a random intercept for cluster. ICC is 
  # calculated as the fraction of group variance, with residual variance = pi^2/6.
  q <- quantile(data$eventtime, probs = seq(0.2, 0.8, by = 0.2), na.rm = TRUE)
  q <- unique(q)
  tmin <- min(data$eventtime, na.rm = TRUE)
  tmax <- max(data$eventtime, na.rm = TRUE)
  q <- q[q > tmin & q < tmax]
  
  icc_glmm <- NA
  discretized_db <- NULL
  
  if (length(q) >= 1) {
    discretized_db <- tryCatch({
      survSplit(Surv(eventtime, status) ~ ., data = data[, c("eventtime", "status", "hospital")],
                cut = q, start = "tstart", end = "tstop")
    }, error = function(e) NULL)
    
    if (!is.null(discretized_db)) {
      discretized_db$interval <- findInterval(discretized_db$tstart, c(0, q)) + 1
      mod_glmm <- tryCatch({
        glmer(status ~ as.factor(interval) + (1|hospital), data = discretized_db, 
              family = binomial(link = "cloglog"), nAGQ = 5)
      }, error = function(e) NULL)
      if (!is.null(mod_glmm) && !isSingular(mod_glmm)) {
        var_glmm <- mod_glmm@theta^2
        icc_glmm <- var_glmm / (var_glmm + pi^2 / 6)
      }
    }
  }
  
  # --- Method 10: CoxME on discretized time (McCune, 2025)
  # Survival times are discretized (expanded dataset), and a CoxME 
  # model is fitted with a Gaussian random intercept for cluster. ICC is 
  # the fraction of group variance, with residual variance approximated as $\pi^2 / 6$.
  icc_cox_exploded <- NA
  if (!is.null(discretized_db)) {
    mod_cox_exploded <- tryCatch({
      coxme(Surv(tstart, tstop, status) ~ 1 + (1|hospital), data = discretized_db)
    }, error = function(e) NULL)
    if (!is.null(mod_cox_exploded)) {
      var_exploded <- VarCorr(mod_cox_exploded)$hospital[[1]]
      icc_cox_exploded <- var_exploded / (var_exploded + pi^2 / 6)
    }
  }
  
  # --- Output Result Table
  data.frame(
    Method = paste0(1:10, ". ", c("Weibull", "Proportion of event", "Martingale",
                                  "Censoring indicators", "Observed event times",
                                  "Exponential-Normal", "CoxME Gaussian frailty",
                                  "Cox gamma frailty (NP)", "GLMM logit discretized",
                                  "CoxME discretized time")),
    ICC = round(c(icc_weibull, icc_proportion, icc_martingale, icc_censoring, 
                  icc_event, icc_exp_norm, icc_coxme, icc_gamma_np, 
                  icc_glmm, icc_cox_exploded), 3)
  )
}



# =============================================================================
# BOOTSTRAP ICC
# =============================================================================
# This function performs cluster-level bootstrapping by sampling entire 
# clusters (e.g., hospitals) with replacement. For each iteration, all 
# subjects within the selected clusters are included. 
#
# USAGE NOTE:
# This method is essential for assessing the ICC and its Confidence Interval (CI) 
# when working with a single real-world dataset. In contrast, for simulation 
# studies, a Monte Carlo approach was preferred (generating multiple 
# datasets with fixed parameters).

bootstrap_icc <- function(data, B = 100, cluster_var = "hospital") {
  
  boot_results <- vector("list", B)
  
  clusters <- unique(data[[cluster_var]])
  
  for (b in 1:B) {
    sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)

        data_boot <- do.call(rbind, lapply(seq_along(sampled_clusters), function(i) {
      cl <- sampled_clusters[i]
      subset_data <- data[data[[cluster_var]] == cl, ]
      subset_data[[cluster_var]] <- i 
      subset_data
    }))
    
    data_boot[[cluster_var]] <- factor(data_boot[[cluster_var]])
    
    # --- ICC Estimation on Bootstrap Sample ---
    res <- tryCatch({
      icc_estimation(data_boot)$ICC
    }, error = function(e) {
       rep(NA, 10) 
    })
    
    boot_results[[b]] <- res
  }
  
  boot_mat <- do.call(rbind, boot_results)
  colnames(boot_mat) <- icc_estimation(data)$Method
  
  # --- Statistical Summary ---
  boot_summary <- data.frame(
    Method  = colnames(boot_mat),
    Mean    = apply(boot_mat, 2, mean, na.rm = TRUE),
    SD      = apply(boot_mat, 2, sd, na.rm = TRUE),
    CI_low  = apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    CI_high = apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    Median  = apply(boot_mat, 2, quantile, probs = 0.500, na.rm = TRUE),
    Q1      = apply(boot_mat, 2, quantile, probs = 0.250, na.rm = TRUE),
    Q3      = apply(boot_mat, 2, quantile, probs = 0.750, na.rm = TRUE)
  )
  
  return(list(
    summary    = boot_summary,
    replicates = boot_mat
  ))
}



# # EXAMPLE 1
# df_example1 <- read.csv("dataset/Pbc3.csv")%>%
#   mutate(fac_unit = as.factor(fac_unit))%>%
#   rename(eventtime = time,
#          status = event,
#          hospital = fac_unit)
# 
# results_icc1 <- icc_estimation(df_example1)
# boot_out1 <- bootstrap_icc(df_example1, B = 100)
# print(boot_out1$summary)
# 
# 
# # EXAMPLE 2
# load("dataset/survival_panitumumab_chemio.rda")
# df_example2 <- survival_data_filter_hosp_plus2 %>%
#   rename(eventtime = time,
#          status = dth,
#          hospital = hospital) %>%
#   droplevels()
# 
# df_example2$hospital <- droplevels(df_example2$hospital)
# 
# 
# results_icc2 <- icc_estimation(df_example2)
# boot_out2 <- bootstrap_icc(df_example2, B = 100)
# print(boot_out2$summary)