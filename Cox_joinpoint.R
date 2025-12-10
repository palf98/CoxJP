library(dplyr)
library(survival)
library(lubridate) # used for the generation of the example dataset

# download periodR package from: https://krebsregister.saarland.de/daten-auswertungen-veroeffentlichungen/software/periodr-english/
# paste the downloaded 'periodR' directory to your R installation's library path. You can check this path running ".libPaths()" in R console.
library(periodR) # used for the estimations of survival by the period approach


###################
#### FUNCTIONS ####
###################

# this function gets the largest time point in a vector before the 'goal' time. Used over survfit() outputs
retrieve_time <- function (goal, vector){
  vector[which.min(goal - vector[which(goal - vector >= 0)])]
} 

## DEFAULT: 
# minimum 3 years before first joinpoint (startbound = 3)
# minimum 5 years after last joinpoint (endbound = 5)
# minimum 2 years between joinpoints (inbound = 3)
# bic_threshold = 0: by default, best = T returns model with lowest BIC. If bic_threshold is greater than 0, a more complex model has to have a lower BIC surpassing such threshold with respect to a simpler model.
# best = F: function returns a list with best models with 0, 1 and 2 change points

CoxJP <- function(data, Surv_var, year_var, agegr_var, n_jp = 2, startbound = 3,
                  endbound = 5, inbound = 3, bic_threshold = 0, best = T){
  if (n_jp > 2) return(print("Function is not generalized to more than 2 joinpoints."))
  if(nlevels(data[[agegr_var]]) > 1){
    form <- paste0(Surv_var, " ~ ", year_var, " + strata(", agegr_var, ")")
  } else{
    form <- paste0(Surv_var, " ~ ", year_var)
  }
  JP0  <- coxph(
    as.formula(form),
    data = data)
  if (n_jp == 0){
    return(JP0)
  }
  if (n_jp > 0){
    results <- data.frame(BIC = numeric(), k = integer())
    for (i in (min(data[[year_var]])+startbound):(max(data[[year_var]])-endbound)) {
      data$d1 <- pmax(data[[year_var]] - i, 0)
      mod <- coxph(as.formula(paste0(form, " + d1")),
                   data = data)
      results <- rbind(results, data.frame(BIC = BIC(mod), k = i))
    }
    sel <- results$k[which.min(results$BIC)]
    data$d1 <- pmax(data[[year_var]] - sel, 0)
    JP1 <- coxph(as.formula(paste0(form, " + d1")),
                 data = data)
    JP1$k <- sel
    
    if(n_jp == 1){
      if(best == F){
        print(paste0("Model with 1 change point: ", as.character(sel),"."))
        return(list(JP0, JP1))
      }
      if(BIC(JP0) < BIC(JP1) + bic_threshold){
        print("Model with 0 change points.")
        return(JP0)
      } else{
        print(paste0("Model with 1 change point: ", as.character(sel),"."))
        return(JP1)
      }
    }
    
    if (n_jp > 1){
      results2 <- data.frame(BIC = numeric(), k1 = integer(), k2 = integer())
      for (i in (min(data[[year_var]])+startbound):(max(data[[year_var]])-(endbound+inbound+1))) {
        for (j in (i+inbound):(max(data[[year_var]])-endbound)){
          data$d1 <- pmax(data[[year_var]] - i, 0)
          data$d2 <- pmax(data[[year_var]] - j, 0)
          mod <- coxph(as.formula(paste0(form, " + d1 + d2")),
                       data = data)
          results2 <- rbind(results2, data.frame(BIC = BIC(mod), k1 = i, k2 = j))
        }
      }
      sel2 <- results2[which.min(results2$BIC),]
      data$d1 <- pmax(data[[year_var]] - sel2$k1, 0)
      data$d2 <- pmax(data[[year_var]] - sel2$k2, 0)
      JP2 <- coxph(as.formula(paste0(form, " + d1 + d2")),
                   data = data)
      JP2$k <- c(sel2$k1, sel2$k2)
    }
    if(best == F){
      print(paste0("Model with 1 change point: ", as.character(sel),"."))
      print(paste0("Model with 2 change points: ", as.character(sel2$k1), " & ", as.character(sel2$k2), "."))
      return(list(JP0, JP1, JP2))
    }
    if(BIC(JP0) < BIC(JP1) + bic_threshold){
      print("Model with 0 change points.")
      return(JP0)
    } else{
      if(BIC(JP2) < min(BIC(JP0),BIC(JP1)) + bic_threshold){
        print(paste0("Model with 2 change points: ", as.character(sel2$k1), " & ", as.character(sel2$k2), "."))
        return(JP2)
      }else{
        print(paste0("Model with 1 change point: ", as.character(sel),"."))
        return(JP1)
      }
    }
  }
}


# function used within JPplot()
std_surv <- function(fit, weights, yrs, time = 5) {
  n_groups <- length(fit)
  n_times  <- length(yrs)
  surv_mat <- matrix(NA, nrow = n_times, ncol = n_groups)
  se_mat   <- matrix(NA, nrow = n_times, ncol = n_groups)
  for (i in seq_len(n_groups)) {
    s <- summary(fit[[i]], times = time, extend = TRUE)
    surv_mat[, i] <- as.numeric(s$surv)
    se_mat[, i]   <- as.numeric(s$std.err)
  }
  # Weighted survival
  S_std <- surv_mat %*% weights$weights
  # Variance (delta method): sum of (w^2 * se^2)
  Var_std <- se_mat^2 %*% (weights$weights^2)
  # 95% CI
  sd   <- sqrt(Var_std)
  lwr  <- S_std - 1.96 * sd
  upr  <- S_std + 1.96 * sd
  
  data.frame(
    yrs = yrs,
    surv = as.vector(S_std),
    sd   = as.vector(sd),
    lwr  = as.vector(lwr),
    upr  = as.vector(upr)
  )
}

# use the selected model to generate a plot and obtain AACS statistics
CoxJP_plot <- function(model, data, agegr_name, bootstrap = 50,
                   cohyrs, prdyrs, lastyear = T) {
  # get number of JP in the input model
  if ("d2" %in% names(coef(model))) {NJP = 2
  } else {NJP = ifelse("d1" %in% names(coef(model)), 1, 0)}
  # get weights from input case dataset (whole period)
  nagegr <- length(unique(data[[agegr_name]]))
  weights <- data.frame(table(data[[agegr_name]]) / nrow(data))
  colnames(weights) = c(agegr_name, "weights")

  ## Observed survival by year, by age group
  obs <- matrix(NA, nrow = length(cohyrs), ncol = nagegr)
  prd <- matrix(NA, nrow = length(prdyrs), ncol = nagegr)
  fit <- list()
  for (i in 1:nagegr){
    agegr <- sort(unique(data[[agegr_name]]))[i]
    tmp <- subset(data, get(agegr_name) == agegr)
    tmp_fmt <- data.frame(sex = tmp$sex, diagage = tmp$age, dm = tmp$MoD, 
                          dy = tmp$YoD, fm = tmp$MoFU, fy = tmp$YoFU,vitstat = tmp$status,
                          com3 = 3, Grupo3 = 1, SubGrupo = 11,
                          ExtSubGrupo = 111, agegr = tmp[[agegr_name]])
    obsfits <- with(subset(tmp, YoD %in% cohyrs), survfit(SurvObj ~ YoD))
    # if any year contains no cases, interpolate that year's survival as the mean of neighbouring years 
    if(!all(cohyrs %in% tmp$YoD)){
      presentes <- cohyrs[cohyrs %in% tmp$YoD]
      ausentes <- cohyrs[!(cohyrs %in% tmp$YoD)]
      for (yr in ausentes){
        tmp_obs <- summary(obsfits, times =  retrieve_time(5, obsfits$time), extend = T)$surv
        tmp_obs <- c(tmp_obs[1:which(cohyrs == yr)-1], mean(tmp_obs[(which(cohyrs == yr)-1):which(cohyrs == yr)]), tmp_obs[which(cohyrs == yr):length(tmp_obs)])
        obs[,i] <- tmp_obs
      }
    } else{
      obs[,i] <- summary(obsfits, times =  retrieve_time(5, obsfits$time), extend = T)$surv
    }
    for (j in prdyrs){
      jj <- j - min(prdyrs) + 1
      prd[jj,i] <- periodR::period(subset(tmp_fmt, agegr == agegr), 5, 
                                   surv.probs.males, surv.probs.females, perbeg = j, perend = j)$abs.surv[5]
    }
    
    ## Get CoxJP model fitted survival by age group
    # NJP = 0
    ifelse(NJP == 0, newdat <- data.frame(
      YoD = min(cohyrs, prdyrs):max(cohyrs, prdyrs)),
      # NJP = 1
      ifelse(NJP == 1, newdat <- data.frame(
        YoD = min(cohyrs, prdyrs):max(cohyrs, prdyrs),
        d1 = pmax(min(cohyrs, prdyrs):max(cohyrs, prdyrs) - model$k[1], 0)),
        # NJP = 2
        newdat <- data.frame(
          YoD = min(cohyrs, prdyrs):max(cohyrs, prdyrs),
          d1 = pmax(min(cohyrs, prdyrs):max(cohyrs, prdyrs) - model$k[1], 0),
          d2 = pmax(min(cohyrs, prdyrs):max(cohyrs, prdyrs) - model$k[2], 0))))
    
    if(nagegr > 1) newdat[[agegr_name]] <- factor(agegr, levels = levels(data[[agegr_name]]))
    fit[[i]] <- survfit(model, newdata = newdat)
  }
  colnames(obs) <- sort(unique(data[[agegr_name]]))
  colnames(prd) <- sort(unique(data[[agegr_name]]))
  
  ## Age-standardisation of observed survival by year
  obs_st <- c()
  for (yr in 1:length(cohyrs)) obs_st <- c(obs_st, sum(t(t(obs[yr,]) * weights$weights)))
  names(obs_st) <- cohyrs
  
  prd_st <- c()
  for (yr in 1:length(prdyrs)) prd_st <- c(prd_st, sum(t(t(prd[yr,]) * weights$weights)))
  names(prd_st) <- prdyrs
  
  ## Age-standardisation of fitted survival by year
  S_std <- std_surv(fit = fit, weights = weights, yrs = min(cohyrs, prdyrs):max(cohyrs, prdyrs), time = 5) 
  
  ## Bootstrap: the original model is fitted for resampled datasets to obtain a bootstrap distribution of AACS and its 95% CI
  B <- bootstrap
  if(NJP == 0) AACS <- matrix(NA, ncol = 1, nrow = B)
  if(NJP == 1) AACS <- matrix(NA, ncol = 2, nrow = B)
  if(NJP == 2) AACS <- matrix(NA, ncol = 3, nrow = B)
  
  bootstrap_mods <- list()
  if (bootstrap > 0){
    form <- paste0("SurvObj ~ YoD + strata(", agegr_name, ")")
    if(NJP == 1){
      data$d1 <- pmax(data$YoD - model$k[1], 0)
      form <- paste0("SurvObj ~ YoD + d1 + strata(", agegr_name, ")")
    }
    if(NJP == 2){
      data$d1 <- pmax(data$YoD - model$k[1], 0)
      data$d2 <- pmax(data$YoD - model$k[2], 0)
      form <- paste0("SurvObj ~ YoD + d1 + d2 + strata(", agegr_name, ")")
    }
      for(iter in 1:B){
        if (lastyear == F){ # in our study, last year is not used to fit models, due to lacking follow-up for 2022
          fitdata <- subset(data, YoD < max(prdyrs))
        } else fitdata <- data
        bootstrap_idx <- sample(nrow(fitdata),nrow(fitdata),TRUE)
        bootstrap_mods[[iter]] <- coxph(as.formula(form), fitdata[bootstrap_idx,])
        fit_tmp <- list()
        for (i in 1:nagegr){
          agegr <- sort(unique(data[[agegr_name]]))[i]
          if (nagegr > 1){
            newdat[[agegr_name]] <- factor(agegr, levels = levels(data[[agegr_name]]))
          }
          fit_tmp[[i]] <- survfit(bootstrap_mods[[iter]], newdata = newdat)
        }
        S_std_bs <- std_surv(fit = fit_tmp, weights = weights, yrs = min(cohyrs, prdyrs):max(cohyrs, prdyrs), time = 5) 
        if(NJP == 0){
          AACS[iter,1] <- mean(diff(S_std_bs$surv))*100
        }
        if(NJP == 1){
          AACS[iter,1] <- mean(diff(S_std_bs$surv[1 : which(S_std_bs$yrs==model$k[1])]))*100
          AACS[iter,2] <- mean(diff(S_std_bs$surv[which(S_std_bs$yrs==model$k[1]) : nrow(S_std_bs)]))*100
          
        }
        if(NJP == 2){
          AACS[iter,1] <- mean(diff(S_std_bs$surv[1 : which(S_std_bs$yrs==model$k[1])]))*100
          AACS[iter,2] <- mean(diff(S_std_bs$surv[which(S_std_bs$yrs==model$k[1]) : which(S_std_bs$yrs==model$k[2])]))*100
          AACS[iter,3] <- mean(diff(S_std_bs$surv[which(S_std_bs$yrs==model$k[2]) : nrow(S_std_bs)]))*100
        }
      }
  }
  
  obj <- list()
  obj[["obs_st"]] <- obs_st
  obj[["st_weights"]] <- weights
  obj[["prd_st"]] <- prd_st/100
  obj[["fit"]] <- fit
  obj[["S_std"]] <- S_std
  obj[["BIC"]] <- BIC(model)
  obj[["n"]] <- nrow(data)
  obj[["AACS"]] <- AACS

  plot(S_std$yrs, S_std$surv, type = "n",
       xlim = c(min(S_std$yrs), max(S_std$yrs)), ylim = c(0,1), xlab = "Year of diagnosis", 
       ylab = "5y cumulative survival")
  abline(v = 1:3000, col = "gray70", lty = 3)
  abline(h = seq(0,1,by=0.1), lty = 3, col = "gray70")
  lines(S_std$yrs, S_std$surv, type = "l",
        xlim = c(min(S_std$yrs), max(S_std$yrs)), lwd = 2)
  polygon(c(S_std$yrs, rev(S_std$yrs)),
          c(S_std$lwr, rev(S_std$upr)),
          col = adjustcolor("grey60", alpha.f = 0.4),
          border = NA)
  points(cohyrs, pch = 16, obs_st)
  points(prdyrs, pch = 21, bg = "grey70", prd_st/100)
  
  return(obj)  
  
}

#########################
#### EXAMPLE DATASET ####
#########################

set.seed(123)
sim_data <- data.frame(ID = 1:1000, YoD = sample(1000, x = 1999:2021, replace = T),
                       MoD = sample(1000, x = 1:12, replace = T),
                       age = floor(rbeta(1000, 2, 2.5)*15), sex = sample(1:2, 1000, replace = T))
sim_data$age_group <- cut(sim_data$age, breaks = c(-1, 4, 9, 14), labels = c("0-4", "5-9", "10-14"))
sim_data$hazard <- 0.0007 * exp((-0.03 + c(-0.02, 0 , 0.02)[sim_data$age_group] - 0.08 * (sim_data$YoD > 2009)) 
                                * (sim_data$YoD - 1999)) # simulating trend over diagnosis years + additive effects in hazards by age group + additive effect after 2009
sim_data$surv_time <- rexp(nrow(sim_data), rate = sim_data$hazard)
sim_data$surv_time <- pmin(sim_data$surv_time, 1826)
sim_data$status <- ifelse(sim_data$surv_time >= 1826, 1, 2)
sim_data$dx_dates <- mdy(gsub(" ", "0", 
                              apply(sim_data, 1, function(x) paste(x[which(names(x) == "MoD")],
                                                                   15,
                                                                   x[which(names(x) == "YoD")], 
                                                                   sep = "-"))))
sim_data$fup_dates <- sim_data$dx_dates + sim_data$surv_time
sim_data$MoFU <- month(sim_data$fup_dates)
sim_data$YoFU <- year(sim_data$fup_dates)
sim_data$SurvObj <- Surv(sim_data$surv_time / 365.25, sim_data$status) # survival times have to be supplied in years


mod <- CoxJP(sim_data, Surv_var = "SurvObj", "YoD", "age_group", n_jp = 2,
                        startbound = 3, endbound = 5, inbound = 3, bic_threshold = 0, best = F)
BIC(mod[[1]]); BIC(mod[[2]]); BIC(mod[[3]])
result <- CoxJP_plot(mod[[3]], sim_data, "age_group", cohyrs = 1999:2018, prdyrs = 2019:2021, bootstrap = 100)
quantile(result$AACS[,1], probs = c(0.025, 0.5, 0.975))
quantile(result$AACS[,2], probs = c(0.025, 0.5, 0.975))
quantile(result$AACS[,3], probs = c(0.025, 0.5, 0.975))







