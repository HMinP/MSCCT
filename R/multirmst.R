#' Test of RMST for comparing two or more survival curves
#'
#' Performs the test of Restricted Mean Survival Time for two or more survival curves,
#' i.e. compares the difference of areas under survival curves.
#'
#' For one group, the Restricted Mean Survival Time at time tau (RMST(tau)) is the area
#' under the survival curve between time 0 and tau. The test of RMST compares the
#' RMST(tau) from both groups and tests whether the difference is zero or not. If the
#' difference is not null, then the survival curves cannot be equal.
#'
#' For exactly two groups, a single test is performed. For more than two survival curves,
#' it compares each survival curve to every other curves and tests the global null
#' hypothesis "all curves are equal" against the hypothesis "the curves are not all equal".
#'
#' @param df A dataframe with columns :
#'   * `time` : positive numbers, time-to-event;
#'   * `status` : integer of factor. 0 is (right) censoring, 1 is event;
#'   * `arm` : integers from 0 to n-1 or factor with at least 2 levels.
#'     The group the patient belongs to.
#' @param tau The truncation time, default is the lowest max(time) of each groups;
#' @param nboot Number of bootstrap samples;
#' @param method The correction used for the p-values. Must be in [p.adjust.methods]. Unused for exactly two groups.
#'
#' @return an object of class `multirmst` containing :
#'   * `results` A matrix. Each row represents a comparison of two curves and contains the difference
#'     of RMST, its standard deviation, the p-value and the adjusted p-value;
#'   * `p` The p-value of the global test;
#'   * `nb_tests` : The number of performed tests;
#'   * The parameters `tau`, `method` and `nboot`.
#'
#' @references Royston, P., & Parmar, M. K. (2013). Restricted mean survival time:
#'     an alternative to the hazard ratio for the design and analysis of randomized
#'     trials with a time-to-event outcome. BMC medical research methodology, 13, 1-15.
#'
#' @examples
#' multirmst(data_under_PH, 36)
#' multirmst(data_not_PH, tau=36, method="BH", nboot=1000)
#'
#' @export
multirmst = function(df, tau=-1, nboot=500, method="bonferroni"){
  if (!all(c("time", "status", "arm") %in% colnames(df))){
    stop("The dataframe must contain the columns 'time', 'status' and 'arm'.")
  }
  
  df$status = as.numeric(df$status)
  df$status = df$status - min(df$status)
  if (!all(df$status %in% c(0,1))){stop("'status' must be either 0 or 1.")}
  
  df$arm = as.numeric(df$arm)
  df$arm = df$arm - min(df$arm)
  nb_arms = length(unique(df$arm))
  if (nb_arms < 2){stop("Need at least two groups.")}
  if (!setequal(df$arm,0:(nb_arms-1))){
    stop(paste("Incorrect value for 'arm', must range from 0 to ",nb_arms-1, ".", sep=""))
  }
  
  if (tau == -1){tau = min(tapply(X=df$time, INDEX=df$arm, FUN=max))}
  
  nb_tests = nb_arms*(nb_arms-1) / 2
  label_test = rep(NA,nb_tests)
  results = matrix(NA,nrow=nb_tests,ncol=4)
  k=1
  
  for (i in 0:(nb_arms-2)){
    for (j in (i+1):(nb_arms-1)){
      label_test[k] = paste(i,"VS",j)
      
      ind = (df$arm == i) | (df$arm == j)
      df_ij = df[ind,]
      df_ij$arm = (df_ij$arm - i) / (j-i)
      
      X = boot(df_ij, rmstdiff, R=nboot, tau=tau)
      results[k,1] = X$t0
      results[k,2] = sd(X$t)
      results[k,3] = 2*(1-stats::pnorm(abs(results[k,1]/results[k,2])))
      
      k=k+1
    }
  }
  
  if (nb_tests == 1){
    results = results[, -4, drop=FALSE]
    colnames(results) = c("dRMST","sd","p-value")
    p = results[1,3]
  }
  else {
    p_adjusted = p.adjust(results[,3], method=method)
    results[,4] = p_adjusted
    p = min(p_adjusted)
    rownames(results) = label_test
    colnames(results) = c("dRMST","sd","p","p adjusted")
  }
  
  z = list(results=results, tau=tau, p=p, nb_tests=nb_tests,
           method=method, nboot=nboot)
  class(z) = "multirmst"
  return(z)
}




#' Print method for the multiple test of RMST
#'
#' @param x An object of class `multirmst` as returned by [multirmst()];
#' @param ... For compatibility with the `print` method, unused and to be ignored.
#'
#' @export
print.multirmst = function(x, ...){
  nb_tests=x$nb_tests
  
  cat("(Multiple) test of RMST \n")
  cat("Truncation time :",x$tau," \n")
  if (nb_tests > 1) {cat("Correction :",x$method,"\n\n")}
  print(x$results)
  cat("",end="\n")
  if (nb_tests > 1) {cat("p=",x$p,sep="")}
}




#' Difference of RMST between two groups
#'
#' Computes the difference of RMST between two treatment groups, ie the difference
#' of area under the survival curves.
#'
#' @param df A dataframe with columns `time`, `status` and `arm`
#' @param ind A vector of integer for shuffling the group of patients. Used by boot
#'   to estimate the standard deviation under null hypothesis. Default is non-shuffling.
#' @param tau Truncation time
#'
#' @return The value of the difference of RMST
#' @noRd
rmstdiff = function(df, ind=numeric(), tau=-1){
  time = df$time
  status = df$status
  if (length(ind) == 0){ind = 1:length(time)}
  arm = df$arm[ind] # shuffling the arm variable
  
  time0 = time[arm == 0]
  status0 = status[arm == 0]
  time1 = time[arm == 1]
  status1 = status[arm == 1]
  
  KM1 = survival::survfit(survival::Surv(time1,status1)~1)
  KM0 = survival::survfit(survival::Surv(time0,status0)~1)
  S1 = stats::stepfun(KM1$time, c(1,KM1$surv), right=FALSE)
  S0 = stats::stepfun(KM0$time, c(1,KM0$surv), right=FALSE)
  
  time_evt = unique(sort(c(0, KM0$time[KM0$time < tau], KM1$time[KM1$time < tau], tau)))
  time_int = time_evt[-1]-time_evt[-length(time_evt)]
  
  surv_diff = S1(time_evt) - S0(time_evt)
  surv_diff = surv_diff[-length(surv_diff)]
  drmst = sum(surv_diff*time_int)
  return(drmst)
}

