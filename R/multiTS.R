#' Two-staged test for comparison of two or more survival curves.
#'
#' Performs a Two-Stage test for each pair of survival curves and apply a correction
#' in case of several comparisons.
#'
#' The first stage is a log-rank test. If the first test is significant, then the whole
#' procedure stops and we conclude that the survival curves are different. If it is not
#' significant, then the survival curves are either equal or crossing each other and
#' the log-rank test can't conclude the difference. A second test is performed to distinguish
#' these two cases.
#'
#' @param df A dataframe with columns :
#'   * `time` : positive numbers, time-to-event;
#'   * `status` : integer of factor. 0 is (right) censoring, 1 is event;
#'   * `arm` : integers from 0 to n-1 or factor with at least 2 levels.
#'     The group the patient belongs to.
#' @param eps A number from 0 to 0.5. See reference for interpretation;
#' @param nboot A positive integer, number of bootstrap sample for the second stage;
#' @param method The correction used for the p-values. Must be in [p.adjust.methods]. Unused for exactly two groups;
#'
#' @return An object of class `multiTS` containing:
#'   * `results` : A matrix. Each row represents a comparison of two curves and contains
#'     the p-values for both stage, the global p-value and the adjusted global p-value;
#'   * `p` : The global p-value for the global test;
#'   * `nb_tests` : The number of performed Two-Stage tests;
#'   * the parameters `eps`, `method` and `nboot`.
#'
#' @references
#'    * Qiu, P., & Sheng, J. (2008). A two-stage procedure for comparing
#'      hazard rate functions. Journal of the Royal Statistical Society Series
#'      B: Statistical Methodology, 70(1), 191-208. Chen, Zhongxue & Huang, Hanwen & Qiu, Peihua. (2017).
#'    * An improved two-stage procedure to compare hazard curves. Journal of Statistical
#'      Computation and Simulation. 87. 1-10. 10.1080/00949655.2017.1292276.
#'
#' @export
#'
#' @examples
#' multiTS(data_not_PH, eps=0.05, nboot=200, method="BH")
multiTS = function(df, eps=0.1, nboot=500, method="bonferroni"){
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
  
  
  nb_tests = nb_arms * (nb_arms-1) / 2
  label_test = rep(NA,nb_tests)
  results = matrix(NA,nrow=nb_tests,ncol=4)
  k=1
  
  for (i in 0:(nb_arms-2)){
    for (j in (i+1):(nb_arms-1)){
      label_test[k] = paste(i,"VS",j)
      
      ind = (df$arm == i) | (df$arm == j)
      df_ij = df[ind,]
      df_ij$arm = (df_ij$arm - i) / (j-i)
      
      test1 = multiLR(df_ij)
      U = test1$U
      p1 = test1$p
      V = boot::boot(df_ij, calcul_V, nboot, eps=eps)
      p2 = mean(V$t > V$t0)
      p = 1 - stats::pchisq(-2*log(p1*p2 + 1e-16), 4)
      results[k,1] = p1
      results[k,2] = p2
      results[k,3] = p
      
      k=k+1
    }
  }
  
  if (nb_tests == 1){
    results = as.vector(results[,-4])
    results = setNames(results, c("p1","p2","p"))
    p = unname(results)[3]
  }
  
  else {
    p_adjusted = p.adjust(results[,3], method=method)
    results[,4] = p_adjusted
    p = min(p_adjusted)
    rownames(results) = label_test
    colnames(results) = c("p1","p2","p","adj_p")
  }
  
  
  z = list(results=results, p=p, nb_tests=nb_tests, eps=eps, method=method, nboot=nboot)
  class(z) = "multiTS"
  return(z)
}


#' Print method for the multiple Two-Stage test
#'
#' @param x An object of class `multiTS` as returned by [multiTS()];
#' @param ... For compatibility with the `print` method, unused and to be ignored.
#'
#' @export
print.multiTS = function(x, ...){
  nb_tests = x$nb_tests
  
  cat("(Multiple) Two-Staged test \n")
  if (nb_tests > 1) {cat("Correction :",x$method," \n\n")}
  print(x$results)
  cat("",end="\n")
  if (nb_tests > 1) {cat("p=",x$p,sep="")}
}




poids_S2 = function(df){
  time = df$time
  status = df$status
  arm = df$arm
  
  evt_time = unique(time[status == 1])
  evt_time_ordered = sort(evt_time)
  D = length(evt_time_ordered)
  
  time0 = time[arm == 0]
  status0 = status[arm == 0]
  n0 = length(time0)
  time1 = time[arm == 1]
  status1 = status[arm == 1]
  n1 = length(time1)
  n = n0+n1
  
  KM_censor0 = survival::survfit(survival::Surv(time0,1-status0)~1)
  L0 = stats::stepfun(KM_censor0$time, c(1,KM_censor0$surv), right=FALSE)
  KM_censor1 = survival::survfit(survival::Surv(time1,1-status1)~1)
  L1 = stats::stepfun(KM_censor1$time, c(1,KM_censor1$surv), right=FALSE)
  
  
  KM_global = survfit(Surv(time,status)~1)
  missing = setdiff(evt_time_ordered, KM_global$time)
  missing_ind = evt_time_ordered %in% missing
  delta_S = rep(0,length(evt_time_ordered))
  diff_S = function(t){
    # variation of S for a non-missing observation
    i = which(KM_global$time == t)
    dS = ifelse(i == 1, 1-KM_global$surv[i],
                KM_global$surv[i-1]-KM_global$surv[i])
    return(dS)
  }
  delta_S[!(missing_ind)] = as.numeric(sapply(evt_time_ordered[!missing_ind], diff_S))
  
  
  alpha = (L0(evt_time_ordered)*L1(evt_time_ordered)*delta_S)
  alpha = alpha / (n0/n*L0(evt_time_ordered) + n1/n*L1(evt_time_ordered))
  sums = cumsum(alpha)
  c = sums / (sums[D]-sums)
  mat_c = matrix(c, nrow=D, ncol=D, byrow=TRUE)
  
  weights = ifelse(lower.tri(mat_c), mat_c, -1)
  
  return(weights)
}



calcul_V = function(df, ind=numeric(), eps=0.1){
  if (length(ind) == 0) {ind = 1:(nrow(df))}
  df_shuffled = df
  df_shuffled$arm = df$arm[ind]
  
  D = length(unique(df$time[df$status == 1]))
  De = floor(D*eps)
  weights = poids_S2(df_shuffled)
  
  V_vect = multiLR(df_shuffled, weights=weights[,De:(D-De)])$U
  V = max(as.numeric(V_vect))
  
  return(V)
}

