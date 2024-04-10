########## FUNCTION : Gehan  #######################
##' Perform Gehan test.
#'
#' @param site A dataframe for site measurements. Includes one column of measurements and one column of detect/nondetect designation. 
#' @param background A dataframe for site measurements. Includes one column of measurements and one column of detect/nondetect designation.
#' @param nd.s Indicator variable for \code{site} measurements that determines if the corresponding measurement (row) is a nondetect. A value of \code{1} indicates that the row is a nondetect (\code{0} is a detected measurement). Default is a column of site dataframe named nondetect (\code{site$nondetect}).
#' @param background A dataframe for background measurements.
#' @param measure.b Variable for \code{background} measurements that provides the corresponding concentration. Default is a column of background dataframe named measurement (\code{background$measurement}).
#' @param nd.b Indicator variable for \code{background} measurements that determines if the corresponding measurement (row) is a nondetect. A value of \code{1} indicates that the row is a nondetect (\code{0} is a detected measurement). Default is a column of background dataframe named nondetect (\code{background$nondetect}).
#' @param epsilon  A proportion of the site that has concentrations greater than the background. Refer to Table 4.3 of referenced literature for possible values.
#' @param alpha Type I error rate. Options are \code{0.05} (default), \code{0.01}, or \code{0.10}.
#' @param power Statistical power. Options are \code{0.75} (default), \code{0.90}, or \code{0.95}.
#' @param plot Logical. Will display a histogram of the site and background measurements with other details. (default is \code{FALSE})
#' @param print Logical. Will display Gehan test results. (default is \code{TRUE})
#' @return A list containing two objects:
#' \itemize{
#'   \item \code{List of parameters and resulting statistics}
#'      \item \code{alpha} - Type I error rate. 
#'      \item \code{power} - Statistical power to detect differences. 
#'      \item \code{deltaS} - Magnitude of the difference in median site and background concentrations.
#'      \item \code{G} - Gehan test statistic (analogous to a Z statistic). 
#'      \item \code{p_val} - P-value of the test, computed from a standard normal distribution. 
#'    \item \code{Data frame used to calculate ranks and scores of each observation}
#'      \item \code{Samp Measurement} - Measured chemcical concentration of each sample.
#'      \item \code{Detect (Y/N)} - Inicator measurement detection status where  detect = 0 and nondetect = 1.
#'      \item \code{Samp Location} - Indicator of location sampled where background = 0 and site = 1.
#'      \item \code{d} - Rank determining variable 1. 
#'      \item \code{e} - Rank determining variable 2. 
#'      \item \code{R} - Calculated rank. 
#'   \item \code{aR} - Calculated rank score. 
#' @return A list of the test parameters and statistics: 
#'  
#' 
#' }
#' @details The Gehan test is a nonparametric test that is useful for detecting differences in the median chemical concentrations of 
#' site and background measurements. It is functionally is similar to the more commonly used Wilcoxon-Rank Sum (WRS) test, but it is 
#' more powerful in situations where sample size is small, there is a significant number of nondetect measurements, and the measurements 
#' have different reporting limits. The test evaluates the difference in median site and background chemical concentration by iteratively 
#' ranking the data by size and detection status. The resulting ranks are used to assign each measurement a score that is used to calculate
#' the Gehan statistic (G). G follows a standard normal distribution, so a p-value can be computed to assess significance of the test results.
#' The \code{gehan()} function automatically assesses statistical assumptions of the test based on the user's desired power,
#' significance level, and sample size, but additional details regarding implementation can be found in 
#' "Guidance for Environmental Background Analysis Volume III: Groundwater" (2003). 
#' 
#' The function allows users to store the test results as an object. The output is a list containing a list of parameters used in the test, 
#' the resulting statistics, and the data frame used to calculate the ranks and scores of each observation. 
#'
#' @examples
#' site_ex1 <- data.frame(samples=c(2, 4, 8, 17, 20, 25, 34, 35, 40, 43), nondetect = c(0, 1, rep(0, 5), 1, rep(0,2)))
#' bck_ex1 <- data.frame(samples=c(1, 4, 5, 7, 12, 15, 18, 21, 25, 27), nondetect = c(0, 1, 0, 0, 1, 0, 0, 1, 1, 0))

#' G_result <- gehan(alpha = 0.05, deltaS = 2.0, power = 0.9, site = site_ex1, background = bck_ex1) 
#' 
#' @export
#' @references Naval Facilities Engineering Command. (2003, October).  \emph{Guidance for Environmental Background Analysis Volume III: Groundwater.} https://vsp.pnnl.gov/docs/Draft_Guidance_for_Review.pdf.
#' 
gehan <- function(site, background, alpha, deltaS, power, print=TRUE, plot=TRUE){
  
  ### Verify test assumptions ####
  
  # code in table 4.6 
  
  table_4.6 <- data.frame(deltaS.choices = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.25, 2.5, 2.75, 4.0, 4.5, 5.0), 
                          power0.95_alpha0.01 = c(3972, 998, 448, 255, 166, 117, 88, 69, 58, 47, 40, 35, 31,28, 25, 23, 22, 20, 19, 18, 16, 15, 15, 14, 13, 13),
                          power0.9_alpha0.01 = c(3278, 824, 370, 211, 137, 97, 73, 57, 47, 39, 33, 29, 26, 23, 21, 19, 18, 17, 15, 15, 14, 13, 12, 12, 11, 11),
                          power0.75_alpha0.01 = c(2268, 570, 256, 148, 95, 67, 51, 40, 32, 27, 23, 20, 18, 16, 15, 14, 13, 12, 11, 11, 10, 9, 9, 8, 8, 8),
                          power0.95_alpha0.05 = c(2726, 685, 307, 175, 114, 81, 61, 48, 39, 32, 28, 24, 22, 19, 18, 16, 15, 14, 13, 13, 11, 11, 10, 10, 9, 9),
                          power0.9_alpha0.05 = c(2157, 542, 243, 139, 90, 64, 48, 38, 31, 26, 22, 19, 17, 15, 14, 13, 12, 11, 11, 10, 9, 9, 8, 8, 8, 7),
                          power0.75_alpha0.05 = c(1355, 341, 153, 87, 57, 40, 30, 24, 20, 16, 14, 12, 11, 10, 9, 8, 8, 7, 7, 7, 6, 6, rep(5, 4)),
                          power0.95_alpha0.10 = c(2157, 542, 243, 139, 90, 64, 48, 38, 31, 26, 22, 19, 17, 15, 14, 13, 12, 11, 11, 10, 9, 9, 8, 8, 8, 7),
                          power0.90_alpha0.10 = c(1655, 416, 187, 106, 69, 19, 37, 29, 24, 20, 17, 15, 13, 12, 11, 10, 9, 9, 8, 8, 7, 7, rep(6, 4)),
                          power0.75_alpha0.10 = c(964, 243, 109, 62, 41, 29, 22, 17, 14, 12, 10, 9, 8, 7, 7, 6, 6, 5, 5, 5, rep(4, 6))) 
  
  
  samples.required <- function(deltaS, power, alpha, table_4.6){
    #filter table based off epsilon, power, and alpha
    table_chosen_deltaS <- table_4.6[table_4.6$deltaS.choices == deltaS,]
    # Check if any rows match the filter criteria
    if (nrow(table_chosen_deltaS) > 0) {
      # power = 0.95
      if (power == 0.95) {
        if(alpha == 0.01){
          return(table_chosen_deltaS$power0.95_alpha0.01) 
        }
        else if(alpha == 0.05){
          return(table_chosen_deltaS$power0.95_alpha0.05)
        }
        else if(alpha == 0.10){
          return(table_chosen_deltaS$power0.95_alpha0.10)
        }
      } 
      # power = 0.90
      else if (power == 0.90){
        if(alpha == 0.01){
          return(table_chosen_deltaS$power0.9_alpha0.01) 
        }
        else if(alpha == 0.05){
          return(table_chosen_deltaS$power0.9_alpha0.05)
        }
        else if(alpha == 0.10){
          return(table_chosen_deltaS$power0.9_alpha0.10)
        }       
      } 
      # power = 0.75
      else if (power == 0.75){
        if(alpha == 0.01){
          return(table_chosen_deltaS$power0.75_alpha0.01) 
        }
        else if(alpha == 0.05){
          return(table_chosen_deltaS$power0.75_alpha0.05)
        }
        else if(alpha == 0.10){
          return(table_chosen_deltaS$power0.75_alpha0.10)
        }       
      }
      else{
        stop("invalid alpha") # add condition for invalid power 
      }
    }
  } # samples.required end 
  
  nreq <- samples.required(deltaS = deltaS, 
                           power = power, 
                           alpha = alpha, 
                           table_4.6 = table_4.6) 
  print(nreq) 
  
  condition.check <- if (nrow(site)>=nreq*0.5 & nrow(background)>=nreq*0.5 & nrow(site) == nrow(background)) {
    print("minimum sample sizes for background and site is satisfied")
    print(paste0("required rows = ", nreq))
    print(paste0("rows of site = ", nrow(site)))
    print(paste0("rows of background = ", nrow(background))) 
  }
  else if(nrow(site) != nrow(background)){
    warning("Number of site measurements should be equal to the number of background measurements")
  }
  else{
    warning("sample size requirements not met for desired power and significance level")
    print(paste0("required rows = ", nreq))
    print(paste0("rows of site = ", nrow(site)))
    print(paste0("rows of background = ", nrow(background)))
    return(condition.check)
  } # condition check end 
  
  
  ##### Perform Gehan test #####
  
  # 0) assign binary variable to determine site ID (0 = background, 1 = site)
  # non-detects should already be in input df (detect = 0, non-detect = 1)
  
  site$h <- 1
  background$h <- 0
  
  ### 1) combine lists m and n, keep ID of site/background (P/B)
  
  combined_meas <- rbind(site, background)
  # print(combined_meas)
  
  ### 2) order list smallest to largest 
  sorted_meas <- combined_meas[order(combined_meas$samples, combined_meas$h),]
  print("Data Frame - Combined and Ranked")
  print(sorted_meas)
  
  ### 3) Determine values of D and E 
  
  d <- 0
  d_vec <- c()
  d_rank <- c()
  
  e <- 0 
  e_vec <- c()
  e_rank <- c()
  
  for (i in  1:length(sorted_meas$nondetect)){
    if (sorted_meas$nondetect[i] == 0){ # if meas is detect 
      d_vec <- c(d_vec, 1)
      d = d + 1
      d_rank = c(d_rank, d)
      e_vec <- c(e_vec, 0)
      e = e 
      e_rank <- c(e_rank, e)
    } else { # if meas is non-detect 
      d_vec <- c(d_vec, 0)
      d = d 
      d_rank = c(d_rank, d)
      e_vec <- c(e_vec, 1)
      e = e + 1
      e_rank <- c(e_rank, e)
    }
  }
  # print(d_vec)
  print("D scores")
  print(d_rank)
  # print(e_vec) 
  print("E scores")
  print(e_rank)
  
  sorted_meas$d = d_rank
  sorted_meas$e = e_rank
  
  ### 4) Calculate R for each observation 
  # Find T (total number of non-detects) 
  Tot = length(which(sorted_meas$nondetect == 1))
  # print(Tot)
  # I = 0 (detect): R = d + (T + e)/2 
  # I = 1 (nondetect): R = (T + 1 + d)/2
  R <- 0 
  R_vec <- c()
  for (i in  1:length(sorted_meas$nondetect)){
    if (sorted_meas$nondetect[i] == 0){ 
      R = sorted_meas$d[i] + ((Tot + sorted_meas$e[i])/2)
      R_vec <- c(R_vec, R)
    } else { 
      R = (1 + sorted_meas$d[i] + Tot)/2
      R_vec <- c(R_vec, R)
    }
  }
  print("Rank")
  print(R_vec) 
  
  sorted_meas$R = R_vec
  
  ### 5) Calculate scores a(Ri)
  # a(Ri) = 2Ri - N - 1 
  N = nrow(sorted_meas)
  aR <- 0 
  aR_vec <- c()
  
  for (i in  1:length(sorted_meas$nondetect)){
    aR = 2*sorted_meas$R[i] - N - 1
    aR_vec <- c(aR_vec, aR)
  }
  print("Score")
  print(aR_vec) # all correct values what a moment 
  
  sorted_meas$aR = aR_vec 
  
  ### 6) Calculate Gehan statistic 
  
  m = nrow(bck_ex1)
  n = nrow(site_ex1)
  
  sum_site = with(sorted_meas, sum(aR[h == 1]))
  
  aR_squared = as.numeric(lapply(sorted_meas$aR, function(x) x^2))
  sum_all_squared = sum(aR_squared)
  
  G = sum_site/sqrt((m*n*sum_all_squared)/(N*(N-1)))
  # print(G)
  
  ### 7) Compare to normal distribution 
  p_val <- pnorm(G, lower.tail = FALSE)
  
  ### 8) Print conclusions 
  
  if (print) {
    text.c <- "Gehan Test Results:"
    cat(paste("\n", text.c, "\n", sep = ""))
    row <- paste(rep("=", 87), collapse = "")
    cat(row, "\n")
    cat(
      paste("SAMPLE SIZE:", "\t",
            "Required (power):"),
      paste0(nreq, " (", power, ")"),
      paste(
        "\t" ,
        "Site (n):",
        nrow(site),
        "\t",
        "Background (m):",
        nrow(background),
        "\n"
      )
    )
    if (p_val < alpha) {
      cat(row, "\n")
      cat(
        paste(
          "Gehan Statistic:",
          round(G, 2),
          "\n",
          "\n", 
          "p-value:", 
          round(p_val, 3), 
          "\t"
        )
      )
      cat(paste("DECISION:", "\t","Reject null hypothesis", "\n"))
      cat(row, "\n")
    } else {
      cat(row, "\n")
      cat(
        paste(
          "Gehan Statistic:",
          round(G, 2), 
          "\n",
          "\n", 
          "p-value:", 
          round(p_val, 3), 
          "\t"
        )
      )
      cat(paste("DECISION:", "\t","Fail to reject null hypothesis", "\n"))
      cat(row, "\n")
    }
  }
  Gehan_results <- as.data.frame(sorted_meas)
  colnames(Gehan_results) = c("Samp Measurement", "Detect (Y/N)", "Samp Location", 
                              "d", "e", "R", "aR")
  return(Gehan_results)
} # Gehan function end 

