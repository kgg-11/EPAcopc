########## FUNCTION : Gehan  #######################
##' Perform Gehan test.
#'
#' @param site A dataframe for site measurements. Includes one column of measurements and one column of detect/nondetect designation. 
#' @param background A dataframe for background measurements. Includes one column of measurements and one column of detect/nondetect designation..
#' @param deltaS  Magnitude of the difference in median site and background concentrations.
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
#' }
#' 
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
#' Actual Data (ex 4-12 from EPA)
#' site_ex <- data.frame(samples=c(2, 4, 8, 17, 20, 25, 34, 35, 40, 43), nondetect = c(0, 1, rep(0, 5), 1, rep(0,2)))
#' bck_ex <- data.frame(samples=c(1, 4, 5, 7, 12, 15, 18, 21, 25, 27), nondetect = c(0, 1, 0, 0, 1, 0, 0, 1, 1, 0))
#' G_result <- gehan(alpha = 0.05, deltaS = 2.0, power = 0.9, 
#'                   site = site_ex, measure.s = site_ex$samples, nd.s = site_ex$nondetect, 
#'                   background = bck_ex, measure.b = bck_ex$sample, nd.b = bck_ex$nondetect) 
#' 
#' Randomized Data 
#' site_ex_random <- data.frame(samples = c(floor(runif(15, 1, 50))), nondetect = c(rep(0, 12), rep(1, 3)))
#' bck_ex_random <- data.frame(samples = c(floor(runif(15, 1, 40))), nondetect = c(rep(0, 12), rep(1, 3)))
#' G_result_random <- gehan(alpha = 0.05, deltaS = 2.0, power = 0.9, 
#'                    site = site_ex_random, measure.s = site_ex_random$samples, nd.s = site_ex_random$nondetect,
#'                    background = bck_ex_random, measure.b = bck_ex_random$sample, nd.b = bck_ex_random$nondetect)
#' 
#' Uneven Number of Site and Background Measurements - Not Recommended by EPA 
#' site_ex_uneven <- data.frame(samples = c(floor(runif(15, 1, 50))), nondetect = c(rep(0, 12), rep(1, 3)))
#' bck_ex_uneven <- data.frame(samples = c(floor(runif(20, 1, 50))), nondetect = c(rep(0, 17), rep(1, 3)))
#' G_result_uneven <- gehan(alpha = 0.05, deltaS = 2.0, power = 0.9,
#'                         site = site_ex_uneven, measure.s = site_ex_uneven$samples, nd.s = site_ex_uneven$nondetect, 
#'                         background = bck_ex_uneven, measure.b = bck_ex_uneven$sample, nd.b = bck_ex_uneven$nondetect)
#' 
#' Flex naming 
#' site_ex_flex<- data.frame(meas = c(floor(runif(15, 1, 50))), nd = c(rep(0, 12), rep(1, 3)))
#' bck_ex_flex <- data.frame(blah = c(floor(runif(15, 1, 40))), scoop = c(rep(0, 12), rep(1, 3)))
#' G_result_flex <- gehan(alpha = 0.05, deltaS = 2.0, power = 0.9, 
#'                        site = site_ex_flex, measure.s = site_ex_flex$meas, nd.s = site_ex_flex$nd, 
#'                        background = bck_ex_flex, measure.b = bck_ex_flex$blah, nd.b = bck_ex_flex$scoop)
#' @export
#' @references Naval Facilities Engineering Command. (2003, October).  \emph{Guidance for Environmental Background Analysis Volume III: Groundwater.} https://vsp.pnnl.gov/docs/Draft_Guidance_for_Review.pdf.
#' 
gehan <- function(site, 
                  measure.s, 
                  nd.s, 
                  background, 
                  measure.b, 
                  nd.b, 
                  alpha, 
                  deltaS, 
                  power, 
                  print=TRUE, 
                  plot=TRUE
                  )
  
  {
  
  ### Verify test assumptions ####
  
  # code in table 4.6 
  
  site_df <- data.frame(measure.s, nd.s)
  colnames(site_df) = c("samples", "nondetect") 

  background_df <- data.frame(measure.b, nd.b)
  colnames(background_df) = c("samples", "nondetect") 
  
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
  
  condition.check <- if (nrow(site_df)>=nreq*0.5 & 
                         nrow(background_df)>=nreq*0.5 & 
                         nrow(site_df) == nrow(background_df) &
                         length(which(site_df$nondetect == 1)) <= 0.4*nrow(site_df) &
                         length(which(background_df$nondetect == 1)) <= 0.4*nrow(background_df)) {
    print("minimum sample sizes for background and site is satisfied")
    print(paste0("required rows = ", nreq))
    print(paste0("rows of site = ", nrow(site_df)))
    print(paste0("rows of background = ", nrow(background_df))) 
  }
  else if(nrow(site_df) != nrow(background_df)){
    warning("Number of site measurements should be equal to the number of background measurements, 
            refer to draft guidance for sample replication guidelines")
  }
  else if(length(which(site_df$nondetect == 1)) > 0.4*nrow(site_df)){
    stop("More than 40% of site measurements are nondetects")
  }
  else if(length(which(background_df$nondetect == 1)) > 0.4*nrow(background_df)){
    stop("More than 40% of site measurements are nondetects")
  }
  else{
    warning("sample size requirements not met for desired power and significance level")
    print(paste0("required rows = ", nreq))
    print(paste0("rows of site = ", nrow(site_df)))
    print(paste0("rows of background = ", nrow(background_df)))
  } # condition check end 
  
  
  ##### Perform Gehan test #####
  
  # 0) assign binary variable to determine site ID (0 = background, 1 = site)
  # non-detects should already be in input df (detect = 0, non-detect = 1)
  
  site_df$h <- 1
  background_df$h <- 0
  
  ### 1) combine lists m and n, keep ID of site/background (P/B)
  
  combined_meas <- rbind(site_df, background_df) 

  
  ### 2) order list smallest to largest 
  sorted_meas <- combined_meas[order(combined_meas$nondetect, combined_meas$h),]
  print("Data Frame - Combined and Sorted")
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
  
  sorted_meas$aR = aR_vec 
  
  ### 6) Calculate Gehan statistic 
  
  m = nrow(background_df)
  n = nrow(site_df)
  
  sum_site = with(sorted_meas, sum(aR[h == 1]))
  
  aR_squared = as.numeric(lapply(sorted_meas$aR, function(x) x^2))
  sum_all_squared = sum(aR_squared)
  
  G = sum_site/sqrt((m*n*sum_all_squared)/(N*(N-1)))
  # print(G)
  
  ### 7) Compare to normal distribution 
  p_val <- pnorm(G, lower.tail = FALSE)
  
  
  ### 8) Plots 
  if(plot == TRUE) { 
    binw3 <- round(max(range(site_df$samples),range(background_df$samples))/5)
    
    my_histogram <- ggplot2::ggplot() + 
      ggplot2::geom_histogram(data = background_df, 
                              ggplot2::aes(x = samples, col = I("black"), fill = "b"), 
                              alpha = 0.8, bins = binw3) +
      ggplot2::geom_histogram(data = site_df, 
                              ggplot2::aes(x = samples, col = I("black"), fill = "r"), 
                              alpha = 0.8, bins = binw3) +
      ggplot2::scale_fill_manual(name ="Measurement ID", 
                                 values = c("r" = "#1F968BFF", "b" = "#440154FF"), 
                                 labels=c("b" = "background", "r" = "site")) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = median(site_df$samples)),
                          col='#00CFAC', size=1, linetype = "dashed") + 
      ggplot2::geom_vline(ggplot2::aes(xintercept = median(background_df$samples)), 
                          col='#BD00FF', size=1, linetype = "dashed") + 
      ggplot2::theme_classic() +
      ggplot2::xlab("COPC Concentration") + 
      ggplot2::ylab("Number of Measurements")
    
    print(my_histogram)
  } 
  
  ### 9) Print conclusions 
  
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
        nrow(site_df),
        "\t",
        "Background (m):",
        nrow(background_df),
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
          "\n", 
          "\n"
        )
      )
      cat(paste("DECISION:", "Reject null hypothesis, sufficient evidence of median difference", "\n"))
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
          "\n", 
          "\n"
        )
      )
      cat(paste("DECISION:", "Fail to reject null hypothesis, insufficient evidence of median difference", "\n"))
      cat(row, "\n")
    }
  }
  
  ### 10) Output results if an object has been defined 
  Gehan_scores <- as.data.frame(sorted_meas)
  colnames(Gehan_scores) = c("Samp Measurement", "Detect (Y/N)", "Samp Location", 
                             "d", "e", "R", "aR")
  Gehan_param <- list(alpha = alpha, power = power, deltaS = deltaS, G = G, p_val = p_val) 
  
  Gehan_results <- list(Gehan_param, Gehan_scores)
  return(Gehan_results)
} # Gehan function end








