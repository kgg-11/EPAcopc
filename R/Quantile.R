########## FUNCTION : Quantile  #######################
##' Perform quantile test.
#'
#' @param site Dataframe for site measurements.
#' @param measure.s Variable for \code{site} measurements that provides the corresponding concentration.
#' @param nd.s Indicator variable for \code{site} measurements that determines if the corresponding measurement (row) is a nondetect. A value of \code{1} indicates that the row is a nondetect (\code{0} is a detected measurement).
#' @param background Dataframe for background measurements.
#' @param measure.b Variable for \code{background} measurements that provides the corresponding concentration.
#' @param nd.b Indicator variable for \code{background} measurements that determines if the corresponding measurement (row) is a nondetect. A value of \code{1} indicates that the row is a nondetect (\code{0} is a detected measurement).
#' @param epsilon A proportion of the site that has concentrations greater than the background. Refer to Table 4.4 of "Guidance for Environmental Background Analysis Volume III: Groundwater" (2003) for possible values.
#' @param alpha Type I error rate. Options are only \code{0.01}, \code{0.025}, \code{0.05}, or \code{0.10}. Must be manually input. (default is \code{0.05})
#' @param power Statistical power. Options are only \code{0.80} or \code{0.90}. Must be manually input. (default is \code{0.80})
#' @param plot Logical. Will display a histogram of the site and background measurements with other details. (default is \code{FALSE})
#' @param print Logical. Will display quantile test results. (default is \code{TRUE})
#' @details Based off the protocols described in "Guidance for Environmental Background Analysis Volume III: Groundwater" (2003), the Quantile test is used to help answer the question: are concentrations of the target chemical within the site greater than the background concentration?
#' The test declares a chemical a COPC if k or more of the r measurements are from the site.
#' While the quantile test does not make assumptions of underlying distributions for the site and background distributions, some assumptions must still be met.
#' In particular, all nondetects must be smaller than the largest "r" measurements in the pooled data set.
#' The referenced literature provides the other assumptions with the Quantile test in Table 4-1, as well as advantages and disadvantages to this test.
#' 
#' @examples
#' # example code
#' example_site <-
#'   data.frame(
#'     samples = c(15, 15, 17, 23, 16, 30, 60, 89, 90, 100),
#'     nondetect = c(1, rep(0, 9))
#'   )
#' example_background <-
#'   data.frame(
#'     samples = c(23, 36, 30, 37, 44, 57, 60, 61, 61, 79),
#'     nondetect = c(1, 1, rep(0, 8))
#'   )
#' Quantile(
#'   alpha = 0.05,
#'   epsilon = 0.5,
#'   power = 0.8,
#'   site = example_site,
#'   measure.s = example_site$samples,
#'   nd.s = example_site$nondetect,
#'   background = example_background,
#'   measure.b = example_background$samples,
#'   nd.b = example_background$nondetect,
#'   plot = TRUE
#' )
#' @export
#' @references Naval Facilities Engineering Command. (2003, October).  \emph{Guidance for Environmental Background Analysis Volume III: Groundwater.} https://vsp.pnnl.gov/docs/Draft_Guidance_for_Review.pdf.



Quantile <- function(site, measure.s = site$measurement, nd.s = site$nondetect, background, measure.b = background$measurement, nd.b = background$nondetect, alpha=0.05, epsilon, power=0.80, print=TRUE, plot=FALSE){
  # immediately renaming measure.s and measure,b to "samples"  allow users to have some flexibility while using this current version of the Quantile test
  
  site$samples <- measure.s
  background$samples <- measure.b
  site$nondetect <- nd.s
  background$nondetect <- nd.b

  ###step 1:####
  ## determine sample size needed for both m (background) and n (site) ##
  # based off epsilon anc power input#
  table_4.4 <- data.frame (epsilon.choices = c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0),
                           power0.8_alpha0.01 = c(100, 55, 25, 20, 15, rep(10, 5)),
                           power0.9_alpha0.01 = c(100, 60, 30, 25, 20, 15, rep(10, 4)),
                           power0.8_alpha0.025 = c(100, 40, 20, 15, 15, rep(10, 5)),
                           power0.9_alpha0.025 = c(100, 40, 25, 20, 15, rep(10, 5)),
                           power0.8_alpha0.05 = c(80, 35, 20, 15, rep(10, 6)),
                           power0.9_alpha0.05 = c(100, 40, 20, 15, rep(10, 6)),
                           power0.8_alpha0.10 = c(55, 25, 15, rep(10, 7)),
                           power0.9_alpha0.10 = c(70, 35, 15, 15, rep(10, 6)))
  samples.required <- function(epsilon, power, alpha, table_4.4){
    #filter table based off epsilon, power, and alpha
    filtered_table <- table_4.4[table_4.4$epsilon.choices == epsilon,]
    # Check if any rows match the filter criteria
    if (nrow(filtered_table) > 0) {
      # Return the value from the correct column based on the power
      #power of 0.80
      if (power == 0.8) {
        if(alpha == 0.01){
          return(filtered_table$power0.8_alpha0.01)
        }
        else if(alpha == 0.025){
          return(filtered_table$power0.8_alpha0.025)
        }
        else if(alpha == 0.05){
          return(filtered_table$power0.8_alpha0.05)
        }
        else if(alpha == 0.10){
          return(filtered_table$power0.8_alpha0.10)
        }
      } # end of power = 0.8 statements
      # power 0f 0.90
      else if (power == 0.90){
        if(alpha == 0.01){
          return(filtered_table$power0.9_alpha0.01)
        }
        else if(alpha == 0.025){
          return(filtered_table$power0.9_alpha0.025)
        }
        else if(alpha == 0.05){
          return(filtered_table$power0.9_alpha0.05)
        }
        else if(alpha == 0.10){
          return(filtered_table$power0.9_alpha0.10)
        }
      } #end of power = 0.90
      else{
        return("invalid alpha")
      }
    }
  }
  nreq <- samples.required(epsilon = epsilon,
                           power = power,
                           alpha = alpha,
                           table_4.4 = table_4.4)
  print(c("minimum sample sizes for both site and background samples is: " , nreq))
  
  ## Step 2: ##
  # Check if provided sample sizes are enough for desired power, alpha, and epsilon
  #   Include nondetects
  condition.check <- if (nrow(site)>=nreq & nrow(background)>=nreq) {
    print("minimum sample sizes for background and site is satisfied")
    print(c("required rows", nreq))
    print(c("rows of site", nrow(site)))
    print(c("rows of background", nrow(background)))
  }
  else{
    warning("sample size is too small for power, epsilon, and alpha desired")
    print(c("required rows", nreq, "rows of site", nrow(site), "rows of background", nrow(background)))
    return()
  }
  ## Step 3:##
  # List pooled site and background from smallest to largest #
  site$source <- 1 #deciding that site is a 1
  background$source <- 0 # deciding that background is a 0
  pooled_samples <- plyr::rbind.fill(site, background)
  sorted_pooled <- dplyr::arrange(pooled_samples, samples)
  
  ## Step 4: ##
  # Find values of r and k needed#
  # based off n (site measurements the columns of arrays) and m (background measurments rows)
  m <- nrow(background) # nondetects included
  m_divided <- (m / 5)
  m_ceil <- ceiling(m_divided)
  m_rounded <- (m_ceil * 5)
  m_rounded <- as.character(m_rounded)
  #n transformation
  n <- nrow(site) #nondetects included
  n_divided <- (n / 5)
  n_ceil <- ceiling(n_divided)
  n_rounded <- (n_ceil * 5)
  n_rounded <- as.character(n_rounded)
  
  ################## ARRAY FOR TABLES ##################
  ################## ALPHA 0.01 ##################
  
  ################ CHANGE O TO NA ##################
  V_r <- c(NA ,NA , 11, 13, 16, 19, 22, 25, 28, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, #row =5
           NA, 6,7, 9, 11, 13, 14, 16, 18, 19, 21, 23, 25, 26, 28, 30, NA, NA, NA, NA, # row =10
           3, 7, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, #row =15
           6, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 19, 20, 21,  #row = 20
           4, 7, 4, 5, 6, 7, 8, 9, 9, 10, 11, 12, 12, 13, 14, 15, 16, 16, 17, 18, #row = 25
           4, 3, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12, 13, 14, 14, 15, 15, #row = 30
           2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11, 11, 12, 13, 13, 14, #row = 35
           2, 3, 7, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, #row = 40
           2, 6, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 10, 11, 11, #row = 45
           NA, 4, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9, 10, 10, 10, # row = 50
           NA, 4, 3, 7, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 9, 9, 9, 10, # row = 55
           NA, 4, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, # row = 60
           NA, 4, 3, 3, 6, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, # row = 65
           NA, 2, 6, 3, 7, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, # row = 70
           NA, 2, 4, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, # row = 75
           NA, 2, 4, 3, 3, 6, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, # row = 80
           NA, 2, 4, 3, 3, 7, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, # row = 85
           NA, NA, 4, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, #row = 90
           NA, NA, 4, 6, 3, 3, 6, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, # row = 95
           NA, NA, 4, 4, 3, 3, 7, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6 #row = 100
  )
  
  M_r <- matrix(V_r, nrow=20, ncol = 20, byrow=TRUE)
  V_k <- c(NA ,NA , 1, 13, 16, 19, 22, 25, 28, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, #row =5
           NA, 6,7, 9, 11, 13, 14, 16, 18, 19, 21, 23, 25, 26, 28, 30, NA, NA, NA, NA, # row =10
           3, 6, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18, 19, 21, 22, 23, 24, 26, 27, #row =15
           4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 19, 20, 21,  #row = 20
           3, 5, 4, 5, 6, 7, 8, 9, 9, 10, 11, 12, 12, 13, 14, 15, 16, 16, 17, 18, #row = 25
           3, 3, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12, 13, 14, 14, 15, 15, #row = 30
           2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11, 11, 12, 13, 13, 14, #row = 35
           2, 3, 5, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, #row = 40
           2, 4, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 10, 11, 11, #row = 45
           NA, 3, 3, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9, 10, 10, 10, # row = 50
           NA, 3, 3, 5, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 9, 9, 9, 10, # row = 55
           NA, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, # row = 60
           NA, 3, 3, 3, 5, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, # row = 65
           NA, 2, 4, 3, 5, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, # row = 70
           NA, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, # row = 75
           NA, 2, 3, 3, 3, 5, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, # row = 80
           NA, 2, 3, 3, 3, 5, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, # row = 85
           NA, NA, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7, #row = 90
           NA, NA, 3, 4, 3, 3, 5, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, # row = 95
           NA, NA, 3, 3, 3, 3, 5, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6 #row = 100
           
  )
  M_k <- matrix(V_k, nrow=20, ncol = 20, byrow=TRUE)
  V_a <-c(
    0, 0, 0.008, 0.015, 0.014, 0.013, 0.013, 0.013, 0.012, 0,0,0,0,0,0,0,0,0,0,0, # row 5
    0, 0.005, 0.013, 0.012, 0.011, 0.010, 0.014, 0.013, 0.012, 0.015, 0.014, 0.013, 0.012, 0.015, 0.014, 0.013, 0, 0,0,0, #row =10
    0.009, 0.007, 0.008, 0.012, 0.014, 0.009, 0.011, 0.013, 0.014, 0.011, 0.012, 0.013, 0.014, 0.015, 0.012, 0.013,0.014,  0.015, 0.013, 0.013, # row = 15
    0.009, 0.007, 0.008, 0.012, 0.014, 0.009, 0.011, 0.013, 0.014, 0.011, 0.012, 0.013, 0.014, 0.015, 0.012,0.013, 0.014, 0.015, 0.013, 0.013, #row = 20
    0.009, 0.012, 0.015, 0.013, 0.011, 0.010, 0.009, 0.009, 0.014, 0.012, 0.011, 0.011, 0.015, 0.014, 0.013, 0.012, 0.011, 0.014, 0.014, 0.013, #row = 25
    0.006, 0.012, 0.009, 0.007, 0.006, 0.012, 0.010, 0.008, 0.013, 0.011, 0.009, 0.013, 0.011, 0.010, 0.013, 0.012, 0.011, 0.014, 0.012, 0.015, #row = 30
    0.013, 0.008, 0.006, 0.014, 0.010, 0.007, 0.012, 0.009, 0.014, 0.011, 0.009, 0.013, 0.010, 0.014, 0.011, 0.015, 0.012, 0.011, 0.013, 0.012, #row = 35
    0.008, 0.008, 0.013, 0.007, 0.006, 0.012, 0.008, 0.013, 0.009, 0.013, 0.010, 0.014, 0.011, 0.014, 0.011, 0.014, 0.012, 0.014, 0.012, 0.014, #rorw =40
    0.008, 0.008, 0.013, 0.007, 0.014, 0.008, 0.014, 0.009, 0.013, 0.009, 0.013, 0.009, 0.012, 0.009, 0.012,0.009, 0.012, 0.015, 0.012, 0.014, # row = 45
    0, 0.013, 0.010, 0.005, 0.010, 0.006, 0.010, 0.015, 0.009, 0.013, 0.009, 0.012, 0.009, 0.011, 0.014, 0.011, 0.013, 0.010, 0.012, 0.015, #row = 50
    0, 0.010, 0.008, 0.013, 0.008,0.014,  0.007, 0.011, 0.007, 0.010, 0.014, 0.009, 0.012, 0.008, 0.010, 0.013, 0.009, 0.012, 0.014, 0.011, #row = 55
    0, 0.008, 0.007, 0.014, 0.006, 0.011, 0.006,0.009,  0.013, 0.007, 0.010, 0.014, 0.009, 0.011, 0.014, 0.010, 0.012, 0.015, 0.010, 0.013, #row = 60
    0, 0.007, 0.006, 0.012, 0.006, 0.009, 0.013, 0.007, 0.010, 0.014, 0.008, 0.011, 0.014, 0.009, 0.011, 0.014, 0.009, 0.011, 0.014, 0.010, #row = 65
    0, 0.014, 0.008, 0.010, 0.013, 0.007, 0.011, 0.005, 0.008, 0.011, 0.015, 0.008, 0.011, 0.014, 0.009, 0.011, 0.013, 0.009, 0.011, 0.013, #row = 70
    0, 0.013, 0.014, 0.008, 0.014, 0.006, 0.009, 0.013, 0.006, 0.009, 0.012, 0.007, 0.009, 0.011, 0.014, 0.009, 0.011, 0.013, 0.008, 0.010, #row = 75
    0, 0.011, 0.012, 0.007, 0.012, 0.006, 0.008, 0.011, 0.005, 0.007, 0.010, 0.013, 0.007, 0.009, 0.012, 0.014, 0.008, 0.010, 0.013, 0.015, # row =80
    0, 0.010, 0.010, 0.006, 0.011, 0.013, 0.006, 0.009, 0.013, 0.006, 0.008, 0.011, 0.014, 0.008, 0.010, 0.012, 0.014, 0.008, 0.010, 0.012, # row = 85
    0, 0, 0.009, 0.005, 0.009, 0.014, 0.005, 0.008, 0.011, 0.005, 0.007, 0.009, 0.012, 0.015, 0.008, 0.010, 0.012, 0.014, 0.008, 0.010, # row = 90
    0,0, 0.008, 0.008, 0.008, 0.013, 0.005, 0.007, 0.010, 0.013, 0.006, 0.008, 0.010, 0.013, 0.007, 0.008, 0.010, 0.012, 0.014, 0.008, #row = 95
    0,0, 0.007, 0.014, 0.007, 0.011, 0.013, 0.006, 0.008, 0.011, 0.015, 0.007, 0.009, 0.011, 0.013, 0.007, 0.009, 0.010, 0.012, 0.014 #row 100
    
    
  )
  
  M_a <- matrix(V_a, nrow=20, ncol = 20, byrow=TRUE)
  names <- seq(5, 100, by = 5)
  Array.01 <- array(c(M_r, M_k, M_a), dim=c(20, 20, 3), dimnames=list(names, names))
  
  ################## END FOR ALPHA 0.01 ##################
  
  ################## START FOR ALPHA 0.025 ##################
  
  V_r.025 <- c( NA, NA, 9, 12, 15, 17, 20, 22, 25, rep(NA,11), # row=5
                NA, 7, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26, 27, rep(NA,3), # row=10
                11, 6, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 21, 21, 22, 23, # row=15
                8, 3, 4, 5, 6, 7, 12, 13, 9, 10, 11, 12, 13, 13, 14, 15, 16, 17, 17, 18, # row=20
                2, 8, 6, 7, 5, 6, 10, 7, 8, 13, 9, 10, 11, 11, 12, 13, 13, 14, 15, 15, # row=25
                6, 6, 9, 4, 7, 5, 9, 6, 7, 12, 8, 9, 9, 10, 10, 11, 11, 12, 13, 13, #row=30
                7, 4, 3, 6, 4, 10, 5, 9, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, # row=35
                3, 4, 8, 11, 6, 4, 10, 5, 9, 6, 10, 7, 12, 8, 8, 9, 9, 10, 10, 11, # row=40
                3, 8, 6, 3, 8, 4, 7, 5, 5, 9, 6, 10, 7, 7, 8, 8, 8, 9, 9, 10, # row=45
                NA, 2, 6, 3, 11, 6, 4, 7, 5, 5, 9, 6, 6, 7, 7, 12, 8, 8, 13, 9, # row=50
                NA, 2, 4, 8, 3, 8, 4, 4, 10, 5, 5, 9, 6, 6, 10, 7, 7, 12, 8, 8, # row=55
                NA, 14, 4, 8, 3, 11, 6, 4, 7, 10, 5, 5, 9, 6, 6, 10, 7, 7, 7, 8, #row=60
                NA, 6, 7, 6, 10, 3, 8, 6, 4, 7, 10, 5, 5, 9, 6, 6, 10, 7, 7, 7, # row=65
                NA, 6, 2, 6, 8, 3, 13, 6, 4, 4, 7, 10, 5, 5, 9, 6, 6, 6, 10, 7, #row=70
                NA, 11, 2, 4, 8, 3, 9, 8, 6, 4, 7, 7, 10, 5, 5, 9, 6, 6, 6, 10, # row=75
                NA, 7, 2, 4, 6, 10, 3, 13, 6, 4, 4, 7, 10,5, 5, 5, 9, 6, 6, 6,
                NA, 3, 2, 4, 6, 8, 3, 9, 6, 8, 6, 4, 4, 7, 10, 5, 5, 9, 6, 6,
                NA, NA, 5, 11, 9, 8, 3, 3, 13, 6, 6, 4, 4, 7, 10, 5, 5, 5, 9, 9,
                NA, NA, 10, 2, 4, 6, 10, 3, 11, 8, 6, 4, 4, 7, 7, 10, 5, 5, 5, 9,
                NA, NA, 6, 2, 4, 6, 8, 3, 3, 13, 6, 6, 4, 4, 7, 10, 10, 5, 5, 5)
  M_r.025 <- matrix(V_r.025, nrow=20, ncol=20, byrow=TRUE)
  V_k_.025 <- c(
    NA, NA, 9, 12, 15, 17, 20, 22, 25, rep(NA,11), # row=5
    NA, 6, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26, 27, rep(NA,3), # row=10
    5, 5, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 21, 21, 22, 23, # row=15
    4, 3, 4, 5, 6, 7, 11, 12, 9, 10, 11, 12, 13, 13, 14, 15, 16, 17, 17, 18, # row=20
    2, 5, 5, 6, 5, 6, 9, 7, 8, 12, 9, 10, 11, 11, 12, 13, 13, 14, 15, 15, # row=25
    3, 4, 6, 4, 6, 5, 8, 6, 7, 11, 8, 9, 9, 10, 10, 11, 11, 12, 13, 13, #row=30
    3, 3, 3, 5, 4, 8, 5, 8, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, # row=35
    2, 3, 5, 7, 5, 4, 8, 5, 8, 6, 9, 7, 11, 8, 8, 9, 9, 10, 10, 11, # row=40
    2, 4, 4, 3, 6, 4, 6, 5, 5, 8, 6, 9, 7, 7, 8, 8, 8, 9, 9, 10, # row=45
    NA, 2, 4, 3, 7, 5, 4, 6, 5, 5, 8, 6, 6, 7, 7, 11, 8, 8, 12, 9, # row=50
    NA, 2, 3, 5, 3, 6, 4, 4, 8, 5, 5, 8, 6, 6, 9, 7, 7, 11, 8, 8, # row=55
    NA, 5, 3, 5, 3, 7, 5, 4, 6, 8, 5, 5, 8, 6, 6, 9, 7, 7, 7, 8, #row=60
    NA, 3, 4, 4, 6, 3, 6, 5, 4, 6, 8, 5, 5, 8, 6, 6, 9, 7, 7, 7, # row=65
    NA, 3, 2, 4, 5, 3, 8, 5, 4, 4, 6, 8, 5, 5, 8, 6, 6, 6, 9, 7, #row=70
    NA, 4, 2, 3, 5, 3, 6, 6, 6, 4, 6, 6, 8, 5, 5, 8, 6, 6, 6, 9, # row=75
    NA, 3, 2, 3, 4, 6, 3, 8, 5, 4, 4, 6, 8, 5, 5, 5, 8, 6, 6, 6,
    NA, 2, 2, 3, 4, 5, 3, 6, 6, 6, 5, 4, 4, 6, 8,5, 5, 5, 8, 6, 6,
    NA, NA, 3, 5, 5, 5, 3, 3, 8, 5, 5, 4, 4, 6, 8, 5, 5, 5, 8, 8,
    NA, NA, 4, 2, 3, 4, 6, 3, 7, 6, 5, 4, 4, 6, 6, 8, 5, 5, 5, 8,
    NA, NA, 3, 2, 3, 4, 5, 3, 3, 8, 5, 5, 4, 4, 6, 8, 8, 5, 5, 5
  )
  
  M_k.025 <- matrix(V_r.025, nrow=20, ncol=20, byrow=TRUE)
  
  V_a.025 <- c(
    NA, NA, 0.030, 0.024, 0.021, 0.026, 0.024, 0.028, 0.025, rep(NA,11),
    NA, 0.029, 0.028, 0.022, 0.029, 0.024, 0.029,0.025, 0.029,0.025, 0.029, 0.026, 0.029, 0.026, 0.029, 0.026, 0.029, NA, NA, NA,
    0.030, 0.023,0.021, 0.024, 0.026, 0.027, 0.028, 0.029, 0.030, 0.022, 0.023, 0.023,0.024, 0.025, 0.025, 0.026, 0.021, 0.027, 0.027, 0.027,
    0.023, 0.030, 0.026, 0.024, 0.022, 0.020, 0.021, 0.024, 0.028, 0.026, 0.024, 0.023, 0.022, 0.029, 0.027, 0.026, 0.025, 0.024, 0.029, 0.028,
    0.023, 0.027, 0.021, 0.023, 0.025, 0.020, 0.026, 0.027, 0.023, 0.027, 0.027, 0.024, 0.022, 0.028, 0.025, 0.023, 0.028, 0.025, 0.023, 0.028,
    0.026, 0.026, 0.026, 0.021, 0.029, 0.026, 0.024, 0.029, 0.023, 0.021, 0.025, 0.021, 0.027, 0.023, 0.029, 0.025, 0.030, 0.026, 0.023, 0.027,
    0.030, 0.030, 0.023, 0.020, 0.026, 0.022, 0.027, 0.024, 0.027, 0.020, 0.027, 0.021, 0.027, 0.022, 0.027, 0.022, 0.027, 0.022, 0.027, 0.023,
    0.029, 0.022, 0.028, 0.025, 0.028, 0.030, 0.026, 0.027, 0.023, 0.026, 0.028, 0.024, 0.020, 0.023, 0.029, 0.022, 0.027, 0.021, 0.026, 0.021,
    0.023, 0.029, 0.030, 0.026, 0.021, 0.023, 0.025, 0.020, 0.028, 0.023, 0.024,0.026, 0.022, 0.027, 0.020, 0.025, 0.030, 0.023, 0.027, 0.021,
    NA, 0.025, 0.022, 0.021, 0.027, 0.026, 0.026, 0.028, 0.021, 0.028, 0.022, 0.023, 0.029, 0.020, 0.025, 0.020, 0.022, 0.026, 0.027, 0.023,
    NA, 0.022, 0.029, 0.028, 0.028, 0.021, 0.020, 0.029, 0.021, 0.022, 0.028, 0.022, 0.023, 0.028, 0.029, 0.023, 0.027, 0.023, 0.023, 0.027,
    NA, 0.022, 0.024, 0.021, 0.023, 0.029, 0.024, 0.023, 0.023, 0.024, 0.023, 0.029, 0.022, 0.022, 0.027, 0.027, 0.021, 0.025, 0.030, 0.021,
    NA, 0.028, 0.021, 0.025, 0.025, 0.029, 0.021, 0.029, 0.026, 0.026, 0.026, 0.023, 0.029, 0.022, 0.021, 0.026, 0.026, 0.020, 0.024, 0.028,
    NA, 0.024, 0.029, 0.021, 0.028, 0.025, 0.026, 0.023, 0.022, 0.028, 0.028, 0.027, 0.024, 0.029, 0.022, 0.021, 0.025, 0.029, 0.030, 0.022,
    NA, 0.022, 0.026, 0.028, 0.022, 0.022, 0.028, 0.021, 0.027, 0.024, 0.023, 0.030, 0.029, 0.024, 0.029, 0.021, 0.021, 0.024, 0.028, 0.028,
    NA, 0.028, 0.024, 0.024, 0.028, 0.024, 0.027, 0.027, 0.023, 0.020, 0.026, 0.024, 0.023, 0.020, 0.025, 0.029, 0.021, 0.020, 0.024, 0.027,
    NA, 0.029, 0.021, 0.021, 0.023, 0.028, 0.023, 0.030, 0.020, 0.026, 0.022, 0.028, 0.026, 0.024, 0.021, 0.025, 0.029, 0.021, 0.020, 0.023,
    NA, NA, 0.020, 0.027, 0.023, 0.023, 0.021, 0.028, 0.028, 0.022, 0.029, 0.024, 0.029, 0.028, 0.026, 0.022, 0.025, 0.030, 0.021, 0.025,
    NA, NA, 0.029, 0.029, 0.028, 0.029, 0.023, 0.025, 0.026, 0.020, 0.025, 0.021, 0.026, 0.024, 0.029, 0.027, 0.022, 0.026, 0.030, 0.021,
    NA, NA, 0.029, 0.027, 0.025, 0.025, 0.028, 0.022, 0.029, 0.028, 0.022, 0.028, 0.023, 0.027, 0.025, 0.022, 0.028, 0.022, 0.026, 0.030
  )
  M_a.025 <- matrix(V_a.025, nrow=20, ncol = 20, byrow=TRUE)
  
  names <- seq(5, 100, by = 5)
  Array.025 <- array(c(M_r.025, M_k.025, M_a.025), dim=c(20, 20, 3), dimnames=list(names, names))
  Array.025["50", "50",]
  
  ################## END FOR ALPHA 0.025 ##################
  
  ################## ARRAY FOR ALPHA 0.05 ##################
  
  V_r.05 <- c( NA, NA, 8, 10, 13, 15, 17, 19, 21, rep(NA,11),
               NA, 4, 5, 14, 8,9,10,12,13,14,15,17,18,19,20,21,23,rep(NA,3),
               2,3,4,5,6,7,8,9,9,10,11,12,13,14,15,16,16,17,18,19,
               9,8,6,4,5,9,6,7,8,8,9,10,10,11,12,12,13,14,14,15,
               6,6,3,6,4,5,5,6,11,7,8,8,9,9,19,11,11,11,12,12,
               3,2,10,3,11,4,8,5,6,6,7,7,8,8,9,9,9,10,10,11,
               8,2,6,3,6,4,4,8,5,9,6,6,7,7,8,8,8,9,9,10,
               4,5,4,10,3,6,4,4,8,5,9,6,6,11,7,7,8,8,8,9,
               4,9,2,8,3,8,6,4,4,8,5,5,9,6,6,11,7,7,8,8,
               NA,6,2,6,12,3,8,6,4,4,8,5,5,9,6,6,6,7,7,7,
               NA,3,2,4,8,3,5,6,9,4,4,8,5,5,9,6,6,6,11,7,
               NA,3,5,4,6,3,3,8,6,9,4,4,13,5,5,5,9,6,6,6,
               NA,3,5,2,6,10,3,3,6,6,4,4,4,13,5,5,5,9,6,6,
               NA,8,9,2,4,8,5,3,3,6,6,4,4,4,13,5,5,5,9,9,
               NA,8,6,2,4,6,10,3,3,8,6,9,4,4,5,13,8,5,5,5,
               NA,4,6,5,2,6,8,5,3,3,6,6,9,4,4,7,13,8,5,5,
               NA,4,3,5,2,4,4,10,5,3,3,6,6,9,4,4,7,10,8,5,
               NA,NA,3,5,2,6,6,8,5,3,3,8,6,6,4,4,4,7,10,8,
               NA,NA,3,9,2,2,4,8,10,5,3,3,6,6,9,4,4,4,7,10,
               NA,NA,3,6,5,2,4,6,10,5,3,3,3,6,6,9,4,4,4,7
               
  )
  M_r.05 <-matrix(V_r.05, nrow=20, ncol=20, byrow=TRUE)
  V_k.05 <- c(NA, NA, 8, 10, 13, 15, 17, 19, 21, rep(NA,11),
              NA, 4, 5, 12, 8,9,10,12,13,14,15,17,18,19,20,21,23,rep(NA,3),
              2,3,4,5,6,7,8,9,9,10,11,12,13,14,15,16,16,17,18,19, #row 3
              4,5,5,4,5,8,6,7,8,8,9,10,10,11,12,12,13,14,14,15,
              3,4,3,5,4,5,5,6,10,7,8,8,9,9,19,11,11,11,12,12,
              2,2,6,3,8,4,7,5,6,6,7,7,8,8,9,9,9,10,10,11,
              3,2,4,3,5,4,4,7,5,8,6,6,7,7,8,8,8,9,9,10,
              2,3,3,6,3,5,4,4,7,5,8,6,6,10,7,7,8,8,8,9,
              2,4,2,5,3,6,5,4,4,7,5,5,8,6,6,10,7,7,8,8,
              NA,3,2,4,7,3,6,5,4,4,7,5,5,8,6,6,6,7,7,7,
              NA,2,2,3,5,3,4,5,7,4,4,7,5,5,8,6,6,6,10,7,
              NA,2,3,3,4,3,3,6,5,7,4,4,10,5,5,5,8,6,6,6,
              NA,2,3,2,4,6,3,3,5,5,4,4,4,10,5,5,5,8,6,6,
              NA,3,4,2,3,5,4,3,3,5,5,4,4,4,10,5,5,5,8,8,
              NA,3,3,2,3,4,6,3,3,6,5,7,4,4,5,10,7,5,5,5,
              NA,2,3,3,2,4,5,4,3,3,5,5,7,4,4,6,10,7,5,5,
              NA,2,2,3,2,3,3,6,4,3,3,5,5,7,4,4,6,8,7,5,
              NA,NA,2,3,2,4,4,5,4,3,3,6,5,5,4,4,4,6,8,7,
              NA,NA,2,4,2,2,3,5,6,4,3,3,5,5,7,4,4,4,6,8,
              NA,NA,2,3,3,2,3,4,6,4,3,3,3,5,5,7,4,4,4,6
              
  )
  M_k.05 <- matrix(V_k.05, nrow=20, ncol = 20, byrow=TRUE)
  
  V_a.05 <- V_a.05 <- c(
    NA, NA, 0.051, 0.057, 0.043, 0.048, 0.051, 0.054, 0.056, rep(NA, 11),
    NA, 0.043, 0.057, 0.045, 0.046, 0.052, 0.058, 0.046, 0.050, 0.054, 0.057, 0.049, 0.052, 0.055, 0.057, 0.059, 0.053, NA, NA, NA,
    0.053, 0.052, 0.050, 0.048, 0.046, 0.045, 0.044, 0.043, 0.060, 0.057, 0.055, 0.054, 0.052, 0.051, 0.050, 0.049, 0.058, 0.057, 0.056, 0.055,
    0.040, 0.056, 0.040, 0.053, 0.043, 0.052, 0.056, 0.048, 0.043, 0.057, 0.051, 0.046, 0.057, 0.052, 0.048, 0.057, 0.053, 0.049, 0.057, 0.054,
    0.041, 0.043, 0.046, 0.052, 0.055, 0.041, 0.059, 0.046, 0.042, 0.050, 0.042, 0.053, 0.045, 0.055, 0.048, 0.042, 0.050, 0.058, 0.052, 0.060,
    0.047, 0.058, 0.052, 0.058, 0.045, 0.056, 0.045, 0.054, 0.040, 0.053, 0.041, 0.052, 0.042, 0.051, 0.042, 0.050, 0.059, 0.049, 0.057, 0.049,
    0.046, 0.045, 0.058, 0.043, 0.041, 0.040, 0.057, 0.043, 0.051, 0.052, 0.047, 0.058, 0.043, 0.053, 0.041, 0.049, 0.057, 0.046, 0.053, 0.044,
    0.055, 0.048, 0.057, 0.059, 0.053, 0.048, 0.043, 0.058, 0.042, 0.048, 0.047, 0.042, 0.051, 0.042, 0.045, 0.053, 0.041, 0.048, 0.055, 0.043,
    0.045, 0.047, 0.059, 0.052, 0.042, 0.041, 0.054, 0.045, 0.058, 0.041, 0.046, 0.057, 0.056, 0.047, 0.055, 0.046, 0.047, 0.054, 0.041, 0.047,
    NA, 0.052, 0.050, 0.051, 0.050, 0.049, 0.049, 0.059, 0.047, 0.059, 0.041, 0.045, 0.054, 0.051, 0.043, 0.050, 0.058, 0.042, 0.048, 0.054,
    NA, 0.059, 0.043, 0.056, 0.058, 0.041, 0.041, 0.046, 0.042, 0.048, 0.059, 0.040, 0.043, 0.052, 0.048, 0.040, 0.047, 0.054, 0.043, 0.043,
    NA, 0.052, 0.052, 0.046, 0.059, 0.035, 0.047, 0.043, 0.051, 0.046, 0.049,0.059, 0.052, 0.042, 0.050, 0.058, 0.054,0.044, 0.050, 0.056,
    NA, 0.045, 0.043, 0.053, 0.048, 0.050, 0.040, 0.053, 0.041, 0.055, 0.042, 0.050, 0.060, 0.052, 0.041, 0.048, 0.055, 0.051, 0.041, 0.047,
    NA, 0.057, 0.048, 0.047, 0.055, 0.050, 0.041, 0.046, 0.057, 0.045, 0.058, 0.043, 0.051, 0.060, 0.051, 0.041, 0.047, 0.054, 0.048, 0.057,
    NA, 0.049, 0.056, 0.043, 0.047, 0.054, 0.053, 0.040, 0.051, 0.044, 0.049, 0.041, 0.044, 0.052, 0.060, 0.051, 0.047, 0.046, 0.052, 0.058,
    NA, 0.059, 0.048, 0.053, 0.055, 0.046, 0.055, 0.042, 0.045, 0.055, 0.041, 0.052, 0.043, 0.045, 0.053, 0.058, 0.051, 0.046, 0.045, 0.051,
    NA, 0.054, 0.058, 0.047, 0.050, 0.054, 0.048, 0.056, 0.049, 0.049, 0.059, 0.044, 0.055, 0.046, 0.046, 0.053, 0.059, 0.060, 0.045, 0.044,
    NA, NA, 0.053, 0.041, 0.046, 0.059, 0.051, 0.058, 0.042, 0.044, 0.053, 0.045, 0.047, 0.058, 0.041, 0.047, 0.054, 0.059, 0.060, 0.045,
    NA, NA, 0.048, 0.048, 0.042, 0.056, 0.059, 0.050, 0.058, 0.048, 0.048, 0.056, 0.041, 0.050, 0.040, 0.042, 0.048, 0.054, 0.059, 0.059,
    NA, NA, 0.044, 0.057, 0.054, 0.052, 0.053, 0.056, 0.049, 0.043, 0.043, 0.051, 0.059, 0.044, 0.053, 0.042, 0.043, 0.049, 0.055, 0.059
  )
  M_a.05 <- matrix(V_a.05, nrow=20, ncol = 20, byrow=TRUE)
  names <- seq(5, 100, by = 5)
  Array.05 <- array(c(M_r.05, M_k.05, M_a.05), dim=c(20, 20, 3), dimnames=list(names, names))
  #Array.05["10","10",]
  
  ################## END FOR ALPHA 0.05 ##################
  ################## START FOR ALPHA 0.10 ##################
  V_r.1 <- c( NA,NA,7,8,10,12,14,15,17,rep(NA,11),
              NA,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,rep(NA,3),
              9,10,3,4,5,5,6,7,7,8,9,9,10,11,11,12,13,13,14,15,
              3,2,5,3,4,4,5,10,6,7,7,8,8,9,9,10,10,11,11,12,
              4,7,8,3,3,4,4,8,5,10,6,6,7,7,8,8,8,9,9,10,
              4,5,2,14,3,3,9,4,8,5,5,6,6,6,7,7,7,8,8,8,
              5,3,2,6,5,3,3,9,4,4,8,5,5,6,6,6,6,7,7,7,
              5,3,5,2,12,5,3,6,9,4,4,8,5,5,5,6,6,6,6,7,
              6,3,5,2,6,7,5,3,6,9,4,4,4,8,5,5,5,6,6,6, #row=45
              NA,7,9,7,2,10,5,3,3,6,9,4,4,4,8,5,5,5,5,6,
              NA,4,3,5,2,6,14,5,3,3,6,9,4,4,4,4,8,5,5,5,
              NA,4,3,5,2,2,8,5,5,3,3,6,9,4,4,4,4,8,5,5,
              NA,4,3,5,7,2,6,12,5,5,3,3,6,9,7,4,4,4,8,8,
              NA,5,7,9,5,2,2,8,7,5,3,3,3,6,9,7,4,4,4,4,
              NA,5,7,3,5,7,2,2,10,5,5,3,3,3,6,9,7,4,4,4,
              NA,5,4,3,5,7,2,2,8,14,5,5,3,3,3,6,6,9,4,4,
              NA,4,4,3,9,5,2,2,2,10,7,5,5,3,3,3,6,6,9,4,
              NA,NA,4,3,3,5,7,2,2,8,12,5,5,3,3,3,3,6,6,9,
              NA,NA,4,7,3,5,7,2,2,2,10,14,5,5,3,3,3,3,6,6,
              NA,NA,4,7,3,5,5,2,2,2,6,12,7,5,5,3,3,3,3,6
  )
  V_k.1 <- c( NA,NA,7,8,10,12,14,15,17,rep(NA,11),
              NA,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,rep(NA,3),
              4,6,3,4,5,5,6,7,7,8,9,9,10,11,11,12,13,13,14,15,
              2,2,4,3,4,4,5,9,6,7,7,8,8,9,9,10,10,11,11,12,
              2,4,5,3,3,4,4,7,5,9,6,6,7,7,8,8,8,9,9,10,
              2,3,2,8,3,3,7,4,7,5,5,6,6,6,7,7,7,8,8,8,
              2,2,2,4,4,3,3,7,4,4,7,5,5,6,6,6,6,7,7,7,
              2,2,3,2,7,4,3,5,7,4,4,7,5,5,5,6,6,6,6,7,
              2,2,3,2,4,5,4,3,5,7,4,4,4,7,5,5,5,6,6,6, #row=45
              NA,3,4,4,2,6,4,3,3,5,7,4,4,4,7,5,5,5,5,6,
              NA,2,2,3,2,4,8,4,3,3,5,7,4,4,4,4,7,5,5,5,
              NA,2,2,3,2,2,5,4,4,3,3,5,7,4,4,4,4,7,5,5,
              NA,2,2,3,4,2,4,7,4,4,3,3,5,7,6,4,4,4,7,7,
              NA,2,3,4,3,2,2,5,5,4,3,3,3,5,7,6,4,4,4,4,
              NA,2,3,2,3,4,2,2,6,4,4,3,3,3,5,7,6,4,4,4,
              NA,2,2,2,3,4,2,2,5,8,4,4,3,3,3,5,5,7,4,4,
              NA,2,2,2,4,3,2,2,2,6,5,4,4,3,3,3,5,5,7,4,
              NA,NA,2,2,2,3,4,2,2,5,7,4,4,3,3,3,3,5,5,7,
              NA,NA,2,3,2,3,4,2,2,2,6,8,4,4,3,3,3,3,5,5,
              NA,NA,2,3,2,3,3,2,2,2,4,7,5,4,4,3,3,3,3,5
  )
  V_a.1 <- c(
    NA, NA, 0.083, 0.116, 0.109, 0.104, 0.100, 0.117, 0.112, rep(NA,11),
    NA, 0.105, 0.108, 0.109, 0.109, 0.109, rep(0.109,11), NA,NA,NA,
    0.098, 0.106, 0.112, 0.093, 0.081, 0.117, 0.102, 0.092, 0.118, 0.106, 0.098, 0.118, 0.109, 0.101, 0.118, 0.110, 0.104, 0.118, 0.111, 0.106,
    0.091, 0.103, 0.093, 0.115, 0.085, 0.119, 0.093, 0.084, 0.099, 0.083, 0.102, 0.088, 0.105, 0.092, 0.107, 0.095, 0.108, 0.098, 0.110, 0.100,
    0.119, 0.084, 0.112, 0.080, 0.117, 0.080, 0.107, 0.108, 0.101, 0.088, 0.096, 0.114, 0.093, 0.108, 0.091, 0.104, 0.117, 0.100, 0.112, 0.098,
    0.089, 0.089, 0.106, 0.111, 0.088, 0.119, 0.116, 0.100, 0.093, 0.088, 0.106, 0.080, 0.095, 0.110, 0.087, 0.100, 0.113, 0.092, 0.103, 0.115,
    0.109, 0.119, 0.086, 0.119, 0.091, 0.093, 0.120, 0.112, 0.094, 0.114, 0.107, 0.094, 0.110, 0.081, 0.094, 0.107, 0.120, 0.094, 0.105, 0.116,
    0.087, 0.098, 0.119, 0.107, 0.109, 0.102, 0.097, 0.100, 0.109, 0.090, 0.107, 0.097, 0.086, 0.099, 0.112, 0.082, 0.093, 0.104, 0.116, 0.089,
    0.103, 0.082, 0.094, 0.091, 0.115, 0.086, 0.112, 0.100, 0.101, 0.107, 0.087, 0.102, 0.117, 0.107, 0.091, 0.103, 0.115, 0.083, 0.093, 0.103,
    NA, 0.083, 0.115, 0.097, 0.108, 0.112, 0.090, 0.084, 0.103, 0.102, 0.105, 0.084, 0.098, 0.112, 0.099, 0.084, 0.095, 0.105, 0.116, 0.083,
    NA, 0.109, 0.114, 0.114, 0.095, 0.112, 0.111, 0.098, 0.088, 0.105, 0.103, 0.104, 0.082, 0.095, 0.107, 0.120, 0.107, 0.088, 0.098, 0.108,
    NA, 0.095, 0.100, 0.097, 0.084, 0.109, 0.119, 0.082, 0.105, 0.091, 0.106, 0.103, 0.102, 0.081, 0.092, 0.103, 0.115, 0.100, 0.083, 0.092,
    NA, 0.084, 0.089, 0.082, 0.090, 0.097, 0.110, 0.113, 0.089, 0.111, 0.093, 0.108, 0.104, 0.101, 0.084, 0.090, 0.100, 0.110, 0.094, 0.107,
    NA, 0.115, 0.101, 0.106, 0.112, 0.088, 0.109, 0.114, 0.081, 0.096, 0.083, 0.096, 0.109, 0.104, 0.101, 0.082, 0.088, 0.097, 0.107, 0.117,
    NA, 0.103, 0.088, 0.111, 0.098, 0.101, 0.099, 0.119, 0.117, 0.083, 0.102, 0.085, 0.098, 0.110, 0.105, 0.100, 0.081, 0.086, 0.095, 0.104,
    NA, 0.093, 0.116, 0.101, 0.086, 0.086, 0.091, 0.109, 0.110, 0.110, 0.089, 0.107, 0.088, 0.099, 0.111, 0.105, 0.120, 0.116, 0.084, 0.093,
    NA, 0.106, 0.106, 0.092, 0.117, 0.111, 0.083, 0.101, 0.118, 0.112, 0.084, 0.094, 0.111, 0.090, 0.101, 0.112, 0.105, 0.119, 0.114, 0.083,
    NA, NA, 0.097, 0.085, 0.119, 0.099, 0.095, 0.093, 0.109, 0.108, 0.114, 0.083, 0.099, 0.082, 0.092, 0.102, 0.112, 0.105, 0.119, 0.113,
    NA, NA, 0.089, 0.100, 0.110, 0.089, 0.084, 0.086, 0.102, 0.117, 0.108, 0.117, 0.088, 0.103, 0.084, 0.094, 0.103, 0.113, 0.106, 0.118,
    NA, NA, 0.082, 0.090, 0.102, 0.080, 0.109, 0.080, 0.095, 0.110, 0.118, 0.109, 0.086, 0.093, 0.108, 0.086, 0.095, 0.104, 0.114, 0.106
    
  )
  
  M_r.1 <- matrix(V_r.1, nrow=20, ncol=20, byrow=TRUE)
  M_k.1 <- matrix(V_k.1, nrow=20, ncol = 20, byrow=TRUE)
  M_a.1 <- matrix(V_a.1, nrow=20, ncol = 20, byrow=TRUE)
  names <- seq(5, 100, by = 5)
  Array.1 <- array(c(M_r.1, M_k.1, M_a.1), dim=c(20, 20, 3), dimnames=list(names, names))
  
  
  
  
  
  ################## END FOR ALPHA 0.10 ##################
  ################## END FOR ARRAY FOR TABLES ##################
  
  # getting r, k, and alpha values based off sample sizes m and n
  r_val <- if (alpha == 0.01){
    Array.01[m_rounded, n_rounded, ]
  }
  else if(alpha == 0.025){
    Array.025[m_rounded, n_rounded, ]
  }
  else if(alpha == 0.05){
    Array.05[m_rounded, n_rounded, ]
  }
  else if(alpha == 0.10){
    Array.1[m_rounded, n_rounded, ]
  }
  print(c("r value: ", r_val[1])) #shows user the r value needed
  print(c("k value: ", r_val[2])) #shows the user the k value
  print(c("alpha value: ", r_val[3])) #shows user the alpha
  
  ### Step 5: ###
  # subset sorted_pooled to include only detects, then run the rest if this code based off that detectable only data
  sorted_pooled <- subset(sorted_pooled, nondetect==0) #0 is a detectable location
  ## From sorted list, determine if k >= r ##
  # k is amount of samples from largest r values that are from site ##
  #first, extract "r" largest samples from the pooled data
  r <- r_val[1] #defines r based off the array
  r_samples <- tail(sorted_pooled, r) # this makes a new data frame consisting of the last r objects
  print(r_samples)
  
  #next, determine if sum of these is greater than or equal to k (site =1)
  k <- r_val[2] #defining "k" value based off array
  sum_r_samples <- (sum(r_samples$source))
  
  
  
  
  ########### RESULTS: #############
  if (print) {
    text.a <- "Quantile Test Results:"
    cat(paste("\n", text.a, "\n", sep = ""))
    row <- paste(rep("=", 87), collapse = "")
    cat(row, "\n")
    cat(
      paste("SAMPLE SIZE:", "\t",
            "Recommended (power):"),
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
    if (sum_r_samples >= k) {
      cat(row, "\n")
      cat(paste("Measurements:", "\n"))
      cat(
        paste(
          "Amount of largest values (r):",
          r,
          "\t",
          "# of largest r values that are from site:",
          sum_r_samples,
          "\n"))
      cat(paste0("# of largest r values that need to be from site for COPC with alpha of ", r_val[[3]]," (k): ", # tried to incorporate the exact alpha we are using here ***
                 k,
                 "\n",
                 "\n"
      )
      )
      cat(paste("DECISION:", "\t","Reject null hypothesis. There is evidence that this is a COPC.", "\n"))
      cat(row, "\n")
    } else {
      cat(row, "\n")
      cat(paste("Measurements:", "\n"))
      cat(
        paste(
          "Amount of largest values (r):",
          r,
          "\t",
          "# of largest r values that are from site:",
          sum_r_samples,
          "\n"))
      cat(paste0("# of largest r values that need to be from site for COPC with alpha of ", r_val[[3]]," (k): ", # tried to incorporate the exact alpha we are using here ***
                 k,
                 "\n",
                 "\n"
      )
      )
      cat(paste("DECISION:", "\t","Fail to reject null hypothesis. There is not enough evidence to claim this is a COPC.", "\n"))
      cat(row, "\n")
    }
  }
  ######## Plot #########
  if (plot) {
    
    # data needed for plots
    site_sorted <- subset(sorted_pooled, source==1)
    site_sorted$measure_notNDs <- site_sorted$samples 
    background_sorted <- subset(sorted_pooled, source==0)
    background_sorted$measure_notND <- background_sorted$samples 
    r_samples_site <- subset(r_samples, source==1)
    r_samples_site$measures <- r_samples_site$samples
    
    site$ND_tfs <- ifelse(nd.s==1, TRUE, NA)
    background$ND_tf <- ifelse(nd.b == 1, TRUE, NA) # getting rid of zero's and making them NA's for easier usage in ND bar plot
    
    # 2 ways of dynamically changing the bin size...
    #binw <- round(max(site$measure_max_measure, na.rm=T)/20)
    binw2 <- round(max(range(site_sorted$measure_notNDs, na.rm=T),range(background_sorted$measure_notND, na.rm=T))/15)
    
    # Creating histogram
    my_histogram <- ggplot2::ggplot(color = "black") +
      ggplot2::geom_histogram(
        data = site_sorted,
        ggplot2::aes(x = measure_notNDs, fill = "measure_notNDs"),
        color = "black",
        alpha = 0.8,
        binwidth = binw2,
        na.rm = TRUE
      ) +
      ggplot2::geom_histogram(
        data = background_sorted,
        ggplot2::aes(x = measure_notND, fill = "measure_notND"),
        color = "black",
        alpha = 0.8,
        binwidth = binw2,
        na.rm = TRUE
      ) +
      # # Coloring over the max r values layer --> black eventually
      ggplot2::geom_histogram(
        data=r_samples,
        ggplot2::aes(x=samples, fill="samples"),
        color="black",
        alpha=1,
        binwidth=binw2,
        na.rm=TRUE
      ) +
      # Coloring the site r's a gold color eventually
      ggplot2::geom_histogram(
        data=r_samples_site,
        ggplot2::aes(x=measures, fill="measures"),
        color="black",
        alpha=1,
        binwidth=binw2,
        na.rm=TRUE
      )
    # Customizing histogram more
    my_histogram <- my_histogram +
      ggplot2::theme_classic() +
      ggplot2::xlab("Concentration") +
      ggplot2::ylab(NULL) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      # legend customization
      ggplot2::scale_fill_manual(
        values = c("#1F968BFF", "#440154FF", "#FDE725FF", "black"),
        name = "Measurement",
        labels = c("Background","Site","Site measurements within r","r")
      ) +
      ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.75))
    
    # # Creating bar chart of non-detects
    my_bar <- ggplot2::ggplot() +
      ggplot2::geom_bar(
        data = subset(site,!is.na(ND_tfs)), # Site bars
        ggplot2::aes(ND_tfs, fill = "ND_tfs"),
        alpha = 0.8,
        color = "black"
      ) +
      ggplot2::geom_bar(
        data = subset(background,!is.na(ND_tf)), # Background bars
        ggplot2::aes(ND_tf, fill = "ND_tf"),
        alpha = 0.8,
        color = "black"
      ) +
      ggplot2::theme_classic() +
      ggplot2::xlab("ND") +
      ggplot2::ylab('Number of Measurements') +
      ggplot2::scale_fill_manual(values = c("#1F968BFF", "#440154FF")) +
      ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.75))
    
    # Creating a dynamic y axis that will scale with the largest count in either the histogram or the ND plots.
    data_histogram <- ggplot2::ggplot_build(my_histogram) # pulling count data from histogram
    data_bar <- ggplot2::ggplot_build(my_bar) # pulling count data from bar chart
    max_counts_histogram <- max(data_histogram$data[[1]]$count, data_histogram$data[[2]]$count) # taking max counts from either site/background
    max_counts_bar <- max(data_bar$data[[1]]$count, data_bar$data[[2]]$count) # taking max counts from either site/background
    max_counts <- max(max_counts_histogram, max_counts_bar) # finding absolute max counts
    
    my_bar <- my_bar +
      ggplot2::scale_y_continuous(limits = c(0,max_counts), expand=c(0,0)) # how to create the dynamic y-axis
    my_histogram <- my_histogram +
      ggplot2::scale_y_continuous(limits = c(0,max_counts), expand=c(0,0)) # how to create the dynamic y-axis
    
    # Putting final figure together with bar plot and histogram side-by-side (using `ggpubr` package)
    my_legend <- ggpubr::get_legend(my_histogram) # using legend information from histogram
    
    figure_quantile <-
      ggpubr::ggarrange(
        my_bar + ggpubr::rremove("x.text") + ggpubr::rremove("x.ticks"), # removing parts of axis not needed
        my_histogram + ggpubr::rremove("y.axis") + ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text"), # removing parts of axis for cohesiveness
        ncol = 2,
        nrow = 1,
        common.legend = TRUE,
        legend = "top",
        legend.grob = my_legend,
        widths = c(0.15, 1), # relative widths of bar:histogram... might need to tweak this more to make it dynamic
        align = "h" # align the horizontal axes
      )
    print(figure_quantile)
  }
}
