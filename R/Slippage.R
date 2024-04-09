#' Perform slippage test.
#'
#' @param site A dataframe for site measurements.
#' @param measure.s Variable for \code{site} measurements that provides the corresponding concentration. Default is a column of site dataframe named measurement (\code{site$measurement}).
#' @param nd.s Indicator variable for \code{site} measurements that determines if the corresponding measurement (row) is a nondetect. A value of \code{1} indicates that the row is a nondetect (\code{0} is a detected measurement). Default is a column of site dataframe named nondetect (\code{site$nondetect}).
#' @param background A dataframe for background measurements.
#' @param measure.b Variable for \code{background} measurements that provides the corresponding concentration. Default is a column of background dataframe named measurement (\code{background$measurement}).
#' @param nd.b Indicator variable for \code{background} measurements that determines if the corresponding measurement (row) is a nondetect. A value of \code{1} indicates that the row is a nondetect (\code{0} is a detected measurement). Default is a column of background dataframe named nondetect (\code{background$nondetect}).
#' @param epsilon  A proportion of the site that has concentrations greater than the background. Refer to Table 4.3 of referenced literature for possible values.
#' @param alpha Type I error rate. Options are only \code{0.05} (default) or \code{0.01}.
#' @param power Statistical power. Options are only \code{0.80} (default) or \code{0.90}.
#' @param plot Logical. Will display a histogram of the site and background measurements with other details. (default is \code{FALSE})
#' @param print Logical. Will display slippage test results. (default is \code{TRUE})
#' @return A list of measures used in the Slippage test, including:
#' \itemize{
#'   \item \code{alpha} - Type I error rate.
#'   \item \code{power} - Statistical power.
#'   \item \code{maxbg} - Maximum background concentration measurement.
#'   \item \code{sitemax.n} - The number of detected site measurements greater than the maximum background measurement (\code{maxbg}).
#'   \item \code{sitemax.m} - The corresponding concentration measurements of \code{sitemax.n}.
#'   \item \code{critical.value} - The critical value, determined by sample size of site (n) and background (m) (including nondetects). Based on Table B-2 and Table B-3 from the literature.
#' }
#' @details The slippage test is one of multiple possible nonparametric tests described in "Guidance for Environmental Background Analysis Volume III: Groundwater" (2003) for determining if a chemical is a Chemical of Potential Concern (COPC).
#' Although the slippage test has the advantage of not making any assumptions on the distribution of either site or background measurements, there are still assumptions that need to be met for performing the slippage test. 
#' There are multiple assumptions that need to be upheld during the usage of the \code{slippage()} function.
#' The maximum background measurement must be a detected measurement; otherwise, an error is generated.
#' Other conditions are implemented based on the included values in certain tables of the referenced literature.
#' \itemize{
#' \item The sample size of both \code{site} and \code{background} measurements must be between 1 and 50 due to the constraints of Tables B-2 and B-3 in the referenced literature.
#'  \item The \code{epsilon} and \code{power} arguments must be able to be located in Table 4-3 of the referenced literature; otherwise, an error is generated.
#'  \item  The critical value must be able to be located in Tables B-2 and B-3 of the referenced literature, otherwise, an error is generated.
#'}
#'
#' @examples
#' set.seed(4321)
#' example_site <-
#'   data.frame(measure = round(rchisq(
#'     n = 50, df = 5, ncp = 6
#'   ), 1),
#'   nondetect = rbinom(n = 50, size = 1, p = 0.2))
#' example_background <-
#'   data.frame(measure = round(rchisq(n = 50, df = 5), 1),
#'              nondetect = rbinom(n = 50, size = 1, p = 0.2))
#' example1 <- slippage(
#'   site = example_site,
#'   measure.s = example_site$measure,
#'   nd.s = example_site$nondetect,
#'   background = example_background,
#'   measure.b = example_background$measure,
#'   nd.b = example_background$nondetect,
#'   epsilon = 0.6,
#'   plot = TRUE
#' )
#' example1
#' @export
#' @references Naval Facilities Engineering Command. (2003, October).  *Guidance for Environmental Background Analysis Volume III: Groundwater.* https://vsp.pnnl.gov/docs/Draft_Guidance_for_Review.pdf.

slippage <-
  function(site,
           measure.s = site$measurement,
           nd.s = site$nondetect,
           background,
           measure.b = background$measurement,
           nd.b = background$nondetect,
           epsilon,
           alpha = 0.05,
           power = 0.80,
           print = TRUE,
           plot = FALSE) {
    
  ############ Condition check: if (any row(background|nondetect=False))>max(background|nondetect=True)   ############
  background_detects <- subset(background, nd.b==0)
  background_nondetects <- subset(background, nd.b==1)
  condition_maxBackground <-
    if (nrow(background_nondetects)>0) {
      if (any(background_detects$measure>max(background_nondetects$measure))) {
        TRUE
      } else {
        stop("Assumption violated: Max background measurement is a non-detect")
      }
    } else {TRUE}

  ############ Condition check: Calculate the required sample size for desired power   ############
  n.required <- function(epsilon, power, table_4.3) {
    # recreating Table 4.3 ("Guidance for Environmental Background Analysis Volume III: Groundwater", 2003)
    table_4.3 <- data.frame (epsilon.choices = c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60),
                             power0.8 = c(60, 40, 30, 25, rep(15, 3), rep(10, 3)),
                             power0.9 = c(75, 50, 35, 30, 25, 20 , 20, 15, 10, 10))
    # Filter the dataframe based on epsilon and power inputs
    filtered_df <- table_4.3[table_4.3$epsilon.choices == epsilon,]
    # Check if any rows match the filter criteria
    if (nrow(filtered_df) > 0) {
      # Return epsilon from the correct column based on the power
      if (power == 0.8) {
        return(filtered_df$power0.8)
      } else if (power == 0.9) {
        return(filtered_df$power0.9)
      } else {
        stop("Invalid power value: must be 0.8 or 0.9")
      }
    } else {
      stop("Epsilon value not found")
    }
  }
  nreq <- n.required(epsilon, power, table_4.3) # saves output of n.required to a new value

  # check if required sample size (including non-detects) is enough for power desired
  m <- nrow(background) # nondetects included
  n <- nrow(site) #nondetects included
  condition_sampleSize <- if (n >= nreq & m >= nreq) {
    TRUE
  } else {
    warning("Sample size is small for power desired")
  }

  ############ K Calculation  ############
  # Count the number, K, of detected site measurements that are larger than the largest detected background measurement. In making this determination, ignore all nondetects in the site dataset
  maxbg <- max(background_detects)  # max background measurement that is detected

  site_detects <- subset(site, nd.s==0) # site measurements that are detected
  k <- sum(site_detects>maxbg) # number of site measurements > max background measurement (both detects only)

  ############ CRITICAL VALUE TABLE BELOW: ############  #############################
  # m and n (non-detects included) are used as indices.
  table_B2 <- as.matrix(read.table(text = "/ / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
/ / / / / / / / / / / / 13 14 15 16 17 18 19 20 21 22 23 23 24 25 26 27 28 29 30 31 32 32 33 34 35 36 37 38 39 40 41 41 42 43 44 45 46 47
/ / / / / / 7 8 9 10 11 11 12 13 14 15 15 16 17 18 18 19 20 21 22 22 23 24 25 26 26 27 28 29 30 30 31 32 33 33 34 35 36 37 37 38 39 40 40 41
/ / / / 5 6 7 8 8 9 10 10 11 12 12 13 14 14 15 16 16 17 18 19 19 20 21 21 22 23 23 24 25 25 26 27 27 28 29 30 30 31 32 32 33 34 34 35 35 36
/ / / 4 5 6 6 7 8 8 9 9 10 11 11 12 12 13 14 14 15 15 16 17 17 18 18 19 20 20 21 22 22 23 23 24 25 25 26 26 27 28 28 29 29 30 31 31 32 32
/ / / 4 5 5 6 6 7 8 8 9 9 10 10 11 11 12 12 13 14 14 15 15 16 16 17 17 18 18 19 19 20 21 21 22 22 23 23 24 24 25 25 26 26 27 28 28 29 29
/ / 3 4 5 5 6 6 7 7 8 8 9 9 10 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 21 21 22 22 23 23 24 24 25 25 26 26 26
/ / 3 4 4 5 5 6 6 7 7 8 8 8 9 9 10 10 11 11 12 12 12 13 13 14 14 15 15 16 16 16 17 17 18 18 19 19 19 20 20 21 21 22 22 23 23 23 24 24
/ / 3 4 4 5 5 5 6 6 7 7 8 8 8 9 9 10 10 10 11 11 12 12 12 13 13 14 14 14 15 15 16 16 16 17 17 18 18 18 19 19 20 20 20 21 21 22 22 22
/ / 3 4 4 4 5 5 6 6 6 7 7 7 8 8 9 9 9 10 10 11 11 11 12 12 12 13 13 13 14 14 15 15 15 16 16 16 17 17 18 18 18 19 19 19 20 20 21 21
/ / 3 3 4 4 5 5 5 6 6 6 7 7 7 8 8 9 9 9 10 10 10 11 11 11 12 12 12 13 13 13 14 14 14 15 15 15 16 16 16 17 17 18 18 18 19 19 19 20
/ / 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12 13 13 13 14 14 14 15 15 15 16 16 16 17 17 17 17 18 18 18
/ 2 3 3 4 4 4 5 5 5 6 6 6 7 7 7 7 8 8 8 9 9 9 10 10 10 10 11 11 11 12 12 12 13 13 13 14 14 14 14 15 15 15 16 16 16 17 17 17 17
/ 2 3 3 4 4 4 4 5 5 5 6 6 6 7 7 7 7 8 8 8 9 9 9 9 10 10 10 11 11 11 11 12 12 12 13 13 13 13 14 14 14 15 15 15 15 16 16 16 17
/ 2 3 3 3 4 4 4 5 5 5 5 6 6 6 7 7 7 7 8 8 8 9 9 9 9 10 10 10 10 11 11 11 11 12 12 12 13 13 13 13 14 14 14 14 15 15 15 15 16
/ 2 3 3 3 4 4 4 4 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8 9 9 9 9 10 10 10 10 11 11 11 11 12 12 12 13 13 13 13 14 14 14 14 15 15 15
/ 2 3 3 3 4 4 4 4 5 6 6 6 6 6 6 6 7 7 7 7 8 8 8 8 9 9 9 9 10 10 10 10 11 11 11 11 12 12 12 12 12 13 13 13 13 14 14 14 14
/ 2 3 3 3 3 4 4 4 5 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8 9 9 9 9 9 10 10 10 10 11 11 11 11 12 12 12 12 12 13 13 13 13 14 14
/ 2 3 3 3 3 4 4 4 4 5 5 5 6 6 6 6 6 6 7 7 7 7 8 8 8 8 8 9 9 9 9 10 10 10 10 10 11 11 11 11 12 12 12 12 13 13 13 13 13
/ 2 3 3 3 3 4 4 4 4 5 5 5 5 5 6 6 6 6 7 7 7 7 7 8 8 8 8 8 9 9 9 9 9 10 10 10 10 11 11 11 11 11 12 12 12 12 12 13 13
/ 2 3 3 3 3 4 4 4 4 4 5 5 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 9 9 9 9 9 10 10 10 10 10 11 11 11 11 11 12 12 12 12 12
/ 2 3 3 3 3 3 4 4 4 4 5 5 5 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 9 9 9 9 9 9 10 10 10 10 10 11 11 11 11 11 12 12 12
/ 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 8 9 9 9 9 9 10 10 10 10 10 11 11 11 11 11 11 12
/ 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 8 9 9 9 9 9 9 10 10 10 10 10 11 11 11 11 11
/ 2 2 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 8 9 9 9 9 9 9 10 10 10 10 10 10 11 11 11
/ 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7 8 8 8 8 8 8 9 9 9 9 9 9 10 10 10 10 10 10 11
/ 2 2 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 9 9 9 9 9 9 10 10 10 10 10 10
/ 2 2 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9 9 9 9 10 10 10 10
/ 2 2 3 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 9 9 9 9 9 9 9 10 10 10
/ 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9 9 9 9 9 10
/ 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9 9 9 9 9
/ 2 2 2 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9 9 9 9
/ 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9 9 9
/ 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 9 9 9
/ 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 9 9
/ 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 9
/ 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8
/ 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8
/ 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8
/ 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 8 8 8 8
/ 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 8 8 8
/ 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 8 8
/ 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 8
/ 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
/ 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7
/ 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7
/ 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7
/ 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7
/ 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7
/ 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7",
                                   sep = " ",
                                   na.strings = "/",
                                   stringsAsFactors = F) )

  table_B3 <- as.matrix(read.table(text =  "/ / / / / / / / / / / / / / / / / / / 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 39 40 41 42 43 44 45 46 47 48 49
/ / / / 5 6 7 8 9 9 10 11 12 13 13 14 15 16 16 17 18 19 20 20 21 22 23 23 24 25 26 26 27 28 29 30 30 31 32 33 33 34 35 36 37 37 38 39 40 40
/ / / 4 5 5 6 7 7 8 9 9 10 11 11 12 12 13 14 14 15 16 16 17 18 18 19 19 20 21 21 22 23 23 24 24 25 26 26 27 28 28 29 30 30 31 31 32 33 33
/ / 3 4 4 5 5 6 6 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 17 17 18 18 19 19 20 20 21 21 22 22 23 23 24 24 25 26 26 27 27 28 28
/ 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 9 10 10 11 11 12 12 13 13 14 14 14 15 15 16 16 17 17 18 18 18 19 19 20 20 21 21 22 22 23 23 23 24 24
/ 2 3 3 4 4 4 5 5 6 6 6 7 7 8 8 8 9 9 10 10 10 11 11 12 12 12 13 13 14 14 14 15 15 16 16 16 17 17 18 18 18 19 19 20 20 20 21 21 21
/ 2 3 3 3 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 11 11 11 12 12 12 13 13 13 14 14 14 15 15 15 16 16 16 17 17 18 18 18 19 19 19
/ 2 3 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12 13 13 13 13 14 14 14 15 15 15 16 16 16 17 17 17 17
/ 2 2 3 3 3 4 4 4 5 5 5 5 6 6 6 7 7 7 7 8 8 8 9 9 9 9 10 10 10 11 11 11 11 12 12 12 13 13 13 13 14 14 14 15 15 16 16 16 16
/ 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 8 8 8 8 9 9 9 9 10 10 10 10 11 11 11 11 12 12 12 12 13 13 13 14 14 14 14 15 15
/ 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8 9 9 9 9 9 10 10 10 10 11 11 11 11 12 12 12 12 13 13 13 13 14 14
/ 2 2 3 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 6 7 7 7 7 8 8 8 8 8 9 9 9 9 10 10 10 10 10 11 11 11 11 12 12 12 12 12 13 13
/ 2 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 9 9 9 9 9 10 10 10 10 11 11 11 11 11 12 12 12 12
/ 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 9 9 9 9 9 10 10 10 10 10 11 11 11 11 11 12
/ 2 2 2 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 9 9 9 10 10 10 10 10 11 11 11
/ 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7 8 8 8 8 8 9 9 9 9 9 9 10 10 10 10 10 10
/ 2 2 2 2 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9 9 9 9 10 10 10
/ 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9 9 9 9 10
/ 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9 9 9
1 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9
1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 9
1 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8 8 8
1 2 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8 8 8 8
1 2 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 8 8
1 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 8
1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7
1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7
1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7
1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 7 7 7
1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 7 7
1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6
1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6
1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6
1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6
1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6
1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6
1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 6 6 6
1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 6 6
1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 6
1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5
1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5
1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5
1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5
1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5
1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5
1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5
1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5
1 1 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5
1 1 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5
1 1 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5",
                                   sep = " ",
                                   na.strings = "/",
                                   stringsAsFactors = F))
  Kc <- array(c(table_B2, table_B3), dim = c(50, 50, 2))
  ############
  # Kc: critical value based on table B-2 and table B-3
  critical_value <- if (alpha == 0.01) {
    Kc[[m, n, 1]]
  } else if (alpha == 0.05) {
    Kc[[m, n, 2]]
  } else {
    stop("Critical value not found, check sample size is between 1-50")
  }

  ########### RESULTS: #############
  if (print) {
    text.a <- "Slippage Test Results:"
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
    if (k > critical_value) {
      cat(row, "\n")
      cat(
        paste(
          "Max background (MBC):",
          maxbg,
          "\t",
          "# sites > MBC:",
          k,
          "\t",
          "Critical value:",
          critical_value,
          "\n",
          "\n"
        )
      )
      cat(paste("DECISION:", "\t","Reject null hypothesis", "\n"))
      cat(row, "\n")
    } else {
      cat(row, "\n")
      cat(paste("Measured Concentrations:", "\n"))
      cat(
        paste(
          "Max background (MBC):",
          maxbg,
          "\t",
          "# sites > MBC:",
          k,
          "\t",
          "Critical value:",
          critical_value,
          "\n",
          "\n"
        )
      )
      cat(paste("DECISION:", "\t","Fail to reject null hypothesis", "\n"))
      cat(row, "\n")
    }
  }

  # Data needed for list output and/or ggplot
  site$measure_notNDs <- ifelse(nd.s == 1, NA, measure.s) # removes concentration if it is a nondetect
  site$measure_max <- ifelse(site$measure_notNDs > maxbg, TRUE, NA) # this is only for detects
  site$measure_max_measure <- ifelse(site$measure_max == TRUE, site$measure_notNDs, NA) # "Number of sites greater than background" with their corresponding measurements

  ########### GGPLOT ###########
  if (plot) {

    # Data needed for the plots
    site$ND_tfs <- ifelse(nd.s == 1, TRUE, NA) # getting rid of zero's and making them NA's for easier usage in ND bar plot
    background$measure_notND <- ifelse(nd.b == 1, NA, measure.b) # removes concentration if it is a nondetect
    background$ND_tf <- ifelse(nd.b == 1, TRUE, NA) # getting rid of zero's and making them NA's for easier usage in ND bar plot


    # 2 ways of dynamically changing the bin size...
    #binw <- round(max(site$measure_max_measure, na.rm=T)/20)
    binw2 <- round(max(range(site$measure_notNDs, na.rm=T),range(background$measure_notND, na.rm=T))/15)

    # Creating histogram
    # site measurement layer
    my_histogram <- ggplot2::ggplot(color = "black") +
      ggplot2::geom_histogram(
        data = site,
        ggplot2::aes(x = measure_notNDs, fill = "measure_notNDs"),
        color = "black",
        alpha = 0.8,
        binwidth = binw2,
        na.rm = TRUE
      ) +
      # background measurement layer
      ggplot2::geom_histogram(
        data = background,
        ggplot2::aes(x = measure_notND, fill = "measure_notND"),
        color = "black",
        alpha = 0.8,
        binwidth = binw2,
        na.rm = TRUE
      ) +
      # this layer will change the colors of the site measurements greater than the maximum background measurement to gold.
      ggplot2::geom_histogram(
        data = site,
        ggplot2::aes(x = measure_max_measure, fill = "measure_max_measure"),
        color = "black",
        alpha = 1,
        binwidth = binw2,
        na.rm = TRUE
      )
    # Customizing histogram more
    my_histogram <- my_histogram +
      ggplot2::theme_classic() +
      ggplot2::xlab("Concentration") +
      ggplot2::ylab(NULL) + # we won't need this because we will use the y-axis from the bar plot later on (for our final figure)
      ggplot2::scale_x_continuous(expand = c(0, 0)) + # remove blank space in x axis
      # legend customization
      ggplot2::scale_fill_manual(
        values = c("#FDE725FF", "#1F968BFF", "#440154FF"),
        name = "Measurement",
        labels = c(
          "Site measure > max(background)\n(P/I > B/G)",
          "Background",
          "Site"
        )
      ) +
      ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.75))

    # Creating bar chart of non-detects
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
      ggplot2::scale_fill_manual(values = c("#1F968BFF", "#440154FF")) + # no gold bars because the gold bars are just special site measurements!
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

    figure_slippage <-
      ggpubr::ggarrange(
        my_bar + ggpubr::rremove("x.text") + ggpubr::rremove("x.ticks"), # removing parts of axis not needed
        my_histogram + ggpubr::rremove("y.axis") + ggpubr::rremove("y.ticks") + ggpubr::rremove("y.text"), # removing parts of axis for cohesiveness
        ncol = 2,
        nrow = 1,
        common.legend = TRUE,
        legend = "top",
        legend.grob = my_legend,
        widths = c(0.15, 1), # relative widths of bar + histogram plots
        align = "h" # align the horizontal axes
      )
    print(figure_slippage)
# Output that is saved if slippage is assigned to an object
    my_list<- list(
    alpha= alpha,
    power=power,
    maxbg = maxbg,
    sitemax.n=k,
    sitemax.m= na.omit(site$measure_max_measure),
    critical.value=critical_value
  )
  my_list
  }

}

