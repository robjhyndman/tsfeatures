

compengine <- function(x, ...){
  hctsa_autocorrelation(x, ...)
}

hctsa_autocorrelation <- function(x){
  output <- c(CO_FirstMin = CO_FirstMin_ac(x),
              CO_trevnum = CO_trev_1_num(x),
              CO_Embed2_First0_incircle_1 = CO_Embed2_Basic_tau_incircle(x,1),
              CO_Embed2_First0_incircle_2 = CO_Embed2_Basic_tau_incircle(x,2), 
              SB_MotifTwo_mean_hhh = SB_MotifTwo_mean_hhh(x), 
              PH_Walker_sw_propcross = PH_Walker_prop_01_sw_propcross(x))
 return(output)
}

hctsa_predictability <- function(x){
  output <- c(FC_LocalSimple_mean1_taures = FC_LocalSimple_taures(x, "mean"), 
              FC_LocalSimple_lfit_taures = FC_LocalSimple_taures(x, "lfit"), 
              EN_SampEn_1 = EN_SampEn_5_03_sampen1(x))
  return(output)
}

hctsa_stationarity <- function(x){
  output <- c(SY_StdNthDer_1 = SY_StdNthDer_1(x), 
              SY_SpreadRandomLocal_meantaul_50 = SY_SpreadRandomLocal_100_meantaul(x, 50), 
              SY_SpreadRandomLocal_meantaul_ac2 = SY_SpreadRandomLocal_100_meantaul(x, "ac2"))
  return(output)
}

hctsa_distribution <- function(x){
  output <- c(DN_HistogramMode_10 = DN_HistogramMode(x),
              DN_OutlierInclude_mdrmd = DN_OutlierInclude_abs_001_mdrmd(x))
}


# autocorr ----------------------------------------------------------------


#' Points inside a given circular boundary in a 2-d embedding space from software package \code{hctsa}
#'
#' The time lag is set to the first zero crossing of the autocorrelation function.
#'
#' @param y the input time series
#' @param boundary the given circular boundary, setting to 1 or 2 in CompEngine. Default to 1.
#' @return Points inside a given circular boundary
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
CO_Embed2_Basic_tau_incircle <- function(y, boundary = NULL){
  if(is.null(boundary)){
    warning("`CO_Embed2_Basic_tau_incircle()` using `boundary = 1. Set value with `boundary`.")
    boundary <- 1
  }
  tau <- CO_FirstZero_ac(y)
  xt <- y[1:(length(y)-tau)]# part of the time series
  xtp <- y[(1+tau):length(y)]# time-lagged time series
  N <- length(y) - tau# Length of each time series subsegment
  
  # CIRCLES (points inside a given circular boundary)
  return(sum(xtp^2+xt^2 < boundary)/N)
}

#' The first zero crossing of the autocorrelation function from software package \code{hctsa}
#' 
#' Search up to a maximum of the length of the time series
#' 
#' @param y the input time series
#' @return The first zero crossing of the autocorrelation function
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
CO_FirstZero_ac <- function(y){
  N <- length(y)
  corrs <- acf(y, N-1, plot=FALSE)$acf[-1]
  for(tau in 1:(N-1)){
    if(corrs[tau]<0) return(tau) # we know it starts 1, so first negative will be the zero-crossing
  }
  return(N) # If haven't left yet, set output to sample size
}






#' Time of first minimum in the autocorrelation function from software package \code{hctsa}
#' 
#'
#' @param x the input time series
#' @return The lag of the first minimum
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @examples
#' CO_FirstMin_ac(WWWusage)
#' @export
CO_FirstMin_ac <- function(x){
  # hctsa uses autocorr in MatLab to calculate autocorrelation
  N <- length(x)
  # getting acf for all lags
  # possible delay when sample size is too big 
  autoCorr <- numeric(N-1)
  autoCorr[1:N-1] <- stats::acf(x,lag.max = N-1, plot = FALSE)$acf[-1]
  for(i in 1:length(autoCorr)){
    if(is.na(autoCorr[i])){
      warning("No minimum was found.")
      return(NA)
    }
    if(i==2 && autoCorr[2] > autoCorr[1]) {
      return(1)
    } else if(i>2 && autoCorr[i-2] > autoCorr[i-1] && autoCorr[i-1] < autoCorr[i]){
      return(i-1)
    }
  }
  return(N-1)
}


#' Normalized nonlinear autocorrelation, the numerator of the trev function of a time series from software package \code{hctsa}
#' 
#' Calculates the numerator of the trev function, a normalized nonlinear autocorrelation,
#' The time lag is set to 1.
#' 
#'
#' @param y the input time series
#' @return the numerator of the trev function of a time series
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @examples
#' CO_trev_1_num(WWWusage)
#' @export
CO_trev_1_num <- function(y){
  yn <-  y[1:(length(y)-1)]
  yn1  <-  y[2:length(y)]
  mean((yn1-yn)^3)
}

#' Local motifs in a binary symbolization of the time series
#'
#'
#' Coarse-graining is performed. Time-series values above its mean are given 1, 
#' and those below the mean are 0.
#'
#' @param y the input time series
#' @return Entropie of words in the binary alphabet of length 3.  
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @examples
#' SB_MotifTwo_mean_hhh(WWWusage)
#' @export
#'
SB_MotifTwo_mean_hhh <- function(y){
  yBin <- BF_Binarize_mean(y)
  N <- length(yBin)
  if(N<5) warning('Time series too short')
  
  r1 <- yBin == 1
  r0 <- yBin == 0

  r1 <- r1[1:(length(r1)-1)]
  r0 <- r0[1:(length(r0)-1)]
  
  r00 <- r0 & yBin[2:N] == 0
  r01 <- r0 & yBin[2:N] == 1
  r10 <- r1 & yBin[2:N] == 0
  r11 <- r1 & yBin[2:N] == 1

  r00 <- r00[1:(length(r00)-1)]
  r01 <- r01[1:(length(r01)-1)]
  r10 <- r10[1:(length(r10)-1)]
  r11 <- r11[1:(length(r11)-1)]
  
  r000 <- r00 & yBin[3:N] == 0
  r001 <- r00 & yBin[3:N] == 1
  r010 <- r01 & yBin[3:N] == 0
  r011 <- r01 & yBin[3:N] == 1
  r100 <- r10 & yBin[3:N] == 0
  r101 <- r10 & yBin[3:N] == 1
  r110 <- r11 & yBin[3:N] == 0
  r111 <- r11 & yBin[3:N] == 1
  
  out.ddd <- mean(r000)
  out.ddu <- mean(r001)
  out.dud <- mean(r010)
  out.duu <- mean(r011)
  out.udd <- mean(r100)
  out.udu <- mean(r101)
  out.uud <- mean(r110)
  out.uuu <- mean(r111)
  ppp <- c(out.ddd, out.ddu, out.dud, out.duu, out.udd, out.udu, out.uud, out.uuu)
  out.hhh <- f_entropy(ppp)
  return(out.hhh)
}


#' Converts an input vector into a binarized version
#' 
#' 
#' @param y the input time series
#' @return Time-series values above its mean are given 1, and those below the mean are 0.
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export

BF_Binarize_mean <- function(y){
  y <- y-mean(y)
  Y <-  numeric(length(y))
  Y[y > 0] <-  1
  return(Y)
}

f_entropy <- function(x){
  # entropy of a set of counts, log(0)=0
  -sum(x[x>0]*log(x[x>0]))
}



#' Simulates a hypothetical walker moving through the time domain.
#' 
#' The hypothetical particle (or 'walker') moves in response to values of the
#' time series at each point.
#' The walker narrows the gap between its value and that
#' of the time series by 10%.
#' 
#' 
#' @param y the input time series
#' @return fraction of time series length that walker crosses time series
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
#' 
#' 
PH_Walker_prop_01_sw_propcross <- function(y){
  N <- length(y)
  p <- 0.1
#   walker starts at zero and narrows the gap between its position
#   and the time series value at that point by 0.1, to give the value at the subsequent time step
  w <- numeric(N)
  w[1] <- 0 # start at zero
  for(i in 2:N){
    w[i] = w[i-1] + p*(y[i-1]-w[i-1])
  }
  out.sw_propcross <-  sum((w[1:(N-1)]-y[1:(N-1)])*(w[2:N]-y[2:N]) < 0)/(N-1)
  return(out.sw_propcross)
}


# pred --------------------------------------------------------------------


#' The first zero crossing of the autocorrelation function of the residuals from Simple local time-series forecasting
#' 
#' Simple predictors using the past trainLength values of the time series to
#' predict its next value.
#' 
#' @param y the input time series
#' @param forecastMeth the forecasting method, default to \code{mean}. 
#' \code{mean}: local mean prediction using the past trainLength time-series values.  
#' \code{lfit}: local linear prediction using the past trainLength time-series values.  
#' @param trainLength the number of time-series values to use to forecast the next value. 
#' Default to 1 when using method \code{mean} and 3 when using method \code{lfit}.
#' @return The first zero crossing of the autocorrelation function of the residuals
#' @export
FC_LocalSimple_taures <- function(y, forecastMeth = NULL, trainLength = NULL ){
  if(is.null(forecastMeth)) forecastMeth <- "mean"
  if(!forecastMeth %in% c("mean", "lfit")) stop("`FC_LocalSimple_taures`:Unknown forecasting method")
  if(forecastMeth == "mean" && is.null(trainLength)) trainLength <- 1 
  if(forecastMeth == "lfit" && is.null(trainLength)) trainLength <- 3
  lp <- trainLength
  
  N <- length(y)
  evalr <-  (lp+1):N
  if( length(evalr)==0)  stop('Time series too short for forecasting in `FC_LocalSimple_taures`')
  
  res <- numeric(length(evalr))
  if(forecastMeth == "mean"){
    for(i in 1:length(evalr))
      res[i] <- mean(y[(evalr[i]-lp):(evalr[i]-1)]) - y[evalr[i]]
  }
  if(forecastMeth == "lfit"){
    for(i in 1:length(evalr)){
       # Fit linear
      a <- 1:lp
      b <- y[(evalr[i]-lp):(evalr[i]-1)]
      lm.ab <- lm(b~a, data = data.frame(a,b))
      res[i] <- predict(lm.ab, newdata = data.frame(a=lp+1))-y[evalr[i]]
      # p = polyfit((1:lp)',y(evalr(i)-lp:evalr(i)-1),1)
      #       res(i) = polyval(p,lp+1) - y(evalr(i)); % prediction - value
    }
  }
  out.taures <- CO_FirstZero_ac(res)
  return(out.taures)
}



#' First Sample Entropy of a time series
#' 
#' Modified from the Ben Fulcher's \code{EN_SampEn} which uses code from PhysioNet.
#' The publicly-available PhysioNet Matlab code, sampenc (renamed here to
#' RN_sampenc) is available from:
#' http://www.physionet.org/physiotools/sampen/matlab/1.1/sampenc.m
#' 
#' Embedding dimension is set to 5.
#' The threshold is set to 0.3.
#'
#'
#' @param y the input time series
#' @references cf. "Physiological time-series analysis using approximate entropy and sample
#' entropy", J. S. Richman and J. R. Moorman, Am. J. Physiol. Heart Circ.
#' Physiol., 278(6) H2039 (2000)
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
EN_SampEn_5_03_sampen1 <- function(y){
  M <- 5
  r <- 0.3
  sampEn = PN_sampenc(y,M+1,r)
  return(sampEn)
}




#' Sample Entropy
#' 
#' Modified from the Ben Fulcher version of original code sampenc.m from
#' http://physionet.org/physiotools/sampen/
#' http://www.physionet.org/physiotools/sampen/matlab/1.1/sampenc.m
#' Code by DK Lake (dlake@virginia.edu), JR Moorman and Cao Hanqing.
#' 
#' 
#' @references cf. "Physiological time-series analysis using approximate entropy and sample
#' entropy", J. S. Richman and J. R. Moorman, Am. J. Physiol. Heart Circ.
#' Physiol., 278(6) H2039 (2000)
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
PN_sampenc <- function(y,M,r){
  N <- length(y)
  lastrun <- numeric(N) #zeros(1,N)
  run <- numeric(N) #zeros(1,N)
  A <- numeric(M) #zeros(M,1)
  B <- numeric(M) #zeros(M,1)
    # Get counting:
  for(i in 1:(N-1)){ # go through each point in the time series, counting matches
    y1 <- y[i]
    for(jj in 1:(N-i)){ # compare to points through the rest of the time series
      # Compare to future index, j:
      j <- i + jj
      # This future point, j, matches the time-series value at i:
      if (abs(y[j]-y1) < r){
        run[jj] <- lastrun[jj] + 1 # increase run count for this lag
        M1 <- min(M, run[jj])
        for (m in 1:M1){
          A[m] <- A[m] + 1
          if (j < N){
            B[m] <- B[m] + 1
          }
        }
      } else{
        run[jj] <- 0
      }
    }
    for( j in 1:N-i){
      lastrun[j] <- run[j]
    }
  }
      # Calculate for m <- 1
    NN <- N*(N-1)/2
    p <- A[1]/NN
    e <- -log(p)  
    return(e)
}


# stationarity ------------------------------------------------------------





#' Standard deviation of the first derivative of the time series.
#'
#' Modified from \code{SY_StdNthDer} in \code{kctsa}. Based on an idea by Vladimir Vassilevsky.
#' 
#' @param y the input time series
#' @return Standard deviation of the first derivative of the time series.
#' @references cf. http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
SY_StdNthDer_1 <- function(y){
  if(length(y)<2) stop("Time series is too short to compute differences")
  yd <- diff(y)
  return(sd(yd))
}

#'  Bootstrap-based stationarity measure.
#' 
#' 100 time-series segments of length \code{l} are selected at random from the time series and 
#' the mean of the first zero-crossings of the autocorrelation function in each segment is calculated.
#' 
#' 
#' @param y the input time series
#' @param l the length of local time-series segments to analyze as a positive integer. Can also be a specified character string: "ac2": twice the first zero-crossing of the autocorrelation function
#' @return mean of the first zero-crossings of the autocorrelation function
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
SY_SpreadRandomLocal_100_meantaul <- function(y, l = NULL){
  if(is.null(l)) l <- 50
  if(is.character(l) && "ac2" %in% l) l <- 2*CO_FirstZero_ac(y)
  if(!is.numeric(l)) stop("Unknown specifier `l`")
  numSegs  <-  100
  N <- length(y)
  if(l>0.9*N) stop("This time series is too short. Specify proper segment lengrh in `l`")

  qs <- numeric(numSegs)  
  
  for (j in 1:numSegs){
  # pick a range
  # in this implementation, ranges CAN overlap
  ist <- sample(N-1-l,1) # random start point (not exceeding the endpoint)
  ifh <- ist+l-1 # finish index
  rs <- ist:ifh # sample range (from starting to finishing index)
  ysub <- y[rs] # subsection of the time series
  taul <- CO_FirstZero_ac(ysub)
  qs[j] <- taul
  }
  return(mean(qs,na.rm = TRUE))
  
}


# distribution ------------------------------------------------------------


#' Mode of a data vector
#' 
#' Measures the mode of the data vector using histograms with a given number of bins as suggestion.
#' The value calculated is different from \code{kctsa} and \code{CompEngine} as the histogram edges are calculated differently.
#' 
#' @param y the input data vector
#' @param numBins the number of bins to use in the histogram.
#' @return the mode
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
DN_HistogramMode <- function(y, numBins = 10){

  # Compute the histogram from the data:
    if (is.numeric(numBins)){
      histdata <- hist(y,plot = FALSE)
      binCenters <- histdata$mids
  } else {
    stop('Unknown format for numBins')
  }
  # Compute bin centers from bin edges:
    # binCenters <- mean([binEdges(1:end-1) binEdges(2:end)])
  # Mean position of maximums (if multiple):
    out <- mean(binCenters[which.max(histdata$counts)])
 return(out)
}






#' How median depend on distributional outliers.
#'
#' Measures meidan as more and
#' more outliers are included in the calculation according to a specified rule,
#' of outliers being furthest from the mean, greatest positive, or negative
#' deviations.
#'
#' The threshold for including time-series data points in the analysis increases
#' from zero to the maximum deviation, in increments of 0.01*sigma (by default), 
#' where sigma is the standard deviation of the time series.
#'
#' At each threshold,  proportion of time series points
#' included and median are calculated, and outputs from the
#' algorithm measure how these statistical quantities change as more extreme
#' points are included in the calculation.
#' 
#' Outliers are defined as furthest from the mean.
#' 
#' @param y the input time series (ideally z-scored)
#' @return median  of the median of range indices
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
DN_OutlierInclude_abs_001_mdrmd <- function(y){
  if(length(unique(y))==1) stop("The time series is a constant!")
  if(!BF_iszscored(y)) {
    warning('The input time series should be z-scored')
    isd <- sd(y) # Modified to fit the 0.01*sigma increment in discription
    } else isd <- 1
  N <- length(y)
  inc <- 0.01*isd
  thr <- seq(from = 0, to = max(abs(y)), by = inc)
  tot <- N
  if(length(thr) == 0) stop("peculiar time series")
  
  msDt <- numeric(length(thr))
  msDtp <- numeric(length(thr))
  for (i in 1:length(thr)){
    th <- thr[i] # the threshold
    # Construct a time series consisting of inter-event intervals for parts
    # of the time serie exceeding the threshold, th
    r <- which(abs(y) >= th)
    
    Dt_exc <- diff(r)  # Delta t (interval) time series exceeding threshold
    msDt[i] <- median(r)/(N/2)-1
    msDtp[i] <-  length(Dt_exc)/tot*100 
    # this is just really measuring the distribution: 
    # the proportion of possible values
    # that are actually used in
    # calculation
  }
  
  
   # Trim off where the statistic power is lacking: less than 2% of data
   # included
  trimthr <- 2  # percent
  mj <- which(msDtp > trimthr)[length(which(msDtp > trimthr))]
  if (length(mj) != 0){
    msDt <- msDt[1:mj]
    msDtp <- msDtp[1:mj]
    thr <- thr[1:mj]
  } else stop("the statistic power is lacking: less than 2% of data included")
  
  out.mdrmd <-  median(msDt)
  return(out.mdrmd)
}

#' Crude check for whether a data vector is (eps-close to being) z-scored.
#' 
#' Used for displaying warning messages for functions that require z-scored inputs.
#' 
#' @param x the input time series (or any vector)
#' @return a logical with the verdict.
#' #' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
BF_iszscored <- function(x){
  numericThreshold <- 100*.Machine$double.eps
  iszscored <- ((abs(mean(x)) < numericThreshold) && (abs(sd(x)-1) < numericThreshold))
  return(iszscored)
}

