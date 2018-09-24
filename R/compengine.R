

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
  output <- c()
  return(output)
}



#' Points inside a given circular boundary in a 2-d embedding space from software package \code{hctsa}
#'
#' The time lag is set to the first zero crossing of the autocorrelation function.
#'
#' @param y the input time series
#' @param boundary the given circular boundary, setting to 1 or 2 in CompEngine
#' @return Points inside a given circular boundary
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
CO_Embed2_Basic_tau_incircle <- function(y, boundary){
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


FC_LocalSimple_taures <- function(x, )

