% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thiyanga.R
\name{holt_parameters}
\alias{holt_parameters}
\title{Parameter estimates of Holt's linear trend method}
\usage{
holt_parameters(x)
}
\arguments{
\item{x}{a univariate time series}
}
\value{
A vector of 2 values: alpha, beta
}
\description{
Estimate the smoothing parameter for the level-alpha and
the smoothing parameter for the trend-beta
}
\author{
Thiyanga Talagala
}