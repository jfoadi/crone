#
# This file is part of the package crone.
#
# It includes all auxiliary files used with the main 1D crystallography
# functions and potentially useful to users as well.
#
# It also includes description for all datasets included in the package

# Description of data included in the package

#' Theoretical scattering factors for all atomic species
#' 
#' This dataset is a list where each element is an atomic species.
#' Each element of the list is a dataframe with 3 columns.
#' 
#' @format A list whose elements are dataframes corresponding to atomic
#'         elements. Each dataframe has the following columns:
#' \describe{
#'    \item{lambda}{Specific value of wavelength in angstroms.}
#'    \item{f1}{Real scattering component.}
#'    \item{f2}{Imaginary scattering component.}
#' }
#'
"anomalous_data"

#' Atom names and atomic number
#' 
#' This is a dataframe including a 2-characters string indicating the atomic
#' element name and an integer indicating the atomic number Z.
#' 
#' @format The dataframe has the following columns:
#' \describe{
#'    \item{anames}{2-character string indicating the atomic name.}
#'    \item{Z}{Integer. The atomic number.}
#' }
"atoms"


# Auxiliary functions

#' Find local maxima in a vector of real values
#'
#' @param x  A vector of real numbers 
#' @return A vector of integers corresponding to local maxima positions 
#' in vector x 
#' @examples 
#' t <- seq(-10,10,length=1000)
#' x <- dnorm(t,mean=3)+dnorm(t,mean=7)
#' yM <- local_maxima(x)
#' @export
local_maxima <- function(x)
{
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]])
  {
    y <- y[-1]
  }

  return(y)
}


#' Heaviside function (step function)
#'
#' @param x  A vector of real numbers.
#' @param x0 A real number. The x value at which the function step occurs.
#' @return One of the two numbers 0 or 1.
#' @import graphics
#' @examples 
#' x <- seq(-3,5,length=1000)
#' x0 <- 1
#' y <- heaviside(x,x0)
#' plot(x,y,type="l")
#' # Step up and step down
#' x1 <- seq(-3,5,length=1000)
#' x10 <- 1
#' y1 <- heaviside(x1,x10)
#' x2 <- seq(1,9,length=1000)
#' x20 <- 5
#' y2 <- heaviside(x2,x20)
#' y2 <- 1-y2
#' plot(x1,y1,type="l",xlim=c(-3,9),xlab="x",ylab="y")
#' points(x2,y2,type="l")
#' @export
heaviside <- function(x,x0=0)
{
  value <- (sign(x-x0)+1)/2
  
  return(value)
}


#' Error function for real values
#'
#' @param x  A vector of real numbers. 
#' @return A real number. The integral of the gaussian, centred on zero and
#'    with standard deviation equal to 1, between 0 and x, multiplied by
#'    2/sqrt(pi).
#' @examples 
#' x <- seq(-1,1,length=1000)
#' y <- erf(x)
#' plot(x,y,type="l")
#' 
#' @import stats
#' 
#' @export
erf <- function(x)
{
  value <- 2 * pnorm(x * sqrt(2)) - 1
  
  return(value)
}

# To check whether elements in a vector have numerical value very close
# to the value of others.

isRoughlyEqual <- function(v,eps=0.000001)
{
  n <- length(v)
  lidx <- list()
  for (i in 2:n)
  {
    idx <- which(abs(v[i-1]-v) < eps & !is.na(v))
    if (length(idx) > 0) 
    {
      lidx <- c(lidx,list(idx))
      v[idx] <- NA
    }
  }
  if (!is.na(v[n])) lidx <- c(lidx,list(n))
  
  # First element of all vectors belonging to list are the unique ones
  vidx <- c()
  for (i in 1:length(lidx)) vidx <- c(vidx,lidx[[i]][1])
  
  return(list(vidx=vidx,lidx=lidx))
}