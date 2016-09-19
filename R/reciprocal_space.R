# 
# This file is part of the crone package
#

#' Scattering factor for 1D gaussian atoms
#' 
#' Given unit cell length, atomic number, B factor, occupancy and Miller
#' index h, this function returns the corresponding value of the analytical
#' scattering factor calculated as Fourier transform of the 1D gaussian atom.
#' 
#' @param h Real numeric. One or more 1D Miller indices. h can also have
#'    non-integer values in between integer values. This enables plotting
#'    of scattering curves.
#' @param a Real numeric. Length of 1D unit cell.
#' @param Zj Integer numeric. Atomic number (e.g. Oxygen has Zj <- 8)
#' @param occj Real numeric. Atomic occupancy, a real number between 0 and 1,
#'    where 0 means that the atom is missing in the crystal and 1 mens that
#'    is present in all unit cells of the crystal.
#' @param Bj Real numeric. This is the B factor associated with the thermal
#'    vibration of the atom. It is measured in squared angstroms and it is
#'    equal to 8*pi^2*sigma^2, where sigma is the gaussian atom width.
#' @param k A real number. It controls the standard deviation of the 
#' gaussian function describing the atom and, thus, the shape of the
#' associated peak. The standard deviation sigma is given by:
#'          \code{sigma = k * sqrt(Z)}
#' @return A real numeric. The value of the scattering factor at the sprcified
#'    Miller idex or corresponding real value.
#' @examples 
#' # Values for some Miller indices
#' h <- 0:10
#' a <- 20
#' Zj <- 16
#' Bj <- 18  # Roughly corresponding to sigma 0.23
#' occj <- 1
#' fval <- scafac(h,a,Zj,occj,Bj)
#' plot(h,fval,pch=16,xlab="Miller index",ylab="Scattering factor",
#'      ylim=c(0,16))
#' 
#' # Continuous resolution
#' h <- seq(0,10,length=1000)
#' fval <- scafac(h,a,Zj,occj,Bj)
#' points(h,fval,type="l",col=2)
#' 
#' # Scattering curve for a lighter atom
#' Zj <- 8
#' fval <- scafac(h,a,Zj,occj,Bj)
#' points(h,fval,type="l",col=3)
#' 
#' # Scattering curve for the same atom, just with smaller Bj (colder)
#' Bj <- 10
#' fval <- scafac(h,a,Zj,occj,Bj)
#' points(h,fval,type="l",col=4)
#' @export
scafac <- function(h,a,Zj,occj,Bj=NULL,k=ksigma)
{
  # Contribution for cold atom
  sigma <- k*sqrt(Zj)
  ff <- Zj*occj*exp(-2*pi^2*sigma^2*h^2/a^2)
  
  # Contribution for thermal atom
  if (!is.null(Bj))
  {
    ff <- ff*exp(-0.25*Bj*h^2/a^2)
  }
  
  return(ff)
}