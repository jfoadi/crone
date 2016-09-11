# 
# This file is part of the crone package
#

#' Suggests unit cell side, a, based on atom content
#'
#' @param vZ A vector of atom Z numbers.
#' @param D A real number. The distance between the two furthest atoms
#' in the cell.
#' @param k A real number. It controls the standard deviation of the 
#' gaussian function describing the atom and, thus, the shape of the
#' associated peak. The standard deviation sigma is given by:
#'          \code{sigma = k * sqrt(Z)}
#' @param Ma A real number. Each gaussian atom has tails truncated at a
#' distance of Ma * sigma from its peak.
#' @return A real number that suggests a feasible unit cell side 
#' containing all atoms
#' @examples 
#' # 2 carbon atoms, a sulphur and an oxygen
#' vZ <- c(6,8,16,6)
#' 
#' # Distance between the two carbons is 15 angstroms
#' D <- 15
#' a <- choose_a(vZ,D)
#' @export 
choose_a <- function(vZ,D,k=ksigma,Ma=5)
{
  vsigma <- k*sqrt(vZ)
  a <- D+2*Ma*max(vsigma)
  
  return(a)
}

#' Constant normalizing wrapped gaussian
#' 
#' @param sigma A real number. Is the standard deviation of the wrapped
#' gaussian.
#' @param a A real number. The unit cell side length.
#' @return A real number, the multiplicative constant normalizing the
#' wrapped gaussian atom so that the area under the curve is equal to 1.
#' @examples 
#' Z <- 16  # Sulphur atom
#' sigma <- 0.05*sqrt(Z)
#' #' @export
Kgauss <- function(sigma,a)
{
  K <- (1/(sigma*sqrt(2*pi)))*(1/erf(a/(2*sqrt(2)*sigma)))
  
  return(K)
}

#' Cold gaussian atom
#' 
#' @param x  Point in the 1D cell at which this function is calculated.
#' @param a  A real number. The width of the unit cell in which the gaussian
#' atom is placed.
#' @param  x0  A real number. The point corresponding to the atom's peak.
#' @param Z  An integer number. Z is the atomic number of the atom (Z(H)=1,
#' Z(He)=2,Z(Li)=3,Z(B)=4, etc).
#' @param k  A real number. It controls the standard deviation of the 
#' gaussian function describing the atom and, thus, the shape of the
#' associated peak. The standard deviation sigma is given by:
#'          \code{sigma = k * sqrt(Z)}
#' 
#' @return A vector of length equal to the length of vector x, with
#' values equal to the evaluated gaussian atom.
#' @examples 
#' # Carbon gaussian atom in the middle of a cell
#' a <- 10
#' x0 <- 5
#' Z <- 6
#' x <- seq(0,a,length=1000)
#' rho <- cold_atom_gauss(x,a,x0,Z)
#' plot(x,rho,type="l",xlab="x",ylab=expression(rho))
#' @export
cold_atom_gauss <- function(x,a,x0=0,Z=1,k=ksigma)
{
  # Find sigma
  sigma <- k*sqrt(Z)
  
  # Integration constant
  K <- Kgauss(sigma,a)
  
  # Shift x0 if outside interval [0,a]
  y0 <- x0
  if (x0 < 0 | x0 > a)
  {
    y0 <- x0%%a
  }
  
  # Shift input interval to [0,a]
  shift <- x[1]
  y <- x-shift
  
  # Value (depends on where x0 is and is a pieced function)
  if (y0 <= 0.5*a)
  {
    cuts <- c(0,y0+0.5*a,a)
    idx <- findInterval(x,cuts,rightmost.closed=TRUE) 
    A <- 2-idx
    B <- idx-1
    rho <- A*K*Z*exp(-(x-y0)^2/(2*sigma^2))+
      B*K*Z*exp(-(x-y0-a)^2/(2*sigma^2))
  }
  if (y0 > 0.5*a)
  {
    cuts <- c(0,y0-0.5*a,a)
    idx <- findInterval(x,cuts,rightmost.closed=TRUE)
    A <- 2-idx
    B <- idx-1
    rho <- A*K*Z*exp(-(x-y0+a)^2/(2*sigma^2))+
      B*K*Z*exp(-(x-y0)^2/(2*sigma^2))
  }
  
  return(rho)
}