# 
# This file is part of the crone package
#

#' Suggests unit cell side, a, based on atom content
#'
#' The unit cell side is roughly calculated by adding two times the half-width
#' of the widest gaussian atom to the largest inter-atomic distance. The
#' half-width of the largest gaussian is computed as Ma times the gaussian
#' sigma. If the "P-1" symmetry is present, D is doubled.
#' @param vZ A vector of atom Z numbers.
#' @param D A real number. The distance between the two furthest atoms
#' in the cell.
#' @param SG 2-letters character string. Symmetry. There are only two
#'  symmetries possible when working within 1D crystallography, P1 (no
#'  symmetry)and P-1 (inversion through the origin). SG can be either "P1"
#'  or "P-1" for this function.
#' @param k A real number. It controls the standard deviation of the 
#' gaussian function describing the atom and, thus, the shape of the
#' associated peak. The standard deviation sigma is given by:
#'          \code{sigma = k * sqrt(Z)}
#' @param Ma A real number. Each gaussian atom has tails truncated at a
#' distance of Ma * sigma from its peak.
#' @return A real number that suggests a feasible unit cell side 
#' containing all atoms.
#' @examples 
#' # 2 carbon atoms, a sulphur and an oxygen
#' vZ <- c(6,8,16,6)
#' 
#' # Distance between the two carbons is 15 angstroms
#' D <- 15
#' a <- choose_a(vZ,D)
#' @export 
choose_a <- function(vZ,D,SG="P1",k=ksigma,Ma=5)
{
  if (SG != "P-1" & SG != "P1")
  {
    stop("Wrong Space Group!")
  }
  vsigma <- k*sqrt(vZ)
  if (SG == "P-1") D <- 2*D
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
#' @export
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
  y0 <- x0%%a
  
  # Shift input interval to [0,a]
  y <- x%%a
  
  # Value (depends on where y0 is and is a pieced function)
  if (y0 <= 0.5*a)
  {
    cuts <- c(0,y0+0.5*a,a)
    idx <- findInterval(y,cuts,rightmost.closed=TRUE) 
    A <- 2-idx
    B <- idx-1
    rho <- A*K*Z*exp(-(y-y0)^2/(2*sigma^2))+
      B*K*Z*exp(-(y-y0-a)^2/(2*sigma^2))
  }
  if (y0 > 0.5*a)
  {
    cuts <- c(0,y0-0.5*a,a)
    idx <- findInterval(y,cuts,rightmost.closed=TRUE)
    A <- 2-idx
    B <- idx-1
    rho <- A*K*Z*exp(-(y-y0+a)^2/(2*sigma^2))+
      B*K*Z*exp(-(y-y0)^2/(2*sigma^2))
  }
  
  return(rho)
}


#' Expand content of asymmetric unit to whole unit cell
#' 
#' Atom positions, types, B factors and occupancies are duplicated if
#' input space group (SG) is P-1; otherwise they are left untouched
#' (space group P1). Value of the occupancy for special positions is
#' not checked.
#' 
#' @param a Real numeric. Unit cell length in angstroms.
#' @param vx0 Vector of real numerics. Atom positions in the asymmetric
#'  unit.
#' @param vZ Vector of integers. Atomic numbers of all atoms in the
#'  asymmetric unit.
#' @param vB Vector of real numerics. B factors for all atoms in the
#'  asymmetric unit.
#' @param vocc Vector of real numerics. Occupancies (value between 0
#'  and 1) for all atoms in the asymmetric unit. In this function there
#'  is no mechanism to check whether the occupancy is appropriate in case
#'  the atom is at a special position.
#' @param SG 2-letters character string. There are only two symmetries
#'  possible when working within 1D crystallography, P1 (no symmetry)
#'  and P-1 (inversion through the origin). SG can be either "P1" or
#'  "P-1" for this function.
#' @return A named list with the following elements:
#'  \itemize{
#'    \item{vx0. Vector of real numerics indicating the expanded atomic
#'                positions in the unit cell.}
#'    \item{vZ.  Vector of integers indicating the expanded atomic numbers
#'                for all atoms in the unit cell.}
#'    \item{vB.  Vector of real numerics indicating the expanded B factors
#'                for all atoms in the unit cell.}
#'    \item{vocc. Vector of real numerics indicating the expanded occupancies
#'                  for all atoms in the unit cell.}
#'  }
#' @examples 
#' # Asymmetric unit includes 3 atoms between 0 and a/2
#' a <- 10
#' vx0 <- c(1,2,4)
#' vZ <- c(6,8,6)
#' vB <- c(1,1.2,1.1)
#' vocc <- c(1,1,1)
#' ltmp <- expand_to_cell(a,vx0,vZ,vocc,"P-1")
#' print(ltmp)  # Positions, atomic numbers, etc. have doubled
#' 
#' # Nothing changes if SG is "P1"
#' ltmp <- expand_to_cell(a,vx0,vZ,vocc,"P1")
#' print(ltmp)
#' @export
expand_to_cell <- function(a,vx0,vZ,vB,vocc,SG="P1")
{
  if (SG != "P-1" & SG != "P1")
  {
    stop("Wrong Space Group!")
  }
  if (SG == "P-1")
  {
    # Shift to interval 0-a/2 if needed
    vx0[vx0 > a/2] <- (-vx0[vx0 > a/2])%%a
    
    # Expansion
    svx0 <- (-vx0)%%a
    vx0 <- c(vx0,svx0)
    
    # Expand occupancies, Z and U
    vocc <- c(vocc,vocc)
    vZ <- c(vZ,vZ)
    vB <- c(vB,vB)
  }
  
  # Get rid of duplicated atoms (special positions)
  ltmp <- isRoughlyEqual(vx0)
  idx <- ltmp[[1]]
  if (length(idx) > 0)
  {
    vx0 <- vx0[idx]
    vZ <- vZ[idx]
    vB <- vB[idx]
    vocc <- vocc[idx]
  }
  
  # Final sorting according to atom positions
  idx <- order(vx0)
  vx0 <- vx0[idx]
  vZ <- vZ[idx]
  vB <- vB[idx]
  vocc <- vocc[idx]
  
  return(list(vx0=vx0,vZ=vZ,vB=vB,vocc=vocc))
}


#' Reduce content of unit cell to asymmetric unit.
#' 
#' Atom types, B factors and occupancies are restricted to those corresponding
#' to positions in the asymmetric unit, if the input space group (SG) is P-1.
#' Otherwise they are only sorted according to increasing atom positions. If
#' the starting positions do not correspond to a properly symmetry-expanded
#' set, then a warning is issued. If no symmetry is used (P1) the starting set
#' is forced to have values inside the unit cell and it is sorted according
#' to increasing atom positions.
#' 
#' @param a Real numeric. Unit cell length in angstroms.
#' @param vx0 Vector of real numerics. Atom positions in the unit cell.
#' @param vZ Vector of integers. Atomic numbers of all atoms in the
#'  unit cell.
#' @param vB Vector of real numerics. B factors for all atoms in the
#'  unit cell.
#' @param vocc Vector of real numerics. Occupancies (value between 0
#'  and 1) for all atoms in the unit cell. In this function there
#'  is no mechanism to check whether the occupancy is appropriate in case
#'  the atom is at a special position.
#' @param SG 2-letters character string. There are only two symmetries
#'  possible when working within 1D crystallography, P1 (no symmetry)
#'  and P-1 (inversion through the origin). SG can be either "P1" or
#'  "P-1" for this function.
#' @return A named list with the following elements:
#'  \itemize{
#'    \item{vx0. Vector of real numerics indicating the reduced atomic
#'                positions in the asymmetric unit.}
#'    \item{vZ.  Vector of integers indicating the reduced atomic numbers
#'                for all atoms in the asymmetric unit.}
#'    \item{vB.  Vector of real numerics indicating the reduced B factors
#'                for all atoms in the asymmetric unit.}
#'    \item{vocc. Vector of real numerics indicating the reduced occupancies
#'                  for all atoms in the asymmetric unit.}
#'  }
#' @examples 
#' # Asymmetric unit includes 3 atoms between 0 and a/2
#' a <- 10
#' vx0 <- c(1,2,4)
#' vZ <- c(6,8,6)
#' vB <- c(1,1.2,1.1)
#' vocc <- c(1,1,1)
#' ltmp <- expand_to_cell(a,vx0,vZ,vB,vocc,SG="P-1")
#' print(ltmp)  # Positions, atomic numbers, etc. have doubled
#' 
#' # Now these expanded values are used as input for reduction
#' ltmp2 <- reduce_to_asu(a,ltmp$vx0,ltmp$vZ,ltmp$vB,ltmp$vocc,SG="P-1")
#' 
#' # Compare
#' print(vx0)
#' print(ltmp2$vx0)
#' 
#' # SG is "P1"
#' a <- 16
#' vx0 <- c(1,2,7,9,10,17.5)  # Other arrays are same as before
#' ltmp3 <- reduce_to_asu(a,vx0,ltmp$vZ,ltmp$vB,ltmp$vocc)
#' print(ltmp3) # Now all positions are inside the unit cell
#' @export
reduce_to_asu <- function(a,vx0,vZ,vB,vocc,SG="P1")
{
  # Check on arrays length
  if ((length(vZ) != length(vx0)) |
      (length(vB) != length(vx0)) |
      (length(vocc) != length(vx0)))
    stop("Input arrays must all have same size.")
  if (SG != "P-1" & SG != "P1")
  {
    stop("Wrong Space Group!")
  }
  
  # Move all atoms inside the unit cell
  vx0 <- vx0%%a
  
  # Sort according to increasing atom positions
  idx <- order(vx0)
  vx0 <- vx0[idx]
  vZ <- vZ[idx]
  vB <- vB[idx]
  vocc <- vocc[idx]
  
  # Get rid of duplicated atoms (special positions)
  ltmp <- isRoughlyEqual(vx0)
  idx <- ltmp[[1]]
  if (length(idx) > 0)
  {
   vx0 <- vx0[idx]
   vZ <- vZ[idx]
   vB <- vB[idx]
   vocc <- vocc[idx]
  }
  
  # Save for later check
  inivx0 <- vx0
  
  # Reduction to first half of the cell only for P-1
  if (SG == "P-1")
  {
   # Reduction
   tx0 <- (-vx0)%%a
   idx <- which(tx0 <= 0.5*a)
   nvx0 <- tx0[idx]
   nvZ <- vZ[idx]
   nvB <- vB[idx]
   nvocc <- vocc[idx]
   
   # Sort according to increasing position in asu
   idx <- order(nvx0)
   nvx0 <- nvx0[idx]
   nvZ <- nvZ[idx]
   nvB <- nvB[idx]
   nvocc <- nvocc[idx]
   
   # Expand to P1 to double-check
   ltmp <- expand_to_cell(a,nvx0,nvZ,nvB,nvocc,SG="P-1")
   expobj <- sort(unique(ltmp$vx0))
   if (length(expobj) == length(inivx0))
   {
     dd <- abs(expobj-inivx0)
     idx <- which(dd > 0.000001)
     if (length(idx) > 0)
       warning("Starting set of atoms not in positions corresponding to P-1
                symmetry.")
   }
   if (length(expobj) != length(inivx0))
   {
     warning("Starting set of atoms not in positions corresponding to P-1
              symmetry.")
   }
   
   # Once checked output is OK
   vx0 <- nvx0
   vZ <- nvZ
   vB <- nvB
   vocc <- nvocc
  }
  
  return(list(vx0=vx0,vZ=vZ,vB=vB,vocc=vocc))
}