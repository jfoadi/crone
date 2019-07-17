# 
# This file is part of the crone package
#

#' Suggests unit cell side, a, based on atom content
#'
#' The unit cell side is roughly calculated by adding two times the half-width
#' of the widest gaussian atom to the largest inter-atomic distance. The
#' half-width of the largest gaussian is computed as Ma times the gaussian
#' sigma. If the "P-1" symmetry is present, D is doubled.
#' @param Z A vector of atom Z numbers.
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
#' Z <- c(6,8,16,6)
#' 
#' # Distance between the two carbons is 15 angstroms
#' D <- 15
#' a <- choose_a(Z,D)
#' @export 
choose_a <- function(Z,D,SG="P1",k=ksigma,Ma=5)
{
  if (SG != "P-1" & SG != "P1")
  {
    stop("Wrong Space Group!")
  }
  vsigma <- k*sqrt(Z)
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
#' a <- 15  # Unit cell size
#' kk <- kgauss(sigma,a)
#' print(kk)
#' 
#' @export
kgauss <- function(sigma,a)
{
  K <- (1/(sigma*sqrt(2*pi)))*(1/erf(a/(2*sqrt(2)*sigma)))
  
  return(K)
}

#' Gaussian atom
#' 
#' @param x  Point in the 1D cell at which this function is calculated.
#' @param a  A real number. The width of the unit cell in which the gaussian
#' atom is placed.
#' @param  x0  A real number. The point corresponding to the atom's peak.
#' @param Z  An integer number. Z is the atomic number of the atom (Z(H)=1,
#' Z(He)=2,Z(Li)=3,Z(B)=4, etc).
#' @param B A real number. This is the B factor characterizing the atom's
#'    thermal agitation. It is given as B=8*pi^2*U, where U is the variance
#'    of the position of the atoms' nucleus around the equilibrium position.
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
#' rho <- atom_gauss(x,a,x0,Z)
#' plot(x,rho,type="l",xlab="x",ylab=expression(rho))
#' @export
atom_gauss <- function(x,a,x0=0,Z=1,B=0,k=ksigma)
{
  # Find sigma of cold atom
  sigma0 <- k*sqrt(Z)
  
  # Variance of atomic thermal displacement
  U <- B/(8*pi^2)
  
  # Sigma of thermal atom to take B into account
  sigma <- sqrt(U+sigma0^2)
  
  # Integration constant (including thermal correction factor)
  K <- kgauss(sigma0,a)*sigma0/sigma
  
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
#' barely checked for values outside [0,1] range. Extra-care needed.
#' 
#' @param sdata A named list, normally obtained through the use of
#'  function \code{\link{read_x}}. The list names correspond to 
#'  different object types:
#'  \itemize{
#'    \item{a.     Real numeric. The size of the unit cell.}
#'    \item{SG.    Character string. Space group symbol; either "P1" 
#'                 or "P-1"}
#'    \item{x0.    Vector of real numerics indicating the expanded atomic
#'                 positions in the unit cell.}
#'    \item{Z.     Vector of integers indicating the expanded 
#'                 atomic numbers for all atoms in the unit cell.}
#'    \item{B.    Vector of real numerics indicating the expanded 
#'                B factors for all atoms in the unit cell.}
#'    \item{occ.  Vector of real numerics indicating the expanded 
#'                occupancies for all atoms in the unit cell.}
#'  }
#' @param SG 2-letters character string. There are only two symmetries
#'  possible when working within 1D crystallography, P1 (no symmetry)
#'  and P-1 (inversion through the origin). SG can be either "P1" or
#'  "P-1" for this function. Default is NULL, in which case the space
#'  group is assumed to be equal to the one of the input structure.
#' @return A named list with the following elements:
#'  \itemize{
#'    \item{a.   Real numeric. The size of the unit cell.}
#'    \item{SG.  Character string. Name of the space group, either "P1"
#'               or "P-1".}
#'    \item{x0. Vector of real numerics indicating the expanded atomic
#'                positions in the unit cell.}
#'    \item{Z.  Vector of integers indicating the expanded atomic numbers
#'                for all atoms in the unit cell.}
#'    \item{B.  Vector of real numerics indicating the expanded B factors
#'                for all atoms in the unit cell.}
#'    \item{occ. Vector of real numerics indicating the expanded occupancies
#'                  for all atoms in the unit cell.}
#'  }
#' @examples 
#' # Asymmetric unit includes 3 atoms between 0 and a/2
#' a <- 10
#' SG <- "P-1"
#' x0 <- c(1,2,4)
#' Z <- c(6,8,6)
#' B <- c(1,1.2,1.1)
#' occ <- c(1,1,0.8)
#' sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' ltmp <- expand_to_cell(sdata)
#' print(ltmp)  # Positions, atomic numbers, etc. have doubled
#' 
#' # Nothing changes if imposed SG is "P1" (but you get a warning!)
#' ltmp <- expand_to_cell(sdata,SG="P1")
#' print(ltmp)
#' @export
expand_to_cell <- function(sdata,SG=NULL)
{
  # Check input object is the right one
  in_names <- names(sdata)
  if (length(in_names) != 6) stop("Wrong input object.")
  if (in_names[1] != "a" | in_names[2] != "SG" | in_names[3] != "x0" |
      in_names[4] != "Z" | in_names[5] != "B" | in_names[6] != "occ")
  {
    stop("Wrong input object.")
  }
  
  # Check on arrays length
  if ((length(sdata$Z) != length(sdata$x0)) |
      (length(sdata$B) != length(sdata$x0)) |
      (length(sdata$occ) != length(sdata$x0)))
  {
    stop("Last 4 arrays in input object must all have same size.")
  }
  
  # Data Ok. Carry on
  a <- sdata$a
  vx0 <- sdata$x0
  vZ <- sdata$Z
  vB <- sdata$B
  vocc <- sdata$occ
  if (is.null(SG)) {
    SG <- sdata$SG
  } else {
    if (SG != sdata$SG)
    {
      warning("Chosen space group different from input data one.")
    }
  }
  
  # Only two space groups in 1D
  if (SG != "P-1" & SG != "P1")
  {
    stop("Wrong Space Group!")
  }
  
  # Shift atoms inside unit cell
  vx0 <- vx0%%a
  
  # Eliminate overlapping atoms (no change of occupancy)
  ltmp <- isRoughlyEqual(vx0)
  idx <- ltmp[[1]]
  vx0 <- vx0[idx]
  vZ <- vZ[idx]
  vB <- vB[idx]
  vocc <- vocc[idx]
  
  # Specific to space group P-1
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
  
  # Get rid of duplicated atoms
  ltmp <- isRoughlyEqual(vx0)
  idx <- ltmp[[1]]
  if (length(idx) > 0)
  {
    vx0 <- vx0[idx]
    vZ <- vZ[idx]
    vB <- vB[idx]
    vocc <- vocc[idx]
  }
  
  # Change occupancy for atoms in special position
  idx <- which(abs(vx0) < 0.000001 | abs(vx0-a) < 0.000001 |
                 abs(vx0-0.5*a) < 0.000001)
  vocc[idx] <- vocc[idx]*2
  
  # Check abnormal values of occupancy
  idx <- which(vocc < 0 | vocc > 1)
  if (length(idx) > 0)
    warning("Out-of-scale occupancies produced!")
  
  # Final sorting according to atom positions
  idx <- order(vx0)
  vx0 <- vx0[idx]
  vZ <- vZ[idx]
  vB <- vB[idx]
  vocc <- vocc[idx]
  
  return(list(a=a,SG=SG,x0=vx0,Z=vZ,B=vB,occ=vocc))
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
#' @param adata A named list, normally obtained through the use of
#'  function \code{\link{read_x}}. The list names correspond to 
#'  different object types:
#'  \itemize{
#'    \item{a.     Real numeric. The size of the unit cell.}
#'    \item{SG.    Character string. Space group symbol; either "P1" 
#'                 or "P-1"}
#'    \item{x0.    Vector of real numerics indicating the expanded atomic
#'                 positions in the unit cell.}
#'    \item{Z.     Vector of integers indicating the expanded 
#'                 atomic numbers for all atoms in the unit cell.}
#'    \item{B.    Vector of real numerics indicating the expanded 
#'                B factors for all atoms in the unit cell.}
#'    \item{occ.  Vector of real numerics indicating the expanded 
#'                occupancies for all atoms in the unit cell.}
#'  }
#' @param SG 2-letters character string. There are only two symmetries
#'  possible when working within 1D crystallography, P1 (no symmetry)
#'  and P-1 (inversion through the origin). SG can be either "P1" or
#'  "P-1" for this function. Default is NULL, in which case the space
#'  group is assumed to be equal to the one of the input structure.
#' @return A named list with the following elements:
#'  \itemize{
#'    \item{x0. Vector of real numerics indicating the reduced atomic
#'                positions in the asymmetric unit.}
#'    \item{Z.  Vector of integers indicating the reduced atomic numbers
#'                for all atoms in the asymmetric unit.}
#'    \item{B.  Vector of real numerics indicating the reduced B factors
#'                for all atoms in the asymmetric unit.}
#'    \item{occ. Vector of real numerics indicating the reduced occupancies
#'                  for all atoms in the asymmetric unit.}
#'  }
#' @examples 
#' # Asymmetric unit includes 4 atoms between 0 and a/2, 1 in special position
#' a <- 10
#' SG <- "P-1"
#' x0 <- c(1,2,4,5)
#' Z <- c(6,8,6,16)
#' B <- c(1,1.2,1.1,0.8)
#' occ <- c(1,1,1,0.5)
#' adata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' ltmp <- expand_to_cell(adata)
#' print(ltmp)  # Positions, atomic numbers, etc. have doubled
#' 
#' # Now these expanded values are used as input for reduction
#' ltmp2 <- reduce_to_asu(ltmp)
#' 
#' # Compare
#' print(x0)
#' print(ltmp2$x0)
#' 
#' # SG is "P1"
#' a <- 16
#' SG <- "P1"
#' x0 <- c(1,2,7,9,12,16)
#' Z <- c(6,6,8,8,7,6)
#' B <- c(0,0,0,0,0,0)
#' occ <- c(1,1,1,1,1,1)
#' adata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' ltmp3 <- reduce_to_asu(adata)
#' print(ltmp3) # Now all positions are inside the unit cell
#' @export
reduce_to_asu <- function(adata,SG=NULL)
{
  # Check input object is the right one
  in_names <- names(adata)
  if (length(in_names) != 6) stop("Wrong input object.")
  if (in_names[1] != "a" | in_names[2] != "SG" | in_names[3] != "x0" |
      in_names[4] != "Z" | in_names[5] != "B" | in_names[6] != "occ")
  {
    stop("Wrong input object.")
  }
  
  # Check on arrays length
  if ((length(adata$Z) != length(adata$x0)) |
      (length(adata$B) != length(adata$x0)) |
      (length(adata$occ) != length(adata$x0)))
  {
    stop("Last 4 arrays in input object must all have same size.")
  }
  
  # Data Ok. Carry on
  a <- adata$a
  vx0 <- adata$x0
  vZ <- adata$Z
  vB <- adata$B
  vocc <- adata$occ
  if (is.null(SG)) {
    SG <- adata$SG
  } else {
    if (SG != adata$SG)
    {
      warning("Chosen space group different from input data one.")
    }
  }
  
  # Only two space groups in 1D
  if (SG != "P-1" & SG != "P1")
  {
    stop("Wrong Space Group!")
  }
  
  # Move all atoms inside the unit cell
  vx0 <- vx0%%a
  
  # Eliminate overlapping atoms (no change of occupancy)
  ltmp <- isRoughlyEqual(vx0)
  idx <- ltmp[[1]]
  vx0 <- vx0[idx]
  vZ <- vZ[idx]
  vB <- vB[idx]
  vocc <- vocc[idx]
  
  # Sort according to increasing atom positions
  idx <- order(vx0)
  vx0 <- vx0[idx]
  vZ <- vZ[idx]
  vB <- vB[idx]
  vocc <- vocc[idx]
  
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
   
   # Change occupancy for atoms in special position
   idx <- which(abs(nvx0) < 0.000001 | abs(nvx0-0.5*a) < 0.000001)
   nvocc[idx] <- nvocc[idx]*0.5
   
   # Expand to P1 to double-check
   adata <- standardise_sdata(a,SG="P-1",nvx0,nvZ,nvB,nvocc)
   celltmp <- expand_to_cell(adata)
   ltmp <- isRoughlyEqual(celltmp$x0)
   idx <- ltmp[[1]]
   expobj <- celltmp$x0[idx]
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
  
  return(list(a=a,SG=SG,x0=vx0,Z=vZ,B=vB,occ=vocc))
}


#' Structure of gaussian atoms
#' 
#' Structure formed by all gaussian atoms in the unit cell. Positions, 
#' atomic numbers and thermal factors are given by vectors of a same 
#' length. Each atom forming the structure is also characterised by a given 
#' occupancy (between 0 and 1).
#' 
#' @param sdata A named list, normally obtained through the use of
#'  function \code{\link{read_x}} or \code{\link{standardise_sdata}}. 
#'  The list names correspond to  different object types:
#'  \itemize{
#'    \item{a.     Real numeric. The size of the unit cell.}
#'    \item{SG.    Character string. Space group symbol; either "P1" 
#'                 or "P-1"}
#'    \item{x0.    Vector of real numerics indicating the expanded atomic
#'                 positions in the unit cell.}
#'    \item{Z.     Vector of integers indicating the expanded 
#'                 atomic numbers for all atoms in the unit cell.}
#'    \item{B.    Vector of real numerics indicating the expanded 
#'                B factors for all atoms in the unit cell.}
#'    \item{occ.  Vector of real numerics indicating the expanded 
#'                occupancies for all atoms in the unit cell.}
#'  }
#' @param x  Point in the 1D cell at which this function is calculated.
#'  Default is NULL, in which case a grid is set up internally.
#' @param N Integer. Number of points in the regular grid, if the grid
#'  is not provided directly.
#' @param k  A real number. It controls the standard deviation of the 
#' gaussian function describing the atom and, thus, the shape of the
#' associated peak. The standard deviation sigma is given by:
#'          \code{sigma = k * sqrt(Z)}
#' 
#' @return A named list of length 2: x is the grid (either input by
#'  user or set up internally), rr is a vector of length equal to the 
#'  length of vector x, with values equal to the evaluated gaussian atoms.
#'
#' @examples 
#' # Cell, atom types, positions and B factors
#' a <- 10
#' SG <- "P1"
#' x0 <- c(2,5,7)
#' Z <- c(6,16,8)
#' B <- c(0,0,0)
#' 
#' # All occupancies to 1
#' occ <- c(1,1,1)
#' 
#' # Standard data format
#' sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' 
#' # Grid for unit cell
#' x <- seq(0,a,length=1000)
#' 
#' # Structure density
#' rtmp <- structure_gauss(sdata,x)
#' plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho))
#' 
#' # Now reduce occupancy of sulphur
#' occ[2] <- 0.5
#' sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' rtmp <- structure_gauss(sdata,x)
#' points(rtmp$x,rtmp$rr,type="l",col=2)
#' 
#' # Increase temperature of oxygen
#' B[3] <- 10
#' sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' rtmp <- structure_gauss(sdata,x)
#' points(rtmp$x,rtmp$rr,type="l",col=3)
#' 
#' @export
structure_gauss <- function(sdata,x=NULL,N=NULL,k=ksigma)
{
  # Check input object is the right one
  in_names <- names(sdata)
  if (length(in_names) != 6) stop("Wrong input object.")
  if (in_names[1] != "a" | in_names[2] != "SG" | in_names[3] != "x0" |
      in_names[4] != "Z" | in_names[5] != "B" | in_names[6] != "occ")
  {
    stop("Wrong input object.")
  }
  
  # Check on arrays length
  if ((length(sdata$Z) != length(sdata$x0)) |
      (length(sdata$B) != length(sdata$x0)) |
      (length(sdata$occ) != length(sdata$x0)))
  {
    stop("Last 4 arrays in input object must all have same size.")
  }
  
  # Expand to P1
  sdata <- expand_to_cell(sdata)
  
  # Data Ok. Carry on
  a <- sdata$a
  vx0 <- sdata$x0
  vZ <- sdata$Z
  vB <- sdata$B
  vocc <- sdata$occ
  
  # If x not set up, create it
  if (is.null(x) & !is.null(N))
  {
    x <- seq(0,a,length=N)
  }
  if (is.null(x) & is.null(N))
  {
    N <- length(vx0)*200
    x <- seq(0,a,length=N)
  }
  
  # Initialise density
  rho <- rep(0,times=length(x))
  
  # Cycle through all atoms to increment density
  for (i in 1:length(vx0))
  {
    rho <- rho+vocc[i]*atom_gauss(x,a,vx0[i],vZ[i],vB[i],k)
  }
  
  return(list(x=x,rr=rho))
}