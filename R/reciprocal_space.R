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


#' Load anomalous data for a specific chemical element
#' 
#' Returns a dataframe with f' and f'' at various wavelengths
#' for the specific chemical element.
#' 
#' @param chem_el 1- or 2-letters character string. The chemical symbol of
#'    interest.
#' @return A dataframe with 3 columns, the specific wavelength in angstroms
#'    (lambda), f' (f1) and f'' (f2).
#' @examples 
#' # Load anomalous data for Fe
#' ano_Fe <- load_anomalous_data("Fe")
#' print(ano_Fe[1:10,])
#' @export
load_anomalous_data <- function(chem_el)
{
  # Load specific dataset of f' and f"
  idx <- which(names(anomalous_data) == chem_el)
  ad <- anomalous_data[[idx]]
  
  return(ad)
}


#' Plot of absorption curves
#' 
#' Plot f' and f'' absorption curves for the specified chemical element.
#' Curves can be plotted in specified wavelength regions using parameter
#' "zoom".
#' 
#' @param chem_el 1- or 2-letters character string. The chemical symbol of
#'    interest.
#' @param zoom Real vector of length 2. The two values are the extremes of
#'    the wavelength window inside which to plot the two curves. Default is
#'    for both curves to be plotted across the full available range.
#' @return Nothing, but causes a 2D plot to be displayed in a graphical window.
#' @examples 
#' # No zoom
#' plot_absorption_curves("Fe")
#' 
#' # Zoom
#' plot_absorption_curves("Fe",c(1,3))
#' @export
plot_absorption_curves <- function(chem_el,zoom=NULL)
{
  # Load table related to specific chemical element
  ad <- load_anomalous_data(chem_el)
  
  # Limit plot to input range
  if (!is.null(zoom))
  {
    lset <- which(ad$lambda >= zoom[1])
    rset <- which(ad$lambda <= zoom[2])
    jdx <- na.omit(match(lset,rset))
  }
  if (is.null(zoom)) jdx <- 1:length(ad[,1])
  
  # Range
  tmp <- range(ad$f1[jdx],ad$f2[jdx])
  
  # Plot
  plot(ad$lambda[jdx],ad$f1[jdx],type="l",
       xlab=expression(paste(lambda,"(",ring(A),")")),
       ylab=expression(paste(f,"(",e^{-phantom(0)},")")),
       ylim=c(min(tmp),max(tmp)),main=chem_el)
  points(ad$lambda[jdx],ad$f2[jdx],type="l")
}


#' Calculation of structure factors
#' 
#' Calculates structure factors corresponding to one or more 1D Miller indices,
#' given the atomic content of one unit cell. Anomalous scattering can be
#' included using logical flag "anoflag" (default is FALSE). Crystal structures
#' are always considered with no symetry. Thus a 1D structure with the P-1
#' symmetry will have to be expanded first with function "expand_to_cell".
#' 
#' @param h Real numeric. One or more 1D Miller indices.
#' @param a A real number. The unit cell side length.
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
#' @param anoflag Logical variable. If TRUE it forces scattering factors
#'    to include anomalous contributions. As a consequence, Friedel's
#'    pairs will not be equal.
#' @param aK Real numeric. This is a fudge factor included to decrease the
#'    strength of the anomalous contributions. Without aK the strength is too
#'    high for 1D structures, compared to real 3D structures. So aK helps
#'    bringing down the anomalous contribution within the 5%-10% normally
#'    met with large-size structures. The default value is aK=0.3 (anoK is
#'    included as internal data).
#' @param lbda Real numeric. This is the wavelength in angstroms. Its value
#'    is important in relation to anomalous scattering.
#' @param k A real number. It controls the standard deviation of the 
#'    gaussian function describing the atom and, thus, the shape of the
#'    associated peak. The standard deviation sigma is given by:
#'          \code{sigma = k * sqrt(Z)}
#'    The default value is k=0.05 (ksigma is included as internal data).
#' @param f1f2out Logical variable. This variable controls output of a small
#'    table of f' and f'' values for all chemical elements in the structure.
#'    Default is for the table to be printed.
#' @return A vector of complex numbers, the structure factors corresponding
#'    to the Miller indices in input.
#' @examples 
#' # First create the crystal structure (P1)
#' a <- 10
#' vx0 <- c(1,4,6.5)
#' vZ <- c(8,26,6)
#' vB <- c(18,20,17)
#' vocc <- c(1,1,1)
#' lbda <- 1.7   # Close to Fe's absorption
#' 
#' # Miller indices (including DC (h=0) component)
#' hidx <- 0:10
#' 
#' # Now structure factors without anomalous contribution
#' F <- strufac(hidx,a,vx0,vZ,vB,vocc,lbda=lbda)
#' print(length(F))  # Includes DC component (h=0)
#' print(F[1])       # DC component is real
#' print(Mod(F))     # Amplitudes decrease with increasing
#'                   # resolution (Miller indices)
#'                   
#' # Now try with anomalous contributions
#' F <- strufac(hidx,a,vx0,vZ,vB,vocc,lbda=lbda,anoflag=TRUE)
#' print(F[0])  # DC component is not any longer real
#' 
#' @export
strufac <- function(h,a,vx0,vZ,vB,vocc,anoflag=FALSE,
                    aK=anoK,lbda=1.0,k=ksigma,f1f2out=TRUE)
{
  if (length(vx0) != length(vZ)) 
    stop("Atomic center and atomic number lengths differ!")
  if (length(vx0) != length(vB))
    stop("Atomic center and B-factors lengths differ!")
  if (length(vx0) != length(vocc))
    stop("Atomic center and atomic occupanciess lengths differ!")
  
  # If anomalous flag is on collect f' and f" for each 
  # atomic species
  if (anoflag)
  {
    atomic_species <- unique(vZ)
    idx <- match(atomic_species,atoms$Z)
    AnoFrame <- data.frame()
    for (i in idx)
    {
      chem_el <- names(anomalous_data)[i]
      ad <- load_anomalous_data(chem_el)
      tmp <- abs(ad$lambda-lbda)
      jdx <- which(tmp == min(tmp))
      AnoFrame <- rbind(AnoFrame,data.frame(Z=atoms$Z[i],f1=ad$f1[jdx],
                                            f2=ad$f2[jdx]))
    }
    if (f1f2out) print(AnoFrame)
  }
  
  # Quantities used if the anomalous flag is off
  f1 <- 0
  f2 <- 0
  ph_ano <- 0
  
  # Calculate structure factors
  FRe <- rep(0,times=length(h))
  FIm <- rep(0,times=length(h))
  for (j in 1:length(vZ))
  {
    # Scattering factor
    f <- scafac(h,a,vZ[j],vocc[j],vB[j],k)
    
    # Anomalous contribution
    if (anoflag)
    {
      idx <- match(vZ[j],AnoFrame$Z)
      f1 <- aK*AnoFrame$f1[idx]
      f2 <- aK*AnoFrame$f2[idx]
      ph_ano <- 0.5*pi+2*pi*h*vx0[j]/a
    }
    
    # Complex factors
    ph_reg <- 2*pi*h*vx0[j]/a
    regRe <- cos(ph_reg)
    regIm <- sin(ph_reg)
    anoRe <- cos(ph_ano)
    anoIm <- sin(ph_ano)
    
    # Scattering coefficients
    F1 <- f+f1
    F2 <- f2
    
    # Sum of various components
    FRe <- FRe+(F1*regRe+F2*anoRe)
    FIm <- FIm+(F1*regIm+F2*anoIm)
  }
  F <- complex(real=FRe,imaginary=FIm,length.out=length(h))
  
  return(F)
}


#' From structure factors to density
#' 
#' Given a set of complex structure factors corresponding to a set of 1D
#' Miller indices, the length of the 1D unit cell, the set of Miller indices
#' and the number of grid points used to calculate the 1D density, this
#' function calculates the 1D density corresponding to the given structure
#' factors.
#' 
#' @param a A real number. The unit cell side length.
#' @param F A vector of complex numbers. The structure factors corresponding to
#'    the 1D density to be calculated.
#' @param hidx A vector of integer numbers. The set of 1D Miller indices
#'    corresponding to the set of structure factors F.
#' @param N An integer number. The number of grid points used to calculate the
#'    1D density.
#' @return A vector of N real numbers representing the calculated 1D density at
#'    each of the regular N grid points.
#' @examples 
#' # First create the crystal structure (in P1)
#' a <- 10
#' vx0 <- c(1,4,6.5)
#' vZ <- c(8,26,6)
#' vB <- c(18,20,17)
#' vocc <- c(1,1,1)
#' 
#' # Enough Fourier components (Miller indices)
#' hidx <- 0:20
#' 
#' # Compute the structure factors
#' F <- strufac(hidx,a,vx0,vZ,vB,vocc)
#' 
#' # Number of grid points
#' N <- 100
#' 
#' # Density
#' rtmp <- FToRho(a,F,hidx,N)
#' 
#' # Density plot in the unit cell
#' x <- rtmp$x
#' rho <- rtmp$rr
#' plot(x,rho,type="l",xlab="x",ylab=expression(rho))
#' @export
FToRho <- function(a,F,hidx,N)
{
  # Arrays checks
  if (length(F) != length(hidx))
    stop("Arrays do not have same length")
  
  # Check numbers are compatible
  hmax <- max(hidx)
  if (hmax >= 0.5*N) stop("Grid too coarse. Increase N.")
  
  # Initialise fft array
  G <- rep(0,length=N)
  
  # Fill fft array with structure factors
  G[hidx+1] <- F
  if (hidx[1] == 0) 
  {
    kidx <- N-hidx[2:length(hidx)]+1
    G[kidx] <- Conj(F[2:length(F)])
  }
  if (hidx[1] != 0)
  {
    kidx <- N-hidx+1
    G[kidx] <- Conj(F)
  }
  
  # Density via fft
  rr <- Re(fft(G))/a
  
  # Return x array with same size as rr array
  x <- seq(0,a,length=N+1)[1:N]
  
  return(list(x=x,rr=rr,G=G))
}