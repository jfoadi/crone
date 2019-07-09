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
#' Calculates structure factors corresponding to one or more 1D Miller 
#' indices, given the atomic content of one unit cell. Anomalous scattering 
#' can be included using logical flag "anoflag" (default is FALSE). Crystal 
#' structures are always considered with no symetry. Thus a 1D structure 
#' with the P-1 symmetry will have to be expanded first with function 
#' "expand_to_cell".
#' 
#' @param hidx Real numeric. One or more 1D Miller indices.
#' @param sdata A named list, normally obtained through the use of
#'  functions \code{\link{read_x}} or \code{\link{standardise_sdata}}. 
#'  The list names correspond to different object types:
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
#' @return A named list with two vectors of real numbers, the structure 
#'  factors amplitudes, Fmod, and phases, Fpha, corresponding to the Miller 
#'  indices in input.
#' @examples 
#' # First create the crystal structure (P1)
#' a <- 10
#' SG <- "P1"
#' x0 <- c(1,4,6.5)
#' Z <- c(8,26,6)
#' B <- c(18,20,17)
#' occ <- c(1,1,1)
#' sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' lbda <- 1.7   # Close to Fe's absorption
#' 
#' # Miller indices (including DC (h=0) component)
#' hidx <- 0:10
#' 
#' # Now structure factors without anomalous contribution
#' ftmp <- strufac(hidx,sdata,lbda=lbda)
#' print(length(ftmp$Fmod))  # Includes DC component (h=0)
#' print(ftmp)          # Amplitudes decrease with increasing
#'                      # resolution (Miller indices)
#'                   
#' # Now try with anomalous contributions
#' ftmp <- strufac(hidx,sdata,lbda=lbda,anoflag=TRUE)
#' print(ftmp)  # DC component is not any longer real
#' 
#' @export
strufac <- function(hidx,sdata,anoflag=FALSE,
                    aK=anoK,lbda=1.0,k=ksigma,f1f2out=TRUE)
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
  
  # Data Ok. Carry on.
  
  # Expand to P1
  sdata <- expand_to_cell(sdata)
  
  # Copy arrays
  a <- sdata$a
  vx0 <- sdata$x0
  vZ <- sdata$Z
  vB <- sdata$B
  vocc <- sdata$occ
  
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
  FRe <- rep(0,times=length(hidx))
  FIm <- rep(0,times=length(hidx))
  for (j in 1:length(vZ))
  {
    # Scattering factor
    f <- scafac(hidx,a,vZ[j],vocc[j],vB[j],k)
    
    # Anomalous contribution
    if (anoflag)
    {
      idx <- match(vZ[j],AnoFrame$Z)
      f1 <- aK*AnoFrame$f1[idx]
      f2 <- aK*AnoFrame$f2[idx]
      ph_ano <- 0.5*pi+2*pi*hidx*vx0[j]/a
    }
    
    # Complex factors
    ph_reg <- 2*pi*hidx*vx0[j]/a
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
  
  # Amplitudes and phases
  F <- complex(real=FRe,imaginary=FIm,length.out=length(hidx))
  Fmod <- Mod(F)
  Fpha <- Arg(F)*180/pi
  idx <- which(abs(Im(F)) < 0.000001 & abs(Re(F)) >= 0.000001 & Re(F) > 0)
  Fpha[idx] <- 0.0
  idx <- which(abs(Im(F)) < 0.000001 & abs(Re(F)) >= 0.000001 & Re(F) < 0)
  Fpha[idx] <- 180.0
  
  return(list(Fmod=Fmod,Fpha=Fpha))
}


#' From structure factors to density using Fourier synthesis
#' 
#' Given a set of structure factors, separately as a vector of amplitudes 
#' and a vector of phases in degrees, corresponding to a set of 1D
#' Miller indices, the length of the 1D unit cell, the set of Miller indices
#' and the number of grid points used to calculate the 1D density, this
#' function calculates the 1D density corresponding to the given structure
#' factors.
#' 
#' @param a A real number. The unit cell side length.
#' @param Fmod A vector of real numbers. The structure factors' amplitudes
#'  corresponding to the 1D density to be calculated.
#' @param Fpha A vector of real numbers. The structure factors' phases
#'  (in degrees) corresponding to the 1D density to be calculated.
#' @param hidx A vector of integer numbers. The set of 1D Miller indices
#'    corresponding to the set of structure factors F.
#' @param N An integer number. The number of grid points used to calculate the
#'    1D density.
#' @return A vector of N real numbers representing the calculated 1D density at
#'    each of the regular N grid points.
#' @examples 
#' # First create the crystal structure (in P1)
#' a <- 10
#' SG <- "P1"
#' x0 <- c(1,4,6.5)
#' Z <- c(8,26,6)
#' B <- c(18,20,17)
#' occ <- c(1,1,1)
#' sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' 
#' # Enough Fourier components (Miller indices)
#' hidx <- 0:20
#' 
#' # Compute the structure factors
#' ftmp <- strufac(hidx,sdata)
#' 
#' # Number of grid points
#' N <- 1000
#' 
#' # Density
#' rtmp <- fousynth(a,ftmp$Fmod,ftmp$Fpha,hidx,N)
#' 
#' # Density plot in the unit cell
#' x <- rtmp$x
#' rho <- rtmp$rr
#' plot(x,rho,type="l",xlab="x",ylab=expression(rho))
#' @export
fousynth <- function(a,Fmod,Fpha,hidx,N)
{
  # Arrays checks
  if (length(Fmod) != length(hidx) | length(Fpha) != length(hidx))
    stop("Arrays do not have same length")
  
  # Merge amplitudes and phases in complex values
  F <- complex(modulus=Fmod,argument=Fpha*pi/180,length.out=length(hidx))
  
  # Check numbers are compatible
  hmax <- max(abs(hidx))
  if (hmax >= 0.5*N) stop("Grid too coarse. Increase N.")
  
  # Prepare appropriate Miller indices and structure factors arrays,
  # if needed (averaging)
  nhidx <- 0:hmax
  nF <- rep(as.complex(0),length=(hmax+1))
  idx0 <- which(hidx == 0)
  if (length(idx0) > 0) nF[1] <- F[idx0]
  for (h in 1:hmax)
  {
    k <- h+1
    idxP <- which(hidx == h)
    idxM <- which(hidx == -h)
    if (length(idxP) > 0 & length(idxM) > 0)
    {
      nF[k] <- 0.5*(F[idxP]+Conj(F[idxM]))
    }
    if (length(idxP) > 0 & length(idxM) == 0)
    {
      nF[k] <- F[idxP]
    }
    if (length(idxP) == 0 & length(idxM) > 0)
    {
      nF[k] <- Conj(F[idxM])
    }
  }
  
  # Initialise fft array
  G <- rep(0,length=N)
  
  # Fill fft array with structure factors
  G[nhidx+1] <- nF
  if (nhidx[1] == 0) 
    kidx <- N-nhidx[2:length(nhidx)]+1
  G[kidx] <- Conj(nF[2:length(nF)])
  
  # Density via fft
  rr <- Re(fft(G))/a
  
  # Return x array with same size as rr array
  x <- seq(0,a,length=N+1)[1:N]
  
  return(list(x=x,rr=rr,G=G))
}


#' From density to structure factors using inverse Fourier synthesis
#' 
#' Given a density as vector calculated in N grid points, the unit cell
#' size and an array of Miller indices hidx, this function calculates
#' amplitudes and phases of the structure factors corresponding to this
#' density, via inverse Fourier transform.
#' 
#' @param a A real number. The unit cell side length.
#' @param rho A vector of N real numbers representing the 1D density at
#'    each of the regular N grid point.
#' @param hidx A vector of integer numbers. The set of 1D Miller indices
#'    corresponding to the set of structure factors F, to be calculated.
#' @return A vector of N real numbers representing the calculated 1D 
#'    density at each one of the regular N grid points.
#' @examples 
#' # First create the crystal structure (in P1)
#' a <- 10
#' SG <- "P1"
#' x0 <- c(1,4,6.5)
#' Z <- c(8,26,6)
#' B <- c(18,20,17)
#' occ <- c(1,1,1)
#' sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' 
#' # 10 Miller indices plus DC component
#' hidx <- 0:10
#' 
#' # Compute structure factors
#' ftmp1 <- strufac(hidx,sdata)
#' 
#' # Number of grid points
#' N <- 1000
#' 
#' # Density
#' rtmp <- fousynth(a,ftmp1$Fmod,ftmp1$Fpha,hidx,N)
#' 
#' # Using inverse Fourier to obtain structure factors
#' ftmp2 <- invfousynth(a,rtmp$rr,hidx)
#' 
#' # Comparison
#' print(abs(ftmp1$Fmod-ftmp2$Fmod))
#' print(abs(ftmp1$Fpha-ftmp2$Fpha))
#' 
#' @export
invfousynth <- function(a,rho,hidx)
{
  # Grid length
  N <- length(rho)
  
  # Check numbers are compatible
  hmax <- max(abs(hidx))
  if (hmax >= 0.5*N) 
    stop("Grid too coarse. Decrease max Miller index, or increase grid sampling.")
  
  # Inverse Fourier transform
  G <- a*fft(rho,inverse=TRUE)/N
  
  # Temporary structure factors holder
  #tmpF <- G[1:(hmax+1)]
  
  # Final structure factors
  fullidx <- 0:hmax
  idx <- match(hidx,fullidx)
  F <- G[idx]
  
  # Amplitudes and phases
  Fmod <- Mod(F)
  Fpha <- Arg(F)*180/pi
  idx <- which(abs(Im(F)) < 0.000001 & abs(Re(F)) >= 0.000001 & Re(F) > 0)
  Fpha[idx] <- 0.0
  idx <- which(abs(Im(F)) < 0.000001 & abs(Re(F)) >= 0.000001 & Re(F) < 0)
  Fpha[idx] <- 180.0
  
  return(list(Fmod=Fmod,Fpha=Fpha))
}


#' Find optimal wavelength for anomalous phasing
#' 
#' This function mimics the behaviour of the fluorescent scan operated at
#' synchrotron beamlines to find out the optimal wavelength which maximises
#' the anomalous signal used to phase a structure.
#' 
#' @param chem_el 1- or 2-letters character string. The chemical symbol of
#'    interest.
#' @param lambda_range Real vector of length 2. The two values are the extremes of
#'    the wavelength window inside which to carry out the search. Default is
#'    for the search to be carried out across the full available range.
#' @return A named list with two elements: 1) idx is the integer indicating the
#'    row in the anomalous_data dataframe related to the specific chemical
#'    element, corresponding to the optimal wavelength; 2) the optimal wavelength
#'    in angstroms.
#' @examples 
#' # Optimal wavelength for iron
#' lFe <- fluorescent_scan("Fe")
#' print(lFe$lambda)
#' idx <- lFe$idx
#' 
#' # Load anomalous dataframe for Fe
#' adFe <- load_anomalous_data("Fe")
#' print(adFe[idx,])  # Same wavelength as before!
#' 
#' # Optimal wavelength with window restriction
#' lFe <- fluorescent_scan("Fe",lambda_range=c(6,8))
#' print(lFe$lambda)
#' 
#' @export
fluorescent_scan <- function(chem_el,lambda_range=NULL)
{
  # Load table related to specific chemical element
  ad <- load_anomalous_data(chem_el)
  
  # Restrict search to lambda_range
  if (!is.null(lambda_range))
  {
    lset <- which(ad$lambda >= lambda_range[1])
    rset <- which(ad$lambda <= lambda_range[2])
    jdx <- na.omit(match(lset,rset))
  }
  if (is.null(lambda_range)) jdx <- 1:length(ad[,1])
  
  # Find minimum of f' and select interval around it
  idx <- which(ad[jdx,2] == min(ad[jdx,2]))
  il <- idx-10
  if (il < 1) il <- 1
  ir <- idx+10
  if (ir > length(jdx)) ir <- length(jdx)
  
  # Select wavelength at which the difference between f' 
  # and f'' is largest
  tmp <- abs(ad[jdx[il]:jdx[ir],2]-ad[jdx[il]:jdx[ir],3])
  idx <- which(tmp == max(tmp))
  Clambda <- ad[jdx[il]+idx-1,1]
  idx <- which(abs(ad$lambda-Clambda) < 0.000001)
  
  return(list(idx=idx,lambda=ad[idx,1]))
}


#' Generation of structure factors with errors
#' 
#' This function generates structure factors calculated starting from
#' the given structure and subject to two types of errors: poissonian
#' counting errors due to the statistical nature of the photons hitting
#' the crystal and normal errors due to the slight random shifting of
#' atoms position in all the unit cells forming the lattice.
#' 
#' @param hidx Real numeric. One or more 1D Miller indices.
#' @param sdata A named list, normally obtained through the use of
#'  function \code{\link{read_x}} or \code{\link{standardise_sdata}}. 
#'  The list names correspond to different object types:
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
#' @param vx0err A real number. The standard deviation of the random 
#'  displacement of all atoms composing the structure from their correct
#'  position. Default value is NULL, corresponding to the generation of
#'  structure factors, with no errors, from the correct structure.
#' @param ntrialP Integer. The number of simulated Poisson counts for each
#'    set of structure factor amplitudes. More counts (high ntrialP) return
#'    smaller errors for the structure factor amplitudes. If ntrialP less or 
#'    equal 0, then poissonian counting errors are not generated.
#' @param ntrialG Integer. This is the number of randomly generated shifts of
#'    each atom of the structure from its true position. The shifts follow a
#'    gaussian distribution with mean 0 and standard deviation vx0err.
#' @param anoflag Logical variable. If TRUE it forces scattering factors
#'    to include anomalous contributions. As a consequence, theoretical
#'    Friedel's pairs will not be equal.
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
#' @return  A named list with two elements:
#'  \itemize{
#'    \item{F Array of mean structure factor amplitudes, among all structure
#'             factor arrays simulated with specific errors.}
#'    \item{sF  Array of structure factors errors. These coincide with the
#'               standard deviations of all structure factors arrays simulated
#'               with specific errors.}
#'         }
#' 
#' @examples 
#' 
#' # Load thiocyanate data
#' sdata <- load_structure("thiocyanate")
#' 
#' # Miller indices used
#' hidx <- 1:10
#' 
#' # Correct amplitudes and phases
#' ftmp <- strufac(hidx,sdata)
#' Ftrue <- ftmp$Fmod
#' phitrue <- ftmp$Fpha
#' 
#' # Only poissonian errors
#' ltmp <- sfobs(hidx,sdata,ntrialP=2)
#' print(names(ltmp))
#' Fpois <- ltmp$F
#' 
#' # True density
#' rtmptrue <- fousynth(sdata$a,Fmod=Ftrue,Fpha=phitrue,hidx,1000)
#' plot(rtmptrue$x,rtmptrue$rr,type="l",xlab="x",ylab=expression(rho),
#'  lwd=3)
#' 
#' # Density with poissonian errors
#' rtmppois <- fousynth(sdata$a,Fmod=Fpois,Fpha=phitrue,hidx,1000)
#' points(rtmppois$x,rtmppois$rr,type="l",
#'  lty=2,col=2,lwd=2) # Very small differences
#' 
#' # Only random atomic errors with standard deviation 0.3 angstroms
#' ltmp <- sfobs(hidx,sdata,ntrialP=0,vx0err=0.3)
#' Fcoords <- ltmp$F
#' 
#' # Density with gaussian errors on atom coordinates
#' rtmpcoords <- fousynth(sdata$a,Fmod=Fcoords,Fpha=phitrue,hidx,1000)
#' points(rtmpcoords$x,rtmpcoords$rr,type="l",
#'  lty=3,col=3,lwd=2) # Larger differences
#' @export
sfobs <- function(hidx,sdata,vx0err=NULL,ntrialP=100,ntrialG=100,
                  anoflag=FALSE,aK=anoK,lbda=1,k=ksigma,f1f2out=TRUE)
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
  
  # Data Ok. Carry on.
  
  # Initial sf and errors
  ftmp <- strufac(hidx,sdata,anoflag,aK,lbda,k,f1f2out)
  Fobs <- ftmp$Fmod
  rm(ftmp)
  Ffinal1 <- Fobs
  sFfinal1 <- rep(0,times=length(hidx))
  Ffinal2 <- Fobs
  sFfinal2 <- rep(0,times=length(hidx))
  
  # Find scales for pure-amplitudes noise (if needed)
  if (ntrialP > 0)
  {
    # Scale and shift Fobs interval to increase precision
    # for weak intensities
    M <- max(Fobs,na.rm=TRUE)
    m <- min(Fobs,na.rm=TRUE)
    b <- M-m
    newF <- 1000+1000*(Fobs-m)/b
    
    # Poissonians errors
    F <- rep(0,times=length(Fobs))
    sF <- rep(0,times=length(Fobs))
    for (i in 1:length(hidx))
    {
      tmp <- rpois(n=ntrialP,lambda=round(newF[i],digits=0))
      F[i] <- m+(M-m)*(mean(tmp)-1000)/1000
      sF[i] <- abs((M-m)*sd(tmp)/1000)
    }
    rm(newF,tmp)
    Ffinal1 <- F
    sFfinal1 <- sF
    rm(F,sF)
  }
  
  # If atoms in the structure are affected by errors
  if (!is.null(vx0err))
  {
    # New atomic data
    nx0 <- sdata$x0
    Fmatrix <- matrix(nrow=ntrialG,ncol=length(Fobs))
    for (itrial in 1:ntrialG)
    {
      tmpadata <- sdata
      if (!is.null(vx0err)) nx0 <- sdata$x0+rnorm(n=length(sdata$x0),
                                                  sd=vx0err)
      tmpadata$x0 <- nx0
      ftmp <- strufac(hidx,tmpadata,anoflag,aK,lbda,k,f1f2out)
      Fmatrix[itrial,] <- ftmp$Fmod
    }
    Ffinal2 <- apply(Fmatrix,2,mean,na.rm=TRUE)
    sFfinal2 <- apply(Fmatrix,2,sd,na.rm=TRUE)
  }
  
  # Joining the two contributions into one
  Ffinal <- 0.5*(Ffinal1+Ffinal2)
  sFfinal <- sqrt(sFfinal1^2+sFfinal2^2)
  
  return(list(F=Ffinal,sF=sFfinal))
}
  
  #' Simulation of 1D diffraction pattern
  #' 
  #' Analytic Fourier transform of electron density corresponding to
  #' an array of \code{Ncell} unit cells calculated using numerical 
  #' integration with the trapezoid rule. The diffraction peaks' height is 
  #' proportional to the number of unit cells ( \code{Ncell}). The 
  #' number of diffraction peaks included in the 1D diffraction pattern is 
  #' related to the maximum resolution \code{D} provided in the input. 
  #' The number of grid points for both the simulated electron density and
  #' the resulting diffraction pattern can also be provided as input. A
  #' further input parameter is the radius of the beamstop disc to stop
  #' diffraction close to the incoming beam (as the risulting intensity far
  #' outweight the rest of the diffracted intensities).
  #' 
  #' @param sdata A named list, normally obtained through the use of
  #'  functions \code{\link{read_x}} or \code{\link{standardise_sdata}}. 
  #'  The list names correspond to different object types:
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
  #' @param D Real numeric. Maximum resolution in angstroms.
  #' @param Ncell Positive integer. It is the number of unit cells in 
  #'    the 1D crystal. The default value is \code{Ncell=10}.
  #' @param N Positive integer indicating the number of grid points for the
  #'    electron density. The default value is \code{N=1000}.
  #' @param n Positive integer determining the reciprocal space grid. The
  #'    grid is made of \code{2*n+1} regularly-spaced points from \code{-1/D} 
  #'    to \code{1/D}. The value 0 is always at the centre of the grid. The 
  #'    default value is \code{n=100}.
  #' @param bstop Real numeric. Is the radius of the backstop disc. Intensities
  #'    at points closer to the origin than \code{bstop} are reduced to 
  #'    0. The presence of a backstop improves the contrast for all diffracted
  #'    intensities because it gets rid of the overwhelming intensity corresponding
  #'    to the origin of the reciprocal space. The default is not to include
  #'    any backstop.
  #' @return A named list with two vectors of real numbers, the values of the
  #'    reciprocal space grid points (in 1/angstrom units) \code{xstar} 
  #'    and the intensities \code{Imod}.
  #' @examples 
  #' # Diffraction from just two unit cells of cyanate
  #' sdata <- load_structure("cyanate")
  #' 
  #' # Max resolution is 1 angstroms; no backstop
  #' ltmp <- diffraction(sdata,D=1,Ncell=1)
  #' 
  #' # Plot diffraction pattern
  #' plot(ltmp$xstar,ltmp$Imod,type="l",
  #'  xlab=expression(paste("x"^"*")),ylab="Intensity")
  #' 
  #' # Diffraction from 20 unit cells with backstop of 20 angstroms diametre
  #' ltmp <- diffraction(sdata,D=1,bstop=10)
  #' plot(ltmp$xstar,ltmp$Imod,type="l",
  #'  xlab=expression(paste("x"^"*")),ylab="Intensity")
  #' 
  #' 
  #' @export
diffraction <- function(sdata,D,Ncell=10,N=1000,n=100,bstop=NULL) {
  # Grid
  x <- seq(0,Ncell*sdata$a,length=N)
  
  # Build analytic density
  rtmp <- structure_gauss(sdata,x=x)
  
  # Electron density
  rho <- rtmp$rr
  
  # Grid for the reciprocal space
  xstar <- seq(-1/D,0,length=(n+1))
  xstar <- c(xstar[1:(length(xstar)-1)],seq(0,1/D,length=(n+1)))
  
  # Final real and imaginary parts
  fRe <- rep(0,times=length(xstar))
  fIm <- rep(0,times=length(xstar))
  
  # Outer loop (over values of h)
  for (i in 1:length(xstar)) {
    # Value of reciprocal-space variable
    h <- xstar[i]
    
    # Inner loop (numerical integration over x)
    rhoRe <- rho * cos(2*pi*h*x)
    rhoIm <- rho * sin(2*pi*h*x)
    
    # Numerical integration (trapezoid rule)
    fRe[i] <- sum(diff(x) * (head(rhoRe,-1)+tail(rhoRe,-1)))/2
    fIm[i] <- sum(diff(x) * (head(rhoIm,-1)+tail(rhoIm,-1)))/2
  }
  
  # Diffraction intensities
  #Imod <- sqrt(fRe^2+fIm^2)
  Imod <- fRe^2+fIm^2
  
  # Zero backstop area if requested
  if(!is.null(bstop)) {
    Dstop <- 1/bstop
    idx <- which(abs(xstar)-Dstop <= 0)
    Imod[idx] <- 0
  }
  
  # Output
  return(list(xstar=xstar,Imod=Imod))
}  