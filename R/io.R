# 
# This file is part of the crone package
#

#' Organise atom data in a standard format for later use
#' 
#' Function to put any group consisting of cell size a, space group SG,
#' atom coordinates x0, atomic numbers Z, atomic B factors B and atomic
#' occupancies occ, into a named list with 6 fields. This is then used
#' as standard input/output format throughout the \emph{crone} package.
#' 
#' @param a    Real numeric. Unit cell length in angstroms.
#' @param SG   SG 2-letters character string. There are only two symmetries
#'  possible when working within 1D crystallography, P1 (no 
#'  symmetry) and P-1 (inversion through the origin). SG can be 
#'  either "P1" or "P-1".
#' @param x0    Vector of real numerics indicating the expanded atomic
#'  positions in the unit cell.
#' @param Z     Vector of integers indicating the expanded atomic 
#'  numbers for all atoms in the unit cell.
#' @param B    Vector of real numerics indicating the expanded B factors 
#'  for all atoms in the unit cell.
#' @param occ  Vector of real numerics indicating the expanded occupancies
#'  for all atoms in the unit cell.
#' @return A named list with the following elements:
#'  \itemize{
#'    \item{a.    Real numeric. Unit cell length in angstroms.}
#'    \item{SG.   SG 2-letters character string. There are only two symmetries
#'                possible when working within 1D crystallography, P1 (no 
#'                symmetry) and P-1 (inversion through the origin). SG can be 
#'                either "P1" or "P-1".}
#'    \item{x0.    Vector of real numerics indicating the expanded atomic
#'                 positions in the unit cell.}
#'    \item{Z.     Vector of integers indicating the expanded 
#'                 atomic numbers for all atoms in the unit cell.}
#'    \item{B.    Vector of real numerics indicating the expanded 
#'                B factors for all atoms in the unit cell.}
#'    \item{occ.  Vector of real numerics indicating the expanded 
#'                occupancies for all atoms in the unit cell.}
#'          }
#' @examples 
#' # Create standard format for an arbitrary structure in P1
#' a <- 20
#' SG <- "P1"
#' x0 <- c(2,11,16,19)
#' Z <- c(6,6,16,8)
#' B <- c(13,14,5,10)
#' occ <- c(1,1,1,1)
#' sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' print(sdata)
#' 
#' @export
standardise_sdata <- function(a,SG,x0,Z,B,occ)
{
  # Check last 4 arrays have all same length
  if ((length(Z) != length(x0)) |
      (length(B) != length(x0)) |
      (length(occ) != length(x0)))
  {
    stop("Last 4 arrays in input must all have same size.")
  }
  
  # Create named list with 6 elements
  sdata <- list()
  sdata$a <- a
  sdata$SG <- SG
  sdata$x0 <- x0
  sdata$Z <- Z
  sdata$B <- B
  sdata$occ <- occ
  
  return(sdata)
}


#' Load 1D structure data in workspace.
#' 
#' Function to load data from one of the many 1D structures available
#' within the \emph{crone} package.
#' 
#' @param sname A character string. Name of the structure whose data
#'  are to be loaded in the workspace. It can be one of:
#'  \itemize{
#'    \item{beryllium_fluoride}
#'    \item{carbon_dioxide}
#'    \item{cyanate}
#'    \item{nitronium}
#'    \item{thiocyanate}
#'    \item{xenon_difluoride}
#'  }
#'  Default is NULL, in which case the function returns a list of all
#'  structures available.
#' @return A named list with the following elements:
#'  \itemize{
#'    \item{a    Real numeric. Unit cell length in angstroms.}
#'    \item{SG   2-letters character string. There are only two symmetries
#'                possible when working within 1D crystallography, P1 (no 
#'                symmetry) and P-1 (inversion through the origin). SG can be 
#'                either "P1" or "P-1".}
#'    \item{x0    Vector of real numerics indicating the expanded atomic
#'                 positions in the unit cell.}
#'    \item{Z     Vector of integers indicating the expanded 
#'                 atomic numbers for all atoms in the unit cell.}
#'    \item{B    Vector of real numerics indicating the expanded 
#'                B factors for all atoms in the unit cell.}
#'    \item{occ  Vector of real numerics indicating the expanded 
#'                occupancies for all atoms in the unit cell.}
#'          }
#' @examples 
#' # Load thiocyanate data
#' sdata <- load_structure("thiocyanate")
#' print(sdata)
#' 
#' # Default returns all names of structures included
#' load_structure()
#' 
#' @export
load_structure <- function(sname=NULL)
{
  # All names
  all_names <- c("beryllium_fluoride",
                 "carbon_dioxide","cyanate",
                 "nitronium","thiocyanate",
                 "xenon_difluoride")

  # Load data if they exist
  if (!is.null(sname))
  {
    ans <- (sname %in% all_names)
  
    if (!ans)
    {
      warning("No structure with this name is included in the 
              CRONE package.")
    }
    if (ans)
    {
      stmp <- paste(sname,"_x.dat",sep="")
      datadir <- system.file("extdata",package="crone")
      filename <- file.path(datadir,stmp)
      sdata <- read_x(filename)
    }
    
    return(sdata)
  }
  # If no input return all names
  if (is.null(sname))
  {
    cat("1D structures available for loading:\n\n")
    for (stmp in all_names)
    {
      cat(paste("   ",stmp,"\n",sep=""))
    }
  }
}

#' Write atomic coordinates to a file.
#' 
#' Function to export all information concerning a given structure to a
#' so-called coordinates file of type *_x.dat.
#' 
#' @param filename A character string. Prefix of the output ASCII file to
#'    include all structural information. The file name will be "[Prefix]_x.dat".
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
#' @return This function does not return anything, but will create an ASCII 
#'  file of name *_x.dat which contains all coordinates of the atoms in the 
#'  structure and other type of information.
#' @examples 
#' # Create an arbitrary structure in P1
#' a <- 23
#' SG <- "P1"
#' x0 <- c(2,11,16,19)
#' Z <- c(6,6,16,8)
#' B <- c(13,14,5,10)
#' occ <- c(1,1,1,1)
#' prfx <- "test"
#' sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' write_x(prfx,sdata)
#' 
#' @export
write_x <- function(filename,sdata)
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
  
  # Input OK. Carry on
  a <- sdata$a
  SG <- sdata$SG
  vx0 <- sdata$x0
  vZ <- sdata$Z
  vB <- sdata$B
  vocc <- sdata$occ
  
  # User gives name to file. Then "_x.dat" is appended
  fname <- paste(filename,"_x.dat",sep="")
  
  # First line reports unit cell size. It also creates file
  line <- sprintf("Unit Cell length:    %10.3f\n",a)
  cat(line,file=fname)
  
  # Second line reports space group.
  line <- sprintf("     Space Group:         %s\n",SG)
  cat(line,file=fname,append=TRUE)
  
  # All other lines
  lines <- c()
  for (i in 1:length(vZ))
  {
    # Atomic name
    idx <- match(vZ[i],atoms$Z)
    aname <- as.character(atoms$anames[idx])
    lines <- c(lines,sprintf("%5d  %s  %10.3f  %10.3f  %4.1f",
                             i,aname,vx0[i],vB[i],vocc[i]))
  }
  
  # Write all other lines in file
  cat(lines,file=fname,sep="\n",append=TRUE)
}


#' Read unit cell content (atom and coordinates).
#' 
#' Read unit cell length, space group, atom coordinates and all other
#' parameters from a coordinate file.
#' 
#' @param filename A character string. Existing file that includes all structural
#'    information. The file name in general has the form "[prefix]_x.dat".
#' @return A named list with the following elements:
#'  \itemize{
#'    \item{a.    Real numeric. Unit cell length in angstroms.}
#'    \item{SG.   SG 2-letters character string. There are only two symmetries
#'                possible when working within 1D crystallography, P1 (no 
#'                symmetry) and P-1 (inversion through the origin). SG can be 
#'                either "P1" or "P-1".}
#'    \item{x0.    Vector of real numerics indicating the expanded atomic
#'                 positions in the unit cell.}
#'    \item{Z.     Vector of integers indicating the expanded 
#'                 atomic numbers for all atoms in the unit cell.}
#'    \item{B.    Vector of real numerics indicating the expanded 
#'                B factors for all atoms in the unit cell.}
#'    \item{occ.  Vector of real numerics indicating the expanded 
#'                occupancies for all atoms in the unit cell.}
#'          }
#' @examples
#' datadir <- system.file("extdata",package="crone")
#' filename <- file.path(datadir,"carbon_dioxide_x.dat")
#' sdata <- read_x(filename)
#' print(names(sdata))
#' print(sdata)
#' 
#' @export
read_x <- function(filename)
{
  data <- read.table(filename,skip=2,stringsAsFactors=FALSE)
  chars <- scan(file=filename,what=character(),sep = "\n",quiet=TRUE)
  stmp <- strsplit(chars[1],split=" ")[[1]]
  a <- as.numeric(stmp[length(stmp)])
  stmp <- strsplit(chars[2],split=" ")[[1]]
  SG <- stmp[length(stmp)]
  anames <- as.character(data$V2)
  Z = rep(0,length=length(anames))
  for(j in 1:length(anames))
  {
    if(nchar(anames[j]) == 1)
    {
      char <- paste(" ",anames[j],sep="")
      idx <- match(char,atoms$anames)
      Z[j] <- atoms$Z[idx]
    }
    if(nchar(anames[j]) == 2)
    {
      idx <- match(anames[j],atoms$anames)
      Z[j] <- atoms$Z[idx]
    }
  }
  
  # Prepare output
  sdata <- list()
  sdata$a <- a
  sdata$SG <- SG
  sdata$x0 <- data$V3
  sdata$Z <- Z
  sdata$B <- data$V4
  sdata$occ <- data$V5
  
  return(sdata)
}


#' Standardise reflections data
#' 
#' This function output a list with fields needed by most of the
#' functions dealing with structure factors. It is the equivalent
#' of the function \cite{standardise_sdata}, used to prepare atomic
#' structures data.
#' 
#' @param a    Real numeric. Unit cell length in angstroms.
#' @param SG   SG 2-letters character string. There are only two symmetries
#'  possible when working within 1D crystallography, P1 (no 
#'  symmetry) and P-1 (inversion through the origin). SG can be 
#'  either "P1" or "P-1".
#' @param hidx  Real numeric array. 1D unique (positive in the 1D context)
#'  Miller indices.
#' @param Fp   Real numeric vector. Amplitudes of the positive 
#'  component of Friedel (or Bijvoet) pairs (F+). Default is NULL,
#'  i.e. no Fp included.
#' @param sigFp  Real numeric vector. Errors associated with Fp. Default
#'  is NULL, i.e. no sigFp included.
#' @param Fm   Real numeric vector. Amplitudes of the negative 
#'  component of Friedel (or Bijvoet) pairs (F-). Default is NULL,
#'  i.e. no Fm included.
#' @param sigFm  Real numeric vector. Errors associated with Fm. Default
#'  is NULL, i.e. no sigFm included.
#' @param Fobs  Real numeric array. Amplitudes of structure factors.
#'  If Fp and Fm are not NULL and Fobs is NULL, then Fobs are calculated
#'  as averages of Fp and Fm. If both Fp, Fm and Fobs are included, input
#'  Fobs are used, instead of Fp and Fm averages.
#' @param sigFobs  Real numeric array. Errors associated with Fobs. If
#'  sigFobs = NULL, errors are estimated from Fp and Fm. Default is NULL.
#' @param Phiobs  Real numeric array. Phases (in degrees) of structure 
#' factors obtained with one of the methods used for structure solution.
#' Default is NULL.
#' @param Phicalc  Real numeric array. Phases (in degrees) of structure 
#' factors calculated from the correct 1D structure. They are normally 
#' used to check correctness of Phiobs. Default is NULL.
#' @return A named list with a variable number of elements. Some of them
#'  are always included; others are not:
#'  \itemize{
#'    \item{a    Real numeric. Unit cell length in angstroms. Always 
#'               included.}
#'    \item{SG.   Spacegroup 2-letters character string. There are only two 
#'                symmetries possible when working within 1D 
#'                crystallography, P1 (no symmetry) and P-1 (inversion 
#'                through the origin). SG can be either "P1" or "P-1". 
#'                Always included.}
#'    \item{hidx. Real numeric array. 1D unique (positive in the 1D context) 
#'                Miller indices. Always included.}
#'    \item{Fobs.      Real numeric array. Amplitudes of observed structure 
#'                     factors. Not always included.}
#'    \item{sigFobs.   Real numeric array. Errors associated with Fobs. Not 
#'                     always included.}
#'    \item{Fp.        Real numeric vector. Amplitudes of the positive 
#'                     component of Friedel (or Bijvoet) pairs (F+). Not
#'                     always included.}
#'    \item{sigFp.     Real numeric vector. Errors associated with Fp. 
#'                     Not always included.}
#'    \item{Fm.        Real numeric vector. Amplitudes of the negative 
#'                     component of Friedel (or Bijvoet) pairs (F-). Not always
#'                     included.}
#'    \item{sigFm.     Real numeric vector. Errors associated with Fm. Not
#'                     always included.}
#'    \item{Phiobs.    Real numeric array. Phases (in degrees) of structure 
#'                     factors obtained with one of the methods used for 
#'                     structure solution. Not always included.}
#'    \item{Phicalc.   Real numeric array. Phases (in degrees) of structure 
#'                     factors calculated from the correct 1D structure. 
#'                     They are normally used to check correctness of 
#'                     Phiobs. Not always included.}
#'          }
#'
#' @examples
#' # Create an arbitrary structure with a heavy atom (Fe)
#' a <- 20
#' SG <- "P1"
#' x0 <- c(1,2,6,16)
#' Z <- c(6,8,26,7)
#' B <- c(8,7,5,8)
#' occ <- c(1,1,1,1)
#' sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
#' 
#' # Miller indices, from -5 to 5 (to include negatives for anomalous)
#' hidx <- -5:5
#' 
#' # Experimental structure factors with anomalous contribution
#' # (lambda = 1.74) for creating Fm and Fp. Errors only due to
#' # photons fluctuations.
#' set.seed(9195)   # For demo purposes.
#' ltmp <- sfobs(hidx,sdata,ntrialP=10,anoflag=TRUE,lbda=1.74)
#' 
#' # Fp and sigFp (Miller indices from 1 to 5)
#' isel <- 1:5
#' idx <- match(isel,hidx)
#' Fp <- ltmp$F[idx]
#' sigFp <- ltmp$sF[idx]
#' 
#' # Fm and sigFm
#' isel <- (-1):(-5)
#' idx <- match(isel,hidx)
#' Fm <- ltmp$F[idx]
#' sigFm <- ltmp$sF[idx]
#' 
#' # Now only positive Miller indices
#' hidx <- 1:5
#' 
#' # Create standardised data for reciprocal space
#' fdata <- standardise_fdata(a,SG,hidx,Fp=Fp,sigFp=sigFp,
#'          Fm=Fm,sigFm=sigFm)
#'          
#' # Fp and Fm
#' print(fdata$Fp)
#' print(fdata$sigFp)
#' print(fdata$Fm)
#' print(fdata$sigFm)
#' 
#' # Fobs and sigFobs automatically created
#' print(fdata$Fobs)
#' print(fdata$sigFobs)
#'
#' # Structure factors without anomalous data for the same structure
#' hidx <- 1:5
#' set.seed(9195)   # For demo purposes.
#' ltmp <- sfobs(hidx,sdata,ntrialP=10)
#' Fobs <- ltmp$F
#' sigFobs <- ltmp$sF
#' fdata <- standardise_fdata(a,SG,hidx,Fobs=Fobs,sigFobs=sigFobs)
#' print(fdata)
#' 
#'@export
standardise_fdata <- function(a,SG,hidx,
                              Fobs=NULL,sigFobs=NULL,
                              Fp=NULL,sigFp=NULL,
                              Fm=NULL,sigFm=NULL,
                              Phiobs=NULL,Phicalc=NULL)
{
  # Check Miller indices are all positive or zero
  idx <- which(hidx < 0)
  if (length(idx) > 0) stop("Only unique (positive) 
                            Miller indices allowed.")
  
  # Check data arrays have same length as hidx array
  if (!is.null(Fp) & length(Fp) != length(hidx))
  {
    stop("Length of Fp different from length of Miller indices.")
  }
  if (!is.null(sigFp) & length(sigFp) != length(hidx))
  {
    stop("Length of sigFp different from length of Miller indices.")
  }
  if (!is.null(Fm) & length(Fm) != length(hidx))
  {
    stop("Length of Fm different from length of Miller indices.")
  }
  if (!is.null(sigFm) & length(sigFm) != length(hidx))
  {
    stop("Length of sigFm different from length of Miller indices.")
  }
  if (!is.null(Fobs) & length(Fobs) != length(hidx))
  {
    stop("Length of Fobs different from length of Miller indices.")
  }
  if (!is.null(sigFobs) & length(sigFobs) != length(hidx))
  {
    stop("Length of sigFobs different from length of Miller indices.")
  }
  if (!is.null(Phiobs) & length(Phiobs) != length(hidx))
  {
    stop("Length of Phiobs different from length of Miller indices.")
  }
  if (!is.null(Phicalc) & length(Phicalc) != length(hidx))
  {
    stop("Length of Phicalc different from length of Miller indices.")
  }
  
  # If either Fp or Fm are present, also the other has to be present
  if ((is.null(Fp) & !is.null(Fm)) | (!is.null(Fp) & is.null(Fm)))
  {
    stop("Both Fp and Fm are needed.")
  }
  
  # If Fp and Fm present and Fobs not present, calculates Fobs, sigFobs
  if (is.null(Fobs) & !is.null(Fp) & !is.null(Fm))
  {
    Fobs <- (Fp + Fm)*0.5
  }
  if (is.null(sigFobs) & !is.null(sigFp) & !is.null(sigFm))
  {
    sigFobs <- sqrt(sigFp^2 + sigFm^2)
  }
  
  # Assemble list
  fdata <- list()
  fdata$a <- a
  fdata$SG <- SG
  fdata$hidx <- hidx
  fdata$Fobs <- Fobs
  fdata$sigFobs <- sigFobs
  fdata$Fp <- Fp
  fdata$sigFp <- sigFp
  fdata$Fm <- Fm
  fdata$sigFm <- sigFm
  fdata$Phiobs <- Phiobs
  fdata$Phicalc <- Phicalc
  
  return(fdata)
}

#' Write structure factors to a reflections file
#' 
#' This function writes standardised structure factors data into an ASCII
#' file. The files includes cell size, space group character symbol and
#' Miller indices vector. It can include all of some of observed and/or
#' calculated structure factors amplitudes and phases, either for anomalous
#' or non-anomalous data.
#' 
#' @param filename A character string. Prefix of the structure factors
#' file name. The file name has the form "[prefix]_h.dat".
#' @param fdata  A names list, usually created with functions 
#'  \code{\link{standardise_fdata}} or \code{\link{read_h}}, and consisting of
#'  the following fields:
#'  \itemize{
#'    \item{a    Real numeric. Unit cell length in angstroms. Always 
#'               included.}
#'    \item{SG.   Spacegroup 2-letters character string. There are only two 
#'                symmetries possible when working within 1D 
#'                crystallography, P1 (no symmetry) and P-1 (inversion 
#'                through the origin). SG can be either "P1" or "P-1". 
#'                Always included.}
#'    \item{hidx. Real numeric array. 1D unique (positive in the 1D context) 
#'                Miller indices. Always included.}
#'    \item{Fobs.      Real numeric array. Amplitudes of observed structure 
#'                     factors. Not always included.}
#'    \item{sigFobs.   Real numeric array. Errors associated with Fobs. Not 
#'                     always included.}
#'    \item{Fp.        Real numeric vector. Amplitudes of the positive 
#'                     component of Friedel (or Bijvoet) pairs (F+). Not
#'                     always included.}
#'    \item{sigFp.     Real numeric vector. Errors associated with Fp. 
#'                     Not always included.}
#'    \item{Fm.        Real numeric vector. Amplitudes of the negative 
#'                     component of Friedel (or Bijvoet) pairs (F-). Not always
#'                     included.}
#'    \item{sigFm.     Real numeric vector. Errors associated with Fm. Not
#'                     always included.}
#'    \item{Phiobs.    Real numeric array. Phases (in degrees) of structure 
#'                     factors obtained with one of the methods used for 
#'                     structure solution. Not always included.}
#'    \item{Phicalc.   Real numeric array. Phases (in degrees) of structure 
#'                     factors calculated from the correct 1D structure. 
#'                     They are normally used to check correctness of 
#'                     Phiobs. Not always included.}
#'          }
#'          
#' @return This function does not return anything, but will create an ASCII 
#'  file of name *_h.dat which contains structure factors and other type of 
#'  information.
#'  
#' @examples
#' # Data from thiocyanate structure
#' datadir <- system.file("extdata",package="crone")
#' filename <- file.path(datadir,"thiocyanate_x.dat")
#' sdata <- read_x(filename)
#' 
#' # Miller indices
#' hidx <- 1:10
#' 
#' # Observed structure factors with errors
#' ltmp <- sfobs(hidx,sdata)
#' Fobs <- ltmp$F
#' sigFobs <- ltmp$sF
#' 
#' # Phases from calculated structure factors
#' ftmp <- strufac(hidx,sdata)
#' phicalc <- ftmp$Fpha
#' 
#' # Create standardised fdata structure
#' fdata <- standardise_fdata(sdata$a,sdata$SG,hidx,Fobs=Fobs,
#'  sigFobs=sigFobs,Phicalc=phicalc)
#'  
#' # Name of structure factors file (in home directory)
#' homedir <- Sys.getenv("HOME")
#' fname <- file.path(homedir,"test")
#' 
#' # Write data to file
#' write_h(fname,fdata)
#' 
#' @export
write_h <- function(filename,fdata)
{
  # Initial flag
  anoflag <- FALSE
  
  # Make fdata fields explicit
  a <- fdata$a
  SG <- fdata$SG
  hidx <- fdata$hidx
  Fobs <- fdata$Fobs
  sigFobs <- fdata$sigFobs
  Fp <- fdata$Fp
  sigFp <- fdata$sigFp
  Fm <- fdata$Fm
  sigFm <- fdata$sigFm
  Phiobs <- fdata$Phiobs
  Phicalc <- fdata$Phicalc
  
  # If both Fp and Fm are not NULL Fobs is not expected
  if (length(Fp) != 0 & length(Fm) != 0)
  {
    if (length(Fp) != length(sigFp) | 
        length(Fp) != length(Fm)    |
        length(Fp) != length(sigFm)) 
    {
      msg <- paste("Structure factors information insufficient ",
                   "or structure factors arrays have differing lengths.",
                   sep="")
      stop(msg)
    }
    anoflag <- TRUE
    if (length(Fobs) == 0 & is.null(Fobs))
    {
      Fobs <- (Fp + Fm)*0.5
      sigFobs <- sqrt(sigFp^2 + sigFm^2)
    }
  }
  if (!anoflag & length(Fobs) == 0 & is.null(Fobs) |
      length(sigFobs) == 0 & is.null(sigFobs))
  {
    msg <- "Structure factors information insufficient."
    stop(msg)
  }
  if (length(Fobs) != length(sigFobs))
  {
    msg <- "Structure factors arrays have differing lengths."
    stop(msg)
  }
  
  # Now a reflections file can be written out
  
  # User gives name to file. Then "_h.dat" is appended
  fname <- paste(filename,"_h.dat",sep="")
  
  # First line reports unit cell size. It also creates file
  lines <- sprintf("Unit Cell length:    %10.3f",a)
  lines <- c(lines,paste("Space Group:      ",SG))
  lines <- c(lines," ")
  
  # If F(+) and F(-) are present
  if (anoflag)
  {
    # Column names
    titles <- sprintf("%5s %10s %8s %10s %8s %10s %8s","h",
                      "Fp","sigFp","Fm","sigFm","Fobs","sigFobs")
    if (!is.null(Phiobs) & is.null(Phicalc))
    {
      titles <- sprintf("%5s %10s %8s %10s %8s %10s %8s %8s","h",
                "Fp","sigFp","Fm","sigFm","Fobs","sigFobs","Phiobs")
    }
    if (is.null(Phiobs) & !is.null(Phicalc))
    {
      titles <- sprintf("%5s %10s %8s %10s %8s %10s %8s %8s","h",
                        "Fp","sigFp","Fm","sigFm","Fobs","sigFobs","Phicalc")
    }
    if (!is.null(Phiobs) & !is.null(Phicalc))
    {
      titles <- sprintf("%5s %10s %8s %10s %8s %10s %8s %8s %8s","h",
                        "Fp","sigFp","Fm","sigFm","Fobs","sigFobs",
                        "Phiobs","Phicalc")
    }
    lines <- c(lines,titles)
    for(i in 1:length(hidx))
    {
      if (is.null(Phiobs) & is.null(Phicalc))
      {
        lines <- c(lines,sprintf(
        "%5d %10.3f %8.3f %10.3f %8.3f %10.3f %8.3f",
        hidx[i],Fp[i],sigFp[i],Fm[i],sigFm[i],Fobs[i],sigFobs[i]))
      }
      if (!is.null(Phiobs) & is.null(Phicalc))
      {
        lines <- c(lines,sprintf(
          "%5d %10.3f %8.3f %10.3f %8.3f %10.3f %8.3f %8.1f",
          hidx[i],Fp[i],sigFp[i],Fm[i],sigFm[i],Fobs[i],sigFobs[i],
          Phiobs[i]))
      }
      if (is.null(Phiobs) & !is.null(Phicalc))
      {
        lines <- c(lines,sprintf(
          "%5d %10.3f %8.3f %10.3f %8.3f %10.3f %8.3f %8.1f",
          hidx[i],Fp[i],sigFp[i],Fm[i],sigFm[i],Fobs[i],sigFobs[i],
          Phicalc[i]))
      }
      if (!is.null(Phiobs) & !is.null(Phicalc))
      {
        lines <- c(lines,sprintf(
          "%5d %10.3f %8.3f %10.3f %8.3f %10.3f %8.3f %8.1f %8.1f",
          hidx[i],Fp[i],sigFp[i],Fm[i],sigFm[i],Fobs[i],sigFobs[i],
          Phiobs[i],Phicalc[i]))
      }
    }
  }
  
  # If F(+) F(-) are not present  
  if(!anoflag)
  {
    # Column names
    titles <- sprintf("%5s %10s %8s","h","Fobs","sigFobs")
    if (!is.null(Phiobs) & is.null(Phicalc))
    {
      titles <- sprintf("%5s %10s %8s %8s",
                        "h","Fobs","sigFobs","Phiobs")
    }
    if (is.null(Phiobs) & !is.null(Phicalc))
    {
      titles <- sprintf("%5s %10s %8s %8s",
                        "h","Fobs","sigFobs","Phicalc")
    }
    if (!is.null(Phiobs) & !is.null(Phicalc))
    {
      titles <- sprintf("%5s %10s %8s %8s %8s",
                        "h","Fobs","sigFobs","Phiobs","Phicalc")
    }
    lines <- c(lines,titles)
    for (i in 1:length(hidx))
    {
      if (is.null(Phiobs) & is.null(Phicalc))
      {
        lines <- c(lines,sprintf("%5d %10.3f %8.3f",
                               hidx[i],Fobs[i],sigFobs[i]))
      }
      if (!is.null(Phiobs) & is.null(Phicalc))
      {
        lines <- c(lines,sprintf("%5d %10.3f %8.3f %8.1f",
                                 hidx[i],Fobs[i],sigFobs[i],
                                 Phiobs[i]))
      }
      if (is.null(Phiobs) & !is.null(Phicalc))
      {
        lines <- c(lines,sprintf("%5d %10.3f %8.3f %8.1f",
                                 hidx[i],Fobs[i],sigFobs[i],
                                 Phicalc[i]))
      }
      if (!is.null(Phiobs) & !is.null(Phicalc))
      {
        lines <- c(lines,sprintf("%5d %10.3f %8.3f %8.1f %8.1f",
                                 hidx[i],Fobs[i],sigFobs[i],
                                 Phiobs[i],Phicalc[i]))
      }
    }
  }
  
  # Write all other lines in file
  cat(lines,file=fname,sep="\n")
}


#' Read data from a reflections file
#' 
#' Read data from a *_h.dat-type file containing cell size, spacegroup
#' symbol and amplitudes and/or phases of observed and/or calculated
#' structure factors. This function loads the file data into a standardised
#' named list for structure factors data.
#' 
#' @param filename A character string. Existing file that includes structure 
#'  factors information. The file name in general has the form 
#'  "[prefix]_h.dat".
#' @return A named list with the following elements:
#'  \itemize{
#'    \item{a    Real numeric. Unit cell length in angstroms. Always 
#'               included.}
#'    \item{SG.   Spacegroup 2-letters character string. There are only two 
#'                symmetries possible when working within 1D 
#'                crystallography, P1 (no symmetry) and P-1 (inversion 
#'                through the origin). SG can be either "P1" or "P-1". 
#'                Always included.}
#'    \item{hidx. Real numeric array. 1D unique (positive in the 1D context) 
#'                Miller indices. Always included.}
#'    \item{Fobs.      Real numeric array. Amplitudes of observed structure 
#'                     factors. Not always included.}
#'    \item{sigFobs.   Real numeric array. Errors associated with Fobs. Not 
#'                     always included.}
#'    \item{Fp.        Real numeric vector. Amplitudes of the positive 
#'                     component of Friedel (or Bijvoet) pairs (F+). Not
#'                     always included.}
#'    \item{sigFp.     Real numeric vector. Errors associated with Fp. 
#'                     Not always included.}
#'    \item{Fm.        Real numeric vector. Amplitudes of the negative 
#'                     component of Friedel (or Bijvoet) pairs (F-). Not always
#'                     included.}
#'    \item{sigFm.     Real numeric vector. Errors associated with Fm. Not
#'                     always included.}
#'    \item{Phiobs.    Real numeric array. Phases (in degrees) of structure 
#'                     factors obtained with one of the methods used for 
#'                     structure solution. Not always included.}
#'    \item{Phicalc.   Real numeric array. Phases (in degrees) of structure 
#'                     factors calculated from the correct 1D structure. 
#'                     They are normally used to check correctness of 
#'                     Phiobs. Not always included.}
#'          }
#'  
#' @examples
#' # Observed structure factors amplitudes and calculated phases
#' # from thiocyanate structure
#' datadir <- system.file("extdata",package="crone")
#' filename <- file.path(datadir,"thiocyanate_h.dat")
#' fdata <- read_x(filename)
#' print(names(fdata))
#' print(fdata$Fobs)
#' print(fdata$sigFobs)
#' 
#' @export
read_h <- function(filename)
{
  # Skip lines and read tabular data
  tmp <- read.table(file=filename,header=TRUE,skip=3)
  
  # Read initial lines to extract cell size and space group info
  chars <- scan(file=filename,what=character(),sep = "\n",quiet=TRUE)
  stmp <- strsplit(chars[1],split=" ")[[1]]
  a <- as.numeric(stmp[length(stmp)])
  stmp <- strsplit(chars[2],split=" ")[[1]]
  SG <- stmp[length(stmp)]
  
  # Check column names are ok
  # Data in list structure
  fdata <- list()
  fdata$a <- a
  fdata$SG <- SG
  fdata$hidx <- tmp$h
  fdata$Fobs <- tmp$Fobs
  fdata$sigFobs <- tmp$sigFobs
  fdata$Fp <- tmp$Fp
  fdata$sigFp <- tmp$sigFp
  fdata$Fm <- tmp$Fm
  fdata$sigFm <- tmp$sigFm
  fdata$Phiobs <- tmp$Phiobs
  fdata$Phicalc <- tmp$Phicalc
  
  return(fdata)
}