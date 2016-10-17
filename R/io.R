# 
# This file is part of the crone package
#

#' Write atomic coordinates to a file.
#' 
#' Function to export all information concerning a given structure to a
#' so-called coordinates file of type *_x.dat.
#' 
#' @param filename A character string. Prefix of the output ASCII file to
#'    include all structural information. The file name will be "[Prefix]_x.dat".
#' @param vx0 Vector of real numerics. Atom positions in the asymmetric unit.
#' @param vZ A vector of atom Z numbers.
#' @param vB Vector of real numerics. B factors for all atoms in the
#'  asymmetric unit.
#' @param vocc Vector of real numerics. Occupancies (value between 0
#'  and 1) for all atoms in the asymmetric unit. In this function there
#'  is no mechanism to check whether the occupancy is appropriate in case
#'  the atom is at a special position.
#' @param a Real numeric. Unit cell length in angstroms.
#' @param SG 2-letters character string. There are only two symmetries
#'  possible when working within 1D crystallography, P1 (no symmetry)
#'  and P-1 (inversion through the origin). SG can be either "P1" or
#'  "P-1" for this function.
#' @return This function does not return anything, but will create an ASCII file
#'    of name *_x.dat which contains all coordinates of the atoms in the structure
#'    and other type of information.
#' @examples 
#' # Create an arbitrary structure in P1
#' SG <- "P1"
#' a <- 20
#' vx0 <- c(2,11,16,19)
#' vZ <- c(6,6,16,8)
#' vB <- c(13,14,5,10)
#' vocc <- c(1,1,1,1)
#' prfx <- "test"
#' write_x(prfx,vx0,vZ,vB,vocc,a,SG)
#' 
#' @export
write_x <- function(filename,vx0,vZ,vB,vocc,a,SG)
{
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
#'    \item{data. Data frame with columns having the following names: x0 (atom
#'                coordinates), Z (atomic numbers), B (B-factors), 
#'                occ (occupancies).}
#'          }
#' @examples
#' datadir <- system.file("extdata",package="crone")
#' filename <- file.path(datadir,"carbon_dioxide_x.dat")
#' ltmp <- read_x(filename)
#' print(names(ltmp))
#' print(ltmp$data)
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
  read <- list()
  read$a <- a
  read$SG <- SG
  read$data <- data.frame(x0=data$V3,Z=Z,B=data$V4,occ=data$V5)
  
  return(read)
}