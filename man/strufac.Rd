% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reciprocal_space.R
\name{strufac}
\alias{strufac}
\title{Calculation of structure factors}
\usage{
strufac(hidx, sdata, anoflag = FALSE, aK = anoK, lbda = 1,
  k = ksigma, f1f2out = TRUE)
}
\arguments{
\item{hidx}{Real numeric. One or more 1D Miller indices.}

\item{sdata}{A named list, normally obtained through the use of
functions \code{\link{read_x}} or \code{\link{standardise_sdata}}. 
The list names correspond to different object types:
\itemize{
  \item{a.     Real numeric. The size of the unit cell.}
  \item{SG.    Character string. Space group symbol; either "P1" 
               or "P-1"}
  \item{x0.    Vector of real numerics indicating the expanded atomic
               positions in the unit cell.}
  \item{Z.     Vector of integers indicating the expanded 
               atomic numbers for all atoms in the unit cell.}
  \item{B.    Vector of real numerics indicating the expanded 
              B factors for all atoms in the unit cell.}
  \item{occ.  Vector of real numerics indicating the expanded 
              occupancies for all atoms in the unit cell.}
}}

\item{anoflag}{Logical variable. If TRUE it forces scattering factors
to include anomalous contributions. As a consequence, Friedel's
pairs will not be equal.}

\item{aK}{Real numeric. This is a fudge factor included to decrease the
strength of the anomalous contributions. Without aK the strength is too
high for 1D structures, compared to real 3D structures. So aK helps
bringing down the anomalous contribution within the 5%-10% normally
met with large-size structures. The default value is aK=0.3 (anoK is
included as internal data).}

\item{lbda}{Real numeric. This is the wavelength in angstroms. Its value
is important in relation to anomalous scattering.}

\item{k}{A real number. It controls the standard deviation of the 
gaussian function describing the atom and, thus, the shape of the
associated peak. The standard deviation sigma is given by:
      \code{sigma = k * sqrt(Z)}
The default value is k=0.05 (ksigma is included as internal data).}

\item{f1f2out}{Logical variable. This variable controls output of a small
table of f' and f'' values for all chemical elements in the structure.
Default is for the table to be printed.}
}
\value{
A named list with two vectors of real numbers, the structure 
 factors amplitudes, Fmod, and phases, Fpha, corresponding to the Miller 
 indices in input.
}
\description{
Calculates structure factors corresponding to one or more 1D Miller 
indices, given the atomic content of one unit cell. Anomalous scattering 
can be included using logical flag "anoflag" (default is FALSE). Crystal 
structures are always considered with no symmetry. Thus a 1D structure 
with the P-1 symmetry will have to be expanded first with function 
"expand_to_cell".
}
\examples{
# First create the crystal structure (P1)
a <- 10
SG <- "P1"
x0 <- c(1,4,6.5)
Z <- c(8,26,6)
B <- c(18,20,17)
occ <- c(1,1,1)
sdata <- standardise_sdata(a,SG,x0,Z,B,occ)
lbda <- 1.7   # Close to Fe's absorption

# Miller indices (including DC (h=0) component)
hidx <- 0:10

# Now structure factors without anomalous contribution
ftmp <- strufac(hidx,sdata,lbda=lbda)
print(length(ftmp$Fmod))  # Includes DC component (h=0)
print(ftmp)          # Amplitudes decrease with increasing
                     # resolution (Miller indices)
                  
# Now try with anomalous contributions
ftmp <- strufac(hidx,sdata,lbda=lbda,anoflag=TRUE)
print(ftmp)  # DC component is not any longer real

}
