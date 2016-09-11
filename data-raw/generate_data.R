# 
# This file is part of the crone package
#


# Use this function only once. It returns a list with
# dataframes as elements. Each dataframe corresponds to
# anomalous data for a specific chamical element. The
# list name is anomalous_data. After having used this
# function, and while still in the R workspace, type:
# 
#   save(anomalous_data,file="data/anomalous_data.RData")
#
# The following commented lines are meant to use roxygen2 to
# build up documentation for anomalous_data


  
create_anom_dataframe <- function()
{
  # Vector with atom names
  anames <- c("H","He","Li","Be","B","C","N","O","F","Ne",
              "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
              "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Xe")
  
  anomalous_data <- list()
  for (chem_el in anames)
  {
   # Load specific dataset of f' and f"
   filename <- file.path("data-raw",paste(chem_el,".dat",sep=""))
   ad <- read.table(file=filename)
  
   # Turn energies (in ev) into wavelength (in angstroms)
   ad$V1 <- 12400/ad$V1
  
   # Reverse order (because of inverse relation between energy 
   # and wavelength)
   ad <- ad[order(ad[,1]),]
  
   # Appropriate names
   colnames(ad) <- c("lambda","f1","f2")
   
   anomalous_data <- c(anomalous_data,list(ad))
  }
  names(anomalous_data) <- anames
  
  return(anomalous_data)
}

# Use this function from within main Package directory to generate
# data file "atoms.RData" in data/

generate_atoms_dataframe <- function()
{
  # Vector with atom names
  anames <- c(" H","He","Li","Be"," B"," C"," N"," O"," F","Ne",
              "Na","Mg","Al","Si"," P"," S","Cl","Ar"," K","Ca",
              "Sc","Ti"," V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Xe")

  # Vector with atom numbers (Z)
  anumbers <- c(1:30,54)

  # Dataframe including atom names and numbers
  atoms <- data.frame(anames=anames,Z=anumbers)
  
  # Save as file Rda
  save(atoms,file="data/atoms.RData")
}

# Use this function to generate ksigma and anoK objects
# as internal data (in R/sysdata.rda)
generate_global_vars <- function()
{
  require(devtools)
  
  # k constant determining gaussian sigma
  ksigma <- 0.05
  
  # anoK to soften anomalous scattering effect
  anoK <- 0.3
  
  # Save all data
  use_data(ksigma,anoK,internal=TRUE)
}
