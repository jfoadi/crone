#
# This file is part of the package crone.
#
# It has been added to avoid NOTES like:
#   read_x: no visible binding for global variable 'atoms'
#   strufac: no visible binding for global variable 'anomalous_data'
#
utils::globalVariables(c("atoms", "anomalous_data"))