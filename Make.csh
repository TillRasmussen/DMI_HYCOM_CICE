#!/bin/csh
#
set echo
cd $cwd
#
# --- Usage:  ./Make.com >& Make.log
#
unset echo
module switch PrgEnv-cray PrgEnv-intel
module load cray-netcdf
module list
set echo
make hycom_cice_nuopc
