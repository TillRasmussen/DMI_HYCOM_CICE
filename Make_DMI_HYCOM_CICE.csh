#!/bin/tcsh

#if (! $?CLEANUP) set CLEANUP = 0
set CLEANUP = 0

# compiledir / sourcecode
setenv compiledir   ${cwd}    # compiledir
setenv SRCDIR       ${cwd}    # directory of source code

# location on compiled libraies
setenv ESMF_libdir  ../esmf/
setenv CICE_libdir  ../CICE/compile/
setenv HYCOM_libdir ../HYCOM-src/compile/

# default
set compiler = intel
set DEBUG = false

# Usage / input args
if ( `echo $argv | grep -c "\-h"` ) then
  echo "Usage: $0 {compiler} {DEBUG}"
  echo "\nDefault: $0 ${compiler} ${DEBUG}"
  echo
  echo "ESMF_libdir  : $ESMF_libdir"
  echo "CICE_libdir  : $CICE_libdir"
  echo "HYCOM_libdir : $HYCOM_libdir"
  exit 2
endif

if ($#argv > 0) set compiler = $1
if ($#argv > 1) set DEBUG = $2

echo "_____________________"
echo " ${compiler}"
echo "_____________________"

if ( `hostname | grep -Ec "boreas|aurora|x3000c0s5b0n0"` ) then  # TDS = x3000c0s5b0n0
  source ../../cfg/settings_modules.tcsh $compiler

else if ( `hostname | grep -Ec "thor|freyja"` ) then
  source /opt/modules/default/init/tcsh
  module rm PrgEnv-intel
  module rm PrgEnv-gnu
  module rm PrgEnv-cray
  switch ( ${compiler} )
    case intel:
    case gnu:
    case cray:
      breaksw
    default:
      echo "ERROR: Compiler not supported : "${compiler}
      exit 1
      breaksw
  endsw
  module add PrgEnv-${compiler}
  module add cray-netcdf
endif

module list


#-- ESMF compiler settings --
switch ( ${compiler} )
  case intel:
    set esmf_compiler = intel
    breaksw
  case gnu:
    set esmf_compiler = gfortran
    breaksw
  case cray:
    set esmf_compiler = cce
    breaksw
  default:
    echo "ERROR: Compiler not supported : "$compiler
    exit 1
    breaksw
endsw


#-- Working directory / source directory
if (${DEBUG} == 'true') then
  setenv ESMFMKFILE ${ESMF_libdir}/lib/libg/Unicos.${esmf_compiler}.64.mpi.default/esmf.mk
  setenv CICE_LIB   ${CICE_libdir}/${compiler}_nuopc_dmi_debug/
  setenv HYCOM_LIB  ${HYCOM_libdir}/${compiler}_nuopc_dmi_debug/
  setenv OUTDIR     ${compiledir}/compile/${compiler}_nuopc_dmi_debug/
  setenv CPPFLAGS   -Dlogerror 
  setenv CPPFLAGS   -Dlogmulti 
else
  setenv ESMFMKFILE ${ESMF_libdir}../esmf/lib/libO/Unicos.${esmf_compiler}.64.mpi.default/esmf.mk
  setenv CICE_LIB   ${CICE_libdir}/${compiler}_nuopc_dmi/
  setenv HYCOM_LIB  ${HYCOM_libdir}/${compiler}_nuopc_dmi/
  setenv OUTDIR     ${compiledir}/compile/${compiler}_nuopc_dmi/
  setenv CPPFLAGS   -Dnone
  setenv CPPFLAGS
endif

setenv HC_FILE hycom_${compiler}
setenv CICE_FILE cice_${compiler}

#-- Compile -------------------------
cd ${SRCDIR}
make clean
make hycom_cice_nuopc
if ($status != 0) then
  echo "COMPILE FAILED: ${0}"
  exit 99
else
  echo "*** COMPILE SUCCESSFUL ***"
  mkdir -p $OUTDIR
  \mv *.mod *.o $OUTDIR
  \mv hycom_cice_nuopc $OUTDIR
  echo "Binaries in $OUTDIR"
endif


#if ($CLEANUP == 1) then
#  \rm ${OUTDIR}/{*.o,*.mod}
#endif


