#!/bin/csh
#
# User settings:
# ###############
setenv compiler intel
#setenv compiler gnu
#setenv compiler cray

setenv debug true
#setenv debug false

# For tar:
set ESMFBLDDIR=/home/${user}/esmf/ESMF/esmf_v7_1_0r/
set CICEBLDDIR=/data/${user}/git_cice/build_cice_dmi/
set HYCOMBLDDIR=/data/${user}/git_hycom/build_hycom_dmi/
## For ang:
#set ESMFBLDDIR=/data/ang/sourcecode/coupled/esmf_v7_1-git/
#set CICEBLDDIR=/data/ang/sourcecode/coupled/build_cice_dmi/
#set HYCOMBLDDIR=/data/ang/sourcecode/coupled/build_hycom_dmi/

# End of user settings
# #########################

if (${compiler} == 'cray') then
  setenv esmf_compiler 'cce' 
else if (${compiler} == 'gnu') then
  setenv esmf_compiler 'gfortran'
  module switch PrgEnv-cray PrgEnv-gnu
else if (${compiler} == 'intel') then
  setenv esmf_compiler ${compiler}
  module switch PrgEnv-cray PrgEnv-intel
else
  echo Compiler "${compiler}" not found.
  exit 1
endif
module load cray-netcdf
if (${debug} == 'true') then 
        set echo
	setenv ESMFMKFILE ${ESMFBLDDIR}/esmf/lib/libg/Unicos.${esmf_compiler}.64.mpi.default/esmf.mk
	setenv CICE_LIB ${CICEBLDDIR}/${compiler}_dmi_debug/compile/
        setenv HYCOM_LIB ${HYCOMBLDDIR}/DMI_${compiler}_debug_coupled/
else
        setenv ESMFMKFILE ${ESMFBLDDIR}/esmf/lib/libO/Unicos.${esmf_compiler}.64.mpi.default/esmf.mk
        setenv CICE_LIB ${CICEBLDDIR}/${compiler}_dmi/compile/
        setenv HYCOM_LIB ${HYCOMBLDDIR}/DMI_${compiler}_coupled/
endif
cd $cwd
#
# --- Usage:  ./Make.com >& Make.log
#
unset echo
module list
set echo
make hycom_cice_nuopc
if (${debug} == 'true') then
  mv hycom_cice_nuopc hycom_cice_nuopc_${compiler}_debug
else
  mv hycom_cice_nuopc hycom_cice_nuopc_${compiler}
endif
