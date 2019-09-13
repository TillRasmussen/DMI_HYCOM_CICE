#!/bin/csh
#
# debug must be O (debug off) or g (debug on)
#set echo
setenv compiler intel
#debug should be true or false
setenv debug true
if (${compiler} == 'cray') then
  setenv esmf_compiler 'cce' 
else if (${compiler} == 'gnu') then
  setenv esmf_compiler 'gfortran'
  module switch PrgEnv-cray PrgEnv-gnu
else
  setenv esmf_compiler ${compiler}
  module switch PrgEnv-cray PrgEnv-intel
endif
module load cray-netcdf
if (${debug} == 'true') then 
        set echo
	setenv ESMFMKFILE /home/${user}/esmf/ESMF/esmf_v7_1_0r/esmf/lib/libg/Unicos.${esmf_compiler}.64.mpi.default/esmf.mk
	setenv CICE_LIB /data/${user}/git_cice/build_cice_dmi/${compiler}_dmi_debug/compile/
        setenv HYCOM_LIB /data/${user}/git_hycom/build_hycom_dmi/DMI_${compiler}_debug_coupled/
else
        setenv ESMFMKFILE /home/${user}/esmf/ESMF/esmf_v7_1_0r/esmf/lib/libO/Unicos.${esmf_compiler}.64.mpi.default/esmf.mk
        setenv CICE_LIB /data/${user}/git_cice/build_cice_dmi/${compiler}_dmi/compile/
        setenv HYCOM_LIB /data/${user}/git_hycom/build_hycom_dmi/DMI_${compiler}_coupled/
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
