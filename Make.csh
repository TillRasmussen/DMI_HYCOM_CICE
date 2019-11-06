#!/bin/csh -f
#
# debug must be O (debug off) or g (debug on)
#set echo
set oprdev =
if ( `echo $user | grep dev$` == $user ) set oprdev = dev
if ( `echo $user | grep opr$` == $user ) set oprdev = opr
echo $oprdev
source /opt/modules/default/init/tcsh
setenv src_dir ${cwd}
setenv compiledir /netapp/research/tar/esmf/hycom-cice-src/nuopc_dmi/git_hycom_cice_esmf_nuopc/build_DMI_HYCOM_CICE/
if ($#argv > 0) then
  set compiler = $1
else
  set compiler = intel
  set compiler = cray
  set compiler = gnu
endif
echo $compiler
setenv v_esmf v8_0_0 # alternatively v7_1_0r/esmf
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
echo ${compiler}
if (${debug} == 'true') then 
        set echo
	setenv ESMFMKFILE /home/${user}/esmf/ESMF/esmf_${v_esmf}/lib/libg/Unicos.${esmf_compiler}.64.mpi.default/esmf.mk
	setenv CICE_LIB /data/${user}/git_cice/build_cice_dmi/${compiler}_dmi_debug/compile/
        setenv HYCOM_LIB /data/${user}/git_hycom/build_hycom_dmi/DMI_${compiler}_debug_coupled/
else
        setenv ESMFMKFILE /home/${user}/esmf/ESMF/esmf_${v_esmf}/lib/libO/Unicos.${esmf_compiler}.64.mpi.default/esmf.mk
        setenv CICE_LIB /data/${user}/git_cice/build_cice_dmi/${compiler}_dmi/compile/
        setenv HYCOM_LIB /data/${user}/git_hycom/build_hycom_dmi/DMI_${compiler}_coupled/
endif
setenv HC_FILE hycom_${compiler}
setenv CICE_FILE cice_${compiler}
echo $HC_FILE
echo $CICE_FILE
echo $ESMFMKFILE
echo $CICE_LIB
echo $HYCOM_LIB
#
# --- Usage:  ./Make.com >& Make.log
#
unset echo
module list
set echo
#set pwd_dir ${PWD}
cd $src_dir
make hycom_cice_nuopc
if (${debug} == 'true') then
  mkdir -p ${compiledir}/nuopc_${compiler}_debug
  mv *.mod *.o ${compiledir}/nuopc_${compiler}_debug/
  mv hycom_cice_nuopc ${compiledir}/nuopc_${compiler}_debug
  #mv hycom_cice_nuopc hycom_cice_nuopc_${compiler}_debug
else
  mkdir ${compiledir}/nuopc_${compiler}
  mv *.mod *.o ${compiledir}/nuopc_${compiler}
  mv hycom_cice_nuopc ${compiledir}/nuopc_${compiler}
#  mv hycom_cice_nuopc hycom_cice_nuopc_${compiler}
endif
unset echo

