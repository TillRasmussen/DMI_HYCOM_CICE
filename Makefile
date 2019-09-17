# GNU Makefile template for user ESMF application

################################################################################
################################################################################
## This Makefile must be able to find the "esmf.mk" Makefile fragment in the  ##
## 'include' line below. Following the ESMF User's Guide, a complete ESMF     ##
## installation should ensure that a single environment variable "ESMFMKFILE" ##
## is made available on the system. This variable should point to the         ##
## "esmf.mk" file.                                                            ##
##                                                                            ##
## This example Makefile uses the "ESMFMKFILE" environment variable.          ##
##                                                                            ##
## If you notice that this Makefile cannot find variable ESMFMKFILE then      ##
## please contact the person responsible for the ESMF installation on your    ##
## system.                                                                    ##
## As a work-around you can simply hardcode the path to "esmf.mk" in the      ##
## include line below. However, doing so will render this Makefile a lot less ##
## flexible and non-portable.                                                 ##
################################################################################

#ifneq ($(origin ESMFMKFILE), environment)
#$(error Environment variable ESMFMKFILE was not set.)
#endif

# include /data/ang/sourcecode/coupled/esmf_v7_1-git/esmf/lib/libg/Unicos.intel.64.mpi.default/esmf.mk 
# LANLCICEDIR= /data/ang/sourcecode/coupled/build_cice_dmi/intel_dmi_debug/compile
# HYCOMDIR=/data/ang/sourcecode/coupled/build_hycom_dmi/DMI_intel_debug_coupled
# UTILINCS        = -I$(LANLCICEDIR) -L$(LANLCICEDIR) -lcice_intel -I$(HYCOMDIR) -L$(HYCOMDIR) -lhycom_intel
include $(ESMFMKFILE)
UTILINCS        = -I$(CICE_LIB) -L$(CICE_LIB) -lcice_$(compiler) -I$(HYCOM_LIB) -L$(HYCOM_LIB) -lhycom_$(compiler)
################################################################################
################################################################################

.SUFFIXES: .f90 .F90 .c .C

%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(UTILINCS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) $<

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(UTILINCS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $<
        
%.o : %.c
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<

%.o : %.C
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<


# -----------------------------------------------------------------------------
hycom_cice_nuopc: hycom_cice_nuopc.o esm.o cice_cap.o hycom_cap.o conn.o mod_cb_arrays_nuopc_glue.o mod_hycom_nuopc_glue.o mod_nuopc_options.o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(UTILINCS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS) -lcice_${compiler} -lhycom_${compiler}
# module dependencies:
esmApp.o: esm.o
esm.o: cice_cap.o hycom_cap.o conn.o mod_nuopc_options.o 
hycom_cap.o: mod_cb_arrays_nuopc_glue.o mod_hycom_nuopc_glue.o
mod_hycom_nuopc_glue.o: mod_cb_arrays_nuopc_glue.o
hycom_cice_nuopc.o: esm.o
mod_nuopc_options.o:
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
.PHONY: dust clean distclean info edit
dust:
	rm -f PET*.ESMF_LogFile
clean:
	rm -f esmApp hycom_cice_nuopc *.o *.mod
distclean: dust clean

info:
	@echo ==================================================================
	@echo ESMFMKFILE=$(ESMFMKFILE)
	@echo ==================================================================
	@cat $(ESMFMKFILE)
	@echo ==================================================================

edit:
	nedit esmApp.F90 esm.F90 atm.F90 ocn.F90 conn.F90 &
