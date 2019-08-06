module hycom_cap

  !-----------------------------------------------------------------------------
  ! OCN Component for CESM-BETA; use ESMF and  NUOPC 
  !-----------------------------------------------------------------------------



















  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_Advance   => label_Advance
  
  use MOD_HYCOM, only : HYCOM_Init, HYCOM_Run, HYCOM_Final, &
    end_of_run, end_of_run_cpl
    
  use mod_hycom_nuopc_glue
  use mod_cb_arrays_nuopc_glue

  implicit none
  
  private
  
  ! private internal state to keep instance data
  type InternalStateStruct
    type(hycom_nuopc_glue_type)   :: glue
    integer                       :: slice
  end type

  type InternalState
    type(InternalStateStruct), pointer :: wrap
  end type
  
  public SetServices

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: canonicalUnits
    character(len=64) :: transferOffer
    logical           :: assoc    ! is the farrayPtr associated with internal data
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter :: fldsMax = 50
  integer :: fldsToOcn_num = 0
  type (fld_list_type) :: fldsToOcn(fldsMax)
  integer :: fldsFrOcn_num = 0
  type (fld_list_type) :: fldsFrOcn(fldsMax)
!  integer, parameter, public        :: number_import_fields = 30
!  integer, parameter, public        :: number_export_fields = 8

  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(ocn:SetServices)'
    write(6,*) subname
    rc = ESMF_SUCCESS

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(6,*) 'init done HYCOM'
    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Model Advance step 2
#ifdef TARNOTNEEDED
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
      userRoutine=routine_Run2, phaseLabelList=(/"OCN_RUN_PHASE"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#ENDIF
!!!FINALIZE TILL - HYCOM END
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    call HYCOM_FieldsSetup()
    ! set Component name so it becomes identifiable
    call ESMF_GridCompSet(gcomp, name="HYCOM", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call HYCOM_AdvertiseFields( importState, fldsToOcn_num, fldsToOcn, tag='HYCOM import', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call HYCOM_AdvertiseFields(exportState, fldsFrOcn_num, fldsFrOcn, tag='HYCOM export', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
      
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! This method assumes that the incoming clock indicates the full simulation
    ! interval between startTime and stopTime.
    
    ! local variables    
    type(ESMF_Field)            :: field
    type(ESMF_Grid)             :: gridIn
    type(ESMF_Grid)             :: gridOut
    type(ESMF_array)            :: array
    type(ESMF_VM)               :: vm
    integer                     :: mpiComm
    TYPE(ESMF_Time)             :: startTime, stopTime, hycomRefTime, currTime
    TYPE(ESMF_TimeInterval)     :: interval,timeStep
    real(ESMF_KIND_R8)          :: startTime_r8, stopTime_r8, l_startTime_r8
    type(InternalState)         :: is
    integer                     :: stat
    type(ESMF_CALKIND_FLAG)     :: calkind
    logical                     :: restFlag = .false.        ! initial/restart run (F/T)
    character(len=32)           :: starttype
    character(len=80)           :: pointer_filename          ! restart pointer file !!Alex
    logical                     :: restart_write = .false.   ! write restart
    
    rc = ESMF_SUCCESS
    
    ! Allocate memory for the internal state and set it in the Component.
    allocate(is%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of the internal state memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridCompSetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Prepare to call into HYCOM_Init
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, mpiCommunicator=mpiComm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! Translate currTime and stopTime into HYCOM format
    call ESMF_ClockGet(clock, currTime=currTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Define HYCOM Ref Time
    call ESMF_TimeSet(hycomRefTime, yy=1901, mm=01, dd=01, calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
!TAR    call ESMF_TimeSet(hycomRefTime, yy=0001, mm=01, dd=01, calkindflag=ESMF_CALKIND_NOLEAP, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    interval = currTime - hycomRefTime
    call ESMF_TimeIntervalGet(interval, d_r8=startTime_r8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    interval = stopTime - hycomRefTime
    call ESMF_TimeIntervalGet(interval, d_r8=stopTime_r8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    print *, " HYCOM_INIT -->> startTime_r8=", startTime_r8, "stopTime_r8=", stopTime_r8

    ! get coupling frequency from ocean clock for 1st export
     call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_TimeIntervalGet(timeStep, d_r8=ocn_cpl_frq, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return  ! bail out
       
    ! Get start type run : start-up or continuous run
    ! seq_infodata_mod
!    call ESMF_AttributeGet(exportState, name="start_type", value=starttype, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!    return ! bail out
! This used to be parameters from CESM module seq_infodata_mod
    starttype="continue"
    if (     trim(starttype) == trim("startup")) then
       l_startTime_r8=-startTime_r8
       restFlag = .false.
    else if (trim(starttype) == trim("continue")) then
       l_startTime_r8=startTime_r8
       restFlag = .true.
    else if (trim(starttype) == trim("branch")) then
       l_startTime_r8=startTime_r8
       restFlag = .true.
    else
       call ESMF_LogWrite('hycom_nuopc ERROR: unknown starttype',  &
          ESMF_LOGMSG_ERROR, rc=rc)
       rc = ESMF_RC_OBJ_BAD
       return
    end if

    ! Get the pointer restart name !!Alex
      restart_write = .true.

    call ESMF_LOGWRITE("BEFORE HYCOM_INIT", ESMF_LOGMSG_INFO, rc=rc)
    
    ! Call into the HYCOM initialization  
!    call HYCOM_Init(mpiComm, & ! -->> call into HYCOM <<--
!       hycom_start_dtg=l_startTime_r8, hycom_end_dtg=stopTime_r8, &
!       pointer_filename=pointer_filename,  restart_write=restart_write)
    call HYCOM_Init(mpiComm, & ! -->> call into HYCOM <<--
       hycom_start_dtg=l_startTime_r8, hycom_end_dtg=stopTime_r8, &
        restart_write=restart_write)

    call ESMF_LOGWRITE("AFTER HYCOM_INIT", ESMF_LOGMSG_INFO, rc=rc)
    ! Fill in the glue structure.
    call HYCOM_GlueInitialize(is%wrap%glue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Write some HYCOM distribution info into the Log.
    call HYCOM_TileInfo(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! use the HYCOM Grid that was setup inside of the glue structure
    gridIn  = is%wrap%glue%grid ! for imported Fields
    gridOut = is%wrap%glue%grid ! for exported Fields
    call HYCOM_RealizeFields(importState, gridIn , fldsToOcn_num, fldsToOcn, "Ocean import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call HYCOM_RealizeFields(exportState, gridOut, fldsFrOcn_num, fldsFrOcn, "Ocean export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

#ifdef TARNOTNEEEDED    
    ! conditionally realize or remove Fields in import and export States
    ! also keep track of these Fields on the glue layer
    ! importable fields:
    call HYCOM_GlueFieldsRealize(is%wrap%glue, importState, &
      StandardNames=(/ &
      "sea_ice_area_fraction                  ",    & ! from SEA-ICE
      "downward_x_stress_at_sea_ice_base      ",    & ! from SEA-ICE
      "downward_y_stress_at_sea_ice_base      ",    & ! from SEA-ICE
      "downward_sea_ice_basal_solar_heat_flux ",    & ! from SEA-ICE
      "upward_sea_ice_basal_heat_flux         ",    & ! from SEA-ICE
      "downward_sea_ice_basal_salt_flux       ",    & ! from SEA-ICE
      "downward_sea_ice_basal_water_flux      ",    & ! from SEA-ICE
      "sea_ice_temperature                    ",    & ! from SEA-ICE
      "sea_ice_thickness                      ",    & ! from SEA-ICE
      "sea_ice_x_velocity                     ",    & ! from SEA-ICE
      "sea_ice_y_velocity                     "/),  & ! from SEA-ICE
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !! exportable fields:
    call HYCOM_GlueFieldsRealize(is%wrap%glue, exportState, &
      StandardNames=(/ &
      "sea_surface_temperature                  ",    &
      "upward_sea_ice_basal_available_heat_flux ",    &
      "sea_lev                                  ",    &
      "mixed_layer_depth                        ",    &
      "s_surf                                   ",    &
      "eastward_sea_surface_slope               ",    &
      "northward_sea_surface_slope              ",    &
      "ocn_current_zonal                        ",    &
      "ocn_current_merid                        "/),  &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Import data to HYCOM native structures through glue fields.
   call HYCOM_GlueFieldsDataImport(is%wrap%glue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Export HYCOM native structures to data through glue fields.
    CALL HYCOM_GlueFieldsDataExport(is%wrap%glue, .not. restFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Reset the slice counter
    is%wrap%slice = 1

    ! duplicate the mpi communicator from the current VM 
    call MPI_Comm_dup(mpicomm, mpicom_ocn, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! satisfy some external needs
    call ESMF_AttributeGet(export_state, name="ID", value=OCNID, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridGet(is%wrap%glue%grid, distgrid=distgrid2D, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_DistGridGet(distgrid2D, 0, elementCount=n_elem, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    allocate(fptrSeqIndex(n_elem))
    call ESMF_DistGridGet(distgrid2D, 0, seqIndexList=fptrSeqIndex, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    distgrid = ESMF_DistGridCreate(fptrSeqIndex, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    
    deallocate(fptrSeqIndex)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_DistGridGet(distgrid, delayout=delayout, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_DelayoutGet(delayout, localDeCount=ldeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !print *, 'HYCOM DG DELAYOUT localDECount: ', ldeCount

    lsize = 0
    do lde = 0, ldeCount-1
      call ESMF_DistGridGet(distgrid, lde, elementCount=eleCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      lsize = lsize + eleCount
    enddo

    !print *, 'HYCOM DG DELAYOUT lsize: ', lsize

    !-----------------------------------------
    ! Create dom 
    !-----------------------------------------

    nfields = shr_string_listGetNum(trim(seq_flds_dom_fields))

    dom = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="domain", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(dom, name="mct_names", value=trim(seq_flds_dom_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Set values of dom (needs ocn initialization info)

    call ocn_domain_esmf(dom, is%wrap%glue%grid, rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   
    !----------------------------------------- 
    !  Create o2x 
    !-----------------------------------------

    ! 1d undistributed index of fields, 2d is packed data

    nfields = shr_string_listGetNum(trim(seq_flds_o2x_fields))

    o2x = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="d2x", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(o2x, name="mct_names", value=trim(seq_flds_o2x_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------------------- 
    !  Create x2o 
    !-----------------------------------------

    nfields = shr_string_listGetNum(trim(seq_flds_x2o_fields))

    x2o = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="x2d", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(x2o, name="mct_names", value=trim(seq_flds_x2o_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------------------- 
    ! Add esmf arrays to import and export state 
    !-----------------------------------------

    call ESMF_StateAdd(export_state, (/dom/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(export_state, (/o2x/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    call ESMF_StateAdd(import_state, (/x2o/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    allocate(gsmap_o)
    allocate(dom_o)
   
    call esmf2mct_init(distgrid, OCNID, gsmap_o, mpicom_ocn, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_init(dom, dom_o, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !print *, 'HYCOM GSMAP gsize: ', mct_gsMap_gsize(gsmap_o)

    call ESMF_GridGet(is%wrap%glue%grid, distgrid=distgrid2D, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_DistGridGet(distgrid2D, maxIndexPTile=maxIndex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="gsize", value=mct_gsMap_gsize(gsmap_o), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="ocn_prognostic", value=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="ocnrof_prognostic", value=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="ocn_nx", value=maxIndex(1,1), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="ocn_ny", value=maxIndex(2,1), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !! Create and Realize Importable fields
    call esmfshr_nuopc_create_fields( &
      ocn_import_fields, meshIn, importState, tag='HYCOM import', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    !! Create and Realize Exportable fields
    call esmfshr_nuopc_create_fields( &
      ocn_export_fields, meshOut, exportState, tag='HYCOM export', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

#endif
    
  end subroutine

  subroutine HYCOM_RealizeFields(state, grid, nfields, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Grid), intent(in)                 :: grid
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc

    integer                                     :: i
    type(ESMF_Field)                            :: field
    integer                                     :: npet, nx, ny, pet, elb(2), eub(2), clb(2), cub(2), tlb(2), tub(2)
    type(ESMF_VM)                               :: vm
    character(len=*),parameter  :: subname='(hycom_cap:HYCOM_RealizeFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields
      write(6,*) field_defs(i)%stdname
      write(6,*) field_defs(i)%shortname
      if (field_defs(i)%assoc) then
        write(6, *) subname, tag, ' Field ', field_defs(i)%shortname, ':', &
          lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
          lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2), &
          lbound(field_defs(i)%farrayPtr,3), ubound(field_defs(i)%farrayPtr,3)
        !call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
        
        field = ESMF_FieldCreate(grid=grid, &
          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then
        call NUOPC_Realize(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=rc)
      else
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is not connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=rc)
        call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

    enddo

  end subroutine HYCOM_RealizeFields

  
  !-----------------------------------------------------------------------------

  SUBROUTINE ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep_O
    type(ESMF_Time)             :: hycomRefTime
    type(ESMF_TimeInterval)     :: interval
    real(ESMF_KIND_R8)          :: endTime_r8    ! end of coupling sequence
    type(InternalState)         :: is
    logical                     :: initFlag
    type(ESMF_CALKIND_FLAG)     :: calkind

    integer                     :: fieldCount, i
    character(len=128), allocatable :: fieldNameList(:)
    type(ESMF_Field)            :: field
    character(len=80)           :: pointer_filename     ! restart pointer file !!Alex
    logical                     :: restart_write = .false.
    character(len=32)           :: starttype            ! infodata start type
    logical                     :: restFlag = .false.
    character(len=128)          :: msg
    character*80 :: filenc !!Alex

    rc = ESMF_SUCCESS
    write(6,*) 'HYCOM advance'   
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Get the internal state from Component.
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    ! on the internal Clock object. The NUOPC Layer will update the Clock 
    ! automatically each time before entering ModelAdvance(), but the HYCOM
    ! model must be stepped forward within this method.
    CALL ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep_O,  rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Translate currTime + timeStep into HYCOM format
    call ESMF_TimeSet(hycomRefTime, yy=1901, mm=01, dd=01, calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
!TAR    call ESMF_TimeSet(hycomRefTime, yy=0001, mm=01, dd=01, calkindflag=ESMF_CALKIND_NOLEAP, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! get endtime 
    interval = currTime - hycomRefTime
    call ESMF_TimeIntervalGet(interval, d_r8=endTime_r8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    endTime_r8 = endTime_r8 + ocn_cpl_frq
    
    call ESMF_FieldBundleGet(is%wrap%glue%importFields, fieldCount=fieldCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    write(6,*) 'FIELDBUNDLEPRINT: -> number of field in HYCOM import FieldBundle: ', fieldcount
    call flush(6)
    allocate(fieldNameList(fieldCount))
    call ESMF_FieldBundleGet(is%wrap%glue%importFields, fieldNameList=fieldNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    do i = 1, fieldCount
!!Alex      print *, 'FIELDBUNDLEPRINT: -> ', fieldNameList(i)
      call ESMF_FieldBundleGet(is%wrap%glue%importFields, fieldName=fieldNameList(i), field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      call ESMF_AttributeSet(field, name="HYCOM_CICE_Connected", value=.false., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
    end do
    deallocate(fieldNameList)

    ! Redistribute 1D CESM import Fields to HYCOM Fields stored in glue import fieldbundle
!    do i  = 1, number_import_fields
!      if(cesm2hycom_table(i)%connected) then
!        call HYCOM_RedistCESM2HYCOM(importState, cesm2hycom_table(i)%cesm_stdname, &
!          is%wrap%glue%importFields, cesm2hycom_table(i)%hycom_stdname, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
!      endif
!    enddo
    
    ! Get start type run : start-up or continuous run !!Alex add restFlag
    !call ESMF_AttributeGet(exportState, name="start_type", value=starttype, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out
    starttype='continue'
    if (     trim(starttype) == trim("startup")) then
!       l_startTime_r8=-startTime_r8
       restFlag = .false.
    else if (trim(starttype) == trim("continue")) then
!       l_startTime_r8=startTime_r8
       restFlag = .true.
    else if (trim(starttype) == trim("branch")) then
!       l_startTime_r8=startTime_r8
       restFlag = .true.
    else
       call ESMF_LogWrite('hycom_nuopc ERROR: unknown starttype',  &
          ESMF_LOGMSG_ERROR, rc=rc)
       rc = ESMF_RC_OBJ_BAD
       return
    end if

    !TODO: don't need the additional initialization step once data-dependency
    !TODO: is taken care of during initialize.
    initFlag = .false.
    if (is%wrap%slice==1 .and. (.not. restFlag)) initFlag = .true.
    
    
    ! Import data to HYCOM native structures through glue fields.
    call HYCOM_GlueFieldsDataImport(is%wrap%glue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    
    !call ESMF_VMLogMemInfo('MEMORY Usage BEFORE HYCOM_RUN', rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out
  
    ! Get the pointer restart name 
!    pointer_filename = 'rpointer.ocn' 
!    restart_write = seq_timemgr_RestartAlarmIsOn(ccsm_EClock_o)

    ! Enter the advancing loop over HYCOM_run...
!    do
      ! ...on return the end-of-run flags indicate whether HYCOM has advanced
      ! far enough...
      CALL HYCOM_Run(endtime=endTime_r8,pointer_filename=pointer_filename, restart_write=restart_write) ! -->> call into HYCOM <<--
!      if (end_of_run .or. end_of_run_cpl) exit
!    enddo

    !call ESMF_VMLogMemInfo('MEMORY Usage AFTER HYCOM_RUN', rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    ! Export HYCOM native data through the glue fields.
    call HYCOM_GlueFieldsDataExport(is%wrap%glue, initFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! advance the time slice counter
    is%wrap%slice = is%wrap%slice + 1

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Finalize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    type(InternalState)  :: is
    integer              :: stat

    rc = ESMF_SUCCESS
  
    ! Get the internal state from Component.
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! TODO: Destroy objects inside of internal state.

    ! Deallocate the internal state memory.
    deallocate(is%wrap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of internal state memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


      
  end subroutine

  !-----------------------------------------------------------------------------
#ifdef TARNOTNEEDED
  subroutine routine_Run2(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)              :: compclock
    logical                       :: ocn_present
    logical                       :: ocnrun_alarm
    logical                       :: tight_coupling
    type(ESMF_Array)              :: d2x

    rc = ESMF_SUCCESS
    tight_coupling = .false.

    call ESMF_GridCompGet(gcomp, clock=compclock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(importState, name="ocean_tight_coupling", &
      value=tight_coupling, defaultvalue=.false., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(exportState, &
      name="ocn_present", value=ocn_present, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(exportState, &
      name="ocnrun_alarm", value=ocnrun_alarm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (ocn_present  .and.  ocnrun_alarm  .and.  tight_coupling) then
      call ESMF_LogWrite(trim('HYCOM RUN2 --->'), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      call HYCOM_ModelAdvance(gcomp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
    endif

  end subroutine
  subroutine HYCOM_WriteFieldBundle(grid, fieldbundle, fieldNameList, filePrefix, overwrite, &
    status, timeslice, relaxedflag, rc)
    type(ESMF_Grid),            intent(in)            :: grid
    type(ESMF_FieldBundle),     intent(in)            :: fieldbundle
    character(len=*),           intent(in),  optional :: fieldNameList(:)
    character(len=*),           intent(in),  optional :: filePrefix
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    logical,                    intent(in),  optional :: relaxedflag
    integer,                    intent(out), optional :: rc
  !-----------------------------------------------------------------------------
    ! local variables
    integer                         :: i, itemCount
    type(ESMF_Field)                :: field
    type(ESMF_StateItem_Flag)       :: itemType
    character(len=80)               :: fileName
    character(len=80), allocatable  :: fieldNameList_loc(:)
    type(ESMF_Field)                :: dst2DField
    real(ESMF_KIND_R8), pointer     :: ptr1D(:), ptr2D(:,:)

    if (present(rc)) rc = ESMF_SUCCESS

    if (present(fieldNameList)) then
      allocate(fieldNameList_loc(size(fieldNameList)))
      do i=1, size(fieldNameList)
        fieldNameList_loc(i) = trim(fieldNameList(i))
      enddo
    else
      call ESMF_FieldBundleGet(fieldbundle, fieldCount=itemCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(fieldNameList_loc(itemCount))
      call ESMF_FieldBundleGet(fieldbundle, fieldNameList=fieldNameList_loc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    do i=1, size(fieldNameList_loc)
      call ESMF_FieldBundleGet(fieldbundle, fieldName=fieldNameList_loc(i), field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=FILENAME)) &
        return  ! bail out
      ! -> output to file
      if (present(filePrefix)) then
        write (fileName,"(A)") filePrefix//trim(fieldNameList_loc(i))//".nc"
      else
        write (fileName,"(A)") trim(fieldNameList_loc(i))//".nc"
      endif

      call ESMF_FieldGet(field, farrayPtr=ptr2D, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      
      !call ESMF_LogWrite(trim('HYCOM_WriteFieldBundle: '//fieldNameList_loc(i)), ESMF_LOGMSG_INFO, rc=rc)
      !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      !  line=__LINE__, &
      !  file=__FILE__)) &
      !return ! bail out

!!Alex incompatibility ESMF 7.0.0      call ESMF_FieldWrite(field, file=trim(fileName), variableName=trim(fieldNameList_loc(i)), &
!        overwrite=overwrite, status=status, timeslice=timeslice, rc=rc)
      call ESMF_FieldWrite(field, fileName=trim(fileName), variableName=trim(fieldNameList_loc(i)), &
        overwrite=overwrite, status=status, timeslice=timeslice, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=FILENAME)) &
        return  ! bail out

    enddo

    deallocate(fieldNameList_loc)

  end subroutine

  subroutine RedistAndWriteField(grid, state, fieldNameList, filePrefix, overwrite, &
    status, timeslice, relaxedflag, rc)
    type(ESMF_Grid),            intent(in)            :: grid
    type(ESMF_State),           intent(in)            :: state
    character(len=*),           intent(in),  optional :: fieldNameList(:)
    character(len=*),           intent(in),  optional :: filePrefix
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    logical,                    intent(in),  optional :: relaxedflag
    integer,                    intent(out), optional :: rc
  !-----------------------------------------------------------------------------
    ! local variables
    integer                         :: i, itemCount, elb(2), eub(2), j, k
    type(ESMF_Field)                :: field
    type(ESMF_StateItem_Flag)       :: itemType
    character(len=80)               :: fileName, msg
    character(len=80), allocatable  :: fieldNameList_loc(:)
    type(ESMF_Field)                :: dst2DField
    real(ESMF_KIND_R8), pointer     :: ptr1D(:), ptr2D(:,:)

    if (present(rc)) rc = ESMF_SUCCESS

    dst2DField = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if (present(fieldNameList)) then
      allocate(fieldNameList_loc(size(fieldNameList)))
      do i=1, size(fieldNameList)
        fieldNameList_loc(i) = trim(fieldNameList(i))
      enddo
    else
      call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(fieldNameList_loc(itemCount))
      call ESMF_StateGet(state, itemNameList=fieldNameList_loc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    do i=1, size(fieldNameList_loc)
      call ESMF_StateGet(state, itemName=fieldNameList_loc(i), &
        itemType=itemType, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=FILENAME)) &
        return  ! bail out
      if (itemType == ESMF_STATEITEM_FIELD) then
        ! field is available in the state
        call ESMF_StateGet(state, itemName=fieldNameList_loc(i), field=field, &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=FILENAME)) &
          return  ! bail out
        ! -> output to file
        if (present(filePrefix)) then
          write (fileName,"(A)") filePrefix//trim(fieldNameList_loc(i))//".nc"
        else
          write (fileName,"(A)") trim(fieldNameList_loc(i))//".nc"
        endif

        call ESMF_FieldGet(field, farrayPtr=ptr1D, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
        return ! bail out

        call ESMF_FieldGet(dst2Dfield, farrayPtr=ptr2D, exclusiveLBound=elb, exclusiveUBound=eub, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
        return ! bail out
        ptr2D = 0.0_ESMF_KIND_R8

        write(msg, *) elb, eub, lbound(ptr1D), ubound(ptr1D)
        !call ESMF_LogWrite(trim('RedistAndWriteField: bounds: ')// trim(msg), ESMF_LOGMSG_INFO, rc=rc)
        !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        !  line=__LINE__, &
        !  file=__FILE__)) &
        !return ! bail out

        !call ESMF_FieldRedist(field, dst2DField, routehandle=CESM2HYCOM_RHR8, rc=rc)
!        call copy_1D_to_2D(field, dst2DField, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!        return ! bail out
        
        !call ESMF_LogWrite(trim('RedistAndWriteField: '//fieldNameList_loc(i)), ESMF_LOGMSG_INFO, rc=rc)
        !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        !  line=__LINE__, &
        !  file=__FILE__)) &
        !return ! bail out

!!Alex incompatibility ESMF 7.0.0        call ESMF_FieldWrite(dst2DField, file=trim(fileName), variableName=trim(fieldNameList_loc(i)), &
!          overwrite=overwrite, status=status, timeslice=timeslice, rc=rc)
        call ESMF_FieldWrite(dst2DField, fileName=trim(fileName), variableName=trim(fieldNameList_loc(i)), &
          overwrite=overwrite, status=status, timeslice=timeslice, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=FILENAME)) &
          return  ! bail out

      endif
    enddo

    deallocate(fieldNameList_loc)

    call ESMF_FieldDestroy(dst2DField, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

  end subroutine

!***********************************************************************
!BOP
! !IROUTINE: ocn_domain_esmf
! !INTERFACE:

 subroutine ocn_domain_esmf( dom, grid, rc)

    use iso_c_binding

! !DESCRIPTION:
!  This routine creates the ocean domain and necessary communication routehandles
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

    implicit none
    type(ESMF_Array), intent(inout)     :: dom      ! CESM DOMAIN INFO
    type(ESMF_Grid),  intent(in)        :: grid     ! Native HYCOM 2D Grid
    integer, intent(out)                :: rc

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


    integer ::   &
      i,j, n, iblock

    integer ::   &
      klon,klat,karea,kmask,kfrac ! domain fields

    real(ESMF_KIND_R8),    pointer ::  &
      fptr (:,:)          ! data pointer into ESMF array

    real(ESMF_KIND_R8)  :: &
      frac                ! temporary var to compute frac/mask from KMT

    type(ESMF_DistGrid)  :: distgrid, distgrid2d
    type(ESMF_VM)        :: vm
    character(len=256)   :: msg
    integer              :: n_elem, n_pet, lpet, k
    integer, pointer     :: indexlist(:)
    logical              :: arbIndexFlag
    type(ESMF_Array)     :: lon1d, lat1d, area1d, mask1d
    type(ESMF_Array)     :: plon, plat, area, mask, area2d, mask2d
    integer              :: elb(2,1), eub(2,1), elb1(1,1), eub1(1,1)
    real(ESMF_KIND_R8), pointer  :: tlon(:), tlat(:), tarea(:), fptrLon(:,:), fptrLat(:,:)
    integer(ESMF_KIND_I4), pointer :: tmask(:), fptrSeqIndex(:)
    real(ESMF_KIND_R8)   :: radian, radius, pi
    type(ESMF_TYPEKIND_FLAG)  :: tkf

    type(ESMF_Array)     :: dummy1D, dummy2D
    real(ESMF_KIND_R8), pointer     :: fptr2D(:,:), fptr1D(:), fptr2D_new(:,:)

    type(ESMF_RouteHandle)          :: redist_padding_rh

!-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Retrieve domain data pointer
    call ESMF_ArrayGet(dom, localDe=0, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Retrieve the HYCOM Grid coordinate and mask arrays for reference      
    call ESMF_GridGetCoord(grid, coordDim=1, array=plon, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(grid, coordDim=2, array=plat, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_MASK, array=mask, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(mask, typekind=tkf, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_AREA, array=area, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGet(grid, distgrid=distgrid2d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    area2d = ESMF_ArrayCreate(distgrid2d, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    mask2d = ESMF_ArrayCreate(distgrid2d, typekind=ESMF_TYPEKIND_I4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_VMGet(vm, petCount=n_pet, localPet=lpet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Use the mesh based 1D distgrid to create DOM elements
    call ESMF_ArrayGet(dom, distgrid=distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    lon1D = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(lon1D, farrayPtr = tlon, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    lat1D = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(lat1D, farrayPtr = tlat, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    area1D = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(area1D, farrayPtr = tarea, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    mask1D = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_I4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(mask1D, farrayPtr = tmask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! CESM uses 1 DE per PET
    call ESMF_DistGridGet(distgrid, 0, arbSeqIndexFlag=arbIndexFlag, elementCount=n_elem, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    allocate(indexList(n_elem))
    call ESMF_DistGridGet(distgrid, 0, seqIndexList=indexlist, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

!-------------------------------------------------------------------
!
!  initialize domain type, lat/lon in degrees,
!  area in radians^2, mask is 1 (ocean), 0 (non-ocean)
!  Fill in correct values for domain components
!
!-------------------------------------------------------------------

    klon  = esmfshr_util_ArrayGetIndex(dom,'lon ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    klat  = esmfshr_util_ArrayGetIndex(dom,'lat ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    karea = esmfshr_util_ArrayGetIndex(dom,'area',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kmask = esmfshr_util_ArrayGetIndex(dom,'mask',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kfrac = esmfshr_util_ArrayGetIndex(dom,'frac',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    fptr(:,:) = -9999.0_ESMF_KIND_R8
    n=0

    write(msg, *) 'DUMPING HYCOM INDICES BEGINS:'
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    write(msg, *) 'total number of ocean pet', n_pet, ' local pet number', lpet
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    write(msg, *) 'number of elements on this pet:', n_elem, ' arbflag', arbIndexFlag
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    write(msg, *) 'lpet ', 'n ', 'Index ', 'lon ', 'lat ', &
      'area ', 'frac ', 'mask'
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(plon, exclusiveLBound=elb, exclusiveUBound=eub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(msg, *) 'src shape: ', elb, eub, ' dst shape: ', lbound(fptr), ubound(fptr)
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    write(msg, *) 'plon: ', elb, eub
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(lon1d, exclusiveLBound=elb1, exclusiveUBound=eub1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(msg, *) 'lon1d: ', elb1, eub1
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayRedistStore(plon, lon1d, routehandle=HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedistStore(mask, mask1d, routehandle=HYCOM2CESM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedist(plon, lon1d, routehandle=HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedist(plat, lat1d, routehandle=HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedist(mask, mask1d, routehandle=HYCOM2CESM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedist(area, area1d, routehandle=HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    pi = 3.14159265358
    radian = 180.0_ESMF_KIND_R8/pi
    radius    = 6370.0e3_ESMF_KIND_R8
    do n = elb1(1,1), eub1(1,1)
      fptr(klon , n)          = TLON(n)
      fptr(klat , n)          = TLAT(n)
      fptr(karea, n)          = TAREA(n)/radius/radius
      frac                    = TMASK(n)
      if (frac > 1.0_ESMF_KIND_R8) frac = 1.0_ESMF_KIND_R8
      fptr(kfrac, n)          = frac
      fptr(kmask, n)          = frac
      !write(msg, '(I4,A1,I8,A7,I8,2F10.3,E15.7,2F10.3)') lpet, ' ', n, ' INDEX=',indexlist(n), fptr(klon, n), &
      !  fptr(klat , n), fptr(karea, n), fptr(kfrac, n), fptr(kmask, n)
      !call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
      !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    enddo

    write(msg, *) 'DUMPING HYCOM INDICES ENDS:'
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(plon, farrayPtr=fptrLon, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_ArrayGet(plat, farrayPtr=fptrLat, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_DistGridGet(distgrid2D, 0, arbSeqIndexFlag=arbIndexFlag, elementCount=n_elem, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    allocate(fptrSeqIndex(n_elem))
    call ESMF_DistGridGet(distgrid2D, 0, seqIndexList=fptrSeqIndex, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    write(msg, *) 'HYCOM 2D distribution is arbitrary? ', arbIndexFlag, n_elem
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !do i = elb(1,1), eub(1,1)
    !  do j = elb(2,1), eub(2,1) 
    !    write(msg, '(2I6,2F10.3)') i, j, fptrLon(i,j), fptrLat(i,j)
    !    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    !    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    !  enddo
    !enddo

    !do n = 1, n_elem
    !  write(msg, '(I6,I8)') n, fptrSeqIndex(n)
    !  call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    !  if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    !enddo

    ! Release these routehandles because they have paddings(halo) in hycom memory
    call ESMF_RouteHandleRelease(HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_RouteHandleRelease(HYCOM2CESM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    ! Recompute routehandles that has no paddings. 

    call ESMF_ArrayRedistStore(area1d, area2d, routehandle=CESM2HYCOM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedistStore(mask1D, mask2d, routehandle=CESM2HYCOM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedistStore(area2d, area1d, routehandle=HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedistStore(mask2D, mask1d, routehandle=HYCOM2CESM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !call ESMF_ArrayRedist(area1d, area2d, routehandle=CESM2HYCOM_RHR8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !dummy1D = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out
    !dummy2D = ESMF_ArrayCreate(distgrid2d, typekind=ESMF_TYPEKIND_R8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayGet(dummy2D, farrayPtr=fptr2D, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !fptr2D(:,:) = lpet

    !call ESMF_ArrayWrite(dummy2D, file='scramble.nc', variableName='dummy2D', &
    !  overwrite=.true., timeslice=1, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayRedist(dummy2D, dummy1D, routehandle=HYCOM2CESM_RHR8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayRedist(dummy1D, dummy2D, routehandle=CESM2HYCOM_RHR8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayWrite(dummy2D, file='scramble.nc', variableName='dummy2D', &
    !  overwrite=.true., timeslice=2, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !k = 1
    !do j = 1, ubound(fptr2D, 2)
    !  do i = 1, ubound(fptr2D, 1)
    !    fptr2D(i, j) = fptrSeqIndex(k)
    !    k = k + 1
    !  enddo
    !enddo

    !call ESMF_ArrayWrite(dummy2D, file='scramble.nc', variableName='dummy2D', &
    !  overwrite=.true., timeslice=3, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &

    !call ESMF_ArrayRedist(dummy2D, dummy1D, routehandle=HYCOM2CESM_RHR8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayRedist(dummy1D, dummy2D, routehandle=CESM2HYCOM_RHR8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayWrite(dummy2D, file='scramble.nc', variableName='dummy2D', &
    !  overwrite=.true., timeslice=4, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    deallocate(indexlist)
    deallocate(fptrSeqIndex)

    return

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_domain_esmf

  subroutine dumpRawData(state, filename, distArrayName, fieldlist, rc)

    use shr_string_mod

    type(ESMF_State), intent(in)    :: state
    character(len=*), intent(in)    :: filename
    character(len=*), intent(in)    :: distArrayName
    character(len=*), intent(in)    :: fieldlist
    integer, intent(out)            :: rc

    type(ESMF_Array)                :: d2x
    real(ESMF_KIND_R8), allocatable :: rawdata(:,:)
    integer                         :: gsize, nfields, elb(2,1), eub(2,1), lpet, rec_len
    type(ESMF_VM)                   :: vm

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=distArrayName, array=d2x, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_AttributeGet(state, name="gsize", value=gsize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_ArrayGet(d2x, exclusiveLBound=elb, exclusiveUBound=eub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    nfields = shr_string_listGetNum(trim(fieldlist))
    !nfields = eub(2,1)
!!Alex    print *, 'dumpRawData: nfields: ', nfields

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_VMGet(vm, localPet=lpet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if(lpet == 0) allocate(rawdata(nfields, gsize))

    call ESMF_ArrayGather(d2x, rawdata, 0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if(lpet == 0) then
      inquire (IOLENGTH=rec_len) rawdata
      open(1901,file=filename,status = 'unknown', form='unformatted', access='direct',recl=rec_len)
      write(1901,rec=1) rawdata
      close(1901)
    endif

    if(lpet == 0) deallocate(rawdata)

  end subroutine
#endif TARNOTNEEEDED
! Functions from cesm esmf library
  ! placeholder for the lookup function
      subroutine HYCOM_AdvertiseFields(state,nfields, field_defs, tag, rc)

      type(ESMF_State), intent(inout)               :: state
      integer,intent(in)                            :: nfields
      type(fld_list_type), intent(inout)          :: field_defs(:)
      character(len=*), intent(in)                  :: tag
      integer,          intent(out),     optional   :: rc

      ! local variables
      integer                                       :: i
      character(80)                                 :: shortName
      character(80)                                 :: stdName


      rc = ESMF_SUCCESS

      !nfields = size(fieldList) This should be input
      do i = 1, nfields
      if (.not. NUOPC_FieldDictionaryHasEntry(trim(field_defs(i)%stdname))) then
         call NUOPC_FieldDictionaryAddEntry( &
              standardName=trim(field_defs(i)%stdname), &
              canonicalUnits=trim(field_defs(i)%canonicalUnits), &
              rc=rc)
         if   (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
      endif

         call flush(6)
         call ESMF_LogWrite('Advertise: '//trim(field_defs(i)%stdname), ESMF_LOGMSG_INFO, rc=rc)
         if  (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, &
             file=__FILE__)) &
             return  ! bail out

         call NUOPC_Advertise(state, &
         standardName=field_defs(i)%stdname, &
         name=field_defs(i)%shortname, &
         rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo

      end subroutine

      subroutine HYCOM_FieldsSetup
      character(len=*),parameter  :: subname='(hycom_cap:HYCOM_FieldsSetup)'

!--------- import fields to Sea Ice -------------
      write(6,*) subname
! tcraig, don't point directly into cice data YET (last field is optional in interface)
! instead, create space for the field when it's "realized".
      call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_surface_temperature"         ,"k"  , "will provide")
      call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_surface_salinity"            ,"1"   , "will provide")
      call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_level"                       ,"m"  , "will provide")
      call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_surface_slope_zonal"         ,"1"   , "will provide")
      call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_surface_slope_merid"         ,"1"   , "will provide")
      call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_zonal"               ,"m/s", "will provide")
      call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_merid"               ,"m/s", "will provide")
      call fld_list_add(fldsFrOcn_num, fldsFrOcn, "freezing_melting_potential"      ,"1"   , "will provide")
      call fld_list_add(fldsFrOcn_num, fldsFrOcn, "mixed_layer_depth"               ,"m"  , "will provide")
! fields for export
      call fld_list_add(fldsToOcn_num, fldsToOcn, "sea_ice_fraction"                ,"1"   , "will provide") !
      call fld_list_add(fldsToOcn_num, fldsToOcn, "stress_on_ocn_ice_zonal"         ,"1"   , "will provide") !
      call fld_list_add(fldsToOcn_num, fldsToOCn, "stress_on_ocn_ice_merid"         ,"1"   , "will provide") !
      call fld_list_add(fldsToOcn_num, fldsToOcn, "sea_ice_temperature"             ,"1"   , "will provide") !
      call fld_list_add(fldsToOcn_num, fldsToOcn, "ice_mask"                        ,"1"   , "will provide")
      call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_sw_pen_to_ocn"              ,"1"   , "will provide") !
      call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_fresh_water_to_ocean_rate" ,"1"   , "will provide")!
      call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_salt_rate"                  ,"1"   , "will provide") !
      call fld_list_add(fldsToOcn_num, fldsToOcn, "net_heat_flx_to_ocn"             ,"1"   , "will provide") !
      call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_ice_volume"                 ,"1"   , "will provide")
      call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_snow_volume"                ,"1"   , "will provide")

  end subroutine HYCOM_FieldsSetup

  subroutine fld_list_add(num, fldlist, stdname, canonicalUnits,transferOffer, data, shortname)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    character(len=*),    intent(in)     :: canonicalUnits
    character(len=*),    intent(in)     :: transferOffer
    real(ESMF_KIND_R8), dimension(:,:,:), optional, target :: data
    character(len=*),    intent(in),optional :: shortname

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(hycom_cap:fld_list_add)'
    ! fill in the new entry

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num gt fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=rc)
      return
    endif

    fldlist(num)%stdname        = trim(stdname)
    fldlist(num)%canonicalUnits = trim(canonicalUnits)
    if (present(shortname)) then
       fldlist(num)%shortname   = trim(shortname)
    else
       fldlist(num)%shortname   = trim(stdname)
    endif
    fldlist(num)%transferOffer  = trim(transferOffer)
    if (present(data)) then
      fldlist(num)%assoc        = .true.
      fldlist(num)%farrayPtr    => data
    else
      fldlist(num)%assoc        = .false.
    endif

  end subroutine fld_list_add

  !-----------------------------------------------------------------------------
end module hycom_cap
