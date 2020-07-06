module hycom_cap

  !-----------------------------------------------------------------------------
  ! OCN Component for CESM-BETA; use ESMF and  NUOPC 
  !-----------------------------------------------------------------------------

  use MOD_HYCOM, only : HYCOM_Init, HYCOM_Run, HYCOM_Final, &
    end_of_run, end_of_run_cpl

  use mod_hycom_nuopc_glue
  use mod_dimensions, only : idm,jdm,nbdy
  use mod_cb_arrays_nuopc_glue
  use mod_nuopc_options, only: esmf_write_diagnostics, nuopc_restart, profile_memory, nuopc_tinterval
  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance, &
    model_label_Finalize  => label_Finalize

  implicit none
  
  private
 
  public SetServices
 
  ! private internal state to keep instance data
  type InternalStateStruct
    type(hycom_nuopc_glue_type)   :: glue
  end type

  type InternalState
    type(InternalStateStruct), pointer :: wrap
  end type

  integer  :: import_slice = 0
  integer  :: export_slice = 0  

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: canonicalUnits
    character(len=64) :: transferOffer
    logical           :: assoc    ! is the farrayPtr associated with internal data
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter    :: fldsMax = 50
  integer              :: fldsToOcn_num = 0
  type (fld_list_type) :: fldsToOcn(fldsMax)
  integer              :: fldsFrOcn_num = 0
  type (fld_list_type) :: fldsFrOcn(fldsMax)
  character(len=2048)  :: info
  real(ESMF_KIND_R8)   :: ocn_cpl_frq

  real, save, allocatable, dimension(:,:) :: si_sice

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(ocn:SetServices)'
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

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetClock, &
      specRoutine=SetClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
      specRoutine=hycom_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    ! Local Variables
    type(ESMF_VM)               :: vm
    integer                     :: mpi_comm, me, npes
    TYPE(ESMF_Time)             :: startTime, stopTime, hycomRefTime, currTime
    TYPE(ESMF_TimeInterval)     :: interval,timeStep
    real(ESMF_KIND_R8)          :: startTime_r8, stopTime_r8, l_startTime_r8
!    character(len=32)           :: starttype            ! infodata start type
!    logical                     :: restFlag = .false.        ! initial/restart run (F/T)

    character(len=*),parameter  :: subname='(HYCOM_cap:InitializeAdvertise)'
    rc = ESMF_SUCCESS
    call ESMF_GridCompGet(gcomp, vm=vm, localPet=me,petCount=npes, rc=rc)
!    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, mpiCommunicator=mpi_comm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call HYCOM_FieldsSetup()
    ! Translate currTime and stopTime into HYCOM format
    call ESMF_ClockGet(clock, currTime=currTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
! Day 1 sets start time which is fine. Does it need day2???? TO BE CHECKED
    ! Define HYCOM Ref Time
    call ESMF_TimeSet(hycomRefTime, yy=1900, mm=12, dd=31, calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
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

!    call ESMF_VMGet(vm,localPet=me,petCount=npes)
    if (me==0) then
      print *,"DMI_CPL: HYCOM_INIT -->> startTime_r8=", startTime_r8
      print *,"DMI_CPL:                  stopTime_r8=", stopTime_r8
    endif

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

    ! HYCOM restart or coldstart
    if (nuopc_restart) then
      l_startTime_r8=startTime_r8   ! Restart
    else
      l_startTime_r8=-startTime_r8  ! Coldstart
    endif

    ! Get the pointer restart name !!Alex
      !restart_write = .true. TAR NEEDED for averaged output.

    call ESMF_LOGWRITE("BEFORE HYCOM_INIT", ESMF_LOGMSG_INFO, rc=rc)
    ! -->> call into HYCOM <<--
    call HYCOM_Init(mpi_comm, & 
       hycom_start_dtg=l_startTime_r8, hycom_end_dtg=stopTime_r8)
!tar        restart_write=restart_write)

    call ESMF_LOGWRITE("AFTER HYCOM_INIT", ESMF_LOGMSG_INFO, rc=rc)

    !-- Check HYCOM timesteps + HYCOM steps between ice steps
    if (me==0) then
      write(6,*)'DMI_CPL: HYCOM timestep: baclin = ',baclin
      write(6,*)'DMI_CPL: HYCOM timesteps between ice update: icpfrq = ',icpfrq
    endif
    if ( nuopc_tinterval /= icpfrq*nint(baclin) ) then
      if (me==0) then
        write(6,*)'DMI_CPL: NUOPC timestep: nuopc_tinterval = ',nuopc_tinterval
        write(6,*)'DMI_CPL: ERROR: nuopc_tinterval /= icpfrq*baclin'
      endif
      call ESMF_LogWrite('DMI_CPL: ERROR: HYCOM ice timesteps do not match NUOPC timestep', &
        ESMF_LOGMSG_ERROR, rc=rc)
      rc = ESMF_RC_OBJ_BAD
      return
    endif

    if (me==0) print *,"DMI_CPL: HYCOM_INIT finished"

    ! Ice salinity
    allocate( si_sice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
    si_sice=0.0  ! Initial setting and allocate memory space

    ! set Component name so it becomes identifiable
#ifdef TARNOTNEEDED
    call ESMF_GridCompSet(gcomp, name="HYCOM", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif
    call HYCOM_AdvertiseFields( importState, fldsToOcn_num, fldsToOcn, tag='HYCOM import', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call HYCOM_AdvertiseFields(exportState, fldsFrOcn_num, fldsFrOcn, tag='HYCOM export', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    write(info,*) subname,' --- initialization phase 1 completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)

    if (me==0) print *,"DMI_CPL: HYCOM InitializeAdvertise finished"
      
  end subroutine InitializeAdvertise

  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! This method assumes that the incoming clock indicates the full simulation
    ! interval between startTime and stopTime.
    
    ! local variables
    type(ESMF_VM)               :: vm
    type(ESMF_Field)            :: field
    type(ESMF_Grid)             :: gridIn
    type(ESMF_Grid)             :: gridOut
    type(ESMF_array)            :: array
!    TYPE(ESMF_Time)             :: startTime, stopTime, hycomRefTime, currTime
!    TYPE(ESMF_TimeInterval)     :: interval,timeStep
!    real(ESMF_KIND_R8)          :: startTime_r8, stopTime_r8, l_startTime_r8
    type(InternalState)         :: is
    integer                     :: stat
    type(ESMF_CALKIND_FLAG)     :: calkind
!    logical                     :: restFlag = .false.        ! initial/restart run (F/T)
!    character(len=32)           :: starttype
    character(len=80)           :: pointer_filename          ! restart pointer file !!Alex
    logical                     :: restart_write = .false.   ! write restart
    character(len=*),parameter  :: subname='(HYCOM_cap:InitializeRealize)'   
    rc = ESMF_SUCCESS

    ! Fill in the glue structure.
    call HYCOM_GridInit(GridIn, rc=rc)
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
    !gridIn  = is%wrap%glue%grid ! for imported Fields
    gridOut = GridIn            ! for exported Fields
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

    ! Import data to HYCOM native structures through glue fields.
!   call HYCOM_Import(importState,.true.,rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    ! Export HYCOM native structures to data through glue fields.
    CALL HYCOM_export(exportState, .false., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine InitializeRealize

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

  subroutine ModelAdvance(gcomp, rc)
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
!    character(len=32)           :: starttype            ! infodata start type
!    logical                     :: restFlag = .false.
    character(len=128)          :: msg
    character*80 :: filenc !!Alex
    type(ESMF_VM)               :: vm
    integer                     :: me, npes
    character(len=*),parameter  :: subname='(HYCOM_cap:Model Advance)'

    rc = ESMF_SUCCESS
    call ESMF_GridCompGet(gcomp,vm=vm,localPet=me,petCount=npes, rc=rc)    
!    call ESMF_VMGet(vm,localPet=me,petCount=npes, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (profile_memory) call ESMF_VMLogMemInfo("Entering HYCOM Model_ADVANCE: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info,*) subname,' --- run phase 1 called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
  
    ! Get the internal state from Component.
#ifdef tarnotneeded
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompGetInternalState(gcomp, is, rc    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#else

! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    import_slice=import_slice+1
    export_slice=export_slice+1

#endif
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
    ! Import data to HYCOM native structures through glue fields.
    if (esmf_write_diagnostics >0) then
      if (mod(import_slice,esmf_write_diagnostics)==0) then
        call nuopc_write(state=importState,filenamePrefix='Import_HYCOM', &
          timeslice=import_slice/esmf_write_diagnostics,rc=rc)
      endif
    endif
    call HYCOM_Import(ImportState,.false.,rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !TODO: don't need the additional initialization step once data-dependency
    !TODO: is taken care of during initialize.
    initFlag = .false.
    !if (is%wrap%slice==1 .and. (.not. restFlag)) initFlag = .true.
    
    
    
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
!      CALL HYCOM_Run(endtime=real(endTime_r8,4),pointer_filename=pointer_filename, restart_write=restart_write) ! -->> call into HYCOM <<--
!      if (end_of_run .or. end_of_run_cpl) exit
!    enddo

    if (profile_memory) call ESMF_VMLogMemInfo("Before HYCOM_Run")
    CALL HYCOM_Run(endtime=endTime_r8,restart_write=restart_write)
    if (profile_memory) call ESMF_VMLogMemInfo("After HYCOM_Run")
    ! Export HYCOM native data through the glue fields.
    call HYCOM_export(exportstate,.false.,rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (esmf_write_diagnostics >0) then
      if (mod(export_slice,esmf_write_diagnostics)==0) then
        call nuopc_write(state=exportstate,filenamePrefix='Export_HYCOM',&
          timeslice=export_slice/esmf_write_diagnostics,rc=rc)
      endif
    endif
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (profile_memory) call ESMF_VMLogMemInfo("Leaving HYCOM Model_ADVANCE: ")

  end subroutine ModelAdvance

  !-----------------------------------------------------------------------------

!  subroutine hc_model_finalize(gcomp, rc)
!    type(ESMF_GridComp)  :: gcomp
!    type(ESMF_Clock)     :: clock
!    integer, intent(out) :: rc
!    character(len=*),parameter  :: subname='(HYCOM_cap:hc_model_finalize)'    
!    type(InternalState)  :: is
!    integer              :: stat

!    rc = ESMF_SUCCESS

!   write( nfo,*) subname,' --- finalize called --- '
!    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

!    call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

!    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

!    call HYCOM_Final

!    write(info,*) subname,' --- finalize completed --- '
!    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
!  end subroutine

! HYCOM uses clock as days from
  !-----------------------------------------------------------------------------
  subroutine SetClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_TimeInterval)       :: stabilityTimeStep, timestep
    character(len=*),parameter  :: subname='(hycom_cap:SetClock)'

    rc = ESMF_SUCCESS
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! tcraig: dt is the cice thermodynamic timestep in seconds
    call ESMF_TimeIntervalSet(timestep, s=nint(baclin), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockSet(clock, timestep=timestep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
    call ESMF_TimeIntervalSet(stabilityTimeStep, s=nint(baclin), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetClock(gcomp, clock, stabilityTimeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetClock

  !-----------------------------------------------------------------------------

  subroutine HYCOM_Import(st,initflag,rc)
    use mod_hycom_nuopc_glue, only: jja
    use mod_dimensions, only: ii
!TILL need to change
    type(ESMF_State)     :: st 
    logical              :: initFlag
    integer, intent(out) :: rc
    integer i,j
    real(8) :: xstress, ystress, pang_rev
    real(ESMF_KIND_R8), pointer :: dataPtr_sic(:,:),  dataPtr_sit(:,:),  dataPtr_sitx(:,:), &
                                   dataPtr_sity(:,:), dataPtr_siqs(:,:), dataPtr_sifs(:,:), &
                                   dataPtr_sih(:,:),  dataPtr_sifw(:,:), dataPtr_sifh(:,:), &
                                   dataPtr_sice(:,:)
    call state_getFldPtr(st,"sea_ice_fraction",dataPtr_sic,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call state_getFldPtr(st,"sea_ice_temperature",dataPtr_sit,rc=rc)  
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call state_getFldPtr(st,"mean_ice_volume",dataPtr_sih,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call state_getFldPtr(st,"stress_on_ocn_ice_zonal",dataPtr_sitx,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call state_getFldPtr(st,"stress_on_ocn_ice_merid",dataPtr_sity,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call state_getFldPtr(st,"mean_sw_pen_to_ocn",dataPtr_siqs,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call state_getFldPtr(st,"mean_salt_rate",dataPtr_sifs,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call state_getFldPtr(st,"mean_fresh_water_to_ocean_rate",dataPtr_sifw,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call state_getFldPtr(st,"net_heat_flx_to_ocn",dataPtr_sifh,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call state_getFldPtr(st,"ice_salinity",dataPtr_sice,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return

!NOT SURE ABOUT THESE FOUR. At least siu and siv should be used.
!        public cpl_sifh,    sifh_import   ! Ice Freezing/Melting Heat Flux
!        public cpl_sifw,    sifw_import   ! Ice Net Water Flux
!ALSO FLXICE
    if (iceflg.ge.2 .and. icmflg.ne.3) then
      do j=1,jja
        do i=1,ii
          covice(i,j) = dataPtr_sic(i,j) !Sea Ice Concentration
          si_c  (i,j) = dataPtr_sic(i,j) !Sea Ice Concentration
          if (covice(i,j).gt.0.0) then
            if (frzh(i,j)>0.0) then
              ! --- add energy to move tmxl towards tfrz (only if tmxl < tfrz)
              ! --- include some dependance on sea ice coverage
!              flxice(i,j) = max( dataPtr_sifh(i,j), frzh(i,j)*max(covice(i,j),0.1) )
              flxice(i,j) = frzh(i,j)*covice(i,j)
            else
              flxice(i,j) = dataPtr_sifh(i,j) ! Sea Ice Heat Flux Freezing potential
            endif
            xstress     = -dataPtr_sitx(i,j) ! opposite of what ice sees
            ystress     = -dataPtr_sity(i,j) ! oppostite of what ice sees
            pang_rev    = -pang(i,j)         ! Reverse angle
            si_tx (i,j) =  xstress*cos(pang_rev) - ystress*sin(pang_rev)
            si_ty (i,j) =  xstress*sin(pang_rev) + ystress*cos(pang_rev)
            fswice(i,j) =  dataPtr_siqs(i,j) !Solar Heat Flux thru Ice to Ocean already in swflx
            sflice(i,j) =  dataPtr_sifs(i,j)*1.e3 !Ice Freezing/Melting Salt Flux
            wflice(i,j) =  dataPtr_sifw(i,j) !Ice Water Flux
            temice(i,j) =  dataPtr_sit(i,j)  !Sea Ice Temperature
            si_t  (i,j) =  dataPtr_sit(i,j)  !Sea Ice Temperature
            thkice(i,j) =  dataPtr_sih(i,j)  !Sea Ice Thickness
            si_h  (i,j) =  dataPtr_sih(i,j)  !Sea Ice Thickness
            si_sice(i,j)=  dataPtr_sice(i,j) !Sea Ice Salinity
          else
            si_tx (i,j) =  0.0 !Sea Ice X-Stress into ocean
            si_ty (i,j) =  0.0 !Sea Ice Y-Stress into ocean
            fswice(i,j) =  0.0 !Solar Heat Flux thru Ice to Ocean already in swflx
            flxice(i,j) =  0.0 !freeze/melt potential
            sflice(i,j) =  0.0 !Ice Freezing/Melting Salt Flux
            wflice(i,j) =  0.0 !Ice Water Flux
            temice(i,j) =  0.0 !Sea Ice Temperature
            si_t  (i,j) =  0.0 !Sea Ice Temperature
            thkice(i,j) =  0.0 !Sea Ice Thickness
            si_h  (i,j) =  0.0 !Sea Ice Thickness
            si_sice(i,j)=  0.0 !Sea Ice Salinity
          endif
         
        enddo
      enddo
    elseif (iceflg.ge.2 .and. icmflg.eq.3) then
      do j=1,jja
        do i=1,ii
          si_c(i,j) = dataPtr_sic(i,j) !Sea Ice Concentration
          if (si_c(i,j).gt.0.0) then
            xstress    = -dataPtr_sitx(i,j) ! opposite of what ice sees
            ystress    = -dataPtr_sity(i,j) ! oppostite of what ice sees
            pang_rev   = -pang(i,j)         ! Reverse Angle
            si_tx(i,j) =  xstress*cos(pang_rev) - ystress*sin(pang_rev)
            si_ty(i,j) =  xstress*sin(pang_rev) + ystress*cos(pang_rev)
            si_h (i,j) =  dataPtr_sih(i,j) !Sea Ice Thickness
            si_t (i,j) =  dataPtr_sit(i,j) !Sea Ice Temperature
          else
            si_tx(i,j) = 0.0
            si_ty(i,j) = 0.0
            si_h (i,j) = 0.0
            si_t (i,j) = 0.0
            si_u (i,j) = 0.0
            si_v (i,j) = 0.0
          endif !covice
        enddo
      enddo
    endif
  end subroutine HYCOM_Import

  !-----------------------------------------------------------------------------

  subroutine HYCOM_export(st,initflag,rc)
    use mod_hycom_nuopc_glue, only: jja
    use mod_dimensions, only: ii

    type(ESMF_State)     :: st
    logical              :: initflag
    integer, intent(out) :: rc

    real(kind=ESMF_KIND_R8), pointer  :: dataPtr_sst(:,:)
    real(kind=ESMF_KIND_R8), pointer  :: dataPtr_sss(:,:)
    real(kind=ESMF_KIND_R8), pointer  :: dataPtr_sssz(:,:)
    real(kind=ESMF_KIND_R8), pointer  :: dataPtr_sssm(:,:)
    real(kind=ESMF_KIND_R8), pointer  :: dataPtr_ocncz(:,:) 
    real(kind=ESMF_KIND_R8), pointer  :: dataPtr_ocncm(:,:)
    real(kind=ESMF_KIND_R8), pointer  :: dataPtr_fmpot(:,:)
    real(kind=ESMF_KIND_R8), pointer  :: dataPtr_mld(:,:)

    integer :: i,j
    real(kind=ESMF_KIND_R8) :: hfrz, t2f, tfrz, smxl, tmxl, ssfi, tmlt
    real(kind=ESMF_KIND_R8) :: usur1, usur2, vsur1, vsur2, utot, vtot
    real(kind=ESMF_KIND_R8) :: ssh_n,ssh_s,ssh_e,ssh_w,dhdy
    integer                 :: cplfrq

! do sst and salinity at the same time
    call State_getFldPtr(st,'sea_surface_temperature',dataPtr_sst,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(st,'sea_surface_salinity',dataPtr_sss,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(st,'freezing_melting_potential',dataPtr_fmpot,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return

    cplfrq = nint( ocn_cpl_frq*(86400.0/baclin) )
    if (.not. initFlag) then
      do j=1,jja
        do i=1,ii
          tmxl = 0.5*(temp(i,j,1,2)+temp(i,j,1,1))
          smxl = 0.5*(saln(i,j,1,2)+saln(i,j,1,1))
          dataPtr_sst(i,j) = tmxl ! construct SST [C]
          dataPtr_sss(i,j) = smxl ! construct SSS
          hfrz = min( thkfrz*onem, dpbl(i,j) )
          t2f  = (spcifh*hfrz)/(baclin*dble(icefrq)*dble(icpfrq)*g)
          tfrz = tfrz_0 + smxl*tfrz_s          ! salinity dependent freezing point: HYCOM
          tmlt = tfrz_0 + si_sice(i,j)*tfrz_s  ! salinity dependent melting point:  CICE
          ! Modified melt point of sea ice. Old version only contained the else clause
          if ((tmlt>tfrz) .and. (tmxl>tfrz)) then
            ssfi = t2f*min(tmlt,tmlt-tmxl)
          else
            ssfi = (tfrz-tmxl)*t2f       !W/m^2 into ocean
          endif
          dataPtr_fmpot(i,j) = max(-1000.0,min(1000.0,ssfi)) ! > 0. freezing potential of flxice
          frzh(i,j)=dataPtr_fmpot(i,j)
        enddo
      enddo
    else
      frzh(:,:)          = 0.
      dataPtr_fmpot(:,:) = 0.
      dataPtr_sst(:,:)   = 0.
      dataPtr_sss(:,:)   = 0.
    endif

    call State_getFldPtr(st,'sea_surface_slope_zonal',dataPtr_sssz,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(st,'sea_surface_slope_merid',dataPtr_sssm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(st,'ocn_current_zonal',dataPtr_ocncz,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(st,'ocn_current_merid',dataPtr_ocncm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return


    if (.not. initFlag) then
      do j=1,jja
        do i=1,ii
! Now calculated in mod_hycom.F with accurate loop boundaries to get a bit for bit restart
! Results now differs from when the calculation was done here
            dataPtr_sssz(i,j)  = dhde(i,j)
            dataPtr_sssm(i,j)  = dhdn(i,j)
            dataPtr_ocncz(i,j) = uml(i,j)
            dataPtr_ocncm(i,j) = vml(i,j)
         enddo
       enddo
    else
            dataPtr_sssz(:,:)  = 0.
            dataPtr_sssm(i,j)  = 0.
            dataPtr_ocncz(:,:) = 0.
            dataPtr_ocncm(i,j) = 0.
    endif


    call State_getFldPtr(st,'mixed_layer_depth',dataPtr_mld,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    if (.not. initFlag) then
        do j=1,jja
        do i=1,ii
              dataPtr_mld(i,j) = dpbl(i,j) * qonem    ! convert mld to m
        enddo
        enddo
     else
        do j=1,jja
        do i=1,ii
              dataPtr_mld(i,j) = 0.
         enddo
        enddo
     endif
   end subroutine HYCOM_export
  !-----------------------------------------------------------------------------
  ! placeholder for the lookup function
   subroutine HYCOM_AdvertiseFields(state,nfields, field_defs, tag, rc)

   type(ESMF_State), intent(inout)   :: state
   integer,intent(in)                :: nfields
   type(fld_list_type), intent(inout):: field_defs(:)
   character(len=*), intent(in)      :: tag
   integer, intent(out), optional    :: rc

   ! local variables
   type(ESMF_VM) :: vm
   integer       :: i, me, npes
   character(80) :: shortName
   character(80) :: stdName

   call ESMF_VMGetGlobal(vm=vm, rc=rc)
   call ESMF_VMGet (vm, localPet=me, petCount=npes)
   rc = ESMF_SUCCESS

   !nfields = size(fieldList) This should be input
   if (me==0) write(6,*)'DMI_CPL: Number of HYCOM fields = ',nfields
   do i = 1, nfields
     if (.not. NUOPC_FieldDictionaryHasEntry(trim(field_defs(i)%stdname))) then
       if (me==0) write(6,*) &
         'DMI_CPL: ',trim(field_defs(i)%stdname),' : ',trim(field_defs(i)%canonicalUnits)
       call NUOPC_FieldDictionaryAddEntry( &
         standardName=trim(field_defs(i)%stdname), &
         canonicalUnits=trim(field_defs(i)%canonicalUnits), &
         rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
     endif

     call flush(6)
     call ESMF_LogWrite('Advertise: '//trim(field_defs(i)%stdname), ESMF_LOGMSG_INFO, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
    ! write(6,*) trim(field_defs(i)%stdname), trim(field_defs(i)%shortname) 
     call NUOPC_Advertise(state, &
       standardName=field_defs(i)%stdname, &
       name=field_defs(i)%shortname, &
       rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

   enddo

   end subroutine HYCOM_AdvertiseFields

   subroutine HYCOM_FieldsSetup
   character(len=*),parameter  :: subname='(hycom_cap:HYCOM_FieldsSetup)'

!--------- import fields to Sea Ice -------------
! tcraig, don't point directly into cice data YET (last field is optional in interface)
! instead, create space for the field when it's "realized".
   call fld_list_add(fldsFrOcn_num,fldsFrOcn,"sea_surface_temperature"       ,"C"  ,"will provide")
   call fld_list_add(fldsFrOcn_num,fldsFrOcn,"sea_surface_salinity"          ,"1"  ,"will provide")
   call fld_list_add(fldsFrOcn_num,fldsFrOcn,"sea_surface_slope_zonal"       ,"1"  ,"will provide")
   call fld_list_add(fldsFrOcn_num,fldsFrOcn,"sea_surface_slope_merid"       ,"1"  ,"will provide")
   call fld_list_add(fldsFrOcn_num,fldsFrOcn,"ocn_current_zonal"             ,"m/s","will provide")
   call fld_list_add(fldsFrOcn_num,fldsFrOcn,"ocn_current_merid"             ,"m/s","will provide")
   call fld_list_add(fldsFrOcn_num,fldsFrOcn,"freezing_melting_potential"    ,"1"  ,"will provide")
   call fld_list_add(fldsFrOcn_num,fldsFrOcn,"mixed_layer_depth"             ,"m"  ,"will provide")
! fields for import
   call fld_list_add(fldsToOcn_num,fldsToOcn,"sea_ice_fraction"              ,"1"  ,"will provide") !
   call fld_list_add(fldsToOcn_num,fldsToOcn,"stress_on_ocn_ice_zonal"       ,"1"  ,"will provide") !
   call fld_list_add(fldsToOcn_num,fldsToOCn,"stress_on_ocn_ice_merid"       ,"1"  ,"will provide") !
   call fld_list_add(fldsToOcn_num,fldsToOcn,"sea_ice_temperature"           ,"C"  ,"will provide") !
   call fld_list_add(fldsToOcn_num,fldsToOcn,"mean_sw_pen_to_ocn"            ,"1"  ,"will provide") !
   call fld_list_add(fldsToOcn_num,fldsToOcn,"mean_fresh_water_to_ocean_rate","1"  ,"will provide")!
   call fld_list_add(fldsToOcn_num,fldsToOcn,"mean_salt_rate"                ,"1"  ,"will provide") !
   call fld_list_add(fldsToOcn_num,fldsToOcn,"net_heat_flx_to_ocn"           ,"1"  ,"will provide") !
   call fld_list_add(fldsToOcn_num,fldsToOcn,"mean_ice_volume"               ,"1"  ,"will provide")
   call fld_list_add(fldsToOcn_num,fldsToOcn,"mean_snow_volume"              ,"1"  ,"will provide")
   call fld_list_add(fldsToOcn_num,fldsToOcn,"ice_salinity"                  ,"1"  ,"will provide")

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

  subroutine hycom_model_finalize(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)     :: clock
    type(ESMF_Time)                        :: currTime
    character(len=*),parameter  :: subname='(hycom_cap:hycom_model_finalize)'

    rc = ESMF_SUCCESS

    write(info,*) subname,' --- finalize called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)

    call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call HYCOM_Final()

    write(info,*) subname,' --- finalize completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)

  end subroutine hycom_model_finalize

  subroutine state_GetFldPtr(ST, fldname, fldptr, rc)
    type(ESMF_State), intent(in) :: ST
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:)
    integer, intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(hycom_cap:State_GetFldPtr)'
   
    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (present(rc)) rc = lrc

  end subroutine State_GetFldPtr


  !-----------------------------------------------------------------------------
end module hycom_cap
