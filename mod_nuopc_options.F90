module mod_nuopc_options
  implicit none
  logical,save,public        :: nuopc_restart  ! If true HYCOM will do a restart. Else Coldstart 
  logical,save,public        :: profile_memory ! If true Profiling around {cice,hycom}_run
  logical,save,public        :: reduce_cpl
  
  integer(4), save, public   :: nuopc_tinterval       ! coupling interval 
  integer(4), dimension(4), save, public :: tstart, tend !year, month, day and hour
  integer(4), save, public   :: ocn_petCount, ice_petCount, esmf_write_diagnostics
  !======================
  public nuopc_opt
  contains
  !======================

  !----------------------------------------------------------
  subroutine nuopc_opt()
  !----------------------------------------------------------
    use ESMF
    implicit none
    integer (4), parameter :: funi=504
    character(len = 10) :: nuopc_tstart, nuopc_tend
    integer (4) :: nml_err
    type(ESMF_VM) :: vm
    integer       :: me, npes, rc
    
    namelist /nuopc_nml/ nuopc_tstart, nuopc_tend,nuopc_tinterval, &
                         nuopc_restart, esmf_write_diagnostics,    &
                         ocn_petCount, ice_petCount, profile_memory, &
                         reduce_cpl
    !default values
    !  in operational runs at time=0 where the model has to be restarted
    nuopc_tstart          = '2017031500' ! start time string
    nuopc_tend            = '2017031503' ! end time string
    nuopc_tinterval       = 180 ! *INTEGER* Time step in seconds (cpl interval)
    nuopc_restart         = .true.  ! Is this a cold start or a restart (true is restart)
    ocn_petCount          = 9
    ice_petCount          = 6
    esmf_write_diagnostics = 0  !number of time steps between netcdf dump. 0 No dump.
                                ! *NB*: Set to zero if more than one node (crash).
    reduce_cpl           = .false. ! If true dont couple fields ice thickness and ice termperature as these are not needed
    profile_memory       = .false.
    ! read namelist
    open (funi, file='nuopc_opt', status='old',iostat=nml_err)
    if (nml_err < 0) then
      write(6,'(a)')'ERROR: nuopc_opt namelist not found or not readable'
      call exit(1)
    endif
    nml_err=9
    do while (nml_err > 0)
      read(funi, nml=nuopc_nml,iostat=nml_err)
    end do
    close(funi)
    !-- Print info
    call ESMF_VMGetGlobal(vm=vm, rc=rc)
    call ESMF_VMGet (vm, localPet=me, petCount=npes)
    if (me==0) then
      write(6,*)'DMI_CPL: Reading nuopc_opt'
      write(6,*)'DMI_CPL: nuopc_tstart:   ',nuopc_tstart
      write(6,*)'DMI_CPL: nuopc_tend:     ',nuopc_tend
      write(6,*)'DMI_CPL: nuopc_tinterval:',nuopc_tinterval
      write(6,*)'DMI_CPL: nuopc_restart:  ',nuopc_restart
      write(6,*)'DMI_CPL: ocn_petCount:   ',ocn_petCount
      write(6,*)'DMI_CPL: ice_petCount:   ',ice_petCount
      write(6,*)'DMI_CPL: esmf_write_diagnostics:',esmf_write_diagnostics
      write(6,*)'DMI_CPL: profile_memory: ',profile_memory
      write(6,*)'DMI_CPL: reduce_cpl NOT IMPLEMENTED. SET TO TRUE HAS NO INFLUENCE: ', reduce_cpl
    endif
    !-- Adjust start/end times
    read (nuopc_tstart(1:4),*) tstart(1)
    read (nuopc_tstart(5:6),*) tstart(2)
    read (nuopc_tstart(7:8),*) tstart(3)
    read (nuopc_tstart(9:10),*)tstart(4)
    read (nuopc_tend(1:4),*)   tend(1)
    read (nuopc_tend(5:6),*)   tend(2)
    read (nuopc_tend(7:8),*)   tend(3)
    read (nuopc_tend(9:10),*)  tend(4)
  end subroutine nuopc_opt  
end module mod_nuopc_options
