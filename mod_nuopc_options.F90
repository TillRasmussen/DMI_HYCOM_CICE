module mod_nuopc_options
  implicit none
  logical,save,public        :: nuopc_restart        ! If true model will restart from 
  
  integer (4), save, public     :: nuopc_tinterval       ! coupling interval 
  integer(4), dimension(4), save, public :: tstart, tend !year, month, day and hour
  integer(4), save, public   :: ocn_petCount, ice_petCount, esmf_write_diagnostics
  !======================
  public nuopc_opt
  contains
  !======================

  !----------------------------------------------------------
  subroutine nuopc_opt()
  !----------------------------------------------------------
    implicit none
    integer (4), parameter :: funi=504
    character(len = 10) :: nuopc_tstart, nuopc_tend
    integer (4) :: nml_err
    
    namelist /nuopc/ nuopc_tstart, nuopc_tend,nuopc_tinterval, nuopc_restart, &
                     esmf_write_diagnostics
    !default values
    !  in operational runs at time=0 where the model has to be restarted
    nuopc_tstart          = '2017031500' ! start time string
    nuopc_tend            = '2017031503' ! end time string
    nuopc_tinterval       = 180 ! Time step in seconds (cpl interval
    nuopc_restart         = .true.  ! Is this a cold start or a restart (true is restart)
    ocn_petCount          = 9
    ice_petCount          = 6
    esmf_write_diagnostics = 0
    ! read namelist
    open (funi, file='nuopc_opt', status='old',iostat=nml_err)
    if (nml_err .ne. 0) then
        write (6,'(a)')               &
        'DMI options are not read.'
    endif
    do while (nml_err == 0)
      read(funi, nml=nuopc,iostat=nml_err)
    end do
    close(funi)
    write(6,*) nuopc_tinterval
    read (nuopc_tstart(1:4),*) tstart(1)
    read (nuopc_tstart(5:6),*) tstart(2)
    read (nuopc_tstart(7:8),*) tstart(3)
    read (nuopc_tstart(9:10),*) tstart(4)
    read (nuopc_tend(1:4),*)   tend(1)
    read (nuopc_tend(5:6),*)   tend(2)
    read (nuopc_tend(7:8),*)   tend(3)
    read (nuopc_tend(9:10),*)   tend(4)
  end subroutine nuopc_opt  
end module mod_nuopc_options
