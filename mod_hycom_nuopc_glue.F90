module mod_hycom_nuopc_glue

  !-----------------------------------------------------------------------------
  ! NUOPC glue code for HYCOM
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC

  use mod_xc
  use mod_cb_arrays_nuopc_glue
  use mod_cb_arrays
  use mod_dimensions, only: ii
  implicit none

  integer :: jja, jtdma 
  
  private
  
  type hycom_nuopc_glue_type
    ! native grid:
    type(ESMF_Grid)         :: grid
    ! native mesh:
!    type(ESMF_Mesh)         :: mesh           ! from Grid center stagger 
    ! import fields:
    type(ESMF_FieldBundle)  :: importFields
    ! export fields:
    type(ESMF_FieldBundle)  :: exportFields
  endtype
  
  public hycom_nuopc_glue_type
  public HYCOM_TileInfo, HYCOM_GridInit !HYCOM_GlueInitialize
!  public HYCOM_GlueFieldRealize, HYCOM_GlueFieldsRealize
!tarnotneeded  public HYCOM_GlueFieldsDataImport, HYCOM_GlueFieldsDataExport
  public jja  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine HYCOM_TileInfo(rc)
    integer, intent(out), optional :: rc
    
    character(len=1600)  :: msg
    
    if (present(rc)) rc = ESMF_SUCCESS
    
    call ESMF_LogWrite("<HYCOM Tile Info>", ESMF_LOGMSG_INFO)
    
    write (msg, *) "ipr=", ipr, "   jpr=", jpr, "   ijpr=", ijpr
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
    write (msg, *) "mproc=", mproc, "   nproc=", nproc, "   mnproc=", mnproc
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
    write (msg, *) "i0=", i0, "   ii=", ii, "   j0=", j0, "   jj=", jja
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
!!Alex    write (msg, *) "margin=", margin, "   nreg=", nreg, "   vland=", vland
    write (msg, *) "nbdy=", nbdy, "   nreg=", nreg, "   vland=", vland
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)

    call ESMF_LogWrite("</HYCOM Tile Info>", ESMF_LOGMSG_INFO)
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine HYCOM_GridInit(grid, rc)
!    type(hycom_nuopc_glue_type), intent(inout)  :: glue
    integer, intent(out), optional              :: rc
    integer, allocatable, target  :: deBlockList(:,:,:)
    integer, allocatable, target  :: deBlockList_buf(:,:,:)
    integer, pointer      :: sendData(:), recvData(:)
    type(ESMF_VM)         :: vm
    type(ESMF_Grid), intent(inout)  :: grid
    integer               :: localPet, petCount
    type(ESMF_DistGridConnection), allocatable :: connectionList(:)
    type(ESMF_DistGrid)   :: dg
    type(ESMF_Array)      :: array_plon, array_plat, array_msk, array_area
    
    integer                 :: i, j, dumpUnit
    
    real(kind=ESMF_KIND_R8), pointer  :: farrayPtr(:,:)
    integer, pointer                  :: farrayPtrI(:,:)

    character(len=80)       :: dumpFile
    
    if (present(rc)) rc = ESMF_SUCCESS
    
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (localPet /= mnproc-1) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg="localPet and mnproc must correspond!", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return
    endif

    ! change the domain if HYCOM Arctic patch   
    if (ltripolar) then
       jtdma = jtdm-1
       jja   = min(jj,jtdma-j0)
    else
       jtdma = jtdm
       jja   = jj
    endif
   
    ! prepare the deBlockList needed for DistGrid creation
    
    ! first step: set the local piece of information
    
    allocate(deBlockList(2, 2, ijpr)) ! dimCount, 2, deCount
    allocate(deBlockList_buf(2, 2, ijpr)) ! dimCount, 2, deCount
    
    deBlockList(1, 1, mnproc) = i0+1   ! minIndex 1st dim
    deBlockList(2, 1, mnproc) = j0+1   ! minIndex 2nd dim
    
    deBlockList(1, 2, mnproc) = i0+ii  ! maxIndex 1st dim
    deBlockList(2, 2, mnproc) = j0+jja ! maxIndex 2nd dim

    deBlockList_buf = deBlockList

    ! second step: all-to-all the information so every PET has full deBlockList

    sendData => deBlockList(:,1,mnproc)
    recvData => deBlockList_buf(:,1,1)
    
    call ESMF_VMAllGather(vm, sendData, recvData, 4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    deBlockList = deBlockList_buf
    deallocate(deBlockList_buf)
    ! prepare a single connection for periodicity along i-axis (longitude)
#ifdef tarnotconnectedlon
    allocate(connectionList(1))
    call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
      tileIndexB=1, positionVector=(/itdm, 0/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! ready to create the HYCOM DistGrid from deBlockList with periodic connect.
    dg = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/itdm,jtdma/), &
      deBlockList=deBlockList, connectionList=connectionList, &
      indexflag=ESMF_INDEX_GLOBAL, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    deallocate(connectionList)
#else      
    dg = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/itdm,jtdma/), &
      deBlockList=deBlockList, &
      indexflag=ESMF_INDEX_GLOBAL, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif
    deallocate(deBlockList)
    
    ! dress up "plon" array, considering HYCOM memory layout with halo + padding
    array_plon = ESMF_ArrayCreate(dg, farray=plon, &
      indexflag=ESMF_INDEX_DELOCAL, &
      computationalLWidth=(/nbdy,nbdy/), computationalUWidth=(/nbdy,nbdy/), &
      totalLWidth=(/nbdy,nbdy/), & ! lower corner pinned to memory alloc, float upper
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! dress up "plat" array, considering HYCOM memory layout with halo + padding
    array_plat = ESMF_ArrayCreate(dg, farray=plat, &
      indexflag=ESMF_INDEX_DELOCAL, &
      computationalLWidth=(/nbdy,nbdy/), computationalUWidth=(/nbdy,nbdy/), &
      totalLWidth=(/nbdy,nbdy/), & ! lower corner pinned to memory alloc, float upper
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! dress up "mask" array, considering HYCOM memory layout with halo + padding
    array_msk = ESMF_ArrayCreate(dg, farray=ip, &
      indexflag=ESMF_INDEX_DELOCAL, &
      computationalLWidth=(/nbdy,nbdy/), computationalUWidth=(/nbdy,nbdy/), &
      totalLWidth=(/nbdy,nbdy/), & ! lower corner pinned to memory alloc, float upper
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! dress up "area" array, considering HYCOM memory layout with halo + padding
    array_area = ESMF_ArrayCreate(dg, farray=scp2, &
      indexflag=ESMF_INDEX_DELOCAL, &
      computationalLWidth=(/nbdy,nbdy/), computationalUWidth=(/nbdy,nbdy/), &
      totalLWidth=(/nbdy,nbdy/), & ! lower corner pinned to memory alloc, float upper
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! ready to create the HYCOM Grid from DistGrid and coordinate Arrays
    grid = ESMF_GridCreate(dg, coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger longitude coordinate Array
    call ESMF_GridSetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      coordDim=1, array=array_plon, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger latitude coordinate Array
    call ESMF_GridSetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      coordDim=2, array=array_plat, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the corner stagger latitude coordinate Array
    call ESMF_GridSetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
      coordDim=1, array=array_plon, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger latitude coordinate Array
    call ESMF_GridSetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
      coordDim=2, array=array_plat, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger mask Array
    call ESMF_GridSetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_MASK, array=array_msk, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger area Array
    call ESMF_GridSetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_AREA, array=array_area, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! initialize coupling flags
    cpl_taux      =.false.
    cpl_tauy      =.false.
    cpl_wndspd    =.false.
    cpl_ustara    =.false.
    cpl_airtmp    =.false.
    cpl_vapmix    =.false.
    cpl_swflx     =.false.
    cpl_lwmdnflx  =.false.
    cpl_lwmupflx  =.false.
    cpl_latflx    =.false.
    cpl_sensflx   =.false.
    cpl_precip    =.false.
    cpl_surtmp    =.false.
    cpl_seatmp    =.false.
    cpl_orivers   =.false.
    cpl_irivers   =.false.
   
  end subroutine

  !-----------------------------------------------------------------------------

#ifdef tarnotneeded
!IMPORT
 ! coupled fields from CESM SAVED TO REMEMBER  
      ! identify the exact import field by standard name and copy the data
      twoLevel = .false. ! reset
      if (fieldStdName == "sea_ice_area_fraction") then
        cpl_sic = .true.
        impPtr => sic_import
        !farrayPtr => sic_import
      elseif (fieldStdName == "downward_x_stress_at_sea_ice_base") then
        cpl_sitx = .true.
        impPtr => sitx_import
      elseif (fieldStdName == "downward_y_stress_at_sea_ice_base") then
        cpl_sity = .true.
        impPtr => sity_import
      elseif (fieldStdName == "downward_sea_ice_basal_solar_heat_flux") then
        cpl_siqs = .true.
        impPtr => siqs_import
      elseif (fieldStdName == "upward_sea_ice_basal_heat_flux") then
        cpl_sifh = .true.
        impPtr => sifh_import
      elseif (fieldStdName == "downward_sea_ice_basal_salt_flux") then
        cpl_sifs = .true.
        impPtr => sifs_import
      elseif (fieldStdName == "downward_sea_ice_basal_water_flux") then
        cpl_sifw = .true.
        impPtr => sifw_import
      elseif (fieldStdName == "sea_ice_temperature") then
        cpl_sit = .true.
        impPtr => sit_import
      elseif (fieldStdName == "sea_ice_thickness") then
        cpl_sih = .true.
        impPtr => sih_import
      elseif (fieldStdName == "sea_ice_x_velocity") then
        cpl_siu = .true.
        impPtr => siu_import
      elseif (fieldStdName == "sea_ice_y_velocity") then
        cpl_siv = .true.
        impPtr => siv_import
!To Alex, connect the water flux from river here then 
! convert Kg/m^2/s -> m/s as needed below for these two fields
      endif
      ! copy the data into the right import location
        do j=1,jja
        do i=1,ii
            impPtr(i,j) = farrayPtr(i,j)
        enddo
        enddo
      
    enddo

    ! transfer SEA-ICE imports into native HYCOM variables
    do j=1,jja
    do i=1,ii
      if (iceflg.ge.2 .and. icmflg.ne.3) then
        covice(i,j) = sic_import(i,j) !Sea Ice Concentration
        si_c  (i,j) = sic_import(i,j) !Sea Ice Concentration
        if (covice(i,j).gt.0.0) then
          if (frzh(i,j).gt.0.0) then
             flxice(i,j) = frzh(i,j)        !Sea Ice Heat Flux Freezing potential
          else
             flxice(i,j) = sifh_import(i,j) !Sea Ice Heat Flux Melting potential
          endif
          si_tx (i,j) =  sitx_import(i,j) !Sea Ice X-Stress into ocean
          si_ty (i,j) =  sity_import(i,j) !Sea Ice Y-Stress into ocean
          fswice(i,j) =  siqs_import(i,j) !Solar Heat Flux thru Ice to Ocean already in swflx
          sflice(i,j) =  sifs_import(i,j)*1.e3 !Ice Freezing/Melting Salt Flux
          wflice(i,j) =  sifw_import(i,j) !Ice Water Flux
          temice(i,j) =   sit_import(i,j) !Sea Ice Temperature
          si_t  (i,j) =   sit_import(i,j) !Sea Ice Temperature
          thkice(i,j) =   sih_import(i,j) !Sea Ice Thickness
          si_h  (i,j) =   sih_import(i,j) !Sea Ice Thickness
          si_u  (i,j) =   siu_import(i,j) !Sea Ice X-Velocity
          si_v  (i,j) =   siv_import(i,j) !Sea Ice Y-Velocity
        else
          si_tx (i,j) = 0.0
          si_ty (i,j) = 0.0
          fswice(i,j) = 0.0
          flxice(i,j) = 0.0
          sflice(i,j) = 0.0
          wflice(i,j) = 0.0
          temice(i,j) = 0.0
          si_t  (i,j) = 0.0
          thkice(i,j) = 0.0
          si_h  (i,j) = 0.0
          si_u  (i,j) = 0.0
          si_v  (i,j) = 0.0
        endif !covice
      elseif (iceflg.ge.2 .and. icmflg.eq.3) then
        si_c(i,j) =  sic_import(i,j) !Sea Ice Concentration
        if (si_c(i,j).gt.0.0) then
          si_tx(i,j) = -sitx_import(i,j) !Sea Ice X-Stress into ocean
          si_ty(i,j) = -sity_import(i,j) !Sea Ice Y-Stress into ocean
          si_h (i,j) =   sih_import(i,j) !Sea Ice Thickness
          si_t (i,j) =   sit_import(i,j) !Sea Ice Temperature
          si_u (i,j) =   siu_import(i,j) !Sea Ice X-Velocity
          si_v (i,j) =   siv_import(i,j) !Sea Ice Y-Velocity
        else
          si_tx(i,j) = 0.0
          si_ty(i,j) = 0.0
          si_h (i,j) = 0.0
          si_t (i,j) = 0.0
          si_u (i,j) = 0.0
          si_v (i,j) = 0.0
        endif !covice
      endif !iceflg>=2 (icmflg)
    enddo
    enddo
!END IMPORT
!START EXPORT
  !-----------------------------------------------------------------------------
      ! identify the exact export field by standard name and copy the data
      if (fieldStdName == "sea_surface_temperature") then
        do j=1,jja
        do i=1,ii
          farrayPtr(i,j) = 0.5*(temp(i,j,1,2)+temp(i,j,1,1))  ! construct SST [C]
         farrayPtr(i,j) = farrayPtr(i,j) + 273.15d0            ! [C] -> [K]
        enddo
        enddo
      elseif (fieldStdName == "upward_sea_ice_basal_available_heat_flux") then
      ! get coupling frequency in time steps
        cplfrq = nint( ocn_cpl_frq*(86400.0/baclin) )
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then
! ---     quantities for available freeze/melt heat flux
! ---     relax to tfrz with e-folding time of cplfrq time steps
! ---     assuming the effective surface layer thickness is hfrz
! ---     multiply by dpbl(i,j)/hfrz to get the actual e-folding time
              hfrz = min( thkfrz*onem, dpbl(i,j) )
              t2f  = (spcifh*hfrz)/(baclin*cplfrq*g)
              ! ---     average both available time steps, to avoid time splitting.
              !smxl = 0.5*(saln(i,j,1,2)+saln(i,j,1,1)) !!Alex calculated in mod_hycom.F
              !tmxl = 0.5*(temp(i,j,1,2)+temp(i,j,1,1)) !!Alex calculated in mod_hycom.F 
              smxl = sml(i,j)
              tmxl = tml(i,j)
              tfrz = tfrz_0 + smxl*tfrz_s  !salinity dependent freezing point
              ssfi = (tfrz-tmxl)*t2f       !W/m^2 into ocean
        
              farrayPtr(i,j) = max(-1000.0,min(1000.0,ssfi))
              frzh     (i,j) = max(-1000.0,min(1000.0,ssfi)) ! > 0. freezing potential of flxice
           else
              farrayPtr(i,j) = 0.
              frzh     (i,j) = 0.                            ! > 0. freezing potential of flxice
           endif
        enddo
        enddo        
      elseif (fieldStdName == "sea_lev") then
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then 
!!Alex              farrayPtr(i,j) = srfhgt(i,j)/g   ! convert ssh to m
              farrayPtr(i,j) = sshm(i,j)/g   ! convert ssh to m
           else
              farrayPtr(i,j) = 0.
           endif
        enddo
        enddo
      elseif (fieldStdName == "mixed_layer_depth") then
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then 
              farrayPtr(i,j) = dpbl(i,j) * qonem    ! convert mld to m
           else
              farrayPtr(i,j) = 0.
           endif
         enddo
        enddo
      elseif (fieldStdName == "s_surf") then
        do j=1,jja
        do i=1,ii
          farrayPtr(i,j) = 0.5d0 * (saln(i,j,1,2)+saln(i,j,1,1))
!!Alex          farrayPtr(i,j) = sml(i,j)

        enddo
        enddo
      elseif (fieldStdName == "ocn_current_zonal") then
        do j=1,jja
           do i=1,ii
           if (.not. initFlag) then 
           ! Now calculated in mod_hycom.F with accurate loop boundaries to get a bit for bit restart
           ! Results now differs from when the calculation was done here
              farrayPtr(i,j) = uml(i,j)
           else
              farrayPtr(i,j) = 0.
           endif
        enddo
        enddo
      elseif (fieldStdName == "ocn_current_merid") then
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then 
           ! Now calculated in mod_hycom.F with accurate loop boundaries to get a bit for bit restart
           ! Results now differs from when the calculation was done here
              farrayPtr(i,j) = vml(i,j)
           else
              farrayPtr(i,j) = 0.
           endif
        enddo
        enddo
      elseif (fieldStdName == "eastward_sea_surface_slope") then
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then
           ! Now calculated in mod_hycom.F with accurate loop boundaries to get a bit for bit restart
           ! Results now differs from when the calculation was done here
              farrayPtr(i,j) = dhde(i,j)
           else
              farrayPtr(i,j) = 0.
           endif
         enddo
        enddo
      elseif (fieldStdName == "northward_sea_surface_slope") then
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then 
           ! Now calculated in mod_hycom.F with accurate loop boundaries to get a bit for bit restart
           ! Results now differs from when the calculation was done here
              farrayPtr(i,j) = dhdn(i,j)
           else
              farrayPtr(i,j) = 0.
           endif
         enddo
        enddo
      endif
      
    enddo
#endif

end module
