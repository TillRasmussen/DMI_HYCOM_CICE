      module mod_cb_arrays_nuopc_glue
      
        ! The purpose of this module is to make select common block entities
        ! available as public module variables. Doing this through an extra
        ! Fortran module is cleaner than simply including the entyre common
        ! blocks header in the actual glue module. Further, the common blocks
        ! header is in fixed format, and that means that which ever unit pulls
        ! it in via include also needs to be in fixed format... not something
        ! that would be nice for the actual hycom_nuopc_glue module.
      
        use mod_xc
        use mod_cb_arrays       ! HYCOM saved arrays

        implicit none
    
        private
  
        ! grid
        public plon     ! lon at the p pts
        public plat     ! lat at the p pts
        public qlon     ! lon at the q pts
        public qlat     ! lat at the q pts
        public scp2     ! grid cell area at the p pts
        public scux     ! grid cell area at the p pts
        public scvy     ! grid cell area at the p pts

        ! native exports
        public temp           ! temp at the p pts
        public u, v           ! velocity components
        public ubavg, vbavg   ! barotropic velocities
        public srfhgt         ! sea surface height, g*ssh(m)
        public dpbl           ! turbulent boundary layer depth
        public saln           ! salinity
        public pang           ! ANGLET
        public dhde           ! eastward ssh slope
        public dhdn           ! northward ssh slope
        public uml           ! eastward u velocity
        public vml           ! northward v velocity
        public tml           ! export T average over one coupling sequence
        public sml           ! export S average over one coupling sequence
        public sshm           ! export SSH average over one coupling sequence
        !public frzh           ! export of freezing and melting potential 
        
        ! native imports
        public covice         ! ice coverage (rel.units)
        public thkice         ! grid-cell avg. ice thknss (m)
        public temice         ! ice surface temperature
        public fswice         ! swv  flux under ice
        public flxice         ! heat flux under ice
        public sflice         ! salt flux under ice
        public wflice         ! water flux under ice
        public si_c           ! ice concentration   on p-grid from coupler
        public si_h           ! ice thickness       on p-grid from coupler
        public si_t           ! ice temperature     on p-grid from coupler
        public si_u           ! ice u-velocity      on p-grid from ccoupler
        public si_v           ! ice v-velocity      on p-grid from coupler
        public si_tx          ! x-stress  under ice on p-grid from coupler
        public si_ty          ! y-stress  under ice on p-grid from coupler
        
        ! ocn coupling frequency
        public icefrq         ! in time step in HYCOM
        public icpfrq         ! number of time steps between sea-ice updates
        
        ! scalars and flags
        public g, onem, qonem, thkfrz, baclin, spcifh, pcp_fact
        public tfrz_0, tfrz_s, frzh
        public iceflg, icmflg, sstflg,itest,jtest,ltripolar
        
      end module
