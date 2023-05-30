! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
! Provides non-LTE radiative heating rate for the mesosphere to lower thermosphere (MLT)
!   in the infrared due to CO2 following the approach of Fomichev (1998).
!   Code adapted from the NCAR whole atmosphere code (WACCM) module nlte_fomichev.F90 
!
! -------------------------------------------------------------------------------------------------
module mo_upper_atmosphere_heating_lw
  use mo_rte_kind,         only: wp, wl
  use mo_rte_config,       only: check_values, check_extents
  use mo_rte_util_array,   only: any_vals_less_than, any_vals_outside, extents_are
  use mo_rrtmgp_constants, only: k_boltz, avogad, R_univ_gconst, grav

  implicit none

  private
  save

! Public subroutines
  public :: nlte_fomichev_init, nlte_fomichev_calc

  type, public :: ty_upper_atmosphere_heating_lw
    private

    integer :: nprs67, nprs51, nprs43, nprs35
    integer :: nlev9, nlev6, nco2amt

    ! MLT LW heating LUT coefficients
    real(wp), dimension(:),   allocatable :: xr, ig, co2o
    real(wp), dimension(:,:), allocatable :: ao3
    real(wp), dimension(:,:), allocatable :: a150, b150
    real(wp), dimension(:,:), allocatable :: a360, b360
    real(wp), dimension(:,:), allocatable :: a540, b540
    real(wp), dimension(:,:), allocatable :: a720, b720
    real(wp), dimension(:),   allocatable :: uco2co, uco2ro
    real(wp), dimension(:),   allocatable :: alo
    real(wp), dimension(:),   allocatable :: cor150, cor360, cor540, cor720
    !
  contains
    procedure, public :: load_upper_atmosphere_heating_lw
    procedure, public :: finalize_lw
  end type ty_upper_atmosphere_heating_lw

! Private module data
  type interp_type
     private
     real(wp), pointer :: wgts(:)
     real(wp), pointer :: wgtn(:)
     integer, pointer  :: jjm(:)
     integer, pointer  :: jjp(:)
  end type interp_type
!
! Fomichev IR radiation parameters
!
  integer, parameter :: ncol = 5
  integer, parameter :: nlyrs = 49
  integer, parameter :: nlvls = 50

  integer, parameter :: nrfmc = 59         ! no. of levels of Fomichev parameterization (2-16.5 by 0.25)
  integer, parameter :: nrfmg = 8          ! no. of levels between ground and first calculated level (0-1.75 by 0.25)
  integer, parameter :: nrfm = nrfmc+nrfmg ! total no. of levels of Fomichev paramterization
  integer, parameter :: nrfmnlte = 17      ! no. of levels of NLTE calculation + 1(b.c. at 12.5)
  integer, parameter :: nrfmlte = 43       ! no. of levels of LTE calculation
  integer, parameter :: nrfmlteo3 = 35     ! no. of levels of LTE calculation - O3 ONLY!
  integer, parameter :: nrfmltelv = 9      ! no. of levels used in the LTE integral
  integer, parameter :: nrfmco2 = 4        ! no. of CO2 precalculated profiles

  integer :: i,ix,js

  real(wp), parameter :: rearth = 6.37122e6_wp       ! Earth radius (m)
  real(wp), parameter :: pi = acos(-1._wp)           ! pi

  real(wp), parameter :: co2_mw = 44.009_wp         ! CO2 molecular weight (g/mol)
  real(wp), parameter :: o1_mw = 15.999_wp          ! O molecular weight (g/mol)
  real(wp), parameter :: o2_mw = 31.999_wp          ! O2 molecular weight (g/mol)
  real(wp), parameter :: o3_mw = 47.998_wp          ! O3 molecular weight (g/mol)
  real(wp), parameter :: n2_mw = 28.014_wp          ! N2 molecular weight (g/mol)
  real(wp), parameter :: no_mw = 30.006_wp          ! NO molecular weight (g/mol)
  
  real(wp), parameter :: a10=1.5988_wp            ! reaction constant
  real(wp), parameter :: const=2.63187E11_wp      ! reaction constant
  real(wp), parameter :: constb=9.08795e9_wp      ! reaction constant

  real(wp), parameter :: ptop_co2cool=7.42e-3_wp  ! top pressure level for co2 cool calculation (Pa)

  real(wp) :: o3pxfac(nlyrs)                      ! o3p cooling masking factors on WACCM vertical grids

  logical :: apply_co2_limit = .false.
      
!================================================================================================
contains
!================================================================================================
  function load_upper_atmosphere_heating_lw(this, xr, ig, co2o, ao3, a150, b150, a360, b360, &
                               a540, b540, a720, b720, uco2ro, alo, uco2co, &
                               cor150, cor360, cor540, cor720) result (error_msg)
    !
    class(ty_upper_atmosphere_heating_lw), intent(inout) :: this

    real(wp), dimension(:),   intent(in   ) :: xr, ig, co2o
    real(wp), dimension(:,:), intent(in   ) :: ao3
    real(wp), dimension(:,:), intent(in   ) :: a150, b150
    real(wp), dimension(:,:), intent(in   ) :: a360, b360
    real(wp), dimension(:,:), intent(in   ) :: a540, b540
    real(wp), dimension(:,:), intent(in   ) :: a720, b720
    real(wp), dimension(:),   intent(in   ) :: uco2co, uco2ro
    real(wp), dimension(:),   intent(in   ) :: alo
    real(wp), dimension(:),   intent(in   ) :: cor150, cor360, cor540, cor720
    character(len=128)                      :: error_msg
    !
    ! Local variables
    !
    integer :: nprs67, nprs51, nprs43, nprs35    !
    integer :: nlev9, nlev6, nco2amt
    !
    ! LUT coefficient dimensions
    !
    nprs67 = size(xr)
    nprs51 = size(alo)
    nprs43 = size(a150,dim=1)
    nprs35 = size(ao3,dim=1)
    nlev9  = size(a150,dim=2)
    nlev6  = size(uco2co)
    nco2amt= size(co2o)

    error_msg = ""
    if(.not. extents_are(ao3, nprs35, nlev9)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array ao3 isn't consistently sized"

    if(.not. extents_are(a360, nprs43, nlev9)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array a360 isn't consistently sized"
    if(.not. extents_are(a540, nprs43, nlev9)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array a540 isn't consistently sized"
    if(.not. extents_are(a720, nprs43, nlev9)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array a720 isn't consistently sized"
    if(.not. extents_are(b150, nprs43, nlev9)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array b150 isn't consistently sized"
    if(.not. extents_are(b360, nprs43, nlev9)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array b360 isn't consistently sized"
    if(.not. extents_are(b540, nprs43, nlev9)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array b540 isn't consistently sized"
    if(.not. extents_are(b720, nprs43, nlev9)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array b720 isn't consistently sized"

    if(.not. extents_are(uco2ro, nprs51)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array uco2ro isn't consistently sized"

    if(.not. extents_are(cor150, nlev6)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array cor150 isn't consistently sized"
    if(.not. extents_are(cor360, nlev6)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array cor360 isn't consistently sized"
    if(.not. extents_are(cor540, nlev6)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array cor540 isn't consistently sized"
    if(.not. extents_are(cor720, nlev6)) &
      error_msg = "upper_atmosphere_heating_lw%load_upper_atmosphere_heating_lw(): array cor720 isn't consistently sized"

    if (error_msg /= '') return

    this%nprs67 = nprs67
    this%nprs51 = nprs51
    this%nprs43 = nprs43
    this%nprs35 = nprs35
    this%nlev9  = nlev9
    this%nlev6  = nlev6
    this%nco2amt= nco2amt

    allocate(this%xr(nprs67))
    this%xr = xr

    allocate(this%ig(nlev9))
    this%ig = ig

    allocate(this%co2o(nco2amt))
    this%co2o = co2o

    allocate(this%ao3(nprs35, nlev9))
    this%ao3 = ao3

    allocate(this%a150(nprs43, nlev9), &
             this%a360(nprs43, nlev9), &
             this%a540(nprs43, nlev9), &
             this%a720(nprs43, nlev9), &
             this%b150(nprs43, nlev9), &
             this%b360(nprs43, nlev9), &
             this%b540(nprs43, nlev9), &
             this%b720(nprs43, nlev9))
    this%a150 = a150
    this%a360 = a360
    this%a540 = a540
    this%a720 = a720
    this%b150 = b150
    this%b360 = b360
    this%b540 = b540
    this%b720 = b720
    
    allocate(this%uco2ro(nprs51), &
             this%alo(nprs51), &
             this%uco2co(nlev6))
    this%uco2ro = uco2ro
    this%alo = alo
    this%uco2co = uco2co

    allocate(this%cor150(nlev6), &
             this%cor360(nlev6), &
             this%cor540(nlev6), &
             this%cor720(nlev6))
    this%cor150 = cor150
    this%cor360 = cor360
    this%cor540 = cor540
    this%cor720 = cor720

  end function load_upper_atmosphere_heating_lw

!================================================================================================
  subroutine finalize_lw(this)
  !-------------------------------------------------------------------
  !
  ! Finalize
  !
  !-------------------------------------------------------------------
    class(ty_upper_atmosphere_heating_lw), intent(inout) :: this

    deallocate(this%xr, this%ig, this%co2o, this%ao3, &
               this%a150, this%b150, this%a360, this%b360, &
               this%a540, this%b540, this%a720, this%b720, &
               this%uco2ro, this%alo, this%uco2co, &
               this%cor150, this%cor360, this%cor540, this%cor720)

  end subroutine finalize_lw

!================================================================================================
  subroutine nlte_fomichev_init ( ncol, nlyrs, nlvls, p_lay, p_lev, apply_co2_limit_in, &
                                  k1mb, ktop_co2cool )
  !---------------------------------------------------------------------------------
  !     
  !     Original version from Ray Roble
  !     First adapted to CCM by F. Sassi - November 1999
  !     
  !---------------------------------------------------------------------------------

    ! Input variables
    integer,  intent(in) :: ncol                 ! number of atmospheric columns
    integer,  intent(in) :: nlyrs                ! number of vertical layers
    integer,  intent(in) :: nlvls                ! number of vertical levels

    real(wp), intent(in) :: p_lay(ncol,nlyrs)    ! Pressure layers (Pa)
    real(wp), intent(in) :: p_lev(ncol,nlvls)    ! Pressure levels (Pa)

    logical,  intent(in) :: apply_co2_limit_in   ! .true.  - apply upper limit for co2 defined by co2_limit 
                                                 ! .false. - do not apply upper limit for co2 defined by co2_limit 

    ! Output variables
    integer, dimension(:), allocatable, intent(inout) :: k1mb         ! Highest level index greater than 1 mb (100 Pa)
    integer, dimension(:), allocatable, intent(inout) :: ktop_co2cool ! Layer index defining the top of CO2 cooling 

    ! Local variables
    integer :: i, k
    logical :: apply_co2_limit

    allocate(k1mb(ncol))
    allocate(ktop_co2cool(ncol))

    apply_co2_limit = apply_co2_limit_in
    find_k1mb: do k = nlvls,1,-1
       do i = 1, ncol
        ! Find 1 mbar (or 100 Pa) level.
          if (p_lev(i,k) > 100._wp) then
             k1mb(i) = k
             exit find_k1mb
          endif
       end do
    end do find_k1mb
    
    ktop_co2cool(:) = nlyrs
    do k=nlyrs,1,-1
       do i=1,ncol
          if (p_lay(i,k) < ptop_co2cool) ktop_co2cool(i) = k
       enddo
    enddo

  end subroutine nlte_fomichev_init

!================================================================================================
  subroutine set_matrices( this, ncol, co2sfc, amat, bmat )
  !--------------------------------------------------------------------------------------
  !
  ! Calculate coefficients for the matrix paramerization:
  !
    implicit none

    class(ty_upper_atmosphere_heating_lw), intent(inout) :: this

    integer, intent(in) :: ncol

    real(wp), intent(in) :: co2sfc(ncol)                       ! co2 bottom layer vmr
    real(wp), intent(out) :: amat(ncol,nrfmlte,nrfmltelv)
    real(wp), intent(out) :: bmat(ncol,nrfmlte,nrfmltelv)

  !-----------------------------------------------------------------
  ! Local variables
    real(wp) :: rco2
    real(wp) :: co2int(4), a
    integer :: ii,i,j,isgn

  do ii = 1, ncol
    rco2 = co2sfc(ii)

    amat(ii,1:nrfmlte,1:nrfmltelv)=0.0_wp
    bmat(ii,1:nrfmlte,1:nrfmltelv)=0.0_wp
    do i = 1,nrfmlte
       do j = 1,nrfmltelv

          if((i.le.5).and.(j.eq.2)) goto 1
          isgn = int(sign(1._wp,this%a150(i,j))+sign(1._wp,this%a360(i,j))+ &
               sign(1._wp,this%a540(i,j))+sign(1._wp,this%a720(i,j)))
          co2int(1)=this%a150(i,j)/this%co2o(1)
          co2int(2)=this%a360(i,j)/this%co2o(2)
          co2int(3)=this%a540(i,j)/this%co2o(3)
          co2int(4)=this%a720(i,j)/this%co2o(4)
          if(isgn.eq.-4) then
             co2int(1) = log(-co2int(1))
             co2int(2) = log(-co2int(2))
             co2int(3) = log(-co2int(3))
             co2int(4) = log(-co2int(4))
             a = -exp(a18lin(rco2,this%co2o,co2int,1,4))
          else if (isgn.eq.4) then
             co2int(1) = log(co2int(1))
             co2int(2) = log(co2int(2))
             co2int(3) = log(co2int(3))
             co2int(4) = log(co2int(4))
             a = exp(a18lin(rco2,this%co2o,co2int,1,4))
          else
             call a18int(this%co2o,co2int,rco2,a,4,1)
          end if
          amat(ii,i,j)=a*rco2
          
          isgn = int(sign(1._wp,this%b150(i,j))+sign(1._wp,this%b360(i,j))+ &
               sign(1._wp,this%b540(i,j))+sign(1._wp,this%b720(i,j)))
          co2int(1)=this%b150(i,j)/this%co2o(1)
          co2int(2)=this%b360(i,j)/this%co2o(2)
          co2int(3)=this%b540(i,j)/this%co2o(3)
          co2int(4)=this%b720(i,j)/this%co2o(4)
          if(isgn.eq.-4) then
             co2int(1) = log(-co2int(1))
             co2int(2) = log(-co2int(2))
             co2int(3) = log(-co2int(3))
             co2int(4) = log(-co2int(4))
             a = -exp(a18lin(rco2,this%co2o,co2int,1,4))
          else if (isgn.eq.4) then
             co2int(1) = log(co2int(1))
             co2int(2) = log(co2int(2))
             co2int(3) = log(co2int(3))
             co2int(4) = log(co2int(4))
             a = exp(a18lin(rco2,this%co2o,co2int,1,4))
          else
             call a18int(this%co2o,co2int,rco2,a,4,1)
          end if
          bmat(ii,i,j)=a*rco2
1         continue
       enddo
    enddo
  enddo

  end subroutine set_matrices

!================================================================================================
  subroutine nlte_fomichev_calc (this, ncol, nlyrs, nlvls, k1mb, ktop_co2cool, &
                                 pmid, pint, t, xo2, xo, xo3, xn2, xco2, co2sfc, &
                                 coolf, co2cool_out, o3cool_out, c2scool_out )
!--------------------------------------------------------------------------------
!
!     Original author: F. Sassi (Dec, 1999)
!
!--------------------------------------------------------------------------------
!
!     This routine prepares arrays to be passed to the Fomichev scheme
!
!     RADFMCINTI should have been called at the beginning of the run to
!     create arrays used in this parameterization.
!
!     Concentrations in input arrays are expected in mass mixing ratios 
!     and are converted to volume mixing ratios
!
!     Because the Fomichev scheme has its own grid which runs from 
!     the ground upward, arrays are prepared and inputs are provided
!     with vertical indexing that follows this convention. 
!
!     The pressures at mid-points and interfaces need to be known 
!     Conversion to normalized X coordinate is carried out here.
!
!     Cooling rates are calculated in units of 
!
!                      dT
!           COOLF = Cp ---
!                      dt
!
!     where Cp is the specific heat (ergs g^-1 K^-1), dT/dt is the 
!     temperature tendency (K s^-1). Therefore, units of cooling
!     are 
!     
!         (erg g^-1 K^-1) (K s^-1) = cm^2 s^-3
!
!     COOLF is converted to J/kg/s before output. 
!
!*****************************************************************
!-----------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------

!   Input variables
    class(ty_upper_atmosphere_heating_lw), intent(inout) :: this
    
    integer, intent(in) :: ncol                 ! number of atmospheric columns
    integer, intent(in) :: nlyrs                ! number of model layers
    integer, intent(in) :: nlvls                ! number of model levels
    integer, intent(in) :: k1mb(ncol)           ! Highest level index greater than 1 mb (100 Pa)
    integer, intent(in) :: ktop_co2cool(ncol)   ! Layer index defining the top of CO2 cooling 
                                                ! calculation (below ptop_co2cool)

    real(wp), intent(in) :: pmid(ncol,nlyrs)    ! model pressure at mid-point
    real(wp), intent(in) :: pint(ncol,nlvls)    ! model pressure at interfaces
    real(wp), intent(in) :: t(ncol,nlyrs)       ! Neutral temperature (K)
    real(wp), intent(in) :: xco2(ncol,nlyrs)    ! CO2 profile (mmr)
    real(wp), intent(in) :: xn2(ncol,nlyrs)     ! N2 profile (mmr)
    real(wp), intent(in) :: xo3(ncol,nlyrs)     ! O3 profile (mmr)
    real(wp), intent(in) :: xo(ncol,nlyrs)      ! O profile (mmr)
    real(wp), intent(in) :: xo2(ncol,nlyrs)     ! O2 profile (mmr)
    real(wp), intent(in) :: co2sfc(ncol)        ! co2 bottom layer vmr

!   Output variables
    real(wp), dimension(:,:), allocatable, intent(out) :: coolf       ! Total cooling
    real(wp), dimension(:,:), allocatable, intent(out) :: co2cool_out ! CO2 cooling 
    real(wp), dimension(:,:), allocatable, intent(out) :: o3cool_out  ! O3 cooling
    real(wp), dimension(:,:), allocatable, intent(out) :: c2scool_out ! Cooling to Space

!   Local variables
      
    real(wp) :: rmo2                         ! O2 molecular weight
    real(wp) :: rmo                          ! O molecular weight
    real(wp) :: rmn2                         ! N2 molecular weight
    real(wp) :: rmco2                        ! CO2 molecular weight
    real(wp) :: rmo3                         ! O3 molecular weight
    real(wp) :: xnorm(ncol,nlyrs)            ! normalized X p.s.h. at midpoints
    real(wp) :: xnori(ncol,nlvls)            ! normalized X p.s.h. at interfaces
    real(wp) :: dxnorm(ncol,nlyrs)           ! xnorm(k+1)-xnorm(k)
    real(wp) :: presm(ncol,nlyrs)            ! pressure at midpoint (dyn/cm^2)
    real(wp) :: presi(ncol,nlvls)            ! pressure at interfaces (dyn/cm^2)
    real(wp) :: dpi(ncol,nlyrs)              ! pressure diff. between interfaces (dyn/cm^2)
    real(wp) :: mwair(ncol,nlyrs)            ! mean air molecular weight (g/mole)
    real(wp) :: ndenair(ncol,nlyrs)          ! mean air number density (cm**-3)
    real(wp) :: colco2(ncol,nlyrs)           ! CO2 column number density
    real(wp) :: uco2(ncol,nrfm)              ! column CO2
    real(wp) :: dummyg(nlyrs)                ! dummy
    real(wp) :: dummyx(nlyrs)                ! dummy
    real(wp) :: dummyf(nrfm)                 ! dummy
    real(wp) :: hco2(ncol,nrfm)              ! CO2 cooling in Fomichev grid
    real(wp) :: ho3(ncol,nrfm)               ! O3 cooling in Fomichev grid
    real(wp) :: tf(ncol,nrfm)                ! neutral temp interpolated to Fomichev grid
    real(wp) :: vn2f(ncol,nrfm)              ! N2 vmr interpolated to Fomichev grid
    real(wp) :: vo3f(ncol,nrfm)              ! O3 vmr interpolated to Fomichev grid
    real(wp) :: vof(ncol,nrfm)               ! O vmr interpolated to Fomichev grid
    real(wp) :: vco2f(ncol,nrfm)             ! CO2 vmr interpolated to Fomichev grid
    real(wp) :: vo2f(ncol,nrfm)              ! O2 vmr interpolated to Fomichev grid
    real(wp) :: mwairf(ncol,nrfm)            ! Mean air molecular weight interpolated to Fomichev grid
    real(wp) :: ndenf(ncol,nrfm)             ! Mean air no. density interpolated to Fomichev grid
    real(wp) :: flux(ncol)                   ! Flux boundary condition for cool-to-space
    real(wp) :: vo2(ncol,nlyrs)              ! O2 vmr
    real(wp) :: vo(ncol,nlyrs)               ! O vmr
    real(wp) :: vo3(ncol,nlyrs)              ! O3 vmr
    real(wp) :: vn2(ncol,nlyrs)              ! N2 vmr
    real(wp) :: vco2(ncol,nlyrs)             ! CO2 vmr
    real(wp) :: co2cooln(ncol,nlyrs)         ! CO2 cooling 
    real(wp) :: o3cooln(ncol,nlyrs)          ! O3 cooling
    real(wp) :: hc2s(ncol,nlyrs)             ! cool to space heating
    real(wp) :: ps0                          ! Reference (surface) pressure
    real(wp) :: ti(ncol,nlyrs)               ! T(NLYRS:1:-1)
    real(wp) :: alam(ncol,nrfm)              ! LAMBDA
    real(wp) :: djm(ncol,nrfm)               ! DJM in recurrence formula
    real(wp) :: dj0(ncol,nrfm)               ! DJ0 in recurrence formula
    real(wp) :: aajm(ncol,nrfm)              ! AAJM in recurrrence formula
    real(wp) :: aaj0(ncol,nrfm)              ! AJJ0 in recurrence formula
    real(wp) :: zhgt(ncol,nlyrs)             ! approx. elevation in cm
    real(wp) :: gravhgt(ncol,nlyrs)          ! accelration of gravity in cm/s^2
    real(wp) :: wrk(ncol)
      
    integer :: i,k

    real(wp), parameter :: co2_limit = 720.e-6_wp

    real(wp) :: akbl                        ! Boltzman constant
    real(wp) :: anav                        ! Avogadro Number (molecules/mole)
    real(wp) :: grav0                       ! gravitational constant (cm/s2)
    real(wp) :: ur                          ! universal gas constant (R_star)
    real(wp) :: arad                        ! planet's radius (cm)

    character(len=200) :: errmsg

!----------------------------------------------------------------

! Allocate output cooling arrays
    allocate(coolf(ncol,nlyrs))
    allocate(co2cool_out(ncol,nlyrs))
    allocate(o3cool_out(ncol,nlyrs))
    allocate(c2scool_out(ncol,nlyrs))

! Define molecular weights
    rmco2 = co2_mw
    rmo3  = o3_mw
    rmo2  = o2_mw
    rmo   = o1_mw
    rmn2  = n2_mw

    coolf(1:ncol,1:nlyrs)=0.0_wp

! Initialize planet's radius (cm)
    arad=rearth*1e2_wp

!-----------------------------------------------------------------
!     The pressure at mid point is converted to normalized x coordinate
!                  X = LN (P0 / PRES)
!     where P0 = 1e6 dyn/cm^2 ; PRES is in dyn/cm^2
!     Note: to convert from Pa to dyn/cm^2, multyply by 10
!-----------------------------------------------------------------
    ps0=1.000e6_wp     ! dyn/cm2  (Pa*10)
    do k=1,nlyrs
       do i=1,ncol
          presm(i,k) = pmid(i,k)*10._wp
          xnorm(i,k) = log(ps0/presm(i,k))
! Calculate pressure at interfaces
          presi(i,k) = pint(i,k)*10._wp
          xnori(i,k) = log(ps0/presi(i,k))
       enddo
    enddo
    presi(:ncol,nlvls) = pint(:ncol,nlvls)*10._wp
    xnori(:ncol,nlvls) = log(ps0/presi(:ncol,nlvls))

!-----------------------------------------------------------------
!     Calculate layer thikcness (DPI).
!     For each pressure interface the following is true:
!     
!            Pint(k) = P0 * exp (-Xint(k))
!
!     Thus,
!
!            DPI(k)= Pint(k)-Pint(k+1)= P0 * exp (-Xint(k)) * DX
!
!     where,
!
!            DX = 1 - exp ( - (Xint(k+1)-Xint(k)) )
!
!-----------------------------------------------------------------
    do k=1,nlyrs
       do i=1,ncol
          dxnorm(i,k)=xnori(i,k+1)-xnori(i,k)
       enddo
    enddo

! Pressure difference between interfaces (positive downward)
    do k=1,nlyrs
       do i=1,ncol
          dpi(i,k)=ps0*exp(-xnorm(i,k))*(1._wp-exp(-dxnorm(i,k)))
       enddo
    enddo

!-----------------------------------------------------------------
!     Calculate molecular weight (g/mol) of mean air MWAIR:
!      
!               MWAIR= 1. / [ Sum(i) MMR(i)/MW(i) ] 
!
!     where MMR(i) are  mass mixing ratio of O2,
!     O, N2, and MW(i) are the corresponding
!     molecular weights.
!-----------------------------------------------------------------
    mwair(1:ncol,1:nlyrs)=0.0_wp
    do k=1,nlyrs
       do i=1,ncol
          mwair(i,k)=1._wp/ (                   &
              xo2(i,k)/rmo2                     &
              +xo(i,k)/rmo                      &
              +xn2(i,k)/rmn2                    &
              +xo3(i,k)/rmo3                    &
              +xco2(i,k)/rmco2                  &
              )
       enddo
    enddo

!
!     Elevation is calculated via integration of 
!     the hydrostatic equation 
!
!           Sum(k) g(k) dz = Sum(k) PHI(k)
!
!     where 
!                     g0*a^2
!               g(k)= ------
!                     (a+z)^2
!
!     where g0=980 cm/s^2;a=6.37e8
!
!     and
!
!                                UR * T(i)  DPI(i)
!            PHI(k) = Sum(i=1,k)  ------    ----
!                                  MW(i)    P(i)
!     
!     where UR is the gas universal constant, T is temperature,
!     MW is the mean air molecular weight, DPI is the pressure 
!     layer thickness and P is the mid-point pressure
!     Then,
!
!                   PHI * a
!              z = ---------
!                  (a*g0-PHI)
!
!     do i=1,ncol
!        phi=0.0
!        do k=1,nlyrs
!           kinv=nlyrs-k+1
!           phi=phi+ur*t(i,kinv)/(mwair(i,k))*dpi(i,k)/presm(i,k)
!           zhgt(i,k)=phi*arad/(arad*grav0-phi)
!           grav(i,k)=grav0*(arad/(arad+zhgt(i,k)))**2
!        enddo
!     enddo

    wrk(:ncol) = 0._wp
!   Convert grav from m/s2 to cm/s2
    grav0 = grav * 1.e2_wp
!   Convert Boltzman constant to cgs units
    akbl = k_boltz * 1.e7_wp
!   Convert universal gas constant to cgs units [erg/(mol K)]
    ur = R_univ_gconst * 1.e7_wp
!   Set Avogadro's constant
    anav = avogad

    do k=1,nlyrs
       do i=1,ncol
          wrk(i) = wrk(i) + ur*t(i,k)/(mwair(i,k))*dpi(i,k)/presm(i,k)
          zhgt(i,k) = wrk(i)*arad/(arad*grav0 - wrk(i))
          gravhgt(i,k) = grav0*(arad/(arad + zhgt(i,k)))**2
       enddo
    enddo

    do k = 1,nlyrs
       do i=1,ncol
!-----------------------------------------------------------------
!     Convert mmr to vmr
!-----------------------------------------------------------------
          vo2 (i,k) = xo2 (i,k) *mwair(i,k)/rmo2
          vo  (i,k) = xo  (i,k) *mwair(i,k)/rmo
          vo3 (i,k) = xo3 (i,k) *mwair(i,k)/rmo3
          vn2 (i,k) = xn2 (i,k) *mwair(i,k)/rmn2
          vco2(i,k) = xco2(i,k) *mwair(i,k)/rmco2

! CGB - The Formichev scheme was not designed to support CO2 > 720 ppmv, so
! limit the amount of CO2 used to 720 ppmv. This keeps the model stable, but
! may yield an incorrect scientific result. It would be nice to extend this
! routine to support higher CO2 values. Putting the limiter here means that
! that the other constituents will have their proper mixing ratio caclulated
! (i.e. the mwair is correct), but vco2 will be limited.  Abort the run if CO2
! exceeds the limit at altitudes above 1 mbar unless apply_co2_limit=.true.

          if (vco2(i,k)>co2_limit) then
             if ((.not.apply_co2_limit) .and. (k<k1mb(i))) then
                 write(errmsg,fmt='(a50,f8.2,a4)') &
                 'nlte_fomichev_calc: CO2 has exceeded the limit of ', co2_limit, ' hPa'
                  stop
             endif
          endif

!-----------------------------------------------------------------
!     Calculate mean air number density ( cm^(-3) )
!
!                              P(k)
!                    n  = -------------
!                              Kb*T(k)
!     where:
!           P(k)    is pressure at mid-point
!           AKBL    is Boltzman constant
!           T       is neutral temperature
!-----------------------------------------------------------------
          ndenair(i,k) = presm(i,k)/(akbl*t(i,k))
       enddo

    enddo

    ! apply the CO2 limiter for all levels and columns
    where ( vco2(:ncol,:) > co2_limit )
       vco2(:ncol,:) = co2_limit
    end where
      
!-----------------------------------------------------------------
!     Calculate CO2 vertical column above each level
!     
!     At each mid-point vertical level (j), the following sum is calculated
!
!                                         
!            COLCO2(j) = Sum_(i=ztop:j:-1) NDENAIR(I)*VCO2(I)*DZ = ...
!                      
!                                            ANAV * VCO2(I)
!                      = Sum_(i=nlyrs:j:-1)  ---------------- * DP
!                                            GRAVHGT * MWAIR(I)
!
!     where ANAV is the Avogadro no., VCO2 is CO2 vmr, GRAVHGT is the
!     accelaration of gravity, MWAIR is mean molecular weight, DP
!     is the pressure increment downward.
!     As boundary condition at NLYRS, it is assumed that VCO2 
!     and NDENAIR stay constant above that level.
!-----------------------------------------------------------------
!     do i=1,ncol
!        colco2(i,nlyrs)=anav/gravhgt(i,nlyrs)*vco2(i,nlyrs)/mwair(i,nlyrs)*dpi(i,nlyrs)
!        do k=nlyrs-1,1,-1
!           colco2(i,k)=colco2(i,k+1)                                       &
!               +anav/gravhgt(i,k)*vco2(i,k)/mwair(i,k)*dpi(i,k)
!        enddo
!     enddo

    colco2(:ncol,nlyrs) = anav/gravhgt(:ncol,nlyrs)*vco2(:ncol,nlyrs) &
			   /mwair(:ncol,nlyrs)*dpi(:ncol,nlyrs)
    do k=nlyrs-1,1,-1
       do i=1,ncol
          colco2(i,k)=colco2(i,k+1)                                       &
         +anav/gravhgt(i,k)*vco2(i,k)/mwair(i,k)*dpi(i,k)
       enddo
    enddo

!-----------------------------------------------------------------
!     Linear interpolation from input vertical  grid defined in 
!     XNORM(1:NLYRS) to the Fomichev grid defined in XR(1:NRFM)
!     Interpolation is carried out using mod PRFLINV which
!     is adapted from A18LINV (originally written for TIME/GCM).
!     All fields are interpolated. If XR levels are beyond the
!     limit of the actual input grid, zero values are inserted in
!     the interpolated arrays. Proper handling of the arrays is done
!     in VICCOOLN.
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!     Create TI array which contains T(NLYRS:1:-1)
!-----------------------------------------------------------------
    ti(1:ncol,1:nlyrs)=t(1:ncol,1:nlyrs)

    do i=1,ncol

       dummyx(1:nlyrs)=xnorm(i,1:nlyrs)

!     Temperature
       dummyg(1:nlyrs)=ti(i,1:nlyrs)
       call a18linvne (this%xr,dummyf,dummyx,dummyg,nlyrs,nrfm)
       tf(i,1:nrfm)=dummyf(1:nrfm)

!     O3
       dummyg(1:nlyrs)=vo3(i,1:nlyrs)
       call a18linvne (this%xr,dummyf,dummyx,dummyg,nlyrs,nrfm)
       vo3f(i,1:nrfm)=dummyf(1:nrfm)

!     O2
       dummyg(1:nlyrs)=vo2(i,1:nlyrs)
       call a18linvne (this%xr,dummyf,dummyx,dummyg,nlyrs,nrfm)
       vo2f(i,1:nrfm)=dummyf(1:nrfm)

!     N2
       dummyg(1:nlyrs)=vn2(i,1:nlyrs)
       call a18linvne (this%xr,dummyf,dummyx,dummyg,nlyrs,nrfm)
       vn2f(i,1:nrfm)=dummyf(1:nrfm)

!     O
       dummyg(1:nlyrs)=vo(i,1:nlyrs)
       call a18linvne (this%xr,dummyf,dummyx,dummyg,nlyrs,nrfm)
       vof(i,1:nrfm)=dummyf(1:nrfm)
         
!     CO2
       dummyg(1:nlyrs)=vco2(i,1:nlyrs)
       call a18linvne (this%xr,dummyf,dummyx,dummyg,nlyrs,nrfm)
       vco2f(i,1:nrfm)=dummyf(1:nrfm)

!     COLCO2
       dummyg(1:nlyrs)=colco2(i,1:nlyrs)
       call a18linvne (this%xr,dummyf,dummyx,dummyg,nlyrs,nrfm)
       uco2(i,1:nrfm)=dummyf(1:nrfm)

!     DEN
       dummyg(1:nlyrs)=ndenair(i,1:nlyrs)
       call a18linvne (this%xr,dummyf,dummyx,dummyg,nlyrs,nrfm)
       ndenf(i,1:nrfm)=dummyf(1:nrfm)
         
!     AM
       dummyg(1:nlyrs)=mwair(i,1:nlyrs)
       call a18linvne (this%xr,dummyf,dummyx,dummyg,nlyrs,nrfm)
       mwairf(i,1:nrfm)=dummyf(1:nrfm)
         
    enddo

    ! Use recurrence relation to calculate AL coefficents
    call recur (this,ncol,uco2,tf,vn2f,vo2f,vof,ndenf,alam, &
                djm,dj0,aajm,aaj0)


    ! Do LTE and NLTE parts of cooling
    call viccooln (this,ncol,co2sfc,alam,djm,dj0,aajm,aaj0, &
                   tf,vco2f,vo3f,mwairf,flux,hco2,ho3)

!     Interpolate from Fomichev grid to CCM grid
    do i=1,ncol
       dummyx(1:nlyrs)=xnorm(i,1:nlyrs)

!     HCO2
       dummyf(1:nrfm)=hco2(i,1:nrfm)
       call a18linvne (dummyx,dummyg,this%xr,dummyf,nrfm,nlyrs)
       co2cooln(i,1:nlyrs)=dummyg(1:nlyrs)

!     HO3
       dummyf(1:nrfm)=ho3(i,1:nrfm)
       call a18linvne (dummyx,dummyg,this%xr,dummyf,nrfm,nlyrs)
       o3cooln(i,1:nlyrs)=dummyg(1:nlyrs)

    enddo

!     Do cool-to-space component of cooling
    call cool2space (ncol,nlyrs,ti,mwair,vn2,vo2,vo,vco2,ndenair,xnorm,flux,hc2s)

!     Calculate total cooling

   ! Above ptop_co2cool use cool to space approx.
    do i=1,ncol
       do k=nlyrs,ktop_co2cool(i)+1,-1
          ! Convert to J/kg/s
          coolf(i,k) = (o3cooln(i,k) + hc2s(i,k)) * 1.e-4_wp
       enddo
    enddo

      ! Below ptop_co2cool use nlte calculation
    do  i=1,ncol
       do k=ktop_co2cool(i),1,-1
          ! Convert to J/kg/s
          coolf(i,k) = (co2cooln(i,k) + o3cooln(i,k)) * 1.e-4_wp
       enddo
    enddo

    ! diagnostics ...
    do k=1,nlyrs
       co2cool_out(:ncol,k) = co2cooln(:ncol,k) * 1.e-4_wp
       o3cool_out(:ncol,k) = o3cooln(:ncol,k) * 1.e-4_wp
       c2scool_out(:ncol,k) = hc2s(:ncol,k) * 1.e-4_wp
    enddo

  end subroutine nlte_fomichev_calc

!====================================================================================
  subroutine a18linvne (x,y,xn,yn,n,imax)                               
!--------------------------------------------------------------------------------
!     ****                                                               
!     ****     This procedure performs linear interpolation within the   
!     ****     table defined by the N points (XN(NN),Y(NN)).             
!     ****     Where:                                                    
!     ****                                                               
!     ****       NN = 1,N,1                                              
!     ****                                                               
!     ****       XN(NN) < XN(NN+1) for NN = 1,N-1            
!     ****                                                               
!     ****     Parameters:                                               
!     ****                                                               
!     ****       X(IMAX) = array of IMAX x-values at which linear        
!     ****                 interpolation is required                     
!     ****                                                               
!     ****       XN(N) = array of N abscissae at which function values   
!     ****               are given                                       
!     ****                                                               
!     ****       YN(N) = function values corresponding to abscissae,     
!     ****               XN(N)                                           
!     ****                                                               
!     ****     Output:                                                   
!     ****                                                               
!     ****      Y(IMAX)  The IMAX interpolated values are                
!     ****                     returned in this array                    
!     ****                                                               
!     
!     It has been modified as follows: 
!     if points X are outside the range X(1)..X(N), the last values
!     are assigned. That is
!     IF X(I) > XN(N) THEN Y(I)=YN(N)
!     IF X(I) < XN(1) THEN Y(I)=YN(1)
!     
!--------------------------------------------------------------------------------

    implicit none

!   Arguments
    integer, intent(in) :: imax
    integer, intent(in) :: n
    real(wp), intent(out) :: y(imax)
    real(wp), intent(in) :: x(imax)
    real(wp), intent(in) :: xn(n)
    real(wp), intent(in) :: yn(n)

!   Local variables
    integer kk(imax)                 
    integer i,nn

!     ****                                                               
!     ****     Where:                                                    
!     ****       Y(IMAX) is vector output                                
!     ****                                                               
!     ****       KK is work space                                        
!     ****                                                               
!     ****                                                               
!     ****     Initialize array KK                                       
!     ****                                                               

    do i = 1,imax                                                      
       kk(i) = 0                                                        
    enddo                                                              

!     ****                                                               
!     ****     Locate interval in (XN,YN) in table containing X(I)       
!     ****                                                               

    do nn = 1,n-1                                                      
       do i = 1,imax                                                    
          kk(i) = merge(nn+1,kk(i),(xn(nn+1)-x(i))*(x(i)-xn(nn))>=0._wp)    
       enddo                                                            
    enddo                                                              

!     ****                                                               
!     ****     Check for                                                 
!     ****                                                               
!     ****       X(I) < XN(1),  X(I) > X(N)                              
!     ****                                                               
!     ****       and use linear extrapolation if necessary               
!     ****                                                               

    do i = 1,imax                                                      
       kk(i) = merge(-1,kk(i),xn(1)-x(i)>=0._wp)                            
       kk(i) = merge(-2,kk(i),x(i)-xn(n)>=0._wp)                            
    enddo                                                              

!     ****                                                               
!     ****     Perform interpolation prescribed above                    
!     ****                                                               

    y(:) = 0._wp
      
    do i = 1,imax      

       if (kk(i).gt.0) then
          y(i) = (                            &
               yn(kk(i)-1)*(xn(kk(i))-x(i))   &
               + yn(kk(i))*(x(i)-xn(kk(i)-1)) &
               )/(xn(kk(i))-xn(kk(i)-1))         
       else if (kk(i).eq.-1) then
          y(i)=yn(1)
       else if (kk(i).eq.-2) then
          y(i)=yn(n)
       endif

    enddo                                                              

  end subroutine a18linvne

!=======================================================================================
  subroutine recur (this,ncol,uco2,tf,vn2f,vo2f,vof,ndenf,alam, &
                    djm,dj0,aajm,aaj0)
!-------------------------------------------------------------------------
!
!     Originally written by R. Roble,  modified by F. Sassi (Nov., 1999)

!-------------------------------------------------------------------------
    implicit none

    class(ty_upper_atmosphere_heating_lw), intent(inout) :: this
    
    integer, intent(in) :: ncol                          ! number of atmospheric columns
      
    real(wp), intent(in) :: UCO2(ncol,nrfm)
    real(wp), intent(in) :: tf(ncol,nrfm)
    real(wp), intent(in) :: vn2f(ncol,nrfm)
    real(wp), intent(in) :: vo2f(ncol,nrfm)
    real(wp), intent(in) :: vof(ncol,nrfm)
    real(wp), intent(in) :: ndenf(ncol,nrfm)

    real(wp), intent(out) :: alam(ncol,nrfm)
    real(wp), intent(out) :: djm(ncol,nrfm)
    real(wp), intent(out) :: dj0(ncol,nrfm)
    real(wp), intent(out) :: aajm(ncol,nrfm)
    real(wp), intent(out) :: aaj0(ncol,nrfm)

    real(wp) CO2INT(nrfmco2)
    real(wp) UREF(nrfmco2)
    real(wp) A(ncol)                  
    real(wp) COR(ncol)
    real(wp) UC(ncol)                                           
    real(wp) al(ncol,nrfm)

    real(wp) tt
    real(wp) y
    real(wp) zn2
    real(wp) zo2
    real(wp) zz
    real(wp) rko

    integer i,k,km,ks

!****this constant should be moved to an intialization routine
!                                                                        
!                                                                        
!     ****  UCO2 (CO2 COLUMN AMOUNT) FOR CO2                         
!                                                                        
!     **** CALCULATE COEFICIENTS FOR THE RECCURENCE FORMULA:             
!                                                                        
!     **** BETWEEN X=12.5 AND 13.75 THESE COEFFICIENTS (AL) ARE          
!     **** CALCULATED USING CORRECTIONS TO THE ESCAPE FUNCTION.          
!     **** STARTING FROM X=14.00 AND ABOVE THE PARAMETERIZATION          
!     **** COEFFICIENTS ARE EQUAL TO THE ESCAPE FUNCTION.                
!                        
    al(1:ncol,1:nrfm)=0.0_wp
    ks=0
    do k=1,nrfm                                                         
         
       if (this%xr(k).ge.12.5_wp .and. this%xr(k).le.13.75_wp) then
          ks=ks+1

          co2int(1) = this%cor150(ks)
          co2int(2) = this%cor360(ks)
          co2int(3) = this%cor540(ks)
          co2int(4) = this%cor720(ks)
          uref(1) = this%uco2co(ks)*150._wp/360._wp
          uref(2) = this%uco2co(ks)
          uref(3) = this%uco2co(ks)*540._wp/360._wp
          uref(4) = this%uco2co(ks)*720._wp/360._wp
          do  i=1,ncol
             uc(i) = uco2(i,k)
          enddo
          call a18linv(uc,a,this%uco2ro,this%alo,51,ncol)
          call a18linv(uc,cor,uref,co2int,4,ncol)
          do i=1,ncol
             al(i,k) = exp(cor(i)+a(i))
          enddo

       endif
    enddo

    do k=1,nrfm

       if (this%xr(k).ge.14.00_wp) then

          do i=1,ncol
             uc(i) = uco2(i,k)
          enddo
          call a18linv(uc,a,this%uco2ro,this%alo,51,ncol)
          do  i=1,ncol
             al(i,k) = exp(a(i))
          enddo
            
       endif

    enddo

! 
!     Calculate ALAM
!                     
    alam(1:ncol,1:nrfm)=0.0_wp
    do  k=1,nrfm

!     ALAM is used only for p.s.h. >= 12.5
!     If the current level is below 12.5 s.h., then do nothing
       if (this%xr(k).ge.12.5_wp) then

          do i=1,ncol                                                  
!                                                                     
!     ****  CO2-O2 AND CO2-N2 V-T CONSTANTS                           
!                                                                     
             tt   = tf(i,k) 
             y    = tt**(-1._wp/3._wp)                                          
             zn2  = 5.5e-17_wp*sqrt(tt)+6.7e-10_wp*exp(-83.8_wp*y)            
             zo2  = 1.e-15_wp*exp(23.37_wp-230.9_wp*y+564._wp*y*y)            
             rko  = 3.0e-12_wp                                                    
!                                                                     
!     ****  COLLISIONAL DEACTIVATION RATE:                            
!                                                                     
             zz   = (vn2f(i,k)*zn2 +  vo2f(i,k)*zo2 +  vof (i,k)*rko)*ndenf(i,k)
!
!     ****
!
             alam(i,k) = a10/( a10+zz )

          enddo  ! end-loop in longitude
         
       endif

    enddo        ! end-loop in levels

!     Calculate coefficients of recurrence formula
!     This coefficients are used for 12.75=< p.s.h.=<16.5
!     Outside this range do nothing
!     It uses ALAM (p.s.h. >= 12.5)

    djm(1:ncol,1:nrfm)=0.0_wp
    dj0(1:ncol,1:nrfm)=0.0_wp
    aajm(1:ncol,1:nrfm)=0.0_wp
    aaj0(1:ncol,1:nrfm)=0.0_wp
    do k=1,nrfm

       if (this%xr(k).ge.12.75_wp .and. this%xr(k).le.16.5_wp) then

          km=k-1

          do i=1,ncol

             djm(i,k)  = +.25_wp*(3._wp*al(i,km) +    al(i,k) )
             dj0(i,k)  = +.25_wp*(   al(i,km) + 3._wp*al(i,k) )

             aajm(i,k) = 1._wp-alam(i,km) * ( 1._wp-djm(i,k) )
             aaj0(i,k) = 1._wp-alam(i,k ) * ( 1._wp-dj0(i,k) )

          enddo    ! end-loop in longitude

       endif

    enddo          ! end-loop in levels

  end subroutine recur
                                                                
!======================================================================================
  subroutine a18linv (x,y,xn,yn,n,imax)                               
!-------------------------------------------------------------------------
!
!     ****                                                               
!     ****     This procedure performs linear interpolation within the   
!     ****     table defined by the N points (XN(NN),Y(NN)).             
!     ****     Where:                                                    
!     ****                                                               
!     ****       NN = 1,N,1                                              
!     ****                                                               
!     ****       XN(NN) < XN(NN+1) for NN = 1,N-1            
!     ****                                                               
!     ****     Parameters:                                               
!     ****                                                               
!     ****       X(IMAX) = array of IMAX x-values at which linear        
!     ****                 interpolation is required                     
!     ****                                                               
!     ****       XN(N) = array of N abscissae at which function values   
!     ****               are given                                       
!     ****                                                               
!     ****       YN(N) = function values corresponding to abscissae,     
!     ****               XN(N)                                           
!     ****                                                               
!     ****     Output:                                                   
!     ****                                                               
!     ****      Y(IMAX)  The IMAX interpolated values are                
!     ****                     returned in this array                    
!     ****                                                               
!-------------------------------------------------------------------------

    implicit none

!   Input variables
    integer i, imax, n
    real(wp) y(imax)
    real(wp) x(imax)
    real(wp) xn(n)
    real(wp) yn(n)

!   Local variables
    integer kk(imax)                 
    integer nn

!     ****                                                               
!     ****     Where:                                                    
!     ****       Y(IMAX) is vector output                                
!     ****                                                               
!     ****       KK is work space                                        
!     ****                                                               
!     ****                                                               
!     ****     Initialize array KK                                       
!     ****                                                               

    do i = 1,imax                                                      
       kk(i) = 0                                                        
    enddo                                                              

!     ****                                                               
!     ****     Locate interval in (XN,YN) in table containing X(I)       
!     ****                                                               

    do nn = 1,n-1                                                      
       do i = 1,imax                                                    
          kk(i) = merge(nn+1,kk(i),(xn(nn+1)-x(i))*(x(i)-xn(nn))>=0._wp)    
       enddo                                                            
    enddo                                                              

!     ****                                                               
!     ****     Check for                                                 
!     ****                                                               
!     ****       X(I) < XN(1),  X(I) > X(N)                              
!     ****                                                               
!     ****       and use linear extrapolation if necessary               
!     ****                                                               

    do i = 1,imax                                                      
       kk(i) = merge(2,kk(i),xn(1)-x(i)>=0._wp)                            
       kk(i) = merge(n,kk(i),x(i)-xn(n)>=0._wp)                            
    enddo                                                              

!     ****                                                               
!     ****     Perform interpolation prescribed above                    
!     ****                                                               

    do i = 1,imax                                                      
       y(i) = (yn(kk(i)-1)*(xn(kk(i))-x(i)) + yn(kk(i))*    &
           (x(i)-xn(kk(i)-1)))/(xn(kk(i))-xn(kk(i)-1))         
    enddo                                                              

  end subroutine a18linv

!========================================================================================
  subroutine viccooln (this,ncol,co2sfc,alam,djm,dj0,aajm,aaj0, &
                       tv,co2,o3,am,flux,hco2,ho3)
!-----------------------------------------------------------------
!
!     Original version from Ray Roble.
!     Adapted to CCM by F. Sassi (Nov. 1999)
!
!-----------------------------------------------------------------
!
!     **** This is the mod that calculates LTE and NLTE components
!     **** of the cooling rates.
!
! XL(17)  - the parameters for NLTE region (12.5 <= X <= 16.5)    

    implicit none

    ! Input variables
    class(ty_upper_atmosphere_heating_lw), intent(inout) :: this

    integer, intent(in) :: ncol                          ! number of atmospheric columns

    real(wp), intent(in) :: tv(ncol,nrfm)               ! neutral temp interpolated to Fomichev grid
    real(wp), intent(in) :: o3(ncol,nrfm)               ! O3 vmr interpolated to Fomichev grid
    real(wp), intent(in) :: am(ncol,nrfm)               ! Mean air molecular weight interpolated to Fomichev grid
    real(wp), intent(in) :: co2(ncol,nrfm)              ! CO2 vmr interpolated to Fomichev grid
    real(wp), intent(in) :: djm(ncol,nrfm)              ! DJM coefficient in recurrence formula
    real(wp), intent(in) :: dj0(ncol,nrfm)              ! DJ0 coefficient in recurrence formula
    real(wp), intent(in) :: aajm(ncol,nrfm)             ! AAJM coefficient in recurrence formula
    real(wp), intent(in) :: aaj0(ncol,nrfm)             ! AAJ0 coefficient in recurrence formula
    real(wp), intent(in) :: alam(ncol,nrfm)             ! LAMBDA
    real(wp), intent(in) :: co2sfc(ncol)                ! co2 bottom layer vmr

!   Local variables
    real(wp) fu(ncol,nrfm)
    real(wp) fo3(ncol,nrfm)
    real(wp) h1(ncol)
    real(wp) h2(ncol)
    real(wp) h3(ncol)

    integer i,k,jj,ks,jjs

!   Output variables
    real(wp), intent(out) :: flux(ncol)            ! Flux boundary condition for cool-to-space
    real(wp), intent(out) :: hco2(ncol,nrfm)       ! CO2 cooling in Fomichev grid
    real(wp), intent(out) :: ho3(ncol,nrfm)        ! O3 cooling in Fomichev grid

    real(wp) :: amat(ncol,nrfmlte,nrfmltelv)
    real(wp) :: bmat(ncol,nrfmlte,nrfmltelv)

!--------------------------------------------------------------------
!  update the amat and bmat matrices with time-dependent surf CO2 
!--------------------------------------------------------------------
    call set_matrices(this,ncol,co2sfc,amat,bmat)

    hco2(1:ncol,1:nrfm)=0.0_wp
    ho3(1:ncol,1:nrfm)=0.0_wp
    flux(1:ncol)=0.0_wp
!                                                                    
! grid levels for height integration   (p.s.h. distance = 0.25*IG)   
!
    do k=1,nrfm
       do i=1,ncol                                                 
          fu(i,k)=exp(-960.217_wp/tv(i,k))                                 
          fo3(i,k)=exp(-1500._wp/tv(i,k))                                  
       enddo
    enddo
!                                                                     
! calculate the heating rates for layer below s.h.p. = 12.5           
!   15 um CO2 + 9.6 um O3:                                            
!                                                                     
!     **** COOLING RATE IN BOTH O3 AND CO2 BANDS (X=2-10.5) MATRIX    
!     **** APPROACH                                                   
!
!     Adding KS=K+8 maps NFRMC into NFRM
    do k=1,5                                                      
       ks = k+8                                                     
       do i=1,ncol                                                  
          h2(i) = (amat(i,k,1)+bmat(i,k,1)*fu(i,ks))*fu(i,1)            
          h3(i) = this%ao3(k,1)*fo3(i,1)                                 
       enddo
       do jj=3,nrfmltelv
          jjs = ks+this%ig(jj)                                             
          do i=1,ncol                                              
             h2(i) = h2(i)+(amat(i,k,jj)+bmat(i,k,jj)*fu(i,ks))*fu(i,jjs)    
             h3(i) = h3(i)+this%ao3(k,jj)*fo3(i,jjs)                          
          enddo
       enddo
       do i=1,ncol                                               
          hco2(i,ks) = h2(i)                                           
          ho3(i,ks) = h3(i)*o3(i,ks)                                   
       enddo
    enddo

    do k=6,18                                                      
       ks = k+8                                                     
       do i=1,ncol                                               
          h2(i) = (amat(i,k,1)+bmat(i,k,1)*fu(i,ks))*fu(i,1)               
          h3(i) = this%ao3(k,1)*fo3(i,1)                                    
       enddo
       do jj=2,nrfmltelv
          jjs = ks+this%ig(jj)                                              
          do i=1,ncol                                               
             h2(i) = h2(i)+(amat(i,k,jj)+bmat(i,k,jj)*fu(i,ks))*fu(i,jjs)     
             h3(i) = h3(i)+this%ao3(k,jj)*fo3(i,jjs)                           
          enddo
       enddo
       do i=1,ncol                                            
          hco2(i,ks) = h2(i)                                            
          ho3(i,ks) = h3(i)*o3(i,ks)                                    
       enddo
    enddo

    do k=19,35                                                        
       ks = k+8                                                        
       do i=1,ncol                                                  
          h2(i) = 0._wp                                                      
          h3(i) = 0._wp                                                      
       enddo
       do jj=1,nrfmltelv
          jjs = ks+this%ig(jj)                                                 
          do  i=1,ncol                                                 
             h2(i) = h2(i)+(amat(i,k,jj)+bmat(i,k,jj)*fu(i,ks))*fu(i,jjs)        
             h3(i) = h3(i)+this%ao3(k,jj)*fo3(i,jjs)                              
          enddo
       enddo
       do i=1,ncol                                                  
          hco2(i,ks) = h2(i)                                               
          ho3(i,ks) = h3(i)*o3(i,ks)                                       
       enddo
    enddo
!                                                                         
!     **** COOLING RATE IN CO2 BANDS (X=10.75-12.5, MATRIX APPROACH)      
!                                                                         
    do k=36,43 
                                                  
       ks = k+8                                                        

       do i=1,ncol                                                 
          h2(i) = 0._wp                                                      
       enddo

       do jj=1,nrfmltelv
          jjs = ks+this%ig(jj)
          do i=1,ncol                                                 
             h2(i) = h2(i) + ( amat(i,k,jj) + bmat(i,k,jj)*fu(i,ks) ) * fu(i,jjs)
          enddo
       enddo

       do i=1,ncol
          hco2(i,ks) = h2(i)
          ho3(i,ks) = 0._wp
       enddo

    enddo

!     Define boundary condition at XR=12.5 for recurrence formula
    do k=1,nrfm
       if (this%xr(k).eq.12.5_wp) then
          do i=1,ncol
             h1(i)=hco2(i,k)/(co2(i,k)*(1._wp-alam(i,k))*constb)
          enddo
       endif
    enddo

!   Do the rest of the XR domain 
!   (transition region 12.75 <= p.s.h. <= 16.5)
    do k=1,nrfm

       if (this%xr(k).ge.12.75_wp .and. this%xr(k).le.16.5_wp) then

          do i=1,ncol

             h2(i) = ( aajm(i,k)*h1(i) + djm(i,k)*fu(i,k-1)    &
                      - dj0(i,k)*fu(i,k) ) / aaj0(i,k)

             hco2(i,k) = h2(i)*co2(i,k)*(1._wp-alam(i,k))/am(i,k)*const
             ho3(i,k)  = 0._wp

             h1(i)     = h2(i)

          enddo               ! next longitude
            
          if (this%xr(k).eq.16.5_wp) then

!     Calculate FLUX at the top of the transition region (XR=16.5)
             do  i=1,ncol
                flux(i) = h2(i) + fu(i,nrfm)
             enddo
               
          endif

       endif

    enddo                     ! next NLTE level

  end subroutine viccooln

!============================================================================================
  subroutine cool2space (ncol,nlyrs,t,mwair,vn2,vo2,vo,vco2, &
                         ndenair,xnorm,flux,hc2s)
!-----------------------------------------------------------------
!
!     Adapted from Ray Roble's model by F. Sassi (Nov. 1999)
!     
!     Performs cool-to-space cooling calculations
!     This mod operates on the same vertical grid of the GCM
!
!-----------------------------------------------------------------

    implicit none

    ! Input variables
    integer, intent(in) :: ncol                          ! number of atmospheric columns
    integer, intent(in) :: nlyrs                         ! number of model layers

    real(wp), intent(in) :: t(ncol,nlyrs)                ! neutral temperature
    real(wp), intent(in) :: vn2(ncol,nlyrs)              ! N2 vmr
    real(wp), intent(in) :: vo2(ncol,nlyrs)              ! O2 vmr
    real(wp), intent(in) :: vco2(ncol,nlyrs)             ! CO2 vmr
    real(wp), intent(in) :: vo(ncol,nlyrs)               ! O vmr
    real(wp), intent(in) :: ndenair(ncol,nlyrs)          ! mean air no. density
    real(wp), intent(in) :: mwair(ncol,nlyrs)            ! mean air molecular weight
    real(wp), intent(in) :: flux(ncol)                   ! Radiative flux at the top of the NLTE region
    real(wp), intent(in) :: xnorm(ncol,nlyrs)            ! p.s.h.

    ! Output  variables
    real(wp), intent(out) :: hc2s(ncol,nlyrs)            ! cool-to-space cooling

    ! Local variables
    real(wp) tt
    real(wp) y
    real(wp) zn2
    real(wp) zo2
    real(wp) zz
    real(wp) alam
    real(wp) rko

    integer i,k

!     ****                                        
!     ****   CO2 COOL-TO-SPACE APPROXIMATION      
!     **** 

    hc2s(1:ncol,1:nlyrs)=0.0_wp
    do k = 1,nlyrs

       do i=1,ncol

          if (xnorm(i,k) .gt. 16.5_wp) then

             tt     = t(i,k)                                        
             y      = tt**(-1._wp/3._wp)                                  
             zn2    = 5.5e-17_wp*sqrt(tt) + 6.7e-10_wp * exp(-83.8_wp*y)    
             zo2    = 1.0e-15_wp*exp(23.37_wp - 230.9_wp*y + 564._wp*y*y)
             rko=3.0e-12_wp 
!     
!     ****  COLLISIONAL DEACTIVATION RATE:                       
!
             zz     = (vn2(i,k)*zn2 + vo2(i,k)*zo2 + vo (i,k)*rko)*ndenair(i,k)
!                                                                
!     ****                                                       
!                                                                
             alam   = a10/(a10+zz)
             hc2s(i,k) = const/mwair(i,k)*vco2(i,k)*(1._wp-alam)      &
                         *(flux(i)-exp(-960.217_wp/tt))

          endif
         
       enddo                  ! end-loop in longitude

    enddo                     ! end-loop in levels

  end subroutine cool2space

!=====================================================================
  real(wp) function a18lin (x,xn,yn,m,n)
!-----------------------------------------------------------------
! input:
!  X - argument for which a value of function should be found
!  XN(N),YN(N) - values of function YN(N) at XN(N) grid. X(N) should be
!                ordered so that X(I-1) < X(I).
! output:
!  A18LIN - value of function for X
!-----------------------------------------------------------------

    implicit none

    ! Input:
    integer,intent(in) :: m,n
    real(wp),intent(in) :: x
    real(wp),intent(in) :: xn(n)
    real(wp),intent(in) :: yn(n)

    ! Local:
    integer :: k,i

    k=m-1
    the_loop: do i=m,n
       k=k+1
       if (x-xn(i).le.0._wp) exit the_loop
    enddo the_loop
    if (k.eq.1) k=2

! k has been found so that xn(k).le.x.lt.xn(k+1)

    a18lin=(yn(k)-yn(k-1))/(xn(k)-xn(k-1))*(x-xn(k))+yn(k)

  end function a18lin

!=============================================================================
  subroutine a18int(x1,y1,x2,y2,n1,n2)
!-----------------------------------------------------------------
!
! third order spline interpolation
! input argument and function:  X1(1:N1),Y1(1:N1)
! output argument and function: X2(1:N2)X2(1:N2),Y2(1:N2)
! the necessary conditionts are: X1(I) < X1(I+1), and the same for X2 array.
!
!-----------------------------------------------------------------

    implicit none
!
!   Input variables
    integer,  intent(in) :: n1
    integer,  intent(in) :: n2
    real(wp), intent(in) :: x1(n1)
    real(wp), intent(in) :: y1(n1)
    real(wp), intent(in) :: x2

!   Output variables
    real(wp), intent(out) :: y2
!
!   Local variables
    real(wp) :: a(150),e(150),f(150),h(150),h2,h1,f1,f2,f3,g
    integer :: nvs,k,kr,l
!
    h2=x1(1)
    nvs=n1-1
    do 1 k=1,nvs
      h1=h2
      h2=x1(k+1)
      h(k)=h2-h1
    1 continue
      a(1)=0._wp
      a(n1)=0._wp
      e(n1)=0._wp
      f(n1)=0._wp
      h1=h(n1-1)
      f1=y1(n1-1)
      f2=y1(n1)
      do 2 kr=2,nvs
      k=nvs+2-kr
      h2=h1
      h1=h(k-1)
      f3=f2
      f2=f1
      f1=y1(k-1)
      g=1._wp/(h2*e(k+1)+2._wp*(h1+h2))
      e(k)=-h1*g
      f(k)=(3._wp*((f3-f2)/h2-(f2-f1)/h1)-h2*f(k+1))*g
    2 continue
      g=0._wp
      do 3 k=2,nvs
      g=e(k)*g+f(k)
      a(k)=g
    3 continue
      l=1
      g=x2
      do 6 k=l,nvs
      if(g.gt.x1(k+1)) goto 6
      l=k
      goto 5
    6 continue
      l=nvs
    5 g=g-x1(l)
      h2=h(l)
      f2=y1(l)
      f1=h2**2
      f3=g**2
      y2=f2+g/h2*(y1(l+1)-f2-(a(l+1)*(f1-f3)+         &
                a(l)*(2._wp*f1-3._wp*g*h2+f3))/3._wp)
  end subroutine a18int

!==================================================================================================
  subroutine lininterp_full1d (arrin, yin, nin, arrout, yout, nout)
!-------------------------------------------------------------------------
!
!  Linear interpolation functions (from /cam/tools/interpolate_data.F90)
!
!-------------------------------------------------------------------------

    integer, intent(in) :: nin, nout
    real(wp), intent(in) :: arrin(nin), yin(nin), yout(nout)
    real(wp), intent(out) :: arrout(nout)
    type (interp_type) :: interp_wgts

    call lininterp_init(yin, nin, yout, nout, 1, interp_wgts)
    call lininterp1d(arrin, nin, arrout, nout, interp_wgts)
    call lininterp_finish(interp_wgts)

  end subroutine lininterp_full1d

!==================================================================================================
  subroutine lininterp_init(yin, nin, yout, nout, extrap_method, &
                            interp_wgts, cyclicmin, cyclicmax)
!-------------------------------------------------------------------------
!
! Description:
!   Initialize a variable of type(interp_type) with weights for linear interpolation.
!       this variable can then be used in calls to lininterp1d and lininterp2d.
!   yin is a 1d array of length nin of locations to interpolate from - this array must 
!       be monotonic but can be increasing or decreasing
!   yout is a 1d array of length nout of locations to interpolate to, this array need
!       not be ordered
!   extrap_method determines how to handle yout points beyond the bounds of yin
!       if 0 set values outside output grid to 0 
!       if 1 set to boundary value
!       if 2 set to cyclic boundaries
!         optional values cyclicmin and cyclicmax can be used to set the bounds of the 
!         cyclic mapping - these default to 0 and 360.
!
!-------------------------------------------------------------------------

    ! Input variables
    integer, intent(in) :: nin
    integer, intent(in) :: nout
    real(wp), intent(in) :: yin(nin)           ! input mesh
    real(wp), intent(in) :: yout(nout)         ! output mesh
    integer, intent(in) :: extrap_method       ! if 0 set values outside output grid to 0 
                                               ! if 1 set to boundary value
                                               ! if 2 set to cyclic boundaries
    real(wp), intent(in), optional :: cyclicmin, cyclicmax

    ! Output variables
    type (interp_type), intent(out) :: interp_wgts

    ! Local variables
    real(wp) :: cmin, cmax
    real(wp) :: extrap
    real(wp) :: dyinwrap
    real(wp) :: ratio
    real(wp) :: avgdyin
    integer :: i, j, icount
    integer :: jj
    real(wp), pointer :: wgts(:)
    real(wp), pointer :: wgtn(:)
    integer, pointer :: jjm(:)
    integer, pointer :: jjp(:)
    logical :: increasing
    !
    ! Check validity of input coordinate arrays: must be monotonically increasing,
    ! and have a total of at least 2 elements
    !
    if (nin.lt.2) then
       stop ('LININTERP: Must have at least 2 input points for interpolation')
    end if
    if(present(cyclicmin)) then
       cmin=cyclicmin
    else
       cmin=0_wp
    end if
    if(present(cyclicmax)) then
       cmax=cyclicmax
    else
       cmax=360_wp
    end if
    if(cmax<=cmin) then
       stop ('LININTERP: cyclic min value must be < max value')
    end if
    increasing=.true.
    icount = 0
    do j=1,nin-1
       if (yin(j).gt.yin(j+1)) icount = icount + 1
    end do
    if(icount.eq.nin-1) then
       increasing = .false.
       icount=0
    endif
    if (icount.gt.0) then
       stop ('LININTERP: Non-monotonic input coordinate array found')
    end if
    allocate(interp_wgts%jjm(nout), &
         interp_wgts%jjp(nout), &
         interp_wgts%wgts(nout), &
         interp_wgts%wgtn(nout))

    jjm => interp_wgts%jjm
    jjp => interp_wgts%jjp
    wgts =>  interp_wgts%wgts
    wgtn =>  interp_wgts%wgtn

    !
    ! Initialize index arrays for later checking
    !
    jjm = 0
    jjp = 0

    extrap = 0.
    if(extrap_method.eq.0) then
       !
       ! For values which extend beyond boundaries, set weights
       ! such that values will be 0.
       !
       do j=1,nout
          if(increasing) then
             if (yout(j).lt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 0.
                wgtn(j) = 0.
                extrap = extrap + 1.
             else if (yout(j).gt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 0.
                wgtn(j) = 0.
                extrap = extrap + 1.
             end if
          else
             if (yout(j).gt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 0.
                wgtn(j) = 0.
                extrap = extrap + 1.
             else if (yout(j).lt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 0.
                wgtn(j) = 0.
                extrap = extrap + 1.
             end if
          end if
       end do
    else if(extrap_method.eq.1) then
       !
       ! For values which extend beyond boundaries, set weights
       ! such that values will just be copied.
       !
       do j=1,nout
          if(increasing) then
             if (yout(j).le.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 1.
                wgtn(j) = 0.
                extrap = extrap + 1.
             else if (yout(j).gt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 1.
                wgtn(j) = 0.
                extrap = extrap + 1.
             end if
          else
             if (yout(j).gt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 1.
                wgtn(j) = 0.
                extrap = extrap + 1.
             else if (yout(j).le.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 1.
                wgtn(j) = 0.
                extrap = extrap + 1.
             end if
          end if
       end do
    else if(extrap_method.eq.2) then
       !
       ! For values which extend beyond boundaries, set weights
       ! for circular boundaries 
       !
       dyinwrap = yin(1) + (cmax-cmin) - yin(nin)
       avgdyin = abs(yin(nin)-yin(1))/(nin-1.)
       ratio = dyinwrap/avgdyin
       if (ratio < 0.9 .or. ratio > 1.1) then
          stop ('interpolate_data')
       end if

       do j=1,nout
          if(increasing) then
             if (yout(j) <= yin(1)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
             else if (yout(j) > yin(nin)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)-yin(nin))/dyinwrap
             end if
          else
             if (yout(j) > yin(1)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
             else if (yout(j) <= yin(nin)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin)-yin(nin))/dyinwrap
             end if

          endif
       end do
    end if

    !
    ! Loop though output indices finding input indices and weights
    !
    if(increasing) then
       do j=1,nout
          do jj=1,nin-1
             if (yout(j).gt.yin(jj) .and. yout(j).le.yin(jj+1)) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    else
       do j=1,nout
          do jj=1,nin-1
             if (yout(j).le.yin(jj) .and. yout(j).gt.yin(jj+1)) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    end if

    !
    ! Check that interp/extrap points have been found for all outputs
    !
    icount = 0
    do j=1,nout
       if (jjm(j).eq.0 .or. jjp(j).eq.0) icount = icount + 1
       ratio=wgts(j)+wgtn(j)
       if((ratio<0.9.or.ratio>1.1).and.extrap_method.ne.0) then
          stop ('Bad weight computed in LININTERP_init')
       end if
    end do
    if (icount.gt.0) then
       stop ('LININTERP: Point found without interp indices')
    end if

  end subroutine lininterp_init

!==================================================================================================
  subroutine lininterp1d (arrin, n1, arrout, m1, interp_wgts)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Do a linear interpolation from input mesh to output
    !          mesh with weights as set in lininterp_init.
    !
    !
    ! Author: Jim Edwards
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none

    !-----------------------------------------------------------------------
    !
    ! Input variables
    integer, intent(in) :: n1                 ! number of input latitudes
    integer, intent(in) :: m1                ! number of output latitudes

    real(wp), intent(in) :: arrin(n1)        ! input array of values to interpolate
    type(interp_type), intent(in) :: interp_wgts

    ! Output variables
    real(wp), intent(out) :: arrout(m1)      ! interpolated array

    !
    ! Local variables
    integer j                                ! latitude indices
    integer, pointer :: jjm(:)
    integer, pointer :: jjp(:)

    real(wp), pointer :: wgts(:)
    real(wp), pointer :: wgtn(:)

    jjm => interp_wgts%jjm
    jjp => interp_wgts%jjp
    wgts =>  interp_wgts%wgts
    wgtn =>  interp_wgts%wgtn

    !
    ! Do the interpolation
    !
    do j=1,m1
      arrout(j) = arrin(jjm(j))*wgts(j) + arrin(jjp(j))*wgtn(j)
    end do

  end subroutine lininterp1d

!==================================================================================================
  subroutine lininterp_finish(interp_wgts)

    type(interp_type) :: interp_wgts

    deallocate(interp_wgts%jjm, &
         interp_wgts%jjp, &
         interp_wgts%wgts, &
         interp_wgts%wgtn)

    nullify(interp_wgts%jjm, &
         interp_wgts%jjp, &
         interp_wgts%wgts, &
         interp_wgts%wgtn)

  end subroutine lininterp_finish

end module mo_upper_atmosphere_heating_lw
