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
!   in the near-IR shortwave due to CO2 following the approach of Ogibalov and Fomichev (2003).
!   Code adapted from the NCAR whole atmosphere code (WACCM) module mo_heatnirco2.F90 
!
! -------------------------------------------------------------------------------------------------
module mo_upper_atmosphere_heating_sw
  use mo_rte_kind,      only: wp, wl
  use mo_rte_config,    only: check_values, check_extents
  use mo_rte_util_array,only: any_vals_less_than, any_vals_outside, extents_are

  implicit none

  private
  save

  ! Public subroutines
  public :: heatnirco2_init, heatnirco2
  public :: airglow_init, airglow

  type, public :: ty_upper_atmosphere_heating_sw
    private

    integer :: ndpara, ncolgr

    ! MLT LW heating LUT coefficients
    real(wp), dimension(:),   allocatable :: xspara        ! pressure scale height
    real(wp), dimension(:),   allocatable :: zppara        !
    real(wp), dimension(:),   allocatable :: co2stand      ! CO2 colume mixing ratio
    real(wp), dimension(:,:), allocatable :: colmpara      !
    real(wp), dimension(:,:), allocatable :: corrnormpara  ! 
    !
  contains
    procedure, public :: load_upper_atmosphere_heating_sw
    procedure, public :: finalize_sw
  end type ty_upper_atmosphere_heating_sw
  !
  ! Fomichev near-IR heating parameters
  !
  integer, parameter :: ncol = 5
  integer, parameter :: nlyrs = 49
  integer, parameter :: nlvls = 50

  integer, parameter :: ndpara = 62
  integer, parameter :: ncolgr = 10

  real(wp) :: xspara(ndpara)
  real(wp) :: zppara(ndpara)
  real(wp) :: co2stand(ndpara)
  real(wp) :: colmpara(ndpara,ncolgr)
  real(wp) :: corrnormpara(ndpara,ncolgr)

  real(wp), parameter :: pi      = 2._wp * asin(1._wp)    ! pi
  real(wp), parameter :: d2r     = pi / 180._wp           ! Conversion for degrees to radians
  real(wp), parameter :: rearth  = 6.37122e6_wp           ! Radius of the earch (m)
  real(wp), parameter :: smallvalue = 1.0e-20_wp
  !
  ! Airglow heating parameters
  !
  integer , parameter :: nag      = 3                     ! number of airglow components
  real(wp), parameter :: avogad   = 6.02214e26_wp         ! Avogadro's constant (molecules/kmole)
  real(wp), parameter :: secpday  = 86400._wp             ! seconds per day
  real(wp), parameter :: daypsec  = 1._wp/secpday         ! days per second
  real(wp), parameter :: hc       = 6.62608e-34_wp*2.9979e8_wp/1.e-9_wp
  real(wp), parameter :: wc_o2_1s = 1._wp/762._wp         ! inverse of 762 nm
  real(wp), parameter :: wc_o2_1d = 1._wp/1270._wp        ! inverse of 1.27 microns
  real(wp), parameter :: wc_o1d   = 1._wp/630._wp         ! inverse of 630 nm

  integer :: rid_ag1, rid_ag2, rid_ag3, rid_rxn
  logical :: has_airglow

!================================================================================================
contains
!================================================================================================
  function load_upper_atmosphere_heating_sw(this, xspara, zppara, co2stand, colmpara, &
                               corrnormpara) result (error_msg)
    !
    class(ty_upper_atmosphere_heating_sw), intent(inout) :: this

    real(wp), dimension(:),   intent(in   ) :: xspara, zppara
    real(wp), dimension(:),   intent(in   ) :: co2stand
    real(wp), dimension(:,:), intent(in   ) :: colmpara
    real(wp), dimension(:,:), intent(in   ) :: corrnormpara
    character(len=128)  ::  error_msg
    !
    ! Local variables
    !
    integer :: ndpara, ncolgr
    !
    ! LUT coefficient dimensions
    !
    ndpara = size(colmpara,dim=1)
    ncolgr = size(colmpara,dim=2)

    error_msg = ""
    if(.not. extents_are(xspara, ndpara)) &
      error_msg = "upper_atmosphere_heating_sw%load_upper_atmosphere_heating_sw(): array xspara isn't consistently sized"

    if(.not. extents_are(zppara, ndpara)) &
      error_msg = "upper_atmosphere_heating_sw%load_upper_atmosphere_heating_sw(): array zppara isn't consistently sized"

    if(.not. extents_are(co2stand, ndpara)) &
      error_msg = "upper_atmosphere_heating_sw%load_upper_atmosphere_heating_sw(): array co2stand isn't consistently sized"

    if(.not. extents_are(corrnormpara, ndpara, ncolgr)) &
      error_msg = "upper_atmosphere_heating_sw%load_upper_atmosphere_heating_sw(): array corrnormpara isn't consistently sized"

    if (error_msg /= '') return

    this%ndpara = ndpara
    this%ncolgr = ncolgr

    allocate(this%xspara(ndpara))
    this%xspara = xspara

    allocate(this%zppara(ndpara))
    this%zppara = zppara

    allocate(this%co2stand(ndpara))
    this%co2stand = co2stand

    allocate(this%colmpara(ndpara, ncolgr), &
             this%corrnormpara(ndpara, ncolgr))
    this%colmpara = colmpara
    this%corrnormpara = corrnormpara
    
  end function load_upper_atmosphere_heating_sw

!================================================================================================
  subroutine finalize_sw(this)
  !-------------------------------------------------------------------
  !
  ! Finalize
  !
  !-------------------------------------------------------------------
    class(ty_upper_atmosphere_heating_sw), intent(inout) :: this

    deallocate(this%xspara, this%zppara, this%co2stand, &
               this%colmpara, this%corrnormpara)

  end subroutine finalize_sw

!================================================================================================
  subroutine heatnirco2_init (ncol, nlyrs, z_lay, sza, t_lay, p_lay, vmrco2, &
                              scco2)
!-----------------------------------------------------------------------
! 	... initialization for near-ir co2 heating rate
!-----------------------------------------------------------------------
!       Derive co2 slant path for near-ir co2 heating calculation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 	... input variables
!-----------------------------------------------------------------------
  integer, intent(in)  :: ncol                ! number of atmospheric columns
  integer, intent(in)  :: nlyrs               ! number of vertical layers
  real(wp), intent(in) :: z_lay(ncol,nlyrs)   ! altitude (km)
  real(wp), intent(in) :: sza(ncol)           ! solar zenith angle (degrees)
  real(wp), intent(in) :: t_lay(ncol,nlyrs)   ! model midpoint layer temperature (K)
  real(wp), intent(in) :: p_lay(ncol,nlyrs)   ! model midpoint layer pressure (Pa)
  real(wp), intent(in) :: vmrco2(ncol,nlyrs)  ! co2 concentration (mol/mol)

  real(wp), dimension(:,:), allocatable, intent(inout) :: scco2 ! co2 slant column (molec/cm^2)

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
  integer :: nid(ncol,0:nlyrs-1)           ! Number of layers crossed by the direct beam
                                         ! from the top of the atmosphere to layer 1
                                         ! at the bottom; nid(i), i=0,nlyrs-1
  real(wp) :: dsdh(ncol,0:nlyrs-1,nlyrs) ! Slant path of direct beam through each
                                         ! layer crossed from the top of the 
                                         ! atmosphere to layer 1 at the bottom;
                                         ! dsdh(i,j), i=0,nlyrs-1, j=1, nlyrs-1
  real(wp) :: invariants(ncol,nlyrs)     ! invariant density (mol/cm3)
  real(wp) :: co2cc(ncol,nlyrs)          ! co2 concentration (mol/mol)
  real(wp) :: delz(ncol,nlyrs)           ! distance between layer heights (cm)

  real(wp), parameter :: km2cm = 1.e5_wp ! Conversion factor from km to cm

! Derive slant path for spherical atmosphere
  call sphers( ncol, nlyrs, z_lay, sza, dsdh, nid )

! Derive invariants
  call setinv( ncol, nlyrs, t_lay, p_lay, invariants )

! Convert vmrco2 to co2cc
  co2cc(:,:) = vmrco2(:,:) * invariants(:,:)

! Derive layer thickness (cm)
  delz(:,1:nlyrs-1) = km2cm * (z_lay(:,2:nlyrs) - z_lay(:,1:nlyrs-1))

! Derive CO2 slant column
  allocate(scco2(ncol,nlyrs))
  call slant_col( ncol, nlyrs, delz, dsdh, nid, co2cc, scco2 )

  end subroutine heatnirco2_init

!================================================================================================
  subroutine heatnirco2(this, ncol, nlyrs, vmrco2, scco2, p_lay, &
                        htng_nir)
!-----------------------------------------------------------------------
! 	... input variables
!-----------------------------------------------------------------------
  class(ty_upper_atmosphere_heating_sw), intent(inout) :: this

  integer,  intent(in) :: ncol                     ! number of atmospheric columns
  integer,  intent(in) :: nlyrs                    ! number of vertical layers
  real(wp), intent(in) :: vmrco2(ncol,nlyrs)       ! co2 concentration (mol/mol)
  real(wp), intent(in) :: scco2(ncol,nlyrs)        ! co2 slant column (molec/cm^2)
  real(wp), intent(in) :: p_lay(ncol,nlyrs)        ! model midpoint pressure (Pa)

!-----------------------------------------------------------------------
!     	... output variables
!-----------------------------------------------------------------------
  real(wp), dimension(:,:), allocatable, intent(out) :: htng_nir ! co2 near ir heating (K/day)

!-----------------------------------------------------------------------
! 	... local variables and parameters
!-----------------------------------------------------------------------
  real(wp), parameter :: pa2hPa     = 1.e-2_wp

  integer  :: i, icolm, icolmp1
  integer  :: k, kk, kndx
  real(wp) :: reldcolm
  real(wp) :: delp
  real(wp) :: pinterp
  real(wp) :: colzpint
  real(wp) :: co2std
  real(wp) :: colparai(this%ncolgr)
  real(wp) :: corrnorai(this%ncolgr)

!-----------------------------------------------------------------------
! 	... initializAtion
!-----------------------------------------------------------------------
  allocate(htng_nir(ncol,nlyrs))
  htng_nir(:,:nlyrs) = smallvalue

!-----------------------------------------------------------------------
! 	... vertical and column interpolation
!-----------------------------------------------------------------------
column_loop : &
  do i = 1,ncol

level_loop : &
! ordered bottom to top
    do k = 1,nlyrs
!-----------------------------------------------------------------------
! 	... first setup pressure interpolation
!-----------------------------------------------------------------------
       pinterp  = p_lay(i,k) * pa2hPa
       colzpint = log( scco2(i,k) )
       if( pinterp <= this%xspara(this%ndpara) ) then
          colparai(:)  = this%colmpara(this%ndpara,:)
          corrnorai(:) = this%corrnormpara(this%ndpara,:)
          co2std       = this%co2stand(this%ndpara)
       else if( pinterp > this%xspara(1) ) then
          colparai(:)  = this%colmpara(1,:)
          corrnorai(:) = this%corrnormpara(1,:)
          co2std       = this%co2stand(1)
       else
          do kk = this%ndpara-1,1,-1
             if( pinterp <= this%xspara(kk) ) then
                kndx = kk + 1
                delp         = (pinterp - this%xspara(kndx))/(this%xspara(kk) - this%xspara(kndx))
                colparai(:)  = this%colmpara(kndx,:) + delp*(this%colmpara(kk,:) - this%colmpara(kndx,:))
                corrnorai(:) = this%corrnormpara(kndx,:) &
                               + delp*(this%corrnormpara(kk,:) - this%corrnormpara(kndx,:))
                co2std       = this%co2stand(kndx) + delp*(this%co2stand(kk) - this%co2stand(kndx))
                exit
             end if
          end do
       end if

!-----------------------------------------------------------------------
! Linear interpolation over column density for given altitude point
!-----------------------------------------------------------------------
       if( colzpint < colparai(1) ) then
          htng_nir(i,k) = corrnorai(1)
       else if( colzpint >= colparai(this%ncolgr) ) then
          htng_nir(i,k) = corrnorai(this%ncolgr)
       else
loop1:   do icolm = 1,this%ncolgr-1
           icolmp1 = icolm + 1
           if( colzpint >= colparai(icolm) .and. &
               colzpint <  colparai(icolmp1) ) then
             reldcolm = (colzpint - colparai(icolm)) &
                       /(colparai(icolmp1) - colparai(icolm))
             htng_nir(i,k) = corrnorai(icolm) &
                       + (corrnorai(icolmp1) - corrnorai(icolm))*reldcolm
             exit loop1
           end if
         end do loop1
       end if
         
!-----------------------------------------------------------------------
! From normalized value to the one corresponding to the given vmrco2
!-----------------------------------------------------------------------

       htng_nir(i,k) = htng_nir(i,k) * vmrco2(i,k)/co2std

    end do level_loop

  end do column_loop

  end subroutine heatnirco2

!================================================================================================
  subroutine sphers( ncol, nlyrs, z_lay, zenith_angle, dsdh, nid )
!=============================================================================!
!   Subroutine sphers                                                         !
!=============================================================================!
!   PURPOSE:                                                                  !
!   Calculate slant path over vertical depth ds/dh in spherical geometry.     !
!   Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model  !
!   for computing the radiation field available for photolysis and heating    !
!   at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)   !
!=============================================================================!
!   PARAMETERS:                                                               !
!   NZ      - INTEGER, number of specified altitude levels in the working (I) !
!             grid                                                            !
!   Z       - REAL, specified altitude working grid (km)                  (I) !
!   ZEN     - REAL, solar zenith angle (degrees)                          (I) !
!   DSDH    - REAL, slant path of direct beam through each layer crossed  (O) !
!             when travelling from the top of the atmosphere to layer i;      !
!             DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                             !
!   NID     - INTEGER, number of layers crossed by the direct beam when   (O) !
!             travelling from the top of the atmosphere to layer i;           !
!             NID(i), i = 0..NZ-1                                             !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   Original: Taken By Doug Kinnison from Sasha Madronich, TUV Code, V4.1a,   !
!             on 1/1/02                                                       !
!=============================================================================!

  implicit none

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
  integer, intent(in)  :: ncol,nlyrs                   ! number model vertical levels
  integer, intent(out) :: nid(ncol,0:nlyrs-1)          ! see above

  real(wp), intent (in) :: z_lay(ncol,nlyrs)           ! geometric altitude (km)
  real(wp), intent (in) :: zenith_angle(ncol)          ! zenith_angle (degrees)
  real(wp), intent (out) :: dsdh(ncol,0:nlyrs-1,nlyrs) ! see above

!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
  real(wp) :: radius              ! radius of earth (km)
  real(wp) :: re(ncol)            ! radius of earth plus bottom level height (km)
  real(wp) :: zenrad              ! solar zenith angle (radians)
  real(wp) :: rpsinz
  real(wp) :: const0(ncol)        ! sin of zenrad
  real(wp) :: rj
  real(wp) :: rjp1
  real(wp) :: dsj
  real(wp) :: dhj
  real(wp) :: ga
  real(wp) :: gb
  real(wp) :: sm
  real(wp) :: zd(ncol,0:nlyrs-1)

  integer :: i, j, k
  integer :: id
  integer :: nlayer

!------------------------------------------------------------------------------
!       ... set radius of earth (km)
!------------------------------------------------------------------------------
  radius = rearth*1.e-3_wp   

!------------------------------------------------------------------------------
!       ... set number of layers:
!------------------------------------------------------------------------------
  nlayer = nlyrs - 1

!------------------------------------------------------------------------------
!       ... set zenith angle in radians
!------------------------------------------------------------------------------
  do i=1,ncol
    zenrad = zenith_angle(i) * d2r
    const0(i) = sin( zenrad )

!------------------------------------------------------------------------------
!       ... include the elevation above sea level to the radius of the earth:
!------------------------------------------------------------------------------
    re(i) = radius + z_lay(i,1)

!------------------------------------------------------------------------------
!       ... inverse coordinate of z
!------------------------------------------------------------------------------
    do k = nlayer,0,-1
      zd(i,k) = z_lay(i,k+1) - z_lay(i,1)
    end do
  end do

!------------------------------------------------------------------------------
!       ... initialize dsdh(i,j), nid(i)
!------------------------------------------------------------------------------
  nid(:,:) = 0
  do j = 1,nlyrs
    do i = 1,ncol
      dsdh(i,0:nlyrs-1,j) = 0._wp
    end do
  end do

!------------------------------------------------------------------------------
!       ... calculate ds/dh of every layer
!------------------------------------------------------------------------------
  do i = 1,ncol
    do k = nlayer,1,-1
      rpsinz = (re(i) + zd(i,k)) * const0(i)
      if( zenith_angle(i) <= 90._wp .or. rpsinz >= re(i) ) then
!------------------------------------------------------------------------------
! Find index of layer in which the screening height lies
!------------------------------------------------------------------------------
         id = k 
         if( zenith_angle(i) > 90._wp ) then
            do j = nlayer,1,-1
               if( rpsinz < (zd(i,j) + re(i)) .and.  rpsinz >= (zd(i,j-1) + re(i)) ) then
                  id = j
                  exit
               end if
            end do
         end if

         do j = nlayer,id,-1
           sm = 1._wp
           if( j == id .and. id == k .and. zenith_angle(i) > 90._wp ) then
              sm = -1._wp
           end if
           rj   = re(i) + zd(i,j-1)
           rjp1 = re(i) + zd(i,j)
           dhj  = zd(i,j) - zd(i,j-1)
           ga   = max( rj*rj - rpsinz*rpsinz,0._wp )
           gb   = max( rjp1*rjp1 - rpsinz*rpsinz,0._wp )
           if( id > k .and. j == id ) then
              dsj = sqrt( gb )
           else
              dsj = sqrt( gb ) - sm*sqrt( ga )
           end if
           dsdh(i,k,j) = dsj / dhj
         end do
         nid(i,k) = id
      else
         nid(i,k) = -1
      end if
    end do
  end do

  end subroutine sphers

!================================================================================================
  subroutine slant_col( ncol, nlyrs, delz, dsdh, nid, absden, scol )
!=============================================================================!
!   PURPOSE:                                                                  !
!   Derive Column
!=============================================================================!
!   PARAMETERS:                                                               !
!   NLYRS   - INTEGER, number of specified altitude levels in the working  (I) !
!            grid                                                             !
!   DELZ   - REAL, specified altitude working grid (km)                   (I) !
!   DSDH   - REAL, slant path of direct beam through each layer crossed  (O)  !
!             when travelling from the top of the atmosphere to layer i;      !
!             DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                             !
!   NID    - INTEGER, number of layers crossed by the direct beam when   (O)  !
!             travelling from the top of the atmosphere to layer i;           !
!             NID(i), i = 0..NZ-1                                             !
!            specified altitude at each specified wavelength                  !
!   absden - REAL, absorber concentration, molecules cm-3                     !
!   SCOL   - REAL, absorber Slant Column, molecules cm-2                      !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   09/01  Read in profile from an input file, DEK                            !
!   01/02  Taken from Sasha Madronich's TUV code                              !
!=============================================================================!

  implicit none

!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
  integer,  intent(in)  :: ncol                       ! number of columns
  integer,  intent(in)  :: nlyrs                      ! number oflayers
  integer,  intent(in)  :: nid(ncol,0:nlyrs-1)        ! see above
  real(wp), intent(in)  :: delz(ncol,nlyrs)	      ! layer thickness (cm)
  real(wp), intent(in)  :: dsdh(ncol,0:nlyrs-1,nlyrs) ! see above
  real(wp), intent(in)  :: absden(ncol,nlyrs)         ! absorber concentration (molec. cm-3)
  real(wp), intent(out) :: scol(ncol,nlyrs)	      ! absorber Slant Column (molec. cm-2)

!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
  real(wp), parameter :: largest = 1.e+36_wp

  real(wp) :: sum
  real(wp) :: hscale
  real(wp) :: numer, denom
  real(wp) :: cz(ncol,nlyrs)
  real(wp) :: scol_tmp(ncol,nlyrs)

  integer :: i, id, j, k

!------------------------------------------------------------------------------
!     ... compute column increments (logarithmic integrals)
!------------------------------------------------------------------------------
! absden and delz are ordered bottom-up

  do i = 1,ncol
    do k = nlyrs-1,1,-1
       if( absden(i,k) /= 0._wp .and. absden(i,k+1) /= 0._wp ) then
           cz(i,nlyrs-k) = (absden(i,k) - absden(i,k+1))/log( absden(i,k)/absden(i,k+1) ) * delz(i,k)
       else
           cz(i,nlyrs-k) = .5_wp*(absden(i,k) + absden(i,k+1)) * delz(i,k)
       end if
    end do
  end do

!------------------------------------------------------------------------------
!     ... Include exponential tail integral from infinity to model top
!         specify scale height near top of data.For WACCM-X model, scale
!         height needs to be increased for higher model top
!------------------------------------------------------------------------------
  hscale = 10.e5_wp

! cz is ordered top-down
  do i = 1,ncol
     cz(i,nlyrs-1) = cz(i,nlyrs-1) + hscale * absden(i,1)

!------------------------------------------------------------------------------
!       ...  Calculate vertical and slant column from each level:
!            work downward
!------------------------------------------------------------------------------
     do id = 0,nlyrs-1
       sum = 0._wp
       if( nid(i,id) >= 0 ) then
!------------------------------------------------------------------------------
!       ...  Single pass layers:
!------------------------------------------------------------------------------
         do j = 1, min(nid(i,id), id)
            sum = sum + cz(i,nlyrs-j)*dsdh(i,id,j)
         end do
!------------------------------------------------------------------------------
!       ...  Double pass layers:
!------------------------------------------------------------------------------
         do j = min(nid(i,id),id)+1, nid(i,id)
            sum = sum + 2._wp*cz(i,nlyrs-j)*dsdh(i,id,j)
         end do
       else
         sum = largest
       end if
       scol_tmp(i,nlyrs-id) = sum
    end do
! scol is ordered top-down
    scol_tmp(i,nlyrs) = 0.95_wp*scol_tmp(i,nlyrs-1)
! Reorder scol to be bottom-up
    do k=1,nlyrs
       scol(i,nlyrs+1-k) = scol_tmp(i,k)
    enddo

  end do

  end subroutine slant_col

!================================================================================================
  subroutine setinv( ncol, nlyrs, t_lay, p_lay, invariants )
  !-----------------------------------------------------------------
  !        ... set the invariant densities (molecules/cm**3)
  !-----------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------
  !        ... dummy arguments
  !-----------------------------------------------------------------
  integer,  intent(in)  :: ncol                   ! number of columns
  integer,  intent(in)  :: nlyrs                  ! number of layers
  real(wp), intent(in)  :: t_lay(ncol,nlyrs)      ! temperature
  real(wp), intent(in)  :: p_lay(ncol,nlyrs)      ! pressure (Pa)
  real(wp), intent(out) :: invariants(ncol,nlyrs) ! invariant density array (mol/cm3)

  !-----------------------------------------------------------------
  !        .. local variables
  !-----------------------------------------------------------------
  integer :: i, k

  real(wp), parameter :: Pa_xfac = 10._wp            ! Pascals to dyne/cm^2
  real(wp), parameter :: boltz = 1.38065e-23_wp      ! Boltzmann's constant (J/K/molecule)
  real(wp), parameter :: boltz_cgs = boltz * 1.e7_wp ! erg/K

  !-----------------------------------------------------------------
  !        note: invariants are in cgs density units.
  !              the p_lay array is in pascals and must be
  !	       mutiplied by 10. to yield dynes/cm**2.
  !-----------------------------------------------------------------
  invariants(:,:) = 0._wp

  do i = 1,ncol
    do k = 1,nlyrs
       invariants(i,k) = Pa_xfac * p_lay(i,k) / (boltz_cgs*t_lay(i,k))
    end do
  end do

  end subroutine setinv

!================================================================================================
  subroutine airglow_init(ncol, nlyrs, nrxn, vmr_o2_1d, vmr_o2_1s, vmr_o_1d, reaction_rates)
!-----------------------------------------------------------------------
! 	... initialization for near-ir co2 heating rate
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     	... input variables
!-----------------------------------------------------------------------
  integer,  intent(in) :: ncol                            ! number of columns
  integer,  intent(in) :: nlyrs                           ! number of layers
  integer,  intent(in) :: nrxn                            ! number of airglow reaction components
  real(wp), intent(in) :: vmr_o2_1d(ncol,nlyrs)           ! o2_1d vmr airglow concentrations
  real(wp), intent(in) :: vmr_o2_1s(ncol,nlyrs)           ! o2_1s vmr airglow concentrations
  real(wp), intent(in) :: vmr_o_1d(ncol,nlyrs)            ! o_1d vmr airglow concentrations
  real(wp), intent(in) :: reaction_rates(ncol,nlyrs,nrxn) ! airglow reaction rates (1/cm^3/s)

  has_airglow = .false.

  rid_ag1 = 0
  rid_ag2 = 0
  rid_ag3 = 0
  rid_rxn = 0

  if (maxval(vmr_o2_1d(:,:)) > 0.0_wp) rid_ag1 = 1
  if (maxval(vmr_o2_1s(:,:)) > 0.0_wp) rid_ag2 = 1
  if (maxval(vmr_o_1d(:,:)) > 0.0_wp) rid_ag3 = 1
  if (maxval(reaction_rates) > 0.0_wp) rid_rxn = 1

  has_airglow = rid_ag1 > 0 .and. rid_ag2 > 0 .and. rid_ag3 > 0 .and. rid_rxn > 0

  end subroutine airglow_init

!================================================================================================
  subroutine airglow(ncol, nlyrs, nrxn, o2_1d, o2_1s, o1d, rxt, cp_dry, &
                     aghr)
!-----------------------------------------------------------------------
!      	... derive the airglow heating rates
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     	... input variables
!-----------------------------------------------------------------------

  integer,  intent(in) ::  ncol                  ! number of columns
  integer,  intent(in) ::  nlyrs                 ! number of layers
  integer,  intent(in) ::  nrxn                  ! number of airglow reaction components
  real(wp), intent(in) ::  rxt(ncol,nlyrs,nrxn)  ! rxt rates (1/cm^3/s)
  real(wp), intent(in) ::  o2_1d(ncol,nlyrs)     ! o2_1d concentration (mol/mol)
  real(wp), intent(in) ::  o2_1s(ncol,nlyrs)     ! o2_1s concentration (mol/mol)
  real(wp), intent(in) ::  o1d(ncol,nlyrs)       ! o1d concentration (mol/mol)
  real(wp), intent(in) ::  cp_dry                ! specific heat capacity

!-----------------------------------------------------------------------
!     	... output variables
!-----------------------------------------------------------------------
  real(wp), dimension(:,:), allocatable, intent(out) :: aghr(:,:) ! total airglow heating rate (K/s)

!-----------------------------------------------------------------------
!     	... local variables
!-----------------------------------------------------------------------
  integer  :: k
  real(wp) :: tmp(ncol)
  real(wp) :: ag_rate(ncol,nlyrs,nag)

!-----------------------------------------------------------------------
! 	... initializAtion
!-----------------------------------------------------------------------
  allocate(aghr(ncol,nlyrs))
  aghr(:,:nlyrs) = smallvalue

  if (.not. has_airglow) return

  do k = 1,nlyrs
     tmp(:)          = hc * avogad / cp_dry
     ag_rate(:,k,1)  = tmp(:)*rxt(:,k,1)*o2_1d(:,k)*wc_o2_1d       ! o2_1d -> o2 heating at 1.27 microns
     ag_rate(:,k,2)  = tmp(:)*rxt(:,k,2)*o2_1s(:,k)*wc_o2_1s       ! o2_1s -> o2 heating at 762 nm
     ag_rate(:,k,3)  = tmp(:)*rxt(:,k,3)*o1d(:,k)*wc_o1d           ! o1d -> o heating at 630 nm
     aghr(:,k)       = ag_rate(:,k,1) + ag_rate(:,k,2) + ag_rate(:,k,3) ! total airglow heating
  enddo

  end subroutine airglow

  end module mo_upper_atmosphere_heating_sw
