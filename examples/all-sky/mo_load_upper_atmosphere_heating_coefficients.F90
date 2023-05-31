
module mo_load_upper_atmosphere_heating_coeffs

  use mo_rte_kind,       only: wp
  use mo_upper_atmosphere_heating_lw, only: ty_upper_atmosphere_heating_lw
  use mo_upper_atmosphere_heating_sw, only: ty_upper_atmosphere_heating_sw
  use mo_simple_netcdf,  only: read_field, read_string, var_exists, get_dim_size, &
                               write_field, create_dim, create_var
  use netcdf

  implicit none

  private
  public :: load_upper_atmosphere_heating_lw_coeff, load_upper_atmosphere_heating_sw_coeff
  public :: read_upper_atmosphere_lw_state, read_upper_atmosphere_sw_state
  public :: write_upper_atmosphere_lw_heating, write_upper_atmosphere_sw_heating
  public :: is_lw, is_sw
  ! ----------------------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------------------
  !
  ! read MLT upper atmosphere LW input coefficients from NetCDF file
  !
  subroutine load_upper_atmosphere_heating_lw_coeff(upper_atmosphere_heating_lw_spec, &
                                                    upper_atmosphere_heating_lw_coeff_file)
    class(ty_upper_atmosphere_heating_lw),   intent(inout) :: upper_atmosphere_heating_lw_spec
    character(len=*),           intent(in   ) :: upper_atmosphere_heating_lw_coeff_file
    ! -----------------
    ! Local variables
    integer :: ncid
    integer :: nprs67, nprs51, nprs43, nprs35
    integer :: nlev9, nlev6, nco2amt

    ! LUT coefficients
    real(wp), dimension(:),   allocatable :: xr           ! pressure scale height
    real(wp), dimension(:),   allocatable :: ig           ! heat exchange level index
    real(wp), dimension(:),   allocatable :: co2o         ! co2 volume mixing ratios
    real(wp), dimension(:,:), allocatable :: ao3          ! o3 cooling rate coefficients
    real(wp), dimension(:,:), allocatable :: a150         ! co2 cooling rate a-coefficients at 150 ppmv 
    real(wp), dimension(:,:), allocatable :: b150         ! co2 cooling rate b-coefficients at 150 ppmv 
    real(wp), dimension(:,:), allocatable :: a360         ! co2 cooling rate a-coefficients at 360 ppmv 
    real(wp), dimension(:,:), allocatable :: b360         ! co2 cooling rate b-coefficients at 360 ppmv 
    real(wp), dimension(:,:), allocatable :: a540         ! co2 cooling rate a-coefficients at 540 ppmv 
    real(wp), dimension(:,:), allocatable :: b540         ! co2 cooling rate b-coefficients at 540 ppmv 
    real(wp), dimension(:,:), allocatable :: a720         ! co2 cooling rate a-coefficients at 720 ppmv 
    real(wp), dimension(:,:), allocatable :: b720         ! co2 cooling rate b-coefficients at 720 ppmv 
    real(wp), dimension(:),   allocatable :: uco2co       ! co2 column amounts for escape function correction
    real(wp), dimension(:),   allocatable :: uco2ro       ! co2 column amounts 
    real(wp), dimension(:),   allocatable :: alo          ! co2 escape functions
    real(wp), dimension(:),   allocatable :: cor150       ! corrections for escape function at 150 ppmv
    real(wp), dimension(:),   allocatable :: cor360       ! corrections for escape function at 360 ppmv
    real(wp), dimension(:),   allocatable :: cor540       ! corrections for escape function at 540 ppmv
    real(wp), dimension(:),   allocatable :: cor720       ! corrections for escape function at 720 ppmv

    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(upper_atmosphere_heating_lw_coeff_file), NF90_WRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_upper_atmosphere_heating_lw_coeff(): can't open file " // &
                         trim(upper_atmosphere_heating_lw_coeff_file))

    ! Read LUT coefficient dimensions
    nprs67  = get_dim_size(ncid,'nprs67')
    nprs51  = get_dim_size(ncid,'nprs51')
    nprs43  = get_dim_size(ncid,'nprs43')
    nprs35  = get_dim_size(ncid,'nprs35')
    nlev9   = get_dim_size(ncid,'nlev9')
    nlev6   = get_dim_size(ncid,'nlev6')
    nco2amt = get_dim_size(ncid,'nco2amt')

    allocate(xr(nprs67))
    xr = read_field(ncid, 'xr', nprs67)

    allocate(ig(nlev9))
    ig = read_field(ncid, 'ig', nlev9)

    allocate(co2o(nco2amt))
    co2o = read_field(ncid, 'co2o', nco2amt)

    ! Allocate array for ozone cooling coefficients
    allocate(ao3(nprs35, nlev9))
    ao3 = read_field(ncid, 'ao3', nprs35, nlev9)

    ! Allocate and read arrays for LTE CO2 cooling coefficients
    allocate(a150(nprs43, nlev9), &
             b150(nprs43, nlev9), &
             a360(nprs43, nlev9), &
             b360(nprs43, nlev9), &
             a540(nprs43, nlev9), &
             b540(nprs43, nlev9), &
             a720(nprs43, nlev9), &
             b720(nprs43, nlev9))
    a150 = read_field(ncid, 'a150', nprs43, nlev9)
    b150 = read_field(ncid, 'b150', nprs43, nlev9)
    a360 = read_field(ncid, 'a360', nprs43, nlev9)
    b360 = read_field(ncid, 'b360', nprs43, nlev9)
    a540 = read_field(ncid, 'a540', nprs43, nlev9)
    b540 = read_field(ncid, 'b540', nprs43, nlev9)
    a720 = read_field(ncid, 'a720', nprs43, nlev9)
    b720 = read_field(ncid, 'b720', nprs43, nlev9)

    ! Allocate and read arrays for CO2 escape functions (alo) and CO2 column amounts (uco2ro)
    allocate(uco2ro(nprs51), &
             alo(nprs51))
    uco2ro = read_field(ncid, 'uco2ro', nprs51)
    alo    = read_field(ncid, 'alo', nprs51)

    ! Allocate and read arrays for CO2 column amounts for CO2 escape function correction tables
    allocate(uco2co(nlev6))
    uco2co = read_field(ncid, 'uco2co', nlev6)

    ! Allocate and read arrays for CO2 escape function corrections at four CO2 concentrations
    allocate(cor150(nlev6), &
             cor360(nlev6), &
             cor540(nlev6), &
             cor720(nlev6))
    cor150 = read_field(ncid, 'cor150', nlev6)
    cor360 = read_field(ncid, 'cor360', nlev6)
    cor540 = read_field(ncid, 'cor540', nlev6)
    cor720 = read_field(ncid, 'cor720', nlev6)

    ncid = nf90_close(ncid)

    call stop_on_err(upper_atmosphere_heating_lw_spec%load_upper_atmosphere_heating_lw( &
                         xr, ig, co2o, ao3, a150, b150, a360, b360, a540, b540, a720, b720, &
                         uco2ro, alo, uco2co, cor150, cor360, cor540, cor720))

  end subroutine load_upper_atmosphere_heating_lw_coeff
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! read MLT upper atmosphere SW input coefficients from NetCDF file
  !
  subroutine load_upper_atmosphere_heating_sw_coeff(upper_atmosphere_heating_sw_spec, &
                                                    upper_atmosphere_heating_sw_coeff_file)
    class(ty_upper_atmosphere_heating_sw),   intent(inout) :: upper_atmosphere_heating_sw_spec
    character(len=*),           intent(in   ) :: upper_atmosphere_heating_sw_coeff_file
    ! -----------------
    ! Local variables
    integer :: ncid
    integer :: ndpara, ncolgr

    ! LUT coefficients
    real(wp), dimension(:),   allocatable :: xspara       ! pressure scale height
    real(wp), dimension(:),   allocatable :: zppara       ! TBD
    real(wp), dimension(:),   allocatable :: co2stand     ! co2 volume mixing ratios
    real(wp), dimension(:,:), allocatable :: colmpara     ! TBD
    real(wp), dimension(:,:), allocatable :: corrnormpara ! TBD

    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(upper_atmosphere_heating_sw_coeff_file), NF90_WRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_upper_atmosphere_heating_sw_coeff(): can't open file " // &
                         trim(upper_atmosphere_heating_sw_coeff_file))

    ! Read LUT coefficient dimensions
    ndpara  = get_dim_size(ncid,'ndpara')
    ncolgr  = get_dim_size(ncid,'ncolgr')

    allocate(xspara(ndpara))
    xspara = read_field(ncid, 'xspara', ndpara)

    allocate(zppara(ndpara))
    zppara = read_field(ncid, 'zppara', ndpara)

    allocate(co2stand(ndpara))
    co2stand = read_field(ncid, 'co2stand', ndpara)

    allocate(colmpara(ndpara, ncolgr))
    colmpara = read_field(ncid, 'colmpara', ndpara, ncolgr)

    allocate(corrnormpara(ndpara, ncolgr))
    corrnormpara = read_field(ncid, 'corrnormpara', ndpara, ncolgr)

    ncid = nf90_close(ncid)

    call stop_on_err(upper_atmosphere_heating_sw_spec%load_upper_atmosphere_heating_sw(xspara, &
                         zppara, co2stand, colmpara, corrnormpara))

  end subroutine load_upper_atmosphere_heating_sw_coeff
  !---------------------------------------------------------------------------------------
  !
  ! Read atmospheric state variables for upper atmosphere calculations
  !
  subroutine read_upper_atmosphere_lw_state(fileName, p_lay, p_lev, t_lay, &
                                            vmr_co2, vmr_o2, vmr_o3, vmr_o, vmr_n2)

    character(len=*),                      intent(in   ) :: fileName
    real(wp), dimension(:,:), allocatable, intent(inout) :: p_lay
    real(wp), dimension(:,:), allocatable, intent(inout) :: p_lev
    real(wp), dimension(:,:), allocatable, intent(inout) :: t_lay
    real(wp), dimension(:,:), allocatable, intent(inout) :: vmr_co2
    real(wp), dimension(:,:), allocatable, intent(inout) :: vmr_o2
    real(wp), dimension(:,:), allocatable, intent(inout) :: vmr_o3
    real(wp), dimension(:,:), allocatable, intent(inout) :: vmr_o
    real(wp), dimension(:,:), allocatable, intent(inout) :: vmr_n2

    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, nlev

    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_upper_atmosphere_lw_state: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    nlev = get_dim_size(ncid, 'lev')

    allocate(p_lay(ncol,nlay))
    allocate(p_lev(ncol,nlev))
    allocate(t_lay(ncol,nlay))
    allocate(vmr_co2(ncol,nlay))
    allocate(vmr_o2(ncol,nlay))
    allocate(vmr_o3(ncol,nlay))
    allocate(vmr_o(ncol,nlay))
    allocate(vmr_n2(ncol,nlay))

    p_lay   =  read_field(ncid,   'p_lay', ncol, nlay)
    p_lev   =  read_field(ncid,   'p_lev', ncol, nlev)
    t_lay   =  read_field(ncid,   't_lay', ncol, nlay)
    vmr_co2 =  read_field(ncid, 'vmr_co2', ncol, nlay)
    vmr_o2  =  read_field(ncid,  'vmr_o2', ncol, nlay)
    vmr_o3  =  read_field(ncid,  'vmr_o3', ncol, nlay)
    vmr_o   =  read_field(ncid,   'vmr_o', ncol, nlay)
    vmr_n2  =  read_field(ncid,  'vmr_n2', ncol, nlay)

    ncid = nf90_close(ncid)

  end subroutine read_upper_atmosphere_lw_state
  ! --------------------------------------------------------------------------------------
  !
  ! Read atmospheric state variables for upper atmosphere calculations
  !
  subroutine read_upper_atmosphere_sw_state(fileName, p_lay, p_lev, t_lay, &
                                            vmr_h2o, vmr_co2, sza)

    character(len=*),                      intent(in   ) :: fileName
    real(wp), dimension(:,:), allocatable, intent(inout) :: p_lay
    real(wp), dimension(:,:), allocatable, intent(inout) :: p_lev
    real(wp), dimension(:,:), allocatable, intent(inout) :: t_lay
    real(wp), dimension(:,:), allocatable, intent(inout) :: vmr_h2o
    real(wp), dimension(:,:), allocatable, intent(inout) :: vmr_co2
    real(wp), dimension(:),   allocatable, intent(inout) :: sza

    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, nlev

    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_upper_atmosphere_sw_state: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    nlev = get_dim_size(ncid, 'lev')

    allocate(p_lay(ncol,nlay))
    allocate(p_lev(ncol,nlev))
    allocate(t_lay(ncol,nlay))
    allocate(vmr_h2o(ncol,nlay))
    allocate(vmr_co2(ncol,nlay))
    allocate(sza(ncol))

    p_lay   =  read_field(ncid,              'p_lay', ncol, nlay)
    p_lev   =  read_field(ncid,              'p_lev', ncol, nlev)
    t_lay   =  read_field(ncid,              't_lay', ncol, nlay)
    vmr_h2o =  read_field(ncid,            'vmr_h2o', ncol, nlay)
    vmr_co2 =  read_field(ncid,            'vmr_co2', ncol, nlay)
    sza     =  read_field(ncid, 'solar_zenith_angle', ncol)

    ncid = nf90_close(ncid)

  end subroutine read_upper_atmosphere_sw_state
  !---------------------------------------------------------------------------------------
  !
  ! Write upper atmosphere LW heating rate
  !
  subroutine write_upper_atmosphere_lw_heating(fileName, upper_atmosphere_heating_lw_merge)

    character(len=*),                      intent(in   ) :: fileName
    real(wp), dimension(:,:), allocatable, intent(inout) :: upper_atmosphere_heating_lw_merge

    ! -------------------
    integer :: ncid
    integer :: ncol, nlay

    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_upper_atmosphere_lw_heating: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')

    call create_var(ncid,  "heating_rate_upper_atmosphere", ["col", "lay"], [ncol, nlay]) 
    call stop_on_err(write_field(ncid,  "heating_rate_upper_atmosphere", upper_atmosphere_heating_lw_merge))

    ncid = nf90_close(ncid)

  end subroutine write_upper_atmosphere_lw_heating
  !---------------------------------------------------------------------------------------
  !
  ! Write upper atmosphere SW heating rate
  !
  subroutine write_upper_atmosphere_sw_heating(fileName, upper_atmosphere_heating_sw_merge)

    character(len=*),                      intent(in   ) :: fileName
    real(wp), dimension(:,:), allocatable, intent(inout) :: upper_atmosphere_heating_sw_merge

    ! -------------------
    integer :: ncid
    integer :: ncol, nlay

    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_upper_atmosphere_sw_heating: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')

    call create_var(ncid,  "heating_rate_upper_atmosphere", ["col", "lay"], [ncol, nlay]) 
    call stop_on_err(write_field(ncid,  "heating_rate_upper_atmosphere", upper_atmosphere_heating_sw_merge))

    ncid = nf90_close(ncid)

  end subroutine write_upper_atmosphere_sw_heating
  !------------------------------------------------------------------------------------------------------
  !
  ! Does this file contain variables needed to do SW calculations ?
  !
  function is_sw(fileName)
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_sw

    integer :: ncid

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("is_sw: can't find file " // trim(fileName))

    is_sw = var_exists(ncid, 'solar_zenith_angle')
    ncid = nf90_close(ncid)
  end function is_sw

  !------------------------------------------------------------------------------------------------------
  function is_lw(fileName)
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_lw

    is_lw = .not. is_sw(fileName)
  end function is_lw

  ! -----------------------------------------------------------------------------------
    subroutine stop_on_err(msg)
      !
      ! Print error message and stop
      !
      use iso_fortran_env, only : error_unit
      character(len=*), intent(in) :: msg
      if(len_trim(msg) > 0) then
        write (error_unit,*) trim(msg)
        stop
      end if
    end subroutine

end module
