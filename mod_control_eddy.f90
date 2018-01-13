module control_eddy
  use params
  use zeta, only : uc, vc, wc, ssh, advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy, &
                   curlnonl, betav, stretchp, err_cor, curlpgrad, curlhdiff, curlvdiff, res, &
                   advu, advv, advw, advVx, advVy, advVz, curlmet, err_nlsub, err_nldecomp, rr_rev, rr_cha, &
                   tlat, tlong, z_t
  implicit none
  private
  public :: load_params, get_yyyymmdd, init_zetavars, init_outputvars, create_outputfiles, loadave_vars, &
            write_outputfiles, release_vars, close_outputfiles

  !!----------------------------------------------------------------------------
  ! Input namelist file name
  character(len = *), parameter, public :: fn_namelist = 'input_eddy.nml'

  !!----------------------------------------------------------------------------
  ! * Average mode *
  ! a) avenm = "a" or avenmn = "m". For monthly and annual (default mean
  ! b) nda: Every fixed number of days
  character(len = 5), public :: avenm = ""
  integer, public :: nda = 365

  !!----------------------------------------------------------------------------
  ! * Input file info *
  ! Input file directory (suppose in the same directory)
  character(len = 300), public :: fn_in_dir = "/lustre/atlas1/cli115/proj-shared/jritchie/SIO/MODEL_DATA/yelpatch60/IOC/DATA/DAT_IO_all/nday1/"
  ! Input file prefix (before "year")
  character(len = 100), public :: fn_in_pfx = "ia_top_tx0.1_v2_yel_patc_1948_intel.pop.h.nday1."
  ! Input file suffix (subregion)
  character(len = 10), public :: fn_in_sfx = ".IO"
  ! Input file date delimiter (between yyyy, mm, dd)
  character(len = 1), public :: fn_in_dlm = '-'
  ! Domain boundary indices of the input file (relative to the entire globe. Here, from 3600x2400 )
  integer, public :: xl_ref = 1489, xr_ref = 1901, yd_ref = 1101, yu_ref = 1401

  ! Input file directory (suppose in the same directory)
  character(len = 300), public :: fn_zeta_dir = "/lustre/atlas1/cli115/proj-shared/hewang/data/vortbud/"
  ! Input file prefix (before "year")
  character(len = 100), public :: fn_zeta_pfx = "zeta_"
  ! Input file suffix (default is subregion)
  character(len = 10), public :: fn_zeta_sfx = ""
  ! Input file date delimiter (between yyyy, mm, dd. Default is the same as input)
  character(len = 1), public :: fn_zeta_dlm = ""

  ! Starting and ending year
  integer, public :: yrst = 2005, yred = 2009
  ! list of yr
  integer :: yrlist_in(60) = 0
  integer, allocatable, public :: yrlist(:)
  ! list of sections
  integer, allocatable, public :: sec(:,:)
  integer, public :: nyr, nsec ! number of years, number of sections

  !!----------------------------------------------------------------------------
  ! * output file info *
  ! Output file directory
  character(len = 300), public :: fn_out_dir = "/lustre/atlas/proj-shared/cli115/hewang/data/"
  ! Output file prefix
  character(len = 50), public :: fn_out_pfx = "zeta_meaneddy_"
  ! Output file suffix (default is subregion)
  character(len = 10), public :: fn_out_sfx = ""
  ! Output file date delimiter (between yyyy, mm, dd. Default is the same as input)
  character(len = 1), public :: fn_out_dlm = ""

  !!----------------------------------------------------------------------------
  ! Output variables
  real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
      curlnonlm, betavm, stretchpm, errcorm, curlpgradm, curlhdiffm, curlvdiffm, resm, err_nlsubm, &
      advuT, advvT, advwT, advVxT, advVyT, advVzT, curlmetT, err_nldecompT

  ! Output file handles
  integer, public :: ncid_out
  integer, public :: varid_out_curlnonl , varid_out_betav    , varid_out_stretchp , varid_out_errcor , &
                     varid_out_curlpgrad, varid_out_curlhdiff, varid_out_curlvdiff, varid_out_res    , &
                     varid_out_errnlsub , &
                     varid_out_advum    , varid_out_advvm    , varid_out_advwm , &
                     varid_out_advVxm   , varid_out_advVym   , varid_out_advVzm, &
                     varid_out_curlmetm , varid_out_errnldcmpm, &
                     varid_out_advue    , varid_out_advve    , varid_out_advwe , &
                     varid_out_advVxe   , varid_out_advVye   , varid_out_advVze, &
                     varid_out_curlmete , varid_out_errnldcmpe

  namelist /avemode/ avenm, nda
  namelist /bnds/ xl_reg, xr_reg, yd_reg, yu_reg, subreg, yrst, yred, yrlist_in, &
                  xi_dp, yi_dp, zi_dpst, zi_dped, ti_dp
  namelist /grid/ fn_grid, fn_dz, fn_cons
  namelist /input/ xl_ref, xr_ref, yd_ref, yu_ref, fn_in_dir, fn_in_pfx, fn_in_sfx, fn_in_dlm, &
                   fn_zeta_dir, fn_zeta_pfx, fn_zeta_sfx, fn_zeta_dlm
  namelist /output/ fn_out_dir, fn_out_pfx, fn_out_sfx, fn_out_dlm

contains
subroutine load_params()
  use, intrinsic :: IEEE_ARITHMETIC
  implicit none
  integer :: nml_error
  integer :: iyr, isec

  open(101, file=fn_namelist, status="old", iostat=nml_error)
  read(101, nml=avemode, iostat=nml_error)
    print*, "input <avemode> error: ", nml_error
  read(101, nml=bnds, iostat=nml_error)
    print*, "input <bnds> error: ", nml_error
  read(101, nml=grid, iostat=nml_error)
    print*, "input <grid> error: ", nml_error
  read(101, nml=input, iostat=nml_error)
    print*, "input <input> error: ", nml_error
  read(101, nml=output, iostat=nml_error)
    print*, "input <output> error: ", nml_error
  close(101)

  MVALUE = ieee_value(MVALUE, IEEE_QUIET_NAN)
  print*, 'Missing value: ', MVALUE

  ! Relative domain boundaries
  xl = xl_reg - xl_ref + 1
  xr = xr_reg - xl_ref + 1
  yd = yd_reg - yd_ref + 1
  yu = yu_reg - yd_ref + 1

  nx = xr - xl + 1
  ny = yu - yd + 1
  nyr = yred - yrst + 1

  ! Date format
  ! Two options in input.nml for yr
  ! (a): starting and ending (st, ed)
  ! (b): list
  if (count(yrlist_in/=0) /= 0 .and. all(yrlist_in >= 0)) then
    nyr = count(yrlist_in/=0)
    allocate(yrlist(nyr))
    yrlist = yrlist_in(1:nyr)
  elseif (yrst /= 0 .and. yrst /= 0) then
    nyr = yred - yrst + 1
    allocate(yrlist(nyr))
    yrlist = (/ (iyr, iyr = yrst, yred) /)
  else
    print*, "Input missing: either year list or yrnm_clm is needed! Aborting"
    stop
  endif

  if (avenm=="m") then
    nda = 0
    nsec = 12
    allocate(sec(nsec, 2))
    do isec = 1, nsec
      sec(isec, 1) = 1 + sum(eom(1:isec)) - eom(isec)
      sec(isec, 2) = sum(eom(1:isec))
    enddo
  elseif (avenm=="a") then
    nda = 0
    nsec = 1
    allocate(sec(nsec, 2))
    sec(1, 1) = 1
    sec(1, 2) = 365
  else
    nsec = ceiling(365. / nda)
    allocate(sec(nsec, 2))
    do isec = 1, nsec
      sec(isec, 1) = (isec - 1) * nda + 1
      sec(isec, 2) =  min(isec * nda, 365)
    enddo
  endif

  if (len(trim(fn_zeta_sfx)) == 0) then
    write(fn_zeta_sfx, '(A, A)') '_', trim(subreg)
  endif

  if (len(trim(fn_zeta_dlm)) == 0) then
    write(fn_zeta_dlm, '(A)') trim(fn_in_dlm)
  endif

    if (len(trim(fn_out_sfx)) == 0) then
    write(fn_out_sfx, '(A)') trim(fn_zeta_sfx)
  endif

  if (len(trim(fn_out_dlm)) == 0) then
    write(fn_out_dlm, '(A)') trim(fn_zeta_dlm)
  endif

  ! Print basic info
  print*, "Domain name :", subreg
  print*, "Domain size : ", "nx: ", nx, "ny: ", ny
  print*, "Year: ", yrlist
  print*, "Mean mode: ", trim(avenm), nda
  print*, "Section info: ", nsec
  print*, sec(:, 1)
  print*, sec(:, 2)
endsubroutine

subroutine get_yyyymmdd(yr, st, ed, avenm, dlm, yyyymmdd)
  implicit none
  integer, intent(in) :: yr, st, ed
  character(len = *), intent(in) :: avenm, dlm
  character(len = 15), intent(out) :: yyyymmdd
  integer :: mm, dd

  if (len(trim(avenm)) == 0) then
    write(yyyymmdd, '(I0.4, A, I0.3, A, I0.3)') yr, trim(dlm), st, trim(dlm), ed
  elseif (avenm == "m") then
    call doy2date(st, mm, dd)
    write(yyyymmdd, '(I0.4, A, I0.2)') yr, trim(dlm), mm
  elseif (avenm == "a") then
    write(yyyymmdd, '(I0.4, A, A)') yr, trim(dlm), 'ann'
  endif
endsubroutine

subroutine create_outputfiles(yyyymmdd)
  use netcdf
  implicit none
  character(len = *), intent(in) :: yyyymmdd
  character(len = 300) :: fn_out, fn_out_err
  integer :: stat_create, stat_defdim, stat_defvar, stat_putatt, stat_inqvar, &
             stat_getvar, stat_putvar, stat_io
  integer :: dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time
  integer :: varid_out_lon, varid_out_lat, varid_out_dep, varid_out_time

  print*, '  '
  print*, 'Creating output file ...'
  write(fn_out, '(A, A, A, A, A)') &
     trim(fn_out_dir), trim(fn_out_pfx), trim(yyyymmdd), trim(fn_out_sfx), '.nc'
  print*, '  Output file: ', trim(fn_out)

  print*, '  Start netcdf define ...'
  stat_create = nf90_create(fn_out, NF90_CLOBBER, ncid_out)
  stat_defdim = nf90_def_dim(ncid_out, "nlon", nx, dimid_out_lon)
  stat_defdim = nf90_def_dim(ncid_out, "nlat", ny, dimid_out_lat)
  stat_defdim = nf90_def_dim(ncid_out, "z_t" , nz, dimid_out_dep)
  stat_defdim = nf90_def_dim(ncid_out, "time", NF90_UNLIMITED, dimid_out_time)

  stat_defvar = nf90_def_var(ncid_out, "TLONG", NF90_FLOAT, &
     (/dimid_out_lon, dimid_out_lat/), varid_out_lon)
  stat_defvar = nf90_def_var(ncid_out, "TLAT",  NF90_FLOAT, &
     (/dimid_out_lon, dimid_out_lat/), varid_out_lat)
  stat_defvar = nf90_def_var(ncid_out, "z_t" ,  NF90_FLOAT, dimid_out_dep , varid_out_dep )
  stat_defvar = nf90_def_var(ncid_out, "time" , NF90_FLOAT, dimid_out_time, varid_out_time)

  stat_defvar = nf90_def_var(ncid_out, "curlnonl" , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlnonl )
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlnonl , "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlnonl , "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlnonl , "long_name", "Curl of nonlinear term (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlnonl , "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "betav"    , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_betav    )
  stat_putatt = nf90_put_att(ncid_out, varid_out_betav    , "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_betav    , "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_betav    , "long_name", "Advection of planetary vorticity term (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_betav    , "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "stretchp",  nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_stretchp  )
  stat_putatt = nf90_put_att(ncid_out, varid_out_stretchp , "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_stretchp , "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_stretchp , "long_name", "Planetary vorticity stretching term (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_stretchp , "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "errcor"  ,  nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errcor    )
  stat_putatt = nf90_put_att(ncid_out, varid_out_errcor   , "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errcor   , "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errcor   , "long_name", "Error from decomposing curl(-fv, fu) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errcor   , "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "curlpgrad", nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlpgrad)
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlpgrad, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlpgrad, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlpgrad, "long_name", "Curl of pressure gradient term (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlpgrad, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "curlhdiff", nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlhdiff)
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlhdiff, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlhdiff, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlhdiff, "long_name", "Curl of horizontal diffusion (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlhdiff, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "curlvdiff", nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlvdiff)
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlvdiff, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlvdiff, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlvdiff, "long_name", "Curl of vertical diffusion (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlvdiff, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "res"      , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_res)
  stat_putatt = nf90_put_att(ncid_out, varid_out_res      , "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_res      , "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_res      , "long_name", "Residual (lhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_res      , "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advu_m"     , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advum)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advum, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advum, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advum, "long_name", "Mean advection of relative vorticity by zonal velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advum, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advv_m"     , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advvm)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advvm, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advvm, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advvm, "long_name", "Mean advection of relative vorticity by meridional velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advvm, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advw_m"     , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advwm)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advwm, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advwm, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advwm, "long_name", "Mean advection of relative vorticity by vertical velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advwm, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advVx_m" , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVxm)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVxm, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVxm, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVxm, "long_name", "Mean twisting of zonal voricity by vertical velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVxm, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advVy_m"     , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVym)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVym, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVym, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVym, "long_name", "Mean twisting of meridional voricity by vertical velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVym, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advVz_m" , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVzm)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVzm, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVzm, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVzm, "long_name", "Mean twisting of vertical voricity by vertical velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVzm, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "curlmet_m"  , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlmetm)
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlmetm, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlmetm, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlmetm, "long_name", "Mean curl of metric term (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlmetm, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "errnldcmp_m", nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnldcmpm)
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnldcmpm, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnldcmpm, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnldcmpm, "long_name", "Mean error from nonlinear term due to decomposition (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnldcmpm, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advu_e"     , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advue)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advue, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advue, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advue, "long_name", "Eddy advection of relative vorticity by zonal velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advue, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advv_e"     , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advve)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advve, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advve, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advve, "long_name", "Eddy advection of relative vorticity by meridional velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advve, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advw_e"     , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advwe)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advwe, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advwe, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advwe, "long_name", "Eddy advection of relative vorticity by vertical velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advwe, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advVx_e" , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVxe)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVxe, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVxe, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVxe, "long_name", "Eddy twisting of zonal voricity by vertical velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVxe, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advVy_e"     , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVye)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVye, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVye, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVye, "long_name", "Eddy twisting of meridional voricity by vertical velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVye, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "advVz_e" , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVze)
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVze, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVze, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVze, "long_name", "Eddy twisting of vertical voricity by vertical velocity (flux form) (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_advVze, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "curlmet_e"  , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlmete)
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlmete, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlmete, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlmete, "long_name", "Eddy curl of metric term (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_curlmete, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "errnldcmp_e", nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnldcmpe)
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnldcmpe, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnldcmpe, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnldcmpe, "long_name", "Eddy error from nonlinear term due to decomposition (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnldcmpe, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_out, "errnlsub"    , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnlsub)
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnlsub, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnlsub, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnlsub, "long_name", "Error from decomposing nonlinear term (rhs)")
  stat_putatt = nf90_put_att(ncid_out, varid_out_errnlsub, "missing_value", MVALUE)

  stat_create = nf90_enddef(ncid_out)
  print*, "    Finished netcdf define!", stat_create

  stat_putvar = nf90_put_var(ncid_out, varid_out_lat,  tlat , &
     start = (/1, 1/), count = (/nx, ny/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_lon,  tlong, &
     start = (/1, 1/), count = (/nx, ny/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_dep,  z_t)
endsubroutine

subroutine init_zetavars()
  print*, "  "
  print*, "Initializing zeta module variables ..."

  allocate(uc(nx, ny, nz), vc(nx, ny, nz), wc(nx, ny, nz))
  allocate(advu (nx, ny, nz), advv (nx, ny, nz), advw (nx, ny, nz), &
           advVx(nx, ny, nz), advVy(nx, ny, nz), advVz(nx, ny, nz), &
           curlmet(nx, ny, nz), err_nldecomp(nx, ny, nz))
  allocate(rr_rev (nx, ny, nz), rr_cha(nx, ny, nz))

  uc = 0.
  vc = 0.
  wc = 0.

  advu  = 0.
  advv  = 0.
  advw  = 0.
  advVx = 0.
  advVy = 0.
  advVz = 0.
  curlmet = 0.
  err_nldecomp = 0.

  rr_rev = 0.
  rr_cha = 0.
endsubroutine

subroutine init_outputvars()
  print*, "  "
  print*, "Initializing output variables ..."

  ! Initializing output fields for curl of momentum terms mode
  allocate(curlnonlm(nx, ny, nz), betavm  (nx, ny, nz), stretchpm  (nx, ny, nz), &
           errcorm (nx, ny, nz), curlpgradm(nx, ny, nz), curlhdiffm(nx, ny, nz), &
           curlvdiffm(nx, ny, nz), resm(nx, ny, nz), err_nlsubm(nx, ny, nz))

  curlnonlm = 0.
  betavm    = 0.
  stretchpm = 0.
  betavm = 0.
  errcorm = 0.
  curlpgradm = 0.
  curlhdiffm = 0.
  curlvdiffm = 0.
  resm = 0.
  err_nlsubm = 0.

  allocate(advuT (nx, ny, nz), advvT (nx, ny, nz), advwT (nx, ny, nz), &
           advVxT(nx, ny, nz), advVyT(nx, ny, nz), advVzT(nx, ny, nz), &
           curlmetT(nx, ny, nz), err_nldecompT(nx, ny, nz))

  advuT = 0.
  advvT = 0.
  advwT = 0.
  advVxT = 0.
  advVyT = 0.
  advVzT = 0.
  curlmetT = 0.
  err_nldecompT = 0.
endsubroutine

subroutine loadave_vars(yr, da, secl)
  use ncio, only : nc_read
  implicit none
  integer, intent(in):: yr, da, secl
  integer :: mm, dd
  character(len = 15) :: yyyymmdd
  character(len = 300) :: fn_in, fn_zeta
  real(kind=kd_r), dimension(nx, ny, nz, 1) :: WORK

  print*, '  '
  print*, 'Loading velocity from input files, averaging over ', secl , 'days'

  call doy2date(da, mm, dd)
  write(yyyymmdd, '(I0.4, A, I0.2, A, I0.2)') yr, trim(fn_in_dlm), mm, trim(fn_in_dlm), dd
  write(fn_in, '(A, A, A, A, A)') &
      trim(fn_in_dir), trim(fn_in_pfx), trim(yyyymmdd), trim(fn_in_sfx), '.nc'

  write(yyyymmdd, '(I0.4, A, I0.2, A, I0.2)') yr, trim(fn_in_dlm), mm, trim(fn_in_dlm), dd
  write(fn_zeta, '(A, A, A, A, A)') &
      trim(fn_zeta_dir), trim(fn_zeta_pfx), trim(yyyymmdd), trim(fn_zeta_sfx), '.nc'

  print*, 'Reading from file : ', trim(fn_in)
  print*, '  and ', trim(fn_zeta)

  call nc_read(fn_zeta, 'curlnonl' , WORK)
  curlnonlm  = curlnonlm  + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'betav'    , WORK)
  betavm     = betavm     + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'stretchp' , WORK)
  stretchpm  = stretchpm  + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'errcor'   , WORK)
  errcorm    = errcorm    + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'curlpgrad', WORK)
  curlpgradm = curlpgradm + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'curlhdiff', WORK)
  curlhdiffm = curlhdiffm + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'curlvdiff', WORK)
  curlvdiffm = curlvdiffm + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'res'      , WORK)
  resm       = resm       + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'errnlsub' , WORK)
  err_nlsubm = err_nlsubm + sum(WORK, 4) / secl

  call nc_read(fn_zeta, 'advu'       , WORK)
  advuT         = advuT         + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'advv'       , WORK)
  advvT         = advvT         + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'advw'       , WORK)
  advwT         = advwT         + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'advVx'   , WORK)
  advVxT        = advVxT        + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'advVy'   , WORK)
  advVyT        = advVyT        + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'advVz'   , WORK)
  advVzT        = advVzT        + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'curlmet'    , WORK)
  curlmetT      = curlmetT      + sum(WORK, 4) / secl
  call nc_read(fn_zeta, 'errnldcmp', WORK)
  err_nldecompT = err_nldecompT + sum(WORK, 4) / secl

  call nc_read(fn_in, 'UVEL', WORK, (/xl, yd, 1, 1/), (/nx, ny, nz, 1/))
  where(WORK > 1e10) WORK = 0.
  uc = uc + sum(WORK, 4) / secl
  call nc_read(fn_in, 'VVEL', WORK, (/xl, yd, 1, 1/), (/nx, ny, nz, 1/))
  where(WORK > 1e10) WORK = 0.
  vc = vc + sum(WORK, 4) / secl
  call nc_read(fn_in, 'WVEL', WORK, (/xl, yd, 1, 1/), (/nx, ny, nz, 1/))
  where(WORK > 1e10) WORK = 0.
  wc = wc + sum(WORK, 4) / secl

  ! print*, 'curlnol: ', curlnonl_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'betav: ', betav_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'stretchp: ', stretchp_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'errcor: ', errcor_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'curlpgrad: ', curlpgrad_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'curlhdiff: ', curlhdiff_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'curlvdiff: ', curlvdiff_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'res: ', res_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'advu: ', advu_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'advv: ', advv_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'advw: ', advw_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'stretchr: ', stretchr_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'tilt: ', tilt_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'curlmet: ', curlmet_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'err_nlsub: ', err_nlsub_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'err_nldecomp: ', err_nldecomp_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)


  ! print*, 'uc: ', uc_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'vc: ', vc_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
  ! print*, 'wc: ', wc_mn(xi_dp, yi_dp, zi_dpst:zi_dped, ti_dp)
endsubroutine

subroutine write_outputfiles()
  use netcdf
  implicit none
  integer :: stat_putvar, stat_io

  print*, '  '
  print*, '  Start writing zeta eddy mean file ...'

  stat_putvar = nf90_put_var(ncid_out, varid_out_curlnonl , curlnonlm , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_betav    , betavm    , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_stretchp , stretchpm , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_errcor   , errcorm   , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_curlpgrad, curlpgradm, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_curlhdiff, curlhdiffm, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_curlvdiff, curlvdiffm, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_res      , resm      , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_errnlsub , err_nlsubm, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))

  stat_putvar = nf90_put_var(ncid_out, varid_out_advum , advu , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_advvm , advv , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_advwm , advw , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_advVxm, advVx, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_advVym, advVy, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_advVzm, advVz, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_curlmetm, curlmet, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_errnldcmpm, err_nldecomp, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))

  stat_putvar = nf90_put_var(ncid_out, varid_out_advue , advuT  - advu , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_advve , advvT  - advv , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_advwe , advwT  - advw , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_advVxe, advVxT - advVx, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_advVye, advVyT - advVy, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_advVze, advVzT - advVz, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_curlmete, curlmetT - curlmet, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_out, varid_out_errnldcmpe, err_nldecompT - err_nldecomp, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))

  stat_io = nf90_close(ncid_out)
  print*, '  Finished writing zeta mean eddy equation file '
endsubroutine

subroutine release_vars()
  print*, " "
  print*, '  Deallocating zeta and output variables ...'

  deallocate(uc, vc, wc, advu, advv, advw, advVx, advVy, advVz, curlmet, err_nldecomp, rr_rev, rr_cha)
  deallocate(curlnonlm, betavm, stretchpm, errcorm, curlpgradm, curlhdiffm, curlvdiffm, resm, err_nlsubm, &
             advuT, advvT, advwT, advVxT, advVyT, advVzT, curlmetT, err_nldecompT)
endsubroutine

subroutine close_outputfiles()
  use netcdf
  implicit none
  integer :: stat_io

  print*, '  '
  print*, 'Closing output files ...'
  stat_io = nf90_close(ncid_out)
endsubroutine

subroutine doy2date(ida, mm, dd)
  implicit none
  integer, intent(in) :: ida
  integer, intent(out) :: mm, dd
  integer :: imn, cumdays

  do imn = 1, 12
    cumdays = sum(eom(1:imn))
    if (ida <= cumdays) exit
  enddo

  mm = imn
  dd = ida - sum(eom(1:mm)) + eom(mm)
endsubroutine
endmodule
