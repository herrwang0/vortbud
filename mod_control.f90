module control
  use params
  use zeta, only : uc, vc, wc, ssh, advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy, &
                   curlnonl, betav, stretchp, err_cor, curlpgrad, curlhdiff, curlvdiff, res, &
                   advu, advv, advw, advVx, advVy, advVz, curlmet, err_nlsub, err_nldecomp, rr_rev, rr_cha, &
                   tlat, tlong, z_t
  implicit none
  private
  public :: load_params, get_yyyymmdd, init_zetavars, create_outputfiles, load_vars, &
            write_outputfiles, release_zetavars, close_outputfiles

  !!----------------------------------------------------------------------------
  ! Input namelist file name
  character(len = *), parameter, public :: fn_namelist = 'input.nml'

  !!----------------------------------------------------------------------------
  ! * Calculation mode *
  ! cmode = 0: zeta equation terms (nonlinear, curl of pgrad, hdiff, vdiff, residual,
  !            betav from Coriolis, stretching from Coriolis, error with Coriolis decomposition)
  ! cmode = 1 (default): including additional terms from the decomosition of nonlinear term
  ! cmode = 2: additional output file including detailed error info from nonlinear term decomposition
  integer, public :: cmode = 1

  !!----------------------------------------------------------------------------
  ! * Input file info *
  ! Input file directory (suppose in the same directory)
  character(len = 300), public :: fn_in_dir = "/lustre/atlas1/cli115/proj-shared/jritchie/SIO/MODEL_DATA/yelpatch60/IOC/DATA/DAT_IO_all/nday1/"
  ! Input file prefix (before "year")
  character(len = 100), public :: fn_in_pfx = "ia_top_tx0.1_v2_yel_patc_1948_intel.pop.h.nday1."
  ! Input file suffix (after "day", or subregion)
  character(len = 10), public :: fn_in_sfx = ".IO"
  ! Input file date delimiter (between yyyy, mm, dd)
  character(len = 1), public :: fn_in_dlm = '-'
  ! Domain boundary indices of the input file (relative to the entire globe. Here, from 3600x2400 )
  integer, public :: xl_ref = 1489, xr_ref = 1901, yd_ref = 1101, yu_ref = 1401

  ! Starting and ending year
  integer, public :: yrst = 2005, yred = 2009
  ! Starting and ending month
  integer, public :: mnst = 1, mned = 12
  ! Starting and ending day
  integer, public :: dast = 1, daed = 31
  ! list of yr
  integer :: yrlist_in(60) = 0
  integer, allocatable, public :: yrlist(:)
  ! list of month
  integer :: mnlist_in(12) = 0 ! default is zero, = 12 months
  integer, allocatable, public :: mnlist(:)
  ! list of days
  integer :: dalist_in(31) = 0 ! default is zero, = 31 days
  integer, allocatable, public :: dalist(:)
  ! Total number of months and years
  integer, public :: nmn, nyr, nda

  character(len=9), public :: yrnm_clm = ""

  !!----------------------------------------------------------------------------
  ! * output file info *
  ! Output file directory
  character(len = 300), public :: fn_out_dir = "/lustre/atlas/proj-shared/cli115/hewang/data/"
  ! Output file prefix
  character(len = 100), public :: fn_out_pfx = "vort_bud_"
  ! Output file suffix (default is subregion)
  character(len = 10), public :: fn_out_sfx = ""
  ! Output file date delimiter (between yyyy, mm, dd. Default is the same as input)
  character(len = 1), public :: fn_out_dlm = ""

  !!----------------------------------------------------------------------------
  ! Output file handles
  integer, public :: ncid_zeta, ncid_error
  integer, public :: varid_out_curlnonl , varid_out_betav    , varid_out_stretchp , varid_out_errcor , &
                     varid_out_curlpgrad, varid_out_curlhdiff, varid_out_curlvdiff, varid_out_res    , &
                     varid_out_advu     , varid_out_advv     , varid_out_advw     ,                    &
                     varid_out_advVx    , varid_out_advVy    , varid_out_advVz    , varid_out_curlmet, &
                     varid_out_errnlsub , varid_out_errnldcmp, varid_out_errnl_rev, varid_out_errnl_cha

  namelist /calcmode/ cmode
  namelist /bnds/ xl_reg, xr_reg, yd_reg, yu_reg, subreg, yrst, yred, mnst, mned, dast, daed, &
                  yrnm_clm, yrlist_in, mnlist_in, dalist_in, xi_dp, yi_dp, zi_dpst, zi_dped, ti_dp
  namelist /grid/ fn_grid, fn_dz, fn_cons
  namelist /input/ xl_ref, xr_ref, yd_ref, yu_ref, fn_in_dir, fn_in_pfx, fn_in_sfx, fn_in_dlm
  namelist /output/ fn_out_dir, fn_out_pfx, fn_out_sfx, fn_out_dlm

contains
subroutine load_params()
  use, intrinsic :: IEEE_ARITHMETIC
  implicit none
  integer :: nml_error
  integer :: iyr, imn, ida

  open(101, file=fn_namelist, status="old", iostat=nml_error)
  read(101, nml=calcmode, iostat=nml_error)
    print*, "input <calcmode> error: ", nml_error
  read(101, nml=bnds, iostat=nml_error)
    print*, "input <bnds> error: ", nml_error
  read(101, nml=grid, iostat=nml_error)
    print*, "input <grid> error: ", nml_error
  read(101, nml=input, iostat=nml_error)
    print*, "input <input> error: ", nml_error
  read(101, nml=output, iostat=nml_error)
    print*, "input <output> error: ", nml_error
  close(101)

  if (cmode < 0 .or. cmode > 2) then
    print*, 'Invalid calculation mode codename! Using default (cmode = 1)'
  endif

  MVALUE = ieee_value(MVALUE, IEEE_QUIET_NAN)
  print*, 'Missing value: ', MVALUE

  ! Relative domain boundaries
  xl = xl_reg - xl_ref + 1
  xr = xr_reg - xl_ref + 1
  yd = yd_reg - yd_ref + 1
  yu = yu_reg - yd_ref + 1

  nx = xr - xl + 1
  ny = yu - yd + 1

  ! Date format
  ! Two options in input.nml for yr, mn and da
  ! (a): starting and ending (st, ed)
  ! (b): list
  ! For year, additonal option is the name for climatology (yrnm_clm). If none is provide, then fatal error.
  ! For month and days, if neither list or st/ed is provided, then assume 12 months and *31 days
  ! If any negative number appears in yrlist, mnlist, dalist, then assume averaging level.
  !    For example, mnlist_in = -1 indicates input file is annually averaged (therefore, no mnlist AND dalist needed).
  if (len(trim(yrnm_clm)) > 0) then
    nyr = 1
    allocate(yrlist(nyr))
    yrlist = -1
  elseif (count(yrlist_in/=0) /= 0 .and. all(yrlist_in >= 0)) then
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

  if (any(mnlist_in < 0)) then
    nmn = 1
    allocate(mnlist(nmn))
    mnlist = -1
  elseif (count(mnlist_in/=0) /= 0 .and. all(mnlist_in >= 0)) then
    nmn = count(mnlist_in/=0)
    allocate(mnlist(nmn))
    mnlist = mnlist_in(1:nmn)
  else
    nmn = mned - mnst + 1
    allocate(mnlist(nmn))
    mnlist = (/ (imn, imn = mnst, mned) /)
  endif

  if (any(dalist_in < 0) .or. any(mnlist_in < 0)) then
    nda = 1
    allocate(dalist(nda))
    dalist = -1
  elseif (count(dalist_in/=0) /= 0 .and. all(dalist_in >= 0)) then
    nda = count(dalist_in/=0)
    allocate(dalist(nda))
    dalist = dalist_in(1:nda)
  else
    nda = daed - dast + 1
    allocate(dalist(nda))
    dalist = (/ (ida, ida = dast, daed) /)
  endif

  if (len(trim(fn_out_sfx)) == 0) then
    write(fn_out_sfx, '(A, A)') '_', trim(subreg)
  endif

  if (len(trim(fn_out_dlm)) == 0) then
    write(fn_out_dlm, '(A)') trim(fn_in_dlm)
  endif

  ! Print basic info
  print*, "Domain name :", subreg
  print*, "Domain size : ", "nx: ", nx, "ny: ", ny
  print*, "Year: ", yrlist
  print*, "Month: ", mnlist
  print*, "Day: ", dalist
  print*, "Calculation mode : ", cmode
endsubroutine

subroutine get_yyyymmdd(yr, mn, da, yrnm_clm, dlm, yyyymmdd)
  implicit none
  integer, intent(in) :: yr, mn, da
  character(len=*), intent(in) :: yrnm_clm, dlm
  character(len=15), intent(out) :: yyyymmdd
  character(len=9) :: yrnm, mnnm, danm

  if (len(trim(yrnm_clm)) > 0) then
    write(yrnm, '(A)') trim(yrnm_clm)
  else
    write(yrnm, '(I0.4)') yr
  endif

  if (mn < 0) then
    write(mnnm, '(A)') ''
  else
    write(mnnm, '(A, I0.2)') trim(dlm), mn
  endif

  if (da < 0) then
    write(danm, '(A)') ''
  else
    write(danm, '(A, I0.2)') trim(dlm), da
  endif

  write(yyyymmdd, '(A, A, A)') trim(yrnm), trim(mnnm), trim(danm)
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
  stat_create = nf90_create(fn_out, cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid_zeta)
  stat_defdim = nf90_def_dim(ncid_zeta, "nlon", nx, dimid_out_lon)
  stat_defdim = nf90_def_dim(ncid_zeta, "nlat", ny, dimid_out_lat)
  stat_defdim = nf90_def_dim(ncid_zeta, "z_t" , nz, dimid_out_dep)
  stat_defdim = nf90_def_dim(ncid_zeta, "time", NF90_UNLIMITED, dimid_out_time)

  stat_defvar = nf90_def_var(ncid_zeta, "TLONG", NF90_FLOAT, &
     (/dimid_out_lon, dimid_out_lat/), varid_out_lon)
  stat_defvar = nf90_def_var(ncid_zeta, "TLAT",  NF90_FLOAT, &
     (/dimid_out_lon, dimid_out_lat/), varid_out_lat)
  stat_defvar = nf90_def_var(ncid_zeta, "z_t" ,  NF90_FLOAT, dimid_out_dep , varid_out_dep )
  stat_defvar = nf90_def_var(ncid_zeta, "time" , NF90_FLOAT, dimid_out_time, varid_out_time)

  stat_defvar = nf90_def_var(ncid_zeta, "curlnonl" , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlnonl )
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonl , "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonl , "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonl , "long_name", "Curl of nonlinear term (rhs)")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonl , "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_zeta, "betav"    , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_betav    )
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_betav    , "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_betav    , "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_betav    , "long_name", "Advection of planetary vorticity term (rhs)")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_betav    , "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_zeta, "stretchp",  nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_stretchp  )
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_stretchp , "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_stretchp , "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_stretchp , "long_name", "Planetary vorticity stretching term (rhs)")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_stretchp , "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_zeta, "errcor"  ,  nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errcor    )
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_errcor   , "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_errcor   , "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_errcor   , "long_name", "Error from decomposing curl(-fv, fu) (rhs)")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_errcor   , "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_zeta, "curlpgrad", nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlpgrad)
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlpgrad, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlpgrad, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlpgrad, "long_name", "Curl of pressure gradient term (rhs)")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlpgrad, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_zeta, "curlhdiff", nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlhdiff)
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlhdiff, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlhdiff, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlhdiff, "long_name", "Curl of horizontal diffusion (rhs)")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlhdiff, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_zeta, "curlvdiff", nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlvdiff)
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlvdiff, "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlvdiff, "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlvdiff, "long_name", "Curl of vertical diffusion (rhs)")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlvdiff, "missing_value", MVALUE)

  stat_defvar = nf90_def_var(ncid_zeta, "res"      , nc_xtype, &
    (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_res)
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_res      , "Units", "1/s^2")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_res      , "coordinates", "TLONG TLAT z_t time")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_res      , "long_name", "Residual (lhs)")
  stat_putatt = nf90_put_att(ncid_zeta, varid_out_res      , "missing_value", MVALUE)

  if (cmode > 0) then
    stat_defvar = nf90_def_var(ncid_zeta, "advu"     , nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advu)
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advu, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advu, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advu, "long_name", "Advection of relative vorticity by zonal velocity (flux form) (rhs)")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advu, "missing_value", MVALUE)

    stat_defvar = nf90_def_var(ncid_zeta, "advv"     , nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advv)
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advv, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advv, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advv, "long_name", "Advection of relative vorticity by meridional velocity (flux form) (rhs)")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advv, "missing_value", MVALUE)

    stat_defvar = nf90_def_var(ncid_zeta, "advw"     , nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advw)
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advw, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advw, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advw, "long_name", "Advection of relative vorticity by vertical velocity (flux form) (rhs)")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advw, "missing_value", MVALUE)

    stat_defvar = nf90_def_var(ncid_zeta, "advVx" , nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVx)
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVx, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVx, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVx, "long_name", "Twisting of zonal voricity by vertical velocity (flux form) (rhs)")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVx, "missing_value", MVALUE)

    stat_defvar = nf90_def_var(ncid_zeta, "advVy"     , nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVy)
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVy, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVy, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVy, "long_name", "Twisting of meridional voricity by vertical velocity (flux form) (rhs)")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVy, "missing_value", MVALUE)

    stat_defvar = nf90_def_var(ncid_zeta, "advVz"     , nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVz)
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVz, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVz, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVz, "long_name", "Twisting of vertical voricity by vertical velocity (flux form) (rhs)")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVz, "missing_value", MVALUE)

    stat_defvar = nf90_def_var(ncid_zeta, "curlmet"  , nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlmet)
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlmet, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlmet, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlmet, "long_name", "Curl of metric term (rhs)")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlmet, "missing_value", MVALUE)

    stat_defvar = nf90_def_var(ncid_zeta, "errnlsub" , nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnlsub)
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "long_name", "Error from nonlinear term due to calculating offline (rhs)")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "missing_value", MVALUE)

    stat_defvar = nf90_def_var(ncid_zeta, "errnldcmp", nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnldcmp)
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnldcmp, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnldcmp, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnldcmp, "long_name", "Error from nonlinear term due to decomposition (rhs)")
    stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnldcmp, "missing_value", MVALUE)
  endif

  stat_create = nf90_enddef(ncid_zeta)
  print*, "    Finished netcdf define!", stat_create

  stat_putvar = nf90_put_var(ncid_zeta, varid_out_lat,  tlat , &
     start = (/1, 1/), count = (/nx, ny/))
  stat_putvar = nf90_put_var(ncid_zeta, varid_out_lon,  tlong, &
     start = (/1, 1/), count = (/nx, ny/))
  stat_putvar = nf90_put_var(ncid_zeta, varid_out_dep,  z_t)

  if (cmode == 2) then
    print*, '  '
    print*, 'Creating output error file ...'
    write(fn_out_err, '(A, A, A, A, A, A)') &
       trim(fn_out_dir), trim(fn_out_pfx), 'decompErrs_', trim(yyyymmdd), trim(fn_out_sfx), '.nc'
    print*, '  Output error file: ', trim(fn_out_err)

    stat_create = nf90_create(fn_out_err, cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid_error)
    stat_defdim = nf90_def_dim(ncid_error, "nlon", nx, dimid_out_lon)
    stat_defdim = nf90_def_dim(ncid_error, "nlat", ny, dimid_out_lat)
    stat_defdim = nf90_def_dim(ncid_error, "z_t" , nz, dimid_out_dep)
    stat_defdim = nf90_def_dim(ncid_error, "time", NF90_UNLIMITED, dimid_out_time)

    stat_defvar = nf90_def_var(ncid_error, "TLONG", NF90_FLOAT, &
       (/dimid_out_lon, dimid_out_lat/), varid_out_lon)
    stat_defvar = nf90_def_var(ncid_error, "TLAT",  NF90_FLOAT, &
       (/dimid_out_lon, dimid_out_lat/), varid_out_lat)
    stat_defvar = nf90_def_var(ncid_error, "z_t" ,  NF90_FLOAT, dimid_out_dep , varid_out_dep )
    stat_defvar = nf90_def_var(ncid_error, "time" , NF90_FLOAT, dimid_out_time, varid_out_time)

    stat_defvar = nf90_def_var(ncid_error, "rr_rev",   nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnl_rev)
    stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_rev, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_rev, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_rev, "missing_value", MVALUE)

    stat_defvar = nf90_def_var(ncid_error, "rr_cha",   nc_xtype, &
      (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnl_cha)
    stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_cha, "Units", "1/s^2")
    stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_cha, "coordinates", "TLONG TLAT z_t time")
    stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_cha, "missing_value", MVALUE)

    stat_create = nf90_enddef(ncid_error)
    print*, "    Finished netcdf define!", stat_create

    stat_putvar = nf90_put_var(ncid_error, varid_out_lat,  tlat , &
       start = (/1, 1/), count = (/nx, ny/))
    stat_putvar = nf90_put_var(ncid_error, varid_out_lon,  tlong, &
       start = (/1, 1/), count = (/nx, ny/))
    stat_putvar = nf90_put_var(ncid_error, varid_out_dep,  z_t)
  endif
endsubroutine

subroutine init_zetavars()
  print*, "  "
  print*, "Initializing zeta module variables ..."

  ! Initializing input fields for zeta module
  allocate(uc    (nx, ny, nz), vc    (nx, ny, nz), &
           advx  (nx, ny, nz), advy  (nx, ny, nz), &
           gradx (nx, ny, nz), grady (nx, ny, nz), &
           hdiffx(nx, ny, nz), hdiffy(nx, ny, nz), &
           vdiffx(nx, ny, nz), vdiffy(nx, ny, nz), &
           ssh   (nx, ny, 1 ))
  uc = 0.
  vc = 0.
  advx = 0.
  advy = 0.
  gradx = 0.
  grady = 0.
  hdiffx = 0.
  hdiffy = 0.
  vdiffx = 0.
  vdiffy = 0.
  ssh = 0.

  ! Initializing output fields for curl of momentum terms mode
  allocate(curlnonl (nx, ny, nz), curlpgrad(nx, ny, nz), res(nx, ny, nz), &
           curlhdiff(nx, ny, nz), curlvdiff(nx, ny, nz), &
           betav    (nx, ny, nz), stretchp (nx, ny, nz), err_cor(nx, ny, nz))
  curlnonl  = 0.
  curlpgrad = 0.
  curlhdiff = 0.
  curlvdiff = 0.
  res = 0.
  betav    = 0.
  stretchp = 0.
  err_cor  = 0.

  ! Initializing output fields for nonlinear term decomposition mode
  if (cmode > 0) then
    allocate(wc(nx, ny, nz))
    allocate(advu (nx, ny, nz), advv (nx, ny, nz), advw (nx, ny, nz), &
             advVx(nx, ny, nz), advVy(nx, ny, nz), advVz(nx, ny, nz), curlmet(nx, ny, nz), &
             err_nlsub(nx, ny, nz), err_nldecomp(nx, ny, nz))
    wc = 0.
    advu  = 0.
    advv  = 0.
    advw  = 0.
    advVx = 0.
    advVy = 0.
    advVz = 0.
    curlmet = 0.
    err_nlsub    = 0.
    err_nldecomp = 0.

    allocate(rr_rev(nx, ny, nz), rr_cha(nx, ny, nz))
    rr_rev = 0.
    rr_cha = 0.
  endif
endsubroutine

subroutine load_vars(yyyymmdd, iyr, imn, idy)
  use popload, only : load_current_day, find_daily_file
  implicit none
  integer, intent(in) :: iyr, imn, idy
  character(len = *), intent(in) :: yyyymmdd
  character(len = 300) :: fn_in
  real(kind=kd_r), dimension(nx, ny, nz, 1) :: WORK
  real(kind=kd_r), dimension(nx, ny, 1 , 1) :: WORK2

  print*, '  '
  print*, 'Loading variables (velocity, mom. terms) from input files'

  ! write(fn_in, '(A, A, A, A, A)') &
  !     trim(fn_in_dir), trim(fn_in_pfx), trim(yyyymmdd), trim(fn_in_sfx), '.nc'
  call find_daily_file(fn_in_dir, fn_in_pfx, iyr, imn, idy, fn_in)
  print*, 'Reading from file : ', trim(fn_in)

  call load_current_day(iyr, imn, idy, fn_in, 'UVEL'  , WORK)
  uc = WORK(:, :, :, 1)
  call load_current_day(iyr, imn, idy, fn_in, 'VVEL'  , WORK)
  vc = WORK(:, :, :, 1)

  call load_current_day(iyr, imn, idy, fn_in, 'ADVU'  , WORK)
  advx = WORK(:, :, :, 1)
  call load_current_day(iyr, imn, idy, fn_in, 'ADVV'  , WORK)
  advy = WORK(:, :, :, 1)

  call load_current_day(iyr, imn, idy, fn_in, 'GRADX' , WORK)
  gradx = WORK(:, :, :, 1)
  call load_current_day(iyr, imn, idy, fn_in, 'GRADY' , WORK)
  grady = WORK(:, :, :, 1)

  call load_current_day(iyr, imn, idy, fn_in, 'HDIFFU', WORK)
  hdiffx = WORK(:, :, :, 1)
  call load_current_day(iyr, imn, idy, fn_in, 'HDIFFV', WORK)
  hdiffy = WORK(:, :, :, 1)

  call load_current_day(iyr, imn, idy, fn_in, 'VDIFFU', WORK)
  vdiffx = WORK(:, :, :, 1)
  call load_current_day(iyr, imn, idy, fn_in, 'VDIFFV', WORK)
  vdiffy = WORK(:, :, :, 1)

  call load_current_day(iyr, imn, idy, fn_in, 'SSH'   , WORK2)
  ssh(:, :, 1) = WORK2(:, :, 1, 1)

  where(abs(uc    ) > 1e10) uc     = 0.
  where(abs(vc    ) > 1e10) vc     = 0.
  where(abs(advx  ) > 1e10) advx   = 0.
  where(abs(advy  ) > 1e10) advx   = 0.
  where(abs(gradx ) > 1e10) gradx  = 0.
  where(abs(grady ) > 1e10) grady  = 0.
  where(abs(hdiffx) > 1e10) hdiffx = 0.
  where(abs(hdiffy) > 1e10) hdiffy = 0.
  where(abs(vdiffx) > 1e10) vdiffx = 0.
  where(abs(vdiffy) > 1e10) vdiffy = 0.

  if(cmode > 0) then
    call load_current_day(iyr, imn, idy, fn_in, 'WVEL'  , WORK)
    wc = WORK(:, :, :, 1)
    where(abs(wc) > 1e10) wc = 0.
  endif

  print*, 'uc: ', uc(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'vc: ', vc(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'Advx', advx(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'Advy', advy(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'Gradx', gradx(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'Grady', grady(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'Hdiffx', hdiffx(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'Hdiffy', hdiffy(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'Vdiffx', vdiffx(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'Vdiffy', vdiffy(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'SSH', ssh(xi_dp, yi_dp, 1)
endsubroutine

subroutine write_outputfiles()
  use netcdf
  implicit none
  integer :: stat_putvar

  print*, '  '
  print*, '  Start writing zeta equation file ...'

  stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlnonl , curlnonl , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_zeta, varid_out_betav    , betav    , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_zeta, varid_out_stretchp , stretchp , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_zeta, varid_out_errcor   , err_cor  , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlpgrad, curlpgrad, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlhdiff, curlhdiff, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlvdiff, curlvdiff, &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  stat_putvar = nf90_put_var(ncid_zeta, varid_out_res      , res      , &
      start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))

  if (cmode > 0) then
    stat_putvar = nf90_put_var(ncid_zeta, varid_out_advu     , advu        , &
        start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
    stat_putvar = nf90_put_var(ncid_zeta, varid_out_advv     , advv        , &
        start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
    stat_putvar = nf90_put_var(ncid_zeta, varid_out_advw     , advw        , &
        start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
    stat_putvar = nf90_put_var(ncid_zeta, varid_out_advVx    , advVx       , &
        start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
    stat_putvar = nf90_put_var(ncid_zeta, varid_out_advVy    , advVy       , &
        start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
    stat_putvar = nf90_put_var(ncid_zeta, varid_out_advVz    , advVz       , &
        start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
    stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlmet  , curlmet     , &
        start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
    stat_putvar = nf90_put_var(ncid_zeta, varid_out_errnlsub , err_nlsub   , &
        start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
    stat_putvar = nf90_put_var(ncid_zeta, varid_out_errnldcmp, err_nldecomp, &
        start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
  endif
  print*, '  Finished writing zeta equation file '

  if (cmode == 2) then
    print*, '  '
    print*, '  Start writing error file ...'

    stat_putvar = nf90_put_var(ncid_error, varid_out_errnl_rev, rr_rev, &
       start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
    stat_putvar = nf90_put_var(ncid_error, varid_out_errnl_cha, rr_cha, &
       start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
    print*, '  Finished writing error file '
  endif
endsubroutine

subroutine release_zetavars()
  print*, " "
  print*, '  Deallocating zeta input and output variables ...'

  deallocate(curlnonl, curlpgrad, res, curlhdiff, curlvdiff, betav, stretchp, err_cor)
  deallocate(uc, vc, advx, advy, grady, gradx, hdiffx, hdiffy, vdiffx, vdiffy, ssh)

  if (cmode > 0) then
    deallocate(advu, advv, advw, advVx, advVy, advVz, curlmet, err_nlsub, err_nldecomp)
    deallocate(wc)
    deallocate(rr_rev, rr_cha)
  endif
endsubroutine

subroutine close_outputfiles()
  use netcdf
  implicit none
  integer :: stat_io

  print*, '  '
  print*, 'Closing output files ...'
  stat_io = nf90_close(ncid_zeta)
  if (cmode == 2) then
    stat_io = nf90_close(ncid_error)
  endif

  ! print*, 'Releasing input variables ...'
  ! deallocate(uc, vc, advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy, ssh)
  ! if (cmode > 0) then
  !   deallocate(wc)
  ! endif
endsubroutine

endmodule
