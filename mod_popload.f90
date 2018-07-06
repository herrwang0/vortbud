module popload
  use ncio, only : nc_read
  private
  public :: load_current_mon, load_current_day, find_daily_file !, fix_missing, ifmiss, current_fn_in

contains
! subroutine get_fn_mon(iyr, imn, fn_in_dir, fn_in_pfx, fn_in_sfx, fn)
!   implicit none
!   integer, intent(in) :: iyr, imn
!   character(len=*), intent(in) :: fn_in_dir, fn_in_pfx, fn_in_sfx
!   character(len=*), intent(inout) :: fn
!
!   write(fn, '(A, A, I4, A, I0.2, A, A)') &
!     trim(fn_in_dir), trim(fn_in_pfx), iyr, '-', imn, trim(fn_in_sfx), '.nc'
! endsubroutine

subroutine load_current_day(iyr, imn, idy, fn, varname, var)
  use netcdf
  use params, only : kd_r, xl, yd, xi_dp, yi_dp, zi_dpst, zi_dped
  implicit none
  integer, intent(in) :: iyr, imn, idy
  character(len=*), intent(in) :: fn, varname
  integer :: nx, ny, nz
  real(kind=kd_r), dimension(:, :, :, :), intent(inout) :: var
  integer :: stat_io, ncid

  nx = size(var, 1)
  ny = size(var, 2)
  nz = size(var, 3)

  stat_io = nf90_open(trim(fn), NF90_NOWRITE, ncid)
  if (stat_io /= NF90_NOERR) then
    call fix_missing(iyr, imn, idy, fn, varname, var)
  else
    call nc_read(fn, varname, var, (/xl, yd, 1, 1/), (/nx, ny, nz, 1/))
  endif

  contains
    subroutine fix_missing(iyr, imn, idy, fn, varname, var)
      use netcdf
      use params, only : kd_r, xl, yd, xi_dp, yi_dp, zi_dpst, zi_dped
      implicit none

      integer, intent(in) :: iyr, imn, idy
      character (len=*), intent(in) :: fn
      character (len=*), intent(in) :: varname
      integer :: nx, ny, nz
      real(kind=kd_r), intent(inout), dimension(:, :, :, :) :: var
      integer :: idx
      character (len = 300) :: fn1, fn2
      real(kind=kd_r), dimension(:, :, :, :), allocatable :: temp1, temp2

      nx = size(var, 1)
      ny = size(var, 2)
      nz = size(var, 3)
      allocate(temp1(nx, ny, nz, 1), temp2(nx, ny, nz, 1))

      if (iyr == 2007 .and. imn == 1 .and. idy == 31) then
        print*, 'Encountering missing day: Jan 31st, 2007, start interpreting '
        idx = index(fn, '-01-31')
        fn1 = fn
        fn1(idx:idx+5) = '-01-30'
        fn2 = fn
        fn2(idx:idx+5) = '-02-01'
        print*, "Using file ", trim(fn1)
        print*, "  and file ", trim(fn2)

        call nc_read(fn , varname, temp1, (/xl, yd, 1, 1/), (/nx, ny, nz, 1/))
        call nc_read(fn1, varname, temp2, (/xl, yd, 1, 1/), (/nx, ny, nz, 1/))
      endif

      if (iyr == 2005 .and. imn == 11 .and. idy == 26) then
        print*, 'Encountering missing day: Nov 26th, 2005, start interpreting '
        idx = index(fn, '-11-26')
        fn1 = fn
        fn1(idx:idx+5) = '-11-25'
        fn2 = fn
        fn2(idx:idx+5) = '-11-27'
        print*, "Using file ", trim(fn1)
        print*, "  and file ", trim(fn2)

        call nc_read(fn, varname, temp1, (/xl, yd, 1, 1/), (/nx, ny, nz, 1/))
        call nc_read(fn, varname, temp2, (/xl, yd, 1, 1/), (/nx, ny, nz, 1/))
      endif

      var = (temp1 + temp2) / 2.
    endsubroutine
endsubroutine

subroutine load_current_mon(iyr, imn, fn, varname, var)
  use netcdf
  use params, only : kd_r, xl, yd, xi_dp, yi_dp, zi_dpst, zi_dped
  implicit none
  integer, intent(in) :: iyr, imn
  character(len=*), intent(in) :: fn, varname
  real(kind=kd_r), dimension(:, :, :, :), intent(inout) :: var
  integer :: ncid, dimid, stat_io, stat_inqdimid, stat_inqdim
  integer :: nx, ny, nz, nt_r
  real(kind=kd_r), dimension(:, :, :, :), allocatable :: temp    ! for missing days

  stat_io = nf90_open(trim(fn), NF90_NOWRITE, ncid)
  stat_inqdimid = nf90_inq_dimid(ncid, 'time', dimid)
  stat_inqdim   = nf90_inquire_dimension(ncid, dimid, len = nt_r)

  call nc_read(fn, varname, var, (/xl, yd, 1, 1/), (/nx, ny, nz, nt_r/))
  call fix_missing(iyr, imn, fn, varname, var)

  contains
    subroutine fix_missing(iyr, imn, fn, varname, var)
      use netcdf
      use params, only : kd_r, xl, yd, xi_dp, yi_dp, zi_dpst, zi_dped
      implicit none

      integer, intent(in) :: iyr, imn
      character (len=*), intent(in) :: fn
      character (len=*), intent(in) :: varname
      integer :: nx, ny, nz
      real(kind=kd_r), intent(inout), dimension(:, :, :, :) :: var
      integer :: idx
      character (len = 300) :: fn1
      real(kind=kd_r), dimension(:, :, :, :), allocatable :: temp1, temp2

      if (iyr == 2007 .and. imn == 1) then
        print*, 'Encountering missing day: Jan 31st, 2007, start interpreting '

        nx = size(var, 1)
        ny = size(var, 2)
        nz = size(var, 3)
        allocate(temp1(nx, ny, nz, 1), temp2(nx, ny, nz, 1))

        idx = index(fn, '-01')
        fn1 = fn
        fn1(idx:idx+2) = '-02'
        print*, "Using additional file ", trim(fn1)

        call nc_read(fn , varname, temp1, (/xl, yd, 1, 30/), (/nx, ny, nz, 1/))
        call nc_read(fn1, varname, temp2, (/xl, yd, 1,  1/), (/nx, ny, nz, 1/))

        var(:, :, :, 31) = (temp1(:, :, :, 1) + temp2(:, :, :, 1)) / 2.
      endif

      if (iyr == 2005 .and. imn == 11) then
        print*, 'Encountering missing day: Nov 26th, 2005, start interpreting '

        nx = size(var, 1)
        ny = size(var, 2)
        nz = size(var, 3)
        allocate(temp1(nx, ny, nz, 1), temp2(nx, ny, nz, 1))

        call nc_read(fn, varname, temp1, (/xl, yd, 1, 25/), (/nx, ny, nz, 1/))
        call nc_read(fn, varname, temp2, (/xl, yd, 1, 26/), (/nx, ny, nz, 1/))

        var(:, :, :, 27:30) = var (:, :, :, 26:29)
        var(:, :, :, 26   ) = (temp1(:, :, :, 1) + temp2(:, :, :, 1)) / 2.
      endif
    endsubroutine
endsubroutine


subroutine find_daily_file(fn_in_dir, fn_in_pfx, yr, mon, day, fn)
    use netcdf
    implicit none
    character(len = *), intent(in) :: fn_in_dir, fn_in_pfx
    integer, intent(in) :: yr, mon, day
    character(len = 300) :: fn_short
    character(len = 300), intent(inout) :: fn
    integer :: iostat, ncid

    write(fn_short, '(A, I4, A, I0.2, A, I0.2, A)') &
        trim(fn_in_pfx), yr, '-', mon, '-', day, '.nc'

    write(fn, '(A, A)') &
        trim(fn_in_dir), trim(fn_short)
    iostat = nf90_open(trim(fn), NF90_NOWRITE, ncid)
    if (iostat == nf90_noerr) return

    write(fn, '(A, A, I4, A, I0.2, A, A)') &
        trim(fn_in_dir), '/pop_', yr, '-', mon, '/', trim(fn_short)
    iostat = nf90_open(trim(fn), NF90_NOWRITE, ncid)
    if (iostat == nf90_noerr) return

    write(fn, '(A, A, I4, I0.2, A, A)') &
        trim(fn_in_dir), 'z_', yr, mon, '/', trim(fn_short)
    iostat = nf90_open(trim(fn), NF90_NOWRITE, ncid)
    if (iostat == nf90_noerr) return

    select case (day)
        case (1:9)
            write(fn, '(A, A, I4, A, I0.2, A, A)') &
                trim(fn_in_dir), '/pop_', yr, '-', mon, '_1/', trim(fn_short)
        case (10:19)
            write(fn, '(A, A, I4, A, I0.2, A, A)') &
                trim(fn_in_dir), '/pop_', yr, '-', mon, '_2/', trim(fn_short)
        case (20:)
            write(fn, '(A, A, I4, A, I0.2, A, A)') &
                trim(fn_in_dir), '/pop_', yr, '-', mon, '_3/', trim(fn_short)
    endselect
    return
endsubroutine

! function ifmiss(yr, mn, dy)
!   implicit none
!   integer, intent(in) :: yr
!   integer, intent(in) :: mn, dy
!   integer :: missinglist(2, 3)
!   integer :: il
!   logical :: ifmiss
!
!   missinglist(1, :) = (/2007,  1, 31/)
!   missinglist(2, :) = (/2005, 11, 26/)
!   ifmiss = .false.
!
!     do il = 1, size(missinglist, 1)
!       if (missinglist(il, 1) == yr .and. missinglist(il, 2) == mn .and. missinglist(il, 3) == dy) then
!         ifmiss = .true.
!         exit
!       endif
!     enddo
! endfunction
endmodule
