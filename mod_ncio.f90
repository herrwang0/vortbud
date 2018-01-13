module ncio
  use params, only : kd_r
  private
  public :: nc_read

contains

! Reading variable (varname) from netCDF file (fn).
!   var is a 4-D assumed-shape array, and dimst, dimct, dimsd (to be added)
!   are 4-element arrays for indexing.
!   For variables with dimension less than four, var needs to be filled with len=1 dimensions
subroutine nc_read(fn, varname, var, dimst, dimct, dimsd)
  use netcdf
  implicit none

  character (len=*), intent(in) :: fn
  character (len=*), intent(in) :: varname
  real(kind=kd_r), dimension(:,:,:,:), intent(inout) :: var
  integer, dimension(:), optional, intent(in) :: dimst, dimct, dimsd
  integer, dimension(:), allocatable :: dimst_local, dimct_local
  character (len=10) :: dimname

  integer :: stat_io, stat
  integer :: ncid, varid
  integer :: xtype, ndims, dimlen, idim, idim_in
  integer, dimension(NF90_MAX_VAR_DIMS) :: dimID

  dimID = 0

 ! open file
  stat_io = nf90_open(path = trim(fn), mode = NF90_NOWRITE, ncid = ncid)
  if (stat_io /= NF90_NOERR) then
    print*, 'Opening file error: ', trim(fn), stat_io
    print*, 'Abort'
    return
  endif

 ! inquire variable info
  stat = nf90_inq_varid(ncid = ncid, name = trim(varname), varid = varid)
  if (stat /= NF90_NOERR) then
    print*, 'Failed to find variable: ', trim(varname), stat
    print*, 'Abort'
    return
  endif

  stat = nf90_inquire_variable(ncid = ncid, varid = varid, &
     xtype = xtype, ndims = ndims, dimids = dimID)
  if (stat /= NF90_NOERR) then
    print*, 'Inquire variable error: ', trim(varname), stat
  endif

  if(ndims > 4) then
    print*, 'Variable has over four dimensions, aborting'
    print*, 'Dimension ID: ', dimID(1:ndims)
    return
  endif

    ! stat = nf90_inquire_dimension(ncid = ncid, dimid = dimID(idim), name = dimname, len = dimlen(idim))
  if (present(dimst) .and. present(dimct)) then
    allocate(dimst_local(ndims), dimct_local(ndims))
    dimst_local = 1
    dimct_local = 1

    idim_in = 0
    do idim = 1, ndims
      idim_in = idim_in + 1
      if(ndims < 4 .and. dimst(idim_in) == 1 .and. dimct(idim_in) == 1) cycle

      dimst_local(idim) = dimst(idim_in)
      dimct_local(idim) = dimct(idim_in)
    enddo
    stat = nf90_get_var(ncid = ncid, varid = varid, values = var, &
        start = dimst_local, count = dimct_local)
  else
    stat = nf90_get_var(ncid = ncid, varid = varid, values = var)
  endif

  if (stat /= NF90_NOERR) then
    print*, 'Failed to read variable: ', trim(varname)
    print*, 'Abort'
    return
  endif

! start
! A vector of integers specifying the index in the variable from which the first (or only) of the data values will be read. The indices are relative to 1, so for example, the first data value of a variable would have index (1, 1, ..., 1). The elements of start correspond, in order, to the variable's dimensions. Hence, if the variable is a record variable, the last index would correspond to the starting record number for writing the data values.
! By default, start(:) = 1.
!
! count
! A vector of integers specifying the number of indices selected along each dimension. To read a single value, for example, specify count as (1, 1, ..., 1). The elements of count correspond, in order, to the variable's dimensions. Hence, if the variable is a record variable, the last element of count corresponds to a count of the number of records to read.
! By default, count(:numDims) = shape(values) and count(numDims + 1:) = 1, where numDims = size(shape(values)).
!
! stride
! A vector of integers that specifies the sampling interval along each dimension of the netCDF variable. The elements of the stride vector correspond, in order, to the netCDF variable's dimensions (stride(1) gives the sampling interval along the most rapidly varying dimension of the netCDF variable). Sampling intervals are specified in type-independent units of elements (a value of 1 selects consecutive elements of the netCDF variable along the corresponding dimension, a value of 2 selects every other element, etc.).
! By default, stride(:) = 1.
endsubroutine

subroutine nc_read2d(fn, varname, var2d, dimst, dimct, dimsd)
  use netcdf
  implicit none

  integer, parameter :: ndims_varin = 2

  character (len=*), intent(in) :: fn
  character (len=*), intent(in) :: varname
  real(kind=kd_r), dimension(:,:), intent(inout) :: var2d
  integer, dimension(:), optional, intent(in) :: dimst, dimct, dimsd
!  integer, dimension(:), allocatable :: dimst_local, dimct_local, dimsd_local

  integer :: stat_io, stat
  integer :: ncid, varid
  integer :: xtype, ndims, dimlen

  integer :: nx, ny
  nx = size(var2d, 1)
  ny = size(var2d, 2)

! open file
  stat_io = nf90_open(path = trim(fn), mode = NF90_NOWRITE, ncid = ncid)
  if (stat_io /= NF90_NOERR) then
    print*, 'Opening file error: ', trim(fn), stat_io
    print*, 'Abort'
    return
  endif

! inquire variable info
  stat = nf90_inq_varid(ncid = ncid, name = trim(varname), varid = varid)
  if (stat /= NF90_NOERR) then
    print*, 'Failed to locate variable: ', trim(varname), stat
    print*, 'Abort'
    return
  endif

  stat = nf90_inquire_variable(ncid = ncid, varid = varid, &
     xtype = xtype, ndims = ndims)
  if (stat /= NF90_NOERR) then
    print*, 'Inquire variable error: ', trim(varname), stat
  endif

  if (present(dimst) .and. present(dimct)) then
    stat = nf90_get_var(ncid = ncid, varid = varid, values = var2d, &
        start = dimst, count = dimct)
  else
    stat = nf90_get_var(ncid = ncid, varid = varid, values = var2d)
  endif

  if (stat /= NF90_NOERR) then
    print*, 'Failed to read variable: ', trim(varname)
    print*, 'Abort'
    return
  endif

! start
! A vector of integers specifying the index in the variable from which the first (or only) of the data values will be read. The indices are relative to 1, so for example, the first data value of a variable would have index (1, 1, ..., 1). The elements of start correspond, in order, to the variable's dimensions. Hence, if the variable is a record variable, the last index would correspond to the starting record number for writing the data values.
! By default, start(:) = 1.
!
! count
! A vector of integers specifying the number of indices selected along each dimension. To read a single value, for example, specify count as (1, 1, ..., 1). The elements of count correspond, in order, to the variable's dimensions. Hence, if the variable is a record variable, the last element of count corresponds to a count of the number of records to read.
! By default, count(:numDims) = shape(values) and count(numDims + 1:) = 1, where numDims = size(shape(values)).
!
! stride
! A vector of integers that specifies the sampling interval along each dimension of the netCDF variable. The elements of the stride vector correspond, in order, to the netCDF variable's dimensions (stride(1) gives the sampling interval along the most rapidly varying dimension of the netCDF variable). Sampling intervals are specified in type-independent units of elements (a value of 1 selects consecutive elements of the netCDF variable along the corresponding dimension, a value of 2 selects every other element, etc.).
! By default, stride(:) = 1.
endsubroutine

subroutine nc_read3d(fn, varname, var3d, dimst, dimct, dimsd)
  use netcdf
  implicit none

  integer, parameter :: ndims_varin = 3

  character (len=*), intent(in) :: fn
  character (len=*), intent(in) :: varname
  real(kind=kd_r), dimension(:,:,:), intent(inout) :: var3d
  integer, dimension(:), optional, intent(in) :: dimst, dimct, dimsd
!  integer, dimension(:), allocatable :: dimst_local, dimct_local, dimsd_local

  integer :: stat_io, stat
  integer :: ncid, varid
  integer :: xtype, ndims, dimlen

  integer :: nx, ny, nz
  nx = size(var3d, 1)
  ny = size(var3d, 2)
  nz = size(var3d, 3)


! open file
  stat_io = nf90_open(path = trim(fn), mode = NF90_NOWRITE, ncid = ncid)
  if (stat_io /= NF90_NOERR) then
    print*, 'Opening file error: ', trim(fn), stat_io
    print*, 'Abort'
    return
  endif

! inquire variable info
  stat = nf90_inq_varid(ncid = ncid, name = trim(varname), varid = varid)
  if (stat /= NF90_NOERR) then
    print*, 'Failed to locate variable: ', trim(varname), stat
    print*, 'Abort'
    return
  endif

  stat = nf90_inquire_variable(ncid = ncid, varid = varid, &
     xtype = xtype, ndims = ndims)
  if (stat /= NF90_NOERR) then
    print*, 'Inquire variable error: ', trim(varname), stat
  endif

  if (present(dimst) .and. present(dimct)) then
    stat = nf90_get_var(ncid = ncid, varid = varid, values = var3d, &
        start = dimst, count = dimct)
  else
    stat = nf90_get_var(ncid = ncid, varid = varid, values = var3d)
  endif

  if (stat /= NF90_NOERR) then
    print*, 'Failed to read variable: ', trim(varname)
    print*, 'Abort'
    return
  endif

! start
! A vector of integers specifying the index in the variable from which the first (or only) of the data values will be read. The indices are relative to 1, so for example, the first data value of a variable would have index (1, 1, ..., 1). The elements of start correspond, in order, to the variable's dimensions. Hence, if the variable is a record variable, the last index would correspond to the starting record number for writing the data values.
! By default, start(:) = 1.
!
! count
! A vector of integers specifying the number of indices selected along each dimension. To read a single value, for example, specify count as (1, 1, ..., 1). The elements of count correspond, in order, to the variable's dimensions. Hence, if the variable is a record variable, the last element of count corresponds to a count of the number of records to read.
! By default, count(:numDims) = shape(values) and count(numDims + 1:) = 1, where numDims = size(shape(values)).
!
! stride
! A vector of integers that specifies the sampling interval along each dimension of the netCDF variable. The elements of the stride vector correspond, in order, to the netCDF variable's dimensions (stride(1) gives the sampling interval along the most rapidly varying dimension of the netCDF variable). Sampling intervals are specified in type-independent units of elements (a value of 1 selects consecutive elements of the netCDF variable along the corresponding dimension, a value of 2 selects every other element, etc.).
! By default, stride(:) = 1.
endsubroutine
endmodule
