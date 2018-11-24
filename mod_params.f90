module params
    use netcdf
    implicit none
    private

    type boundary_indices
        ! Domain name (used as a suffix of output files)
        character(len = 10), public :: subreg = 'sc'
        ! Domain boundary indices (relative to the entire globe. Here, from 3600x2400 )
        integer :: xl_reg = 1482, xr_reg = 1881, yd_reg = 1032, yu_reg = 1452
        ! Domain boundary indices of the input file (relative to the entire globe. Here, from 3600x2400 )
        integer, public :: xl_ref = 1422, xr_ref = 2150, yd_ref = 978, yu_ref = 1510
        ! Domain boundary indices relative to the input files
        integer :: xl, xr, yd, yu
        ! Number of x, y cells, number of years
        integer :: nx, ny, nz
        ! x, y, z index for display
        integer :: xi_dp, yi_dp, zi_dpst, zi_dped, ti_dp
    endtype

    type grid_files
        ! * Grid files info (with full path included. Domain is global for all files) *
        ! Basic 2-D grid info
        character(len = 300) :: grid = "./ocn_static_grid.nc"
        ! Dzt and Dzu calculated offline
        character(len = 300) :: dz = "./ocn_static_dz.nc"
        ! Model constants (omega, grav etc.)
        character(len = 300) :: cons = "./ocn_constant.nc"
    endtype

    !!--------------------------------------------------------------------------
    ! * Parameters *
    ! Real kind value used throughout
    integer, parameter, public :: kd_r = selected_real_kind(7)
    ! Number of vertical layers
    integer, parameter, public :: nzgl = 42
    ! Number of days in each month
    integer, parameter, dimension(12), public :: eom = (/31, 28, 31, 30, 31, 30, &
                                                         31, 31, 30, 31, 30, 31/)
    ! Missing valeu (NaN)
    real(kind=kd_r), public :: MVALUE
    ! Output file precision (Default = NF90_FLOAT)
    integer, public :: nc_xtype = NF90_FLOAT

    !!--------------------------------------------------------------------------
    type(boundary_indices), public :: B
    type(grid_files), public :: fngrid
    !!--------------------------------------------------------------------------

    character(len=30), public :: fn_namelist = "input.nml"
    character(len=30), public :: fmt_exp = "E12.4E2"
    character(len=30), public :: fmt_flt = "F10.4"
endmodule
