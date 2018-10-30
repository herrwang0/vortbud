module params
    use netcdf
    implicit none
    private
    !!--------------------------------------------------------------------------
    ! * Parameters *
    ! Real kind value used throughout
    integer, parameter, public :: kd_r = selected_real_kind(7)
    ! Number of vertical layers
    integer, parameter, public :: nz = 42
    ! Number of days in each month
    integer, parameter, dimension(12), public :: eom = (/31, 28, 31, 30, 31, 30, &
                                                         31, 31, 30, 31, 30, 31/)
    ! Missing valeu (NaN)
    real(kind=kd_r), public :: MVALUE
    ! Output file precision (Default = NF90_FLOAT)
    integer, public :: nc_xtype = NF90_FLOAT

    !!--------------------------------------------------------------------------
    ! * Domain info *
    ! Domain name (used as a suffix of output files)
    character(len = 10), public :: subreg = 'sc'
    ! Domain boundary indices (relative to the entire globe. Here, from 3600x2400 )
    integer, public :: xl_reg = 1482, xr_reg = 1881, yd_reg = 1032, yu_reg = 1452
    ! Domain boundary indices relative to the input files
    integer, public :: xl, xr, yd, yu
    ! Number of x, y cells, number of years
    integer, public :: nx, ny
    ! x, y, z index for display
    integer, public :: xi_dp, yi_dp, zi_dpst, zi_dped, ti_dp

    !!--------------------------------------------------------------------------
    ! * Grid files info (with full path included. Domain is global for all files) *
    ! Basic 2-D grid info
    character(len = 300), public :: fn_grid = "./ocn_static_grid.nc"
    ! Dzt and Dzu calculated offline
    character(len = 300), public :: fn_dz = "./ocn_static_dz.nc"
    ! Model constants (omega, grav etc.)
    character(len = 300), public :: fn_cons = "./ocn_constant.nc"
endmodule
