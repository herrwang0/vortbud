! Functions from POP model
module popfun
  use params, only : kd_r
  private
  public :: u2t, t2u, zcurl

contains
function u2t(uvar, uarea, tarea)
  implicit none
  real(kind=kd_r), dimension(:, :), intent(in) :: uvar, uarea, tarea
  real(kind=kd_r), dimension(:, :), allocatable :: u2t
  integer :: nx, ny

  nx = size(uvar, 1)
  ny = size(uvar, 2)

  allocate(u2t(nx, ny))
  u2t = 0.

  u2t(2:nx, 2:ny) = uvar(1:nx-1, 1:ny-1) * uarea(1:nx-1, 1:ny-1) + &
                    uvar(1:nx-1, 2:ny  ) * uarea(1:nx-1, 2:ny  ) + &
                    uvar(2:nx  , 1:ny-1) * uarea(2:nx  , 1:ny-1) + &
                    uvar(2:nx  , 2:ny  ) * uarea(2:nx  , 2:ny  )

  u2t = u2t / tarea / 4.
  return
endfunction u2t

function t2u(tvar, tarea, uarea)
  implicit none
  real(kind=kd_r), dimension(:, :), intent(in) :: tvar, tarea, uarea
  real(kind=kd_r), dimension(:, :), allocatable :: t2u
  integer :: nx, ny

  nx = size(tvar, 1)
  ny = size(tvar, 2)

  allocate(t2u(nx, ny))
  t2u = 0.

  t2u(1:nx-1, 1:ny-1) = tvar(1:nx-1, 1:ny-1) * tarea(1:nx-1, 1:ny-1) + &
                        tvar(1:nx-1, 2:ny  ) * tarea(1:nx-1, 2:ny  ) + &
                        tvar(2:nx  , 1:ny-1) * tarea(2:nx  , 1:ny-1) + &
                        tvar(2:nx  , 2:ny  ) * tarea(2:nx  , 2:ny  )

  t2u = t2u / uarea / 4.
  return
endfunction t2u

! This is slightly different from the POP code. Instead of specifiy dxu and dyu as
! the weight function, the more general form of the function here allows dzu to be
! included in the weight function.
function zcurl(varu, varv, wgtu, wgtv, wgt)
  use derives, only : dd_xsw, dd_ysw
  implicit none

  real(kind=kd_r), dimension(:, :), intent(in) :: varu, varv, wgtu, wgtv, wgt
  real(kind=kd_r), dimension(:, :), allocatable :: zcurl
  integer :: nx, ny

  nx = size(varu, 1)
  ny = size(varu, 2)

  allocate(zcurl(nx, ny))

  zcurl = dd_xsw(varv * wgtv, wgt) - dd_ysw(varu * wgtu, wgt)
  return
endfunction zcurl
endmodule
