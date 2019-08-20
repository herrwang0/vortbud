module derives
use params, only : kd_r
private
public :: ddx, ddy, ddx_w, ddy_s, ddx_chain, ddy_chain, &
          ddx_w_chain, ddy_s_chain, meanx_w, meany_s, shiftx_w, shifty_s

contains
! southwest of v (curl_dx)
subroutine ddx(v, wgt1, wgt2, dvdx)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in)   , dimension(:,:) :: v
  real(kind=kd_r), intent(in)   , dimension(:,:) :: wgt1, wgt2
  real(kind=kd_r), intent(inout), dimension(:,:) :: dvdx
  nx = size(v, 1)
  ny = size(v, 2)

  dvdx = 0.
  dvdx(2:nx, 2:ny) = &
    ((v(2:nx  , 2:ny  ) * wgt1(2:nx  , 2:ny  ) + v(2:nx  , 1:ny-1) * wgt1(2:nx  , 1:ny-1)) - &
     (v(1:nx-1, 2:ny  ) * wgt1(1:nx-1, 2:ny  ) + v(1:nx-1, 1:ny-1) * wgt1(1:nx-1, 1:ny-1)))  &
     * 0.5 / wgt2(2:nx, 2:ny)
endsubroutine

! southwest of v (curl_dy)
subroutine ddy(v, wgt1, wgt2, dvdy)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in)   , dimension(:,:) :: v
  real(kind=kd_r), intent(in)   , dimension(:,:) :: wgt1, wgt2
  real(kind=kd_r), intent(inout), dimension(:,:) :: dvdy
  nx = size(v, 1)
  ny = size(v, 2)

  dvdy = 0.
  dvdy(2:nx, 2:ny) = &
    ((v(2:nx  , 2:ny  ) * wgt1(2:nx  , 2:ny  ) + v(1:nx-1, 2:ny  ) * wgt1(1:nx-1, 2:ny  )) - &
     (v(2:nx  , 1:ny-1) * wgt1(2:nx  , 1:ny-1) + v(1:nx-1, 1:ny-1) * wgt1(1:nx-1, 1:ny-1)))  &
     * 0.5 / wgt2(2:nx, 2:ny)
endsubroutine

! west of v (adv dx)
subroutine ddx_w(v, wgt, dvdx)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in)   , dimension(:,:) :: v
  real(kind=kd_r), intent(in)   , dimension(:,:) :: wgt
  real(kind=kd_r), intent(inout), dimension(:,:) :: dvdx
  nx = size(v, 1)
  ny = size(v, 2)

  dvdx = 0.
  dvdx(2:nx,  :  ) = (v(2:nx,  :  ) - v(1:nx-1,  :    )) / wgt(2:nx,  :  )
endsubroutine

! south of v (adv dy)
subroutine ddy_s(v, wgt, dvdy)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in)   , dimension(:,:) :: v
  real(kind=kd_r), intent(in)   , dimension(:,:) :: wgt
  real(kind=kd_r), intent(inout), dimension(:,:) :: dvdy
  nx = size(v, 1)
  ny = size(v, 2)

  dvdy = 0.
  dvdy( :  , 2:ny) = (v( :  , 2:ny) - v( :    , 1:ny-1)) / wgt( :  , 2:ny)
endsubroutine

subroutine ddx_chain(u1, u10, u2, u20, wgt1, wgt2, u10du2dx, u20du1dx, rrc, rr)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in)   , dimension(:,:) :: u1, u10, u2, u20, wgt1, wgt2
  real(kind=kd_r), intent(inout), dimension(:,:) :: u10du2dx, u20du1dx, rrc, rr
  real(kind=kd_r), allocatable  , dimension(:,:) :: du1u2dx, WORK, ONES, r1, r2

  nx = size(u1, 1)
  ny = size(u1, 2)

  allocate(du1u2dx(nx, ny), WORK(nx, ny), ONES(nx, ny), r1(nx, ny), r2(nx, ny))
  ONES = 1.

  call ddx(u1 * u2, wgt1, wgt2, du1u2dx)

  call ddx(u2, wgt1, wgt2, WORK)
  u10du2dx = u10 * WORK

  call ddx(u1, wgt1, wgt2, WORK)
  u20du1dx = u20 * WORK

  call ddx_chainerror(u1, u10, u2, u20, wgt1, wgt2, r1)

  call ddx(ONES, wgt1, wgt2, WORK)
  r2 = - u10 * u20 * WORK

  rrc = r1 + r2
  rr = du1u2dx - u10du2dx - u20du1dx
endsubroutine


subroutine ddy_chain(u1, u10, u2, u20, wgt1, wgt2, &
                     u10du2dy, u20du1dy, rrc, rr)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in)   , dimension(:,:) :: u1, u10, u2, u20, wgt1, wgt2
  real(kind=kd_r), intent(inout), dimension(:,:) :: u10du2dy, u20du1dy, rrc, rr
  real(kind=kd_r), allocatable  , dimension(:,:) :: du1u2dy, WORK, ONES, r1, r2

  nx = size(u1, 1)
  ny = size(u1, 2)

  allocate(du1u2dy(nx, ny), WORK(nx, ny), ONES(nx, ny), r1(nx, ny), r2(nx, ny))
  ONES = 1.

  call ddy(u1 * u2, wgt1, wgt2, du1u2dy)

  call ddy(u2, wgt1, wgt2, WORK)
  u10du2dy = u10 * WORK

  call ddy(u1, wgt1, wgt2, WORK)
  u20du1dy = u20 * WORK

  call ddy_chainerror(u1, u10, u2, u20, wgt1, wgt2, r1)

  call ddy(ONES, wgt1, wgt2, WORK)
  r2 = - u10 * u20 * WORK

  rrc = r1 + r2
  rr = du1u2dy - u10du2dy - u20du1dy
endsubroutine


subroutine ddx_w_chain(u1, u10, u2, u20, wgt, u10du2dx, u20du1dx, rrc, rr)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in)   , dimension(:,:) :: u1, u10, u2, u20, wgt
  real(kind=kd_r), intent(inout), dimension(:,:) :: u10du2dx, u20du1dx, rrc, rr
  real(kind=kd_r), allocatable  , dimension(:,:) :: du1u2dx, WORK

  nx = size(u1, 1)
  ny = size(u1, 2)

  allocate(du1u2dx(nx, ny), WORK(nx, ny))

  call ddx_w(u1 * u2, wgt, du1u2dx)

  call ddx_w(u2, wgt, WORK)
  u10du2dx = u10 * WORK

  call ddx_w(u1, wgt, WORK)
  u20du1dx = u20 * WORK

  call ddx_w_chainerror(u1, u10, u2, u20, wgt, rrc)

  rr = du1u2dx - u10du2dx - u20du1dx
endsubroutine


subroutine ddy_s_chain(u1, u10, u2, u20, wgt, u10du2dy, u20du1dy, rrc, rr)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in)   , dimension(:,:) :: u1, u10, u2, u20, wgt
  real(kind=kd_r), intent(inout), dimension(:,:) :: u10du2dy, u20du1dy, rrc, rr
  real(kind=kd_r), allocatable  , dimension(:,:) :: du1u2dy, WORK

  nx = size(u1, 1)
  ny = size(u1, 2)

  allocate(du1u2dy(nx, ny), WORK(nx, ny))

  call ddy_s(u1 * u2, wgt, du1u2dy)

  call ddy_s(u2, wgt, WORK)
  u10du2dy = u10 * WORK

  call ddy_s(u1, wgt, WORK)
  u20du1dy = u20 * WORK

  call ddy_s_chainerror(u1, u10, u2, u20, wgt, rrc)

  rr = du1u2dy - u10du2dy - u20du1dy
endsubroutine


subroutine ddx_chainerror(u1, u10, u2, u20, wgt1, wgt2, r1)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in), dimension(:,:) :: u1, u10, u2, u20, wgt1, wgt2
  real(kind=kd_r), intent(inout), dimension(:,:) :: r1
  real(kind=kd_r), dimension(:,:), allocatable :: difu1_ne, difu1_nw, difu1_se, difu1_sw, &
                                       difu2_ne, difu2_nw, difu2_se, difu2_sw

  nx = size(u1, 1)
  ny = size(u1, 2)

  allocate(difu1_ne(nx, ny), difu1_nw(nx, ny), difu1_se(nx, ny), difu1_sw(nx, ny), &
           difu2_ne(nx, ny), difu2_nw(nx, ny), difu2_se(nx, ny), difu2_sw(nx, ny))

  difu1_ne             = u1                 - u10
  difu1_nw(2:nx,  :  ) = u1(1:nx-1,  :)     - u10(2:nx,  :  )
  difu1_se( :  , 2:ny) = u1( :    , 1:ny-1) - u10( :  , 2:ny)
  difu1_sw(2:nx, 2:ny) = u1(1:nx-1, 1:ny-1) - u10(2:nx, 2:ny)

  difu2_ne             = u2                 - u20
  difu2_nw(2:nx,  :  ) = u2(1:nx-1,  :)     - u20(2:nx,  :  )
  difu2_se( :  , 2:ny) = u2( :    , 1:ny-1) - u20( :  , 2:ny)
  difu2_sw(2:nx, 2:ny) = u2(1:nx-1, 1:ny-1) - u20(2:nx, 2:ny)

  r1 = 0.
  r1(2:nx, 2:ny) = &
    ((difu1_ne(2:nx, 2:ny) * difu2_ne(2:nx, 2:ny) * wgt1(2:nx  , 2:ny  ) +  &
      difu1_se(2:nx, 2:ny) * difu2_se(2:nx, 2:ny) * wgt1(2:nx  , 1:ny-1)) - &
     (difu1_nw(2:nx, 2:ny) * difu2_nw(2:nx, 2:ny) * wgt1(1:nx-1, 2:ny  ) +  &
      difu1_sw(2:nx, 2:ny) * difu2_sw(2:nx, 2:ny) * wgt1(1:nx-1, 1:ny-1)))  &
    * 0.5 / wgt2(2:nx, 2:ny)

endsubroutine


subroutine ddy_chainerror(u1, u10, u2, u20, wgt1, wgt2, r1)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in), dimension(:,:) :: u1, u10, u2, u20, wgt1, wgt2
  real(kind=kd_r), intent(inout), dimension(:,:) :: r1
  real(kind=kd_r), dimension(:,:), allocatable :: difu1_ne, difu1_nw, difu1_se, difu1_sw, &
                                       difu2_ne, difu2_nw, difu2_se, difu2_sw

  nx = size(u1, 1)
  ny = size(u1, 2)

  allocate(difu1_ne(nx, ny), difu1_nw(nx, ny), difu1_se(nx, ny), difu1_sw(nx, ny), &
           difu2_ne(nx, ny), difu2_nw(nx, ny), difu2_se(nx, ny), difu2_sw(nx, ny))

  difu1_ne             = u1                 - u10
  difu1_nw(2:nx,  :  ) = u1(1:nx-1,  :)     - u10(2:nx,  :  )
  difu1_se( :  , 2:ny) = u1( :    , 1:ny-1) - u10( :  , 2:ny)
  difu1_sw(2:nx, 2:ny) = u1(1:nx-1, 1:ny-1) - u10(2:nx, 2:ny)

  difu2_ne             = u2                 - u20
  difu2_nw(2:nx,  :  ) = u2(1:nx-1,  :)     - u20(2:nx,  :  )
  difu2_se( :  , 2:ny) = u2( :    , 1:ny-1) - u20( :  , 2:ny)
  difu2_sw(2:nx, 2:ny) = u2(1:nx-1, 1:ny-1) - u20(2:nx, 2:ny)

  r1 = 0.
  r1(2:nx, 2:ny) = &
    ((difu1_ne(2:nx, 2:ny) * difu2_ne(2:nx, 2:ny) * wgt1(2:nx  , 2:ny  ) +  &
      difu1_nw(2:nx, 2:ny) * difu2_nw(2:nx, 2:ny) * wgt1(1:nx-1, 2:ny  )) - &
     (difu1_se(2:nx, 2:ny) * difu2_se(2:nx, 2:ny) * wgt1(2:nx  , 1:ny-1) +  &
      difu1_sw(2:nx, 2:ny) * difu2_sw(2:nx, 2:ny) * wgt1(1:nx-1, 1:ny-1)))  &
    * 0.5 / wgt2(2:nx, 2:ny)

endsubroutine

subroutine ddx_w_chainerror(u1, u10, u2, u20, wgt, r1)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in), dimension(:,:) :: u1, u10, u2, u20, wgt
  real(kind=kd_r), intent(inout), dimension(:,:) :: r1
  real(kind=kd_r), dimension(:,:), allocatable :: difu1_e, difu1_w, difu2_e, difu2_w

  nx = size(u1, 1)
  ny = size(u1, 2)

  allocate(difu1_e(nx, ny), difu1_w(nx, ny), difu2_e(nx, ny), difu2_w(nx, ny))

  difu1_e          = u1            - u10
  difu1_w(2:nx, :) = u1(1:nx-1, :) - u10(2:nx, :)

  difu2_e          = u2            - u20
  difu2_w(2:nx, :) = u2(1:nx-1, :) - u20(2:nx, :)

  r1 = 0.
  r1(2:nx, :) = &
    (difu1_e(2:nx, :) * difu2_e(2:nx, :) - difu1_w(2:nx, :) * difu2_w(2:nx, :)) &
     / wgt(2:nx, :)
endsubroutine

subroutine ddy_s_chainerror(u1, u10, u2, u20, wgt, r1)
  implicit none
  integer :: nx, ny
  real(kind=kd_r), intent(in), dimension(:,:) :: u1, u10, u2, u20, wgt
  real(kind=kd_r), intent(inout), dimension(:,:) :: r1
  real(kind=kd_r), dimension(:,:), allocatable:: difu1_n, difu1_s, difu2_n, difu2_s

  nx = size(u1, 1)
  ny = size(u1, 2)

  allocate(difu1_n(nx, ny), difu1_s(nx, ny), difu2_n(nx, ny), difu2_s(nx, ny))

  difu1_n          = u1            - u10
  difu1_s(:, 2:ny) = u1(:, 1:ny-1) - u10(:, 2:ny)

  difu2_n          = u2            - u20
  difu2_s(:, 2:ny) = u2(:, 1:ny-1) - u20(:, 2:ny)

  r1 = 0.
  r1(:, 2:ny) = &
    (difu1_n(:, 2:ny) * difu2_n(:, 2:ny) - difu1_s(:, 2:ny) * difu2_s(:, 2:ny)) &
     / wgt(:, 2:ny)
endsubroutine

function meanx_w(var)
  implicit none
  real(kind=kd_r), dimension(:, :), intent(in) :: var
  real(kind=kd_r), dimension(:, :), allocatable :: meanx_w
  integer :: nx, ny

  nx = size(var, 1)
  ny = size(var, 2)

  allocate(meanx_w(nx, ny))
  meanx_w = 0.

  meanx_w(2:nx, :) = (var(1:nx-1, :) + var(2:nx, :)) / 2
  return
endfunction meanx_w

function meany_s(var)
  implicit none
  real(kind=kd_r), dimension(:, :), intent(in) :: var
  real(kind=kd_r), dimension(:, :), allocatable :: meany_s
  integer :: nx, ny

  nx = size(var, 1)
  ny = size(var, 2)

  allocate(meany_s(nx, ny))
  meany_s = 0.

  meany_s(:, 2:ny) = (var(:, 1:ny-1) + var(:, 2:ny)) / 2
  return
endfunction meany_s

function shiftx_w(var)
  implicit none
  real(kind=kd_r), dimension(:, :), intent(in) :: var
  real(kind=kd_r), dimension(:, :), allocatable :: shiftx_w
  integer :: nx, ny

  nx = size(var, 1)
  ny = size(var, 2)

  allocate(shiftx_w(nx, ny))
  shiftx_w = 0.
  shiftx_w(2:nx,:) = var(1:nx-1,:)
  return
endfunction shiftx_w


function shifty_s(var)
  implicit none
  real(kind=kd_r), dimension(:, :), intent(in) :: var
  real(kind=kd_r), dimension(:, :), allocatable :: shifty_s
  integer :: nx, ny

  nx = size(var, 1)
  ny = size(var, 2)

  allocate(shifty_s(nx, ny))
  shifty_s = 0.
  shifty_s(:,2:ny) = var(:,1:ny-1)
  return
endfunction shifty_s

endmodule
