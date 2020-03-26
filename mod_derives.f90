module derives
    use params, only : kd_r
    private
    public :: dd_xw, dd_ys, dd_xsw, dd_ysw, dd_xw_chain, dd_ys_chain, &
              dd_xsw_chain, dd_ysw_chain, mean_xw, mean_ys, shift_xe, shift_yn

    contains
    ! Decomposing derivatives of two variables' product with error terms
    subroutine dd_xw_chain(u1, u10, u2, u20, wgt, u10du2, u20du1, rrc, rr)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)    :: u1, u10, u2, u20, wgt
        real(kind=kd_r), dimension(:,:), intent(inout) :: u10du2, u20du1
        real(kind=kd_r), dimension(:,:), intent(inout), optional :: rrc, rr
        real(kind=kd_r), dimension(:,:), allocatable   :: du1u2

        nx = size(u1, 1)
        ny = size(u1, 2)

        u10du2 = u10 * dd_xw(u2, wgt)
        u20du1 = u20 * dd_xw(u1, wgt)

        ! calculate error term
        if (present(rrc)) then
            rrc = dd_xw_chainerror(u1, u10, u2, u20, wgt)
        endif
        ! calculate error term 
        if (present(rr)) then
            allocate(du1u2(nx, ny))
            du1u2 = dd_xw(u1 * u2, wgt)
            rr = du1u2 - u10du2 - u20du1
        endif
    endsubroutine dd_xw_chain

    ! Decomposing derivatives of two variables' product with error terms
    subroutine dd_ys_chain(u1, u10, u2, u20, wgt, u10du2, u20du1, rrc, rr)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)    :: u1, u10, u2, u20, wgt
        real(kind=kd_r), dimension(:,:), intent(inout) :: u10du2, u20du1
        real(kind=kd_r), dimension(:,:), intent(inout), optional :: rrc, rr
        real(kind=kd_r), dimension(:,:), allocatable   :: du1u2

        nx = size(u1, 1)
        ny = size(u1, 2)

        u10du2 = u10 * dd_ys(u2, wgt)
        u20du1 = u20 * dd_ys(u1, wgt)

      ! calculate error term
        if (present(rrc)) then
            rrc = dd_ys_chainerror(u1, u10, u2, u20, wgt)
        endif
        ! calculate error term 
        if (present(rr)) then
            allocate(du1u2(nx, ny))
            du1u2 = dd_ys(u1 * u2, wgt)
            rr = du1u2 - u10du2 - u20du1
        endif
    endsubroutine dd_ys_chain

    ! Decomposing derivatives of two variables' product with error terms
    subroutine dd_xsw_chain(u1, u10, u2, u20, wgt, u10du2, u20du1, rrc, rr)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)    :: u1, u10, u2, u20, wgt
        real(kind=kd_r), dimension(:,:), intent(inout) :: u10du2, u20du1
        real(kind=kd_r), dimension(:,:), intent(inout), optional :: rrc, rr
        real(kind=kd_r), dimension(:,:), allocatable   :: du1u2, ONES

        nx = size(u1, 1)
        ny = size(u1, 2)

        u10du2 = u10 * dd_xsw(u2, wgt)
        u20du1 = u20 * dd_xsw(u1, wgt)

        if (present(rrc)) then
            allocate(ONES(nx, ny))
            ONES = 1.
            rrc = dd_xsw_chainerror_wgt(u1, u10, u2, u20, ONES, wgt) - u10 * u20 * dd_xsw(ONES, wgt)
        endif
        if (present(rr)) then
            allocate(du1u2(nx, ny))
            du1u2 = dd_xsw(u1 * u2, wgt)
            rr = du1u2 - u10du2 - u20du1
        endif
    endsubroutine

    ! Decomposing derivatives of two variables' product with error terms
    subroutine dd_ysw_chain(u1, u10, u2, u20, wgt, u10du2, u20du1, rrc, rr)
      implicit none
      integer :: nx, ny
      real(kind=kd_r), dimension(:,:), intent(in)    :: u1, u10, u2, u20, wgt
      real(kind=kd_r), dimension(:,:), intent(inout) :: u10du2, u20du1
      real(kind=kd_r), dimension(:,:), intent(inout), optional :: rrc, rr
      real(kind=kd_r), dimension(:,:), allocatable   :: du1u2, ONES

      nx = size(u1, 1)
      ny = size(u1, 2)

      u10du2 = u10 * dd_ysw(u2, wgt)
      u20du1 = u20 * dd_ysw(u1, wgt)

      if (present(rrc)) then
          allocate(ONES(nx, ny))
          ONES = 1.
          rrc = dd_ysw_chainerror_wgt(u1, u10, u2, u20, ONES, wgt) - u10 * u20 * dd_ysw(ONES, wgt)
      endif
      if (present(rr)) then
          allocate(du1u2(nx, ny))
          du1u2 = dd_ysw(u1 * u2, wgt)
          rr = du1u2 - u10du2 - u20du1
      endif
    endsubroutine

    ! Error term from decompose the derivitive of two variables's product 
    function dd_xw_chainerror(u1, u10, u2, u20, wgt)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)  :: u1, u10, u2, u20, wgt
        real(kind=kd_r), dimension(:,:), allocatable :: difu1_e, difu1_w, difu2_e, difu2_w, dd_xw_chainerror

        nx = size(u1, 1)
        ny = size(u1, 2)

        allocate(difu1_e(nx, ny), difu1_w(nx, ny), difu2_e(nx, ny), difu2_w(nx, ny), dd_xw_chainerror(nx, ny))
        difu1_e          = u1            - u10
        difu1_w(2:nx, :) = u1(1:nx-1, :) - u10(2:nx, :)

        difu2_e          = u2            - u20
        difu2_w(2:nx, :) = u2(1:nx-1, :) - u20(2:nx, :)

        dd_xw_chainerror = 0.
        dd_xw_chainerror(2:nx, :) = &
          (difu1_e(2:nx, :) * difu2_e(2:nx, :) - difu1_w(2:nx, :) * difu2_w(2:nx, :)) &
          / wgt(2:nx, :)
        return
    endfunction dd_xw_chainerror

    ! Error term from decompose the derivitive of two variables's product 
    function dd_ys_chainerror(u1, u10, u2, u20, wgt)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)  :: u1, u10, u2, u20, wgt
        real(kind=kd_r), dimension(:,:), allocatable :: difu1_n, difu1_s, difu2_n, difu2_s, dd_ys_chainerror

        nx = size(u1, 1)
        ny = size(u1, 2)

        allocate(difu1_n(nx, ny), difu1_s(nx, ny), difu2_n(nx, ny), difu2_s(nx, ny),  dd_ys_chainerror(nx, ny))
        difu1_n          = u1            - u10
        difu1_s(:, 2:ny) = u1(:, 1:ny-1) - u10(:, 2:ny)

        difu2_n          = u2            - u20
        difu2_s(:, 2:ny) = u2(:, 1:ny-1) - u20(:, 2:ny)

        dd_ys_chainerror = 0.
        dd_ys_chainerror(:, 2:ny) = &
          (difu1_n(:, 2:ny) * difu2_n(:, 2:ny) - difu1_s(:, 2:ny) * difu2_s(:, 2:ny)) &
          / wgt(:, 2:ny)
        return
    endfunction dd_ys_chainerror

    ! Error term from decompose the derivitive of two variables's product multiplied by an arbitrary wgt1
    function dd_xsw_chainerror_wgt(u1, u10, u2, u20, wgt1, wgt2)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)  :: u1, u10, u2, u20, wgt1, wgt2
        real(kind=kd_r), dimension(:,:), allocatable :: difu1_ne, difu1_nw, difu1_se, difu1_sw, &
                                                        difu2_ne, difu2_nw, difu2_se, difu2_sw, & 
                                                        dd_xsw_chainerror_wgt

        nx = size(u1, 1)
        ny = size(u1, 2)

        allocate(difu1_ne(nx, ny), difu1_nw(nx, ny), difu1_se(nx, ny), difu1_sw(nx, ny), &
                 difu2_ne(nx, ny), difu2_nw(nx, ny), difu2_se(nx, ny), difu2_sw(nx, ny), &
                 dd_xsw_chainerror_wgt(nx, ny))

        difu1_ne             = u1                 - u10
        difu1_nw(2:nx,  :  ) = u1(1:nx-1,  :)     - u10(2:nx,  :  )
        difu1_se( :  , 2:ny) = u1( :    , 1:ny-1) - u10( :  , 2:ny)
        difu1_sw(2:nx, 2:ny) = u1(1:nx-1, 1:ny-1) - u10(2:nx, 2:ny)

        difu2_ne             = u2                 - u20
        difu2_nw(2:nx,  :  ) = u2(1:nx-1,  :)     - u20(2:nx,  :  )
        difu2_se( :  , 2:ny) = u2( :    , 1:ny-1) - u20( :  , 2:ny)
        difu2_sw(2:nx, 2:ny) = u2(1:nx-1, 1:ny-1) - u20(2:nx, 2:ny)

        dd_xsw_chainerror_wgt = 0.
        dd_xsw_chainerror_wgt(2:nx, 2:ny) = &
          ((difu1_ne(2:nx, 2:ny) * difu2_ne(2:nx, 2:ny) * wgt1(2:nx  , 2:ny  ) +  &
            difu1_se(2:nx, 2:ny) * difu2_se(2:nx, 2:ny) * wgt1(2:nx  , 1:ny-1)) - &
           (difu1_nw(2:nx, 2:ny) * difu2_nw(2:nx, 2:ny) * wgt1(1:nx-1, 2:ny  ) +  &
            difu1_sw(2:nx, 2:ny) * difu2_sw(2:nx, 2:ny) * wgt1(1:nx-1, 1:ny-1)))  &
          * 0.5 / wgt2(2:nx, 2:ny)
        return
    endfunction dd_xsw_chainerror_wgt

    ! Error term from decompose the derivitive of two variables's product multiplied by an arbitrary wgt1
    function dd_ysw_chainerror_wgt(u1, u10, u2, u20, wgt1, wgt2)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)  :: u1, u10, u2, u20, wgt1, wgt2
        real(kind=kd_r), dimension(:,:), allocatable :: difu1_ne, difu1_nw, difu1_se, difu1_sw, &
                                                        difu2_ne, difu2_nw, difu2_se, difu2_sw, & 
                                                        dd_ysw_chainerror_wgt

        nx = size(u1, 1)
        ny = size(u1, 2)

        allocate(difu1_ne(nx, ny), difu1_nw(nx, ny), difu1_se(nx, ny), difu1_sw(nx, ny), &
                difu2_ne(nx, ny), difu2_nw(nx, ny), difu2_se(nx, ny), difu2_sw(nx, ny), &
                dd_ysw_chainerror_wgt(nx, ny))

        difu1_ne             = u1                 - u10
        difu1_nw(2:nx,  :  ) = u1(1:nx-1,  :)     - u10(2:nx,  :  )
        difu1_se( :  , 2:ny) = u1( :    , 1:ny-1) - u10( :  , 2:ny)
        difu1_sw(2:nx, 2:ny) = u1(1:nx-1, 1:ny-1) - u10(2:nx, 2:ny)

        difu2_ne             = u2                 - u20
        difu2_nw(2:nx,  :  ) = u2(1:nx-1,  :)     - u20(2:nx,  :  )
        difu2_se( :  , 2:ny) = u2( :    , 1:ny-1) - u20( :  , 2:ny)
        difu2_sw(2:nx, 2:ny) = u2(1:nx-1, 1:ny-1) - u20(2:nx, 2:ny)

        dd_ysw_chainerror_wgt = 0.
        dd_ysw_chainerror_wgt(2:nx, 2:ny) = &
          ((difu1_ne(2:nx, 2:ny) * difu2_ne(2:nx, 2:ny) * wgt1(2:nx  , 2:ny  ) +  &
            difu1_nw(2:nx, 2:ny) * difu2_nw(2:nx, 2:ny) * wgt1(1:nx-1, 2:ny  )) - &
           (difu1_se(2:nx, 2:ny) * difu2_se(2:nx, 2:ny) * wgt1(2:nx  , 1:ny-1) +  &
            difu1_sw(2:nx, 2:ny) * difu2_sw(2:nx, 2:ny) * wgt1(1:nx-1, 1:ny-1)))  &
          * 0.5 / wgt2(2:nx, 2:ny)
        return
    endfunction dd_ysw_chainerror_wgt

    ! DvDx located to west of v (adv dx), devided by wgt
    function dd_xw(v, wgt)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)  :: v
        real(kind=kd_r), dimension(:,:), intent(in)  :: wgt
        real(kind=kd_r), dimension(:,:), allocatable :: dd_xw

        nx = size(v, 1)
        ny = size(v, 2)

        allocate(dd_xw(nx, ny))
        dd_xw = (v - shift_xe(v)) / wgt
        dd_xw(1, :) = 0.
        ! dd_xw = 0.
        ! dd_xw(2:nx,  :  ) = (v(2:nx,  :  ) - v(1:nx-1,  :    )) / wgt(2:nx,  :  )
        return
    endfunction dd_xw

    ! DvDy located to south of v (adv dy), devided by wgt
    function dd_ys(v, wgt)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)  :: v
        real(kind=kd_r), dimension(:,:), intent(in)  :: wgt
        real(kind=kd_r), dimension(:,:), allocatable :: dd_ys
        nx = size(v, 1)
        ny = size(v, 2)

        allocate(dd_ys(nx, ny))
        dd_ys = (v - shift_yn(v)) / wgt
        dd_ys(:, 1) = 0.
        ! dd_ys = 0.
        ! dd_ys( :  , 2:ny) = (v( :  , 2:ny) - v( :    , 1:ny-1)) / wgt( :  , 2:ny)
        return
    endfunction dd_ys

    ! DvDx located to southwest of v (curl dx), devided by wgt
    function dd_xsw(v, wgt)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)  :: v
        real(kind=kd_r), dimension(:,:), intent(in)  :: wgt
        real(kind=kd_r), dimension(:,:), allocatable :: dd_xsw
        nx = size(v, 1)
        ny = size(v, 2)

        allocate(dd_xsw(nx, ny))
        dd_xsw = mean_ys(v - shift_xe(v)) / wgt
        dd_xsw(1, :) = 0.
        dd_xsw(:, 1) = 0. 
        ! dd_xsw = 0.
        ! dd_xsw(2:nx, 2:ny) = &
        !   ((v(2:nx  , 2:ny  ) + v(2:nx  , 1:ny-1)) - (v(1:nx-1, 2:ny  ) + v(1:nx-1, 1:ny-1))) &
        !    * 0.5 / wgt(2:nx, 2:ny)
        return
    endfunction dd_xsw

    ! DvDy located to southwest of v (curl dy), devided by wgt
    function dd_ysw(v, wgt)
        implicit none
        integer :: nx, ny
        real(kind=kd_r), dimension(:,:), intent(in)  :: v
        real(kind=kd_r), dimension(:,:), intent(in)  :: wgt
        real(kind=kd_r), dimension(:,:), allocatable :: dd_ysw
        nx = size(v, 1)
        ny = size(v, 2)

        allocate(dd_ysw(nx, ny))
        dd_ysw = mean_xw(v - shift_yn(v)) / wgt
        dd_ysw(1, :) = 0.
        dd_ysw(:, 1) = 0. 
        ! dvdy = 0.
        ! dvdy(2:nx, 2:ny) = &
        !   ((v(2:nx  , 2:ny  ) + v(1:nx-1, 2:ny  )) - (v(2:nx  , 1:ny-1) + v(1:nx-1, 1:ny-1)))  &
        !    * 0.5 / wgt(2:nx, 2:ny)
        return
    endfunction dd_ysw

    ! Arithmetic mean along the x direction (for point i, mean_xw(var(i)) is the interpolated value to the west of it) 
    function mean_xw(var)
        implicit none
        real(kind=kd_r), dimension(:, :), intent(in) :: var
        real(kind=kd_r), dimension(:, :), allocatable :: mean_xw
        integer :: nx, ny

        nx = size(var, 1)
        ny = size(var, 2)

        allocate(mean_xw(nx, ny))
        mean_xw = 0.
        mean_xw(2:nx, :) = (var(1:nx-1, :) + var(2:nx, :)) / 2
        return
    endfunction mean_xw

    ! Arithmetic mean along the y direction (for point j, mean_ys(var(j)) is the interpolated value to the south of it) 
    function mean_ys(var)
        implicit none
        real(kind=kd_r), dimension(:, :), intent(in) :: var
        real(kind=kd_r), dimension(:, :), allocatable :: mean_ys
        integer :: nx, ny

        nx = size(var, 1)
        ny = size(var, 2)

        allocate(mean_ys(nx, ny))
        mean_ys = 0.
        mean_ys(:, 2:ny) = (var(:, 1:ny-1) + var(:, 2:ny)) / 2
        return
    endfunction mean_ys

    ! Make the matrix one index higher in the x direction (push it to the east, i-1 --> i)
    function shift_xe(var)
        implicit none
        real(kind=kd_r), dimension(:, :), intent(in) :: var
        real(kind=kd_r), dimension(:, :), allocatable :: shift_xe
        integer :: nx, ny

        nx = size(var, 1)
        ny = size(var, 2)

        allocate(shift_xe(nx, ny))
        shift_xe = 0.
        shift_xe(2:nx,:) = var(1:nx-1,:)
        return
    endfunction shift_xe

    ! Make the matrix one index higher in the y direction (push it to the north, j-1 --> j)
    function shift_yn(var)
        implicit none
        real(kind=kd_r), dimension(:, :), intent(in) :: var
        real(kind=kd_r), dimension(:, :), allocatable :: shift_yn
        integer :: nx, ny

        nx = size(var, 1)
        ny = size(var, 2)

        allocate(shift_yn(nx, ny))
        shift_yn = 0.
        shift_yn(:,2:ny) = var(:,1:ny-1)
        return
    endfunction shift_yn
endmodule
