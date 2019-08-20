Program main
    use params
    use io
    use zeta, only : load_const, calc_zeta, init_zetavars_input, init_zetavars_output, &
                     release_zetavars_input, release_zetavars_output
    use netcdf
    implicit none
    integer :: iyr, it, ida, isec, slyr
    character(len = 20) :: yyyymmdd
    character(len=10):: func, func_m, func_me
    integer :: ncid, iostat

!---------------------Initialization--------------------------------------------
    print*, 'Loading parameters ...'
    call load_params(func, func_m, func_me)
    func = 'adm'

    print*, 'Loading grids and constants ...'
    call load_const()
    print*, ' '

    call init_zetavars_output(func=func)

    call init_zetavars_input(func=func)

    call loadave_mom_sf(func, (/2009/), (/ 1 /), T%yrnm_clm, T%avnm_clm, fn_mom)

    call calc_zeta(func=func)

    call release_zetavars_input(func=func)

    write(fn_vor%fn, '(A, A, A, A, A)') &
        trim(fn_vor%dir), trim(fn_vor%pfx), 'test', trim(fn_vor%sfx), '.nc'
    print*, 'Output files: ', trim(fn_vor%fn)
    print*, '  and ', trim(fn_error%fn)

    ! call output_sf(cmode, fn_vor%fn, fn_error%fn)
endprogram
