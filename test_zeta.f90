Program main
    use params
    use io
    use zeta, only : load_const, calc_zeta, init_zetavars_input, init_zetavars_output, &
                     release_zetavars_input, release_zetavars_output
    use netcdf
    implicit none
    integer :: iyr, it, ida, isec, slsec
    character(len = 20) :: yyyymmdd
    character(len=10):: func, func_m, func_me
    integer :: ncid, iostat

    call load_params(func, func_m, func_me)
    func = "cadme"

    call load_const()

    iyr = 1
    isec = 1
    slsec = 1
    it = T%seclist(isec, 1)
    ! Initializing input/outut fields used by zeta module'
    call init_zetavars_input(func = func)
    call init_zetavars_output(func = func)
    ! Loading (and average) variable fields from input files
    call loadave_mom_sf(func, T%yrlist(iyr:iyr+1-1), &
                            (/ (ida, ida = it, it+slsec-1) /), &
                            T%yrnm_clm, T%avnm_clm, fn_mom)

    ! Calculation
    call calc_zeta(func = func)
    call release_zetavars_input(func = func)
    call release_zetavars_input(func=func)

    ! Output file name and creatation
    call get_yyyymmdd(T%yrlist(iyr), it, &
                        T%yrnm_clm, T%avnm_clm, fn_vor%dlm, yyyymmdd)
    write(fn_vor%fn, '(A, A, A, A, A)') &
            trim(fn_vor%dir), "test_", trim(yyyymmdd), trim(fn_vor%sfx), '.nc'
    write(*, *)
    write(*, '(A)') '  ---------------------------------------------------'
    write(*, '(A, A)') '  Outputing file(s): ', trim(fn_vor%fn)

    ! Output files
    call output_sf(func, fn_vor%fn)
    ! Release working variables
    call release_zetavars_output(func = func)
    write(*, '(A, A)') "  Finished ", yyyymmdd
endprogram
