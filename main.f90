Program main
    use params
    use io
    use zeta, only : load_const, calc_zeta, init_zetavars_input, init_zetavars_output, &
                     release_zetavars_input, release_zetavars_output
    implicit none
    integer :: iyr, it, ida, isec, slyr, slsec
    character(len = 20) :: yyyymmdd
    character(len=10) :: func, func_m, func_me
    logical :: meanfilestat

!---------------------Initialization--------------------------------------------
    ! print*, 'Loading parameters ...'
    ! call load_params(cmode, nyr, ndoy, slyr, sldoy, nsec)
    call load_params(func, func_m, func_me)
    write(*, '(A, A)') "func = ", func
    write(*, '(A, A)') "func_m = ", func_m
    write(*, '(A, A)') "func_me = ", func_me

    ! print*, 'Loading grids and constants ...'
    call load_const()

!--------------------Calculate Mean---------------------------------------------
    write(*, *)
    write(*, '(A)') '-----------------------------------------------------'
    if (trim(func_m) /= "") then
        write(*, '(A)') '  Checking mean files ...'
        if (ifmeanclm) then
            slyr = T%nyr
            if (T%menm_clm=="") then
                write(T%menm_clm, '(I0.4, A, I0.4)') T%yrlist(1), '-', T%yrlist(T%nyr)
            endif
        else
           slyr = 1
        endif
        do iyr = 1, T%nyr, slyr
            do isec = 1, T%nsec, 1
                call get_yyyymmdd(T%yrlist(iyr), 0, T%menm_clm, T%meannm(isec), fn_vorm%dlm, yyyymmdd)
                write(fn_vorm%fn, '(A, A, A, A, A)') &
                    trim(fn_vorm%dir), trim(fn_vorm%pfx), trim(yyyymmdd), trim(fn_vorm%sfx), '.nc'
                write(*, '(A, A)') '  File ', trim(fn_vorm%fn)

                call check_meanfile(func_m, trim(fn_vorm%fn), meanfilestat)
                if (.not. meanfilestat) then
                    write(*, '(A)') '    Mean file not found! Creating ...'
                    call init_zetavars_input(func_m)
                    call loadave_mom_sf(func_m, T%yrlist(iyr:iyr+slyr-1), (/ (ida, ida = T%seclist(isec,1), T%seclist(isec,2)) /), '', '', fn_mom)
                    call init_zetavars_output(func_m)
                    call calc_zeta(func_m)

                    call output_sf(func_m, fn_vorm%fn, '')
                    call release_zetavars_input(func_m)
                    call release_zetavars_output(func_m)
                else
                    write(*, '(A)') '   Status OK!'
                endif
         enddo
        enddo
    endif

!---------------------Main loop-------------------------------------------------
    write(*, '(A)') '-----------------------------------------------------'
    write(*, '(A)') 'Start the main loop'
    do iyr = 1, T%nyr, 1
        do isec = 1, T%nsec, 1
            if (index(func, "d") /=0 ) then
                slsec = 1
            else
                slsec = T%seclist(isec, 2) - T%seclist(isec, 1) + 1
            endif

            if (trim(func) /= "") then
                do it = T%seclist(isec, 1), T%seclist(isec, 2), slsec
                    ! Initializing input/outut fields used by zeta module'
                    call init_zetavars_input(func = func)
                    call init_zetavars_output(func = func)

                    ! Loading variable fields from input files
                    call loadave_mom_sf(func, T%yrlist(iyr:iyr+1-1), &
                                         (/ (ida, ida = it, T%seclist(isec,1)+slsec-1) /), &
                                         T%yrnm_clm, T%avnm_clm, fn_mom)

                    ! Calculation
                    call calc_zeta(func = func)
                    call release_zetavars_input(func = func)

                    if (slsec == 1) then
                        ! Output file name and creatation
                        call get_yyyymmdd(T%yrlist(iyr), it, &
                                          T%yrnm_clm, T%avnm_clm, fn_vor%dlm, yyyymmdd)
                        write(fn_vor%fn, '(A, A, A, A, A)') &
                              trim(fn_vor%dir), trim(fn_vor%pfx), trim(yyyymmdd), trim(fn_vor%sfx), '.nc'
                        write(*, *)
                        write(*, '(A)') '-----------------------------------------------------'
                        write(*, '(A, A)') 'Outputing file(s): ', trim(fn_vor%fn)
                        if (index(func, "-") /= 0) then
                            write(fn_error%fn, '(A, A, A, A, A, A)') &
                                  trim(fn_vor%dir), trim(fn_vor%pfx), 'decompErrs_', trim(yyyymmdd), trim(fn_vor%sfx), '.nc'
                            write(*, '(A, A)') '  and ', trim(fn_error%fn)
                        endif

                        ! Output files
                        call output_sf(func, fn_vor%fn, fn_error%fn)
                          ! call create_output(cmode, fn_zeta, fn_error, ncid_zeta, ncid_error, varids)
                          !
                          ! ! Output file name and creatation
                          ! call write_output(cmode, ncid_zeta, ncid_error, varids)
                          !
                          ! ! Close output file
                          ! call close_output(cmode, ncid_zeta, ncid_error)

                        ! Release working variables
                        call release_zetavars_output(func = func)
                        write(*, '(A, A)') "Finished ", yyyymmdd
                    endif
                 enddo
            endif

            if (trim(func_me) /= "" ) then
                if (slsec == 1 .or. trim(func) == "") then
                    call init_zetavars_output(func = func_me)
                    call loadave_vor_sf(func_me, T%yrlist(iyr:iyr+1-1), &
                                        (/ (ida, ida = T%seclist(isec,1), T%seclist(isec,2)) /), &
                                        '', '', fn_vor)
                endif
                call get_yyyymmdd(T%yrlist(iyr), 0, T%menm_clm, T%meannm(isec), fn_vorm%dlm, yyyymmdd)
                write(fn_vorm%fn, '(A, A, A, A, A)') &
                      trim(fn_vorm%dir), trim(fn_vorm%pfx), trim(yyyymmdd), trim(fn_vorm%sfx), '.nc'
                write(*, '(A, A)') '  Using file: ', trim(fn_vorm%fn)
                call init_zetavars_mean(func_m)
                call load_mean_sf(func_m, fn_vorm%fn)

                call get_yyyymmdd(T%yrlist(iyr), 0, "", T%meannm(isec), fn_vore%dlm, yyyymmdd)
                write(fn_vore%fn, '(A, A, A, A, A)') &
                      trim(fn_vore%dir), trim(fn_vore%pfx), trim(yyyymmdd), trim(fn_vore%sfx), '.nc'
                write(*, '(A, A)') '  Outputing file: ', trim(fn_vore%fn)
                call output_me(func_m, fn_vore%fn)

                call release_zetavars_mean(func_m)
                call release_zetavars_output(func = func_me)
           endif
       enddo
    enddo
endprogram
