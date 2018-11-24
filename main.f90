Program main
    use params
    use io
    use zeta, only : load_const, calc_zeta, init_zetavars_input, init_zetavars_output, &
                     release_zetavars_input, release_zetavars_output
    use netcdf
    implicit none
    integer :: iyr, it, ida, isec, slyr
    character(len = 20) :: yyyymmdd
    integer :: cmode, cmode_m
    integer :: ncid, iostat

!---------------------Initialization--------------------------------------------
    print*, 'Loading parameters ...'
    ! call load_params(cmode, nyr, ndoy, slyr, sldoy, nsec)
    call load_params(cmode, cmode_m)

    print*, 'Loading grids and constants ...'
    call load_const()
    print*, ' '

!--------------------Calculate Mean---------------------------------------------
    if (cmode_m /= -1) then
        write(*, '(A)') '  Calculating mean ...'
        if (ifmeanclm) then
           slyr = T%nyr
           write(T%menm_clm, '(I0.4, A, I0.4)') T%yrlist(1), '-', T%yrlist(T%nyr)
        else
           slyr = 1
        endif
        do iyr = 1, T%nyr, slyr
            do isec = 1, T%nsec, 1
                call get_yyyymmdd(T%yrlist(iyr), 0, T%menm_clm, T%meannm(isec), fn_vorm%dlm, yyyymmdd)
                write(fn_vorm%fn, '(A, A, A, A, A)') &
                    trim(fn_vorm%dir), trim(fn_vorm%pfx), trim(yyyymmdd), trim(fn_vorm%sfx), '.nc'
                write(*, '(A, A)') '  Using file ', trim(fn_vorm%fn)

                iostat = nf90_open(trim(fn_vorm%fn), NF90_NOWRITE, ncid)
                if (iostat /= nf90_noerr) then
                    print*, 'Mean file not found! Creating ...'
                    call init_zetavars_input(cmode_m)
                    call loadave_mom_sf(cmode_m, T%yrlist(iyr:iyr+slyr-1), (/ (ida, ida = T%seclist(isec,1), T%seclist(isec,2)) /), '', '', fn_mom)
                    call init_zetavars_output(cmode_m)
                    call calc_zeta(cmode_m)

                    call output_sf(cmode_m, fn_vorm%fn, '')
                    call release_zetavars_input(cmode_m)
                    call release_zetavars_output(cmode_m)
                endif
           enddo
        enddo
     endif

!---------------------Main loop-------------------------------------------------
    print*, 'Starting looping ...'
    do iyr = 1, T%nyr, 1
        do it = 1, T%nt, 1
            call init_zetavars_output(cmode)
            if (cmode /= -1) then
                ! Initializing input/outut fields used by zeta module'
                call init_zetavars_input(cmode)

                ! Loading variable fields from input files
                call loadave_mom_sf(cmode, T%yrlist(iyr:iyr+1-1), (/ (ida, ida = T%tlist(it,1), T%tlist(it,2)) /), T%yrnm_clm, T%avnm_clm, fn_mom)

                ! Calculation
                call calc_zeta(cmode)
                call release_zetavars_input(cmode)

                if (T%tlist(it, 1) == T%tlist(it, 2)) then
                    ! Output file name and creatation
                    call get_yyyymmdd(T%yrlist(iyr), T%tlist(it,1), T%yrnm_clm, T%avnm_clm, fn_vor%dlm, yyyymmdd)
                    write(fn_vor%fn, '(A, A, A, A, A)') &
                        trim(fn_vor%dir), trim(fn_vor%pfx), trim(yyyymmdd), trim(fn_vor%sfx), '.nc'
                    if (cmode < 0) then
                        write(fn_error%fn, '(A, A, A, A, A, A)') &
                            trim(fn_vor%dir), trim(fn_vor%pfx), 'decompErrs_', trim(yyyymmdd), trim(fn_vor%sfx), '.nc'
                    endif
                    print*, 'Output files: ', trim(fn_vor%fn)
                    print*, '  and ', trim(fn_error%fn)

                    call output_sf(cmode, fn_vor%fn, fn_error%fn)

                    ! call create_output(cmode, fn_zeta, fn_error, ncid_zeta, ncid_error, varids)
                    !
                    ! ! Output file name and creatation
                    ! call write_output(cmode, ncid_zeta, ncid_error, varids)
                    !
                    ! ! Close output file
                    ! call close_output(cmode, ncid_zeta, ncid_error)

                    ! Release working variables
                    call release_zetavars_output(cmode)
                    print*, "Finished ", yyyymmdd
                    cycle
                endif
            else
                call loadave_vor_sf(cmode, T%yrlist(iyr:iyr+1-1), (/ (ida, ida = T%tlist(it,1), T%tlist(it,2)) /), '', '', fn_vor)
            endif

            call init_zetavars_mean(cmode_m)
            call load_mean_sf(cmode_m, fn_vorm%fn)
            call output_me(cmode_m, fn_vore%fn)

            call release_zetavars_mean(cmode_m)
            call release_zetavars_output(cmode)
        enddo
    enddo
endprogram
