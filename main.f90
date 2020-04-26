Program main
    use params
    use io
    use zeta, only : load_const, calc_zeta, init_zetavars_input, init_zetavars_output, &
                     release_zetavars_input, release_zetavars_output
    implicit none
    integer :: iyr, it, ida, isec, slyr, slsec
    character(len = 20) :: yyyymmdd
    character(len=10) :: func, func_m, func_me
    logical :: meanfilestat, ifsavesf

!---------------------Initialization--------------------------------------------
    ! Loading parameters
    call load_params(func, func_m, func_me)
    ! Loading grids and constants
    call load_const()

!--------------------Calculate Mean---------------------------------------------
    if (trim(func_m) /= "") then
        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Checking mean files'
        ! Check if using climatology (ifmeanclm = .True.) or annual average
        if (ifmeanclm) then
            slyr = T%nyr
            ! Deriving "YYYY-YYYY" in the filenames if not provided 
            if (T%menm_clm=="") then
                write(T%menm_clm, '(I0.4, A, I0.4)') T%yrlist(1), '-', T%yrlist(T%nyr)
            endif
        else
           slyr = 1
        endif
        do iyr = 1, T%nyr, slyr
            do isec = 1, T%nsec, 1
                ! Get filename (fn_vorm%fn) for mean fields 
                call get_yyyymmdd(T%yrlist(iyr), 0, T%menm_clm, T%meannm(isec), fn_vorm%dlm, yyyymmdd)
                fn_vorm%fn = trim(fn_vorm%dir) // trim(fn_vorm%pfx) // trim(yyyymmdd) // trim(fn_vorm%sfx) // '.nc'
                write(*, '(A, A)') '  Looking for ', trim(fn_vorm%fn)

                call check_meanfile(func_m, trim(fn_vorm%fn), meanfilestat)
                if (.not. meanfilestat) then
                    write(*, '(A)') '    Not found! Create now'
                    call init_zetavars_input(func_m)
                    call loadave_mom_sf(func_m, T%yrlist(iyr:iyr+slyr-1), &
                                        (/ (ida, ida = T%seclist(isec,1), T%seclist(isec,2)) /), &
                                        T%yrnm_clm, T%avnm_clm, fn_mom)
                    call init_zetavars_output(func_m)
                    call calc_zeta(func_m)
                    call output_sf(fn_vorm%fn)
                    call release_zetavars_input(func_m)
                    call release_zetavars_output()
                else
                    write(*, '(A)') '    Found and status OK!'
                endif
            enddo
        enddo
    endif

!---------------------Main loop-------------------------------------------------
    write(*, '(A)') '-----------------------------------------------------'
    write(*, '(A)') '-----------------------------------------------------'
    write(*, '(A)') 'Start the main loop'
    do iyr = 1, T%nyr, 1
        do isec = 1, T%nsec, 1
            if (trim(func) /= "") then
                ! For climatology input, we average over the section first. Then calculate curl,
                !  calculate nonlinear term offline as the mean, calculate the difference as eddy
                !  and output after each section
                ! For non-clim and decomposing nonlinear, we need to calculate at each timestep,
                !  and save them for m/e decomposition 
                ! For non-clim and no decompose and m/e, we can again first average
                if (trim(T%yrnm_clm) /= "") then
                    slsec = T%seclist(isec, 2) - T%seclist(isec, 1) + 1
                    ifsavesf = .True.
                else
                    if (trim(func_me) /= "" .and. index(func, "d") == 0) then
                        slsec = T%seclist(isec, 2) - T%seclist(isec, 1) + 1
                        ifsavesf = .False.
                    else
                        slsec = 1
                        ifsavesf = .True.
                    endif
                endif
                do it = T%seclist(isec, 1), T%seclist(isec, 2), slsec
                    ! Initializing input/outut fields used by zeta module
                    call init_zetavars_input(func = func)
                    call init_zetavars_output(func = func)
                    ! Loading (and average) variable fields from input files
                    call loadave_mom_sf(func, T%yrlist(iyr:iyr+1-1), &
                                         (/ (ida, ida = it, it+slsec-1) /), &
                                         T%yrnm_clm, T%avnm_clm, fn_mom)
                    ! Calculation
                    call calc_zeta(func = func)
                    call release_zetavars_input(func)

                    ! Output happens in three cases: 
                    !   1) ifdecomp is True. 
                    !   2) Curl only and each section has just one timestep.
                    !   3) Input file is climatology
                    if (ifsavesf) then 
                        ! Get output filename
                        call get_yyyymmdd(T%yrlist(iyr), it, &
                                          T%yrnm_clm, T%avnm_clm, fn_vor%dlm, yyyymmdd)
                        fn_vor%fn = trim(fn_vor%dir) // trim(fn_vor%pfx) // trim(yyyymmdd) // trim(fn_vor%sfx) // '.nc'
                        ! Output fields
                        call output_sf(fn_vor%fn)
                        ! Release working variables
                        call release_zetavars_output()
                    endif
                 enddo
            endif

            if (trim(func_me) /= "" ) then
                ! Load and average previously calculated vorticity fields if needed
                if (slsec == 1 .or. trim(func) == "") then
                    call init_zetavars_output(func = func_me)
                    call loadave_vor_sf(func_me, T%yrlist(iyr:iyr+1-1), &
                                        (/ (ida, ida = T%seclist(isec,1), T%seclist(isec,2)) /), &
                                        '', '', fn_vor)
                endif
                ! Load mean field and Output mean/eddy file
                call get_yyyymmdd(T%yrlist(iyr), 0, T%menm_clm, T%meannm(isec), fn_vorm%dlm, yyyymmdd)
                fn_vorm%fn = trim(fn_vorm%dir) // trim(fn_vorm%pfx) // trim(yyyymmdd) // trim(fn_vorm%sfx) // '.nc'
                call get_yyyymmdd(T%yrlist(iyr), 0, "", T%meannm(isec), fn_vore%dlm, yyyymmdd)
                fn_vore%fn = trim(fn_vore%dir) // trim(fn_vore%pfx) // trim(yyyymmdd) // trim(fn_vore%sfx) // '.nc'
                call output_me(fn_vorm%fn, fn_vore%fn)
                call release_zetavars_output()
           endif
       enddo
    enddo
endprogram
