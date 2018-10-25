Program main
    use params
    use control
    use zeta, only : load_const, calc_zeta, init_zetavars_input, init_zetavars_output, &
                     release_zetavars_input, release_zetavars_output
    use netcdf

    implicit none
    ! integer :: iyr, imn, ida
    integer :: iyr, idoy, ida, isec
    character(len = 300) :: fn_zeta, fn_error, fn_zetam, fn_zetae, yyyymmdd*20
    integer :: cmode, nyr, ndoy, slyr, sldoy, nsec
    integer :: ncid_zeta, ncid_error, ncid_zetam, iostat
      real :: timer_st, timer_ed

!---------------------Initialization--------------------------------------------
    print*, 'Loading parameters ...'
    call load_params(cmode, nyr, ndoy, slyr, sldoy, nsec)
    print*, ' '

    print*, 'Loading grids and constants ...'
    call load_const()
    print*, ' '
!---------------------Main loop-------------------------------------------------
    if (ifcurl .or. ifdecomp) then
        print*, 'Starting looping ...'
        do iyr = 1, nyr, slyr
            do idoy = 1, ndoy, sldoy
                ! Initializing input/outut fields used by zeta module'
                call init_zetavars_input(cmode)
                call init_zetavars_output(cmode)

                ! Loading variable fields from input files
                ! call loadave_vars_sf(cmode, (/ yrlist(iyr) /), (/ dalist(ida) + sum(eom(1:mnlist(imn)-1)) /))
                  call cpu_time(timer_st)
                call loadave_vars_sf(cmode, yrlist(iyr:iyr+slyr-1), doylist(idoy:idoy+sldoy-1))
                  timer_st = timer_ed; call cpu_time(timer_ed)
                  print '("Time = ",f10.2," seconds.")', timer_ed - timer_st

                ! ! Output file name and creatation
                ! call get_yyyymmdd(yrlist(iyr), doylist(idoy), yrnm_clm, avnm_clm, fn_out_dlm, yyyymmdd)
                ! write(fn_zeta, '(A, A, A, A, A)') &
                !     trim(fn_out_dir), trim(fn_out_pfx), trim(yyyymmdd), trim(fn_out_sfx), '.nc'
                ! if (cmode < 0) then
                !     write(fn_error, '(A, A, A, A, A, A)') &
                !         trim(fn_out_dir), trim(fn_out_pfx), 'decompErrs_', trim(yyyymmdd), trim(fn_out_sfx), '.nc'
                ! endif
                ! print*, 'Output files: ', trim(fn_zeta)
                ! print*, '  and ', trim(fn_error)
                !
                ! ! call get_zeta_fn(iyr, imn, ida, yrnm_clm, yyyymmdd)
                !
                ! call create_outputfiles(cmode, fn_zeta, fn_error, ncid_zeta, ncid_error)
                !
                ! Calculation
                call calc_zeta(cmode)
                !
                ! ! Output file name and creatation
                ! call write_outputfiles(cmode, ncid_zeta, ncid_error)
                !
                ! ! Close output file
                ! call close_outputfiles(cmode, ncid_zeta, ncid_error)

                ! Release working variables
                call release_zetavars_input(cmode)
                call release_zetavars_output(cmode)

                print*, "Finished ", yyyymmdd
            enddo
        enddo
    endif

! !--------------------Mean/Eddy decomposition------------------------------------
!     ! check if mean/eddy calculation is needed
!     if (ifmeaneddy) then
!         print*, 'Starting mean eddy decomposition ...'
!         if (ifmeanclm) then
!            slyr = nyr
!            write(yrnm_clm, '(I0.4, A, I0.4)') yrlist(1), '-', yrlist(nyr)
!         endif
!         print*, 'yrnm_clm', yrnm_clm
!         do iyr = 1, nyr, slyr
!             do isec = 1, nsec, 1
!                 call get_mmdd(sec(isec, 1), trim(meanfreq), avnm_clm)
!                 print*, 'avnm_clm', avnm_clm
!
!                 call get_yyyymmdd(yrlist(iyr), 0, yrnm_clm, avnm_clm, fn_out_dlm, yyyymmdd)
!                 write(fn_zetam, '(A, A, A, A, A)') &
!                     trim(fn_outm_dir), trim(fn_out_pfx), trim(yyyymmdd), trim(fn_out_sfx), '_m.nc'
!                 print*, 'Using file ', trim(fn_zetam)
!
!                 iostat = nf90_open(trim(fn_zetam), NF90_NOWRITE, ncid_zetam)
!
!                 call init_zetavars_output(3)
!                 if (iostat /= nf90_noerr) then
!                     print*, 'Mean file not found! Creating ...'
!                     call init_zetavars_input(3)
!                     call loadave_vars_sf(3, yrlist(iyr:iyr+slyr-1), (/ (ida, ida = sec(isec,1), sec(isec,2)) /))
!                     call create_outputfiles(3, fn_zetam, '', ncid_zetam, ncid_error)
!                     call calc_zeta(3)
!                     call write_outputfiles(3, ncid_zetam, 0)
!                     call close_outputfiles(3, ncid_zetam, 0)
!                     call release_zetavars_input(3)
!                 else
!                     print*, 'Loading mean file ...'
!                     call load_mean_sf(fn_zetam)
!                 endif
!
!                 call init_zetavars_input(1)
!                 allocate(advuT(nx, ny, nz), advvT(nx, ny, nz), advwT(nx, ny, nz), &
!                     advVxT(nx, ny, nz), advVyT(nx, ny, nz), advVzT(nx, ny, nz), &
!                     curlmetT(nx, ny, nz), err_nldecompT(nx, ny, nz))
!
!                 call loadave_zeta_sf(yrlist(iyr:iyr+slyr-1), (/ (ida, ida = sec(isec,1), sec(isec,2)) /))
!                 !
!                 !
!                 !
!                 !
!                 ! ! call create_outputfiles(meaneddy, fn_zeta, fn_error)
!                 ! ! call write_outputfiles(3)
!                 ! ! call close_outputfiles(3)
!                 !
!                 call release_zetavars_output(3)
!                 deallocate(advuT, advvT, advwT, advVxT, advVyT, advVzT, curlmetT, err_nldecompT)
!             enddo
!         enddo
!     endif
endprogram
