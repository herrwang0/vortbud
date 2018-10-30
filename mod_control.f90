module control
    use params
    use zeta, only : uc, vc, wc, ssh, advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy, &
                     curlnonl, betav, stretchp, err_cor, curlpgrad, curlhdiff, curlvdiff, res, &
                     curladv, curlmet, err_nlsub, advu, advv, advw, advVx, advVy, advVz, err_nldecomp, &
                     rr_rev, rr_cha, tlat, tlong, z_t, init_zetavars_output
    implicit none
    private
    ! public :: load_params, loadave_vars_sf, get_mmdd, get_yyyymmdd, create_outputfiles, write_outputfiles, close_outputfiles
    public :: load_params, loadave_vars_sf, load_mean_sf, loadave_zeta_sf, get_mmdd, get_yyyymmdd, create_outputfiles, write_outputfiles, close_outputfiles, output_me, init_meaneddy_vars

    !!--------------------------------------------------------------------------
    ! * Input namelist file name *
    character(len = *), parameter, public :: fn_namelist = 'input.nml'

    !!--------------------------------------------------------------------------
    ! * Calculation mode *
    logical, public :: ifcurl, ifdecomp, ifdecompdebug = .False.

    !!--------------------------------------------------------------------------
    ! * Mean field frequency *
    ! Date format
    ! Four options in input.nml (listed by prioity)
    ! (a) meanfreq: period name ("a" for annual, "m" for monthly)
    ! (b) seclist_in_st, seclist_in_ed: list of starting and ending date
    ! (c) secst, seced: starting and ending day (of a year) for a single section
    ! (d) nda_sec: Fixed number of days in each section
    logical, public :: ifmeaneddy, ifmeanclm
    character(len = 5), public :: meanfreq = ""
    integer, dimension(365) :: seclist_in_st = 0, seclist_in_ed = 0
    integer, public :: secst = 0, seced = 0
    integer, public :: nda_sec = 365
    ! list of sections for mean eddy decomposition
    integer, allocatable, public :: seclist(:,:)

    !!--------------------------------------------------------------------------
    ! * Input file info *
    ! Input file directory (suppose in the same directory)
    character(len = 300), public :: fn_in_dir = "/lustre/atlas1/cli115/proj-shared/jritchie/SIO/MODEL_DATA/yelpatch60/IOC/DATA/DAT_IO_all/nday1/"
    ! Input file prefix (before "year")
    character(len = 100), public :: fn_in_pfx = "ia_top_tx0.1_v2_yel_patc_1948_intel.pop.h.nday1."
    ! Input file suffix (after "day", or subregion)
    character(len = 10), public :: fn_in_sfx = ".IO"
    ! Input file date delimiter (between yyyy, mm, dd)
    character(len = 1), public :: fn_in_dlm = '-'
    ! Domain boundary indices of the input file (relative to the entire globe. Here, from 3600x2400 )
    integer, public :: xl_ref = 1489, xr_ref = 1901, yd_ref = 1101, yu_ref = 1401

    !!--------------------------------------------------------------------------
    ! * Date for the main loop *
    ! Date format
    ! Two options in input.nml for yr, mn and da (listed by prioity)
    ! (a): list
    ! (b): starting and ending (st, ed)
    ! For year, additonal option is the name for climatology (yrnm_clm). If none is provide, then fatal error.
    ! For month and days, if neither list or st/ed is provided, then assume 12 months and 31* days
    ! If any negative number appears in mnlist, dalist, then assume average at this level.
    !    For example, mnlist_in = -1 indicates input file is annually averaged (therefore, no mnlist AND dalist needed).
    ! Starting and ending year
    integer, public :: yrst = 2005, yred = 2009
    ! Starting and ending month
    integer, public :: mnst = 1, mned = 12
    ! Starting and ending day
    integer, public :: dast = 1, daed = 31
    ! list of yr
    integer :: yrlist_in(60) = 0
    integer, allocatable, public :: yrlist(:)
    ! list of month
    integer :: mnlist_in(12) = 0 ! default is zero, = 12 months
    integer, allocatable :: mnlist(:)
    ! list of days
    integer :: dalist_in(31) = 0 ! default is zero, = 31 days
    integer, allocatable :: dalist(:)
    ! list of days in year
    integer :: doylist_full(365) = 0
    integer, allocatable, public :: doylist(:)
    ! Year part of the name of climatology file
    character(len=9), public :: yrnm_clm = ""
    ! Generated month/day part of the name if climatology is calculated here
    character(len=9), public :: avnm_clm = ""

    !!----------------------------------------------------------------------------
    ! * output file info *
    ! Output file directory
    character(len = 300), public :: fn_out_dir = "/lustre/atlas/proj-shared/cli115/hewang/data/"
    ! Output file prefix
    character(len = 100), public :: fn_out_pfx = "vort_bud_"
    ! Output file suffix (default is subregion)
    character(len = 10), public :: fn_out_sfx = ""
    ! Output file date delimiter (between yyyy, mm, dd. Default is the same as input)
    character(len = 1), public :: fn_out_dlm = ""

    ! Output file directory (mean)
    character(len = 300), public :: fn_outm_dir = "/lustre/atlas/proj-shared/cli115/hewang/data/"
    ! Output file directory (eddy)
    character(len = 300), public :: fn_oute_dir = "/lustre/atlas/proj-shared/cli115/hewang/data/"

    !!----------------------------------------------------------------------------
    ! netCDF handles
    integer, public :: varid_out_curlnonl , varid_out_betav    , varid_out_stretchp , varid_out_errcor   , &
                       varid_out_curlpgrad, varid_out_curlhdiff, varid_out_curlvdiff, varid_out_res      , &
                       varid_out_curlnonlc, varid_out_errnlsub ,                                           &
                       varid_out_advu     , varid_out_advv     , varid_out_advw     , varid_out_errnldcmp, &
                       varid_out_advVx    , varid_out_advVy    , varid_out_advVz    , varid_out_curlmet  , &
                       varid_out_errnl_rev, varid_out_errnl_cha
    !!----------------------------------------------------------------------------
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
        advuT, advvT, advwT, advVxT, advVyT, advVzT, curlmetT, err_nldecompT
    integer, public :: varid_out_advue    , varid_out_advve    , varid_out_advwe , &
                       varid_out_advVxe   , varid_out_advVye   , varid_out_advVze, &
                       varid_out_curlmete , varid_out_errnldcmpe


    !!----------------------------------------------------------------------------
    namelist /calcmode/ ifcurl, ifdecomp, ifdecompdebug
    namelist /mecmode/ ifmeaneddy, ifmeanclm, meanfreq, nda_sec, seclist_in_st, seclist_in_ed, secst, seced
    namelist /bnds/ xl_reg, xr_reg, yd_reg, yu_reg, subreg, yrst, yred, mnst, mned, dast, daed, &
                    yrnm_clm, yrlist_in, mnlist_in, dalist_in, xi_dp, yi_dp, zi_dpst, zi_dped, ti_dp
    namelist /grid/ fn_grid, fn_dz, fn_cons
    namelist /input/ xl_ref, xr_ref, yd_ref, yu_ref, fn_in_dir, fn_in_pfx, fn_in_sfx, fn_in_dlm
    namelist /output/ fn_out_dir, fn_out_pfx, fn_out_sfx, fn_out_dlm, fn_outm_dir, fn_oute_dir

    contains
    subroutine load_params(cmode, nyr, ndoy, slyr, sldoy, nsec)
        use, intrinsic :: IEEE_ARITHMETIC
        implicit none
        integer, intent(inout) :: cmode, ndoy, nyr, sldoy, slyr, nsec
        integer :: nml_error
        integer :: iyr, imn, nmn, ida, nda, idoy, isec

        open(101, file=fn_namelist, status="old", iostat=nml_error)
        read(101, nml=calcmode, iostat=nml_error)
          print*, "input <calcmode> error: ", nml_error
        read(101, nml=mecmode, iostat=nml_error)
          print*, "input <mecmode> error: ", nml_error
        read(101, nml=bnds, iostat=nml_error)
          print*, "input <bnds> error: ", nml_error
        read(101, nml=grid, iostat=nml_error)
          print*, "input <grid> error: ", nml_error
        read(101, nml=input, iostat=nml_error)
          print*, "input <input> error: ", nml_error
        read(101, nml=output, iostat=nml_error)
          print*, "input <output> error: ", nml_error
        close(101)

        !-------------------------------------------
        ! * Interpreting calculation mode *
        !-------------------------------------------
        print*, '-------------------'
        print*, 'Calculation Modes:'
        print*, "curl", ifcurl
        print*, "decomp", ifdecomp
        print*, "meaneddy", ifmeaneddy

        if     (      ifcurl .and. .not. ifdecomp .and. .not. ifmeaneddy) then
            print*, '-------------------'
            print*, 'Calculating curl of momentum equation only for each input file'
            print*, '(mode 0)'
            print*, '-------------------'
            cmode = 0
        elseif (      ifcurl .and. .not. ifdecomp .and.       ifmeaneddy) then
            print*, '-------------------'
            print*, 'Calculating curl of momentum equation and eddy/mean decomposition of the nonlinear term'
            print*, '(mode 1)'
            if (yrnm_clm == "") then
                print*, "Nonclimatology files are used. Will calculate climatology first" ! ave field then

            else
                print*, "Climatology files are used. ERR_nlsub == curlnonl_eddy"
            endif
            print*, '-------------------'
            cmode = 1
        elseif (      ifcurl .and.       ifdecomp) then
            print*, '-------------------'
            print*, 'Calculating curl of momentum equation and decompositing nonlinear terms'
            print*, '(mode 2)'
            if (yrnm_clm /= "") then
                if (.not. ifmeaneddy) then
                    print*, "Note: A decomposition of nonlinear terms of the climatology &
                             is essentially equivalent to an eddy/mean decomposition."
                else
                    print*, "Climatology files are used. ERR_nlsub == curlnonl_eddy and only mean terms have decompostion."
                endif
                print*, 'ifmeaneddy will be set to F'
                ifmeaneddy = .False.
            endif
            print*, '-------------------'
            cmode = 2
            if (ifdecompdebug) then
               cmode = -2
            endif
        elseif (.not. ifcurl .and.       ifdecomp ) then
            print*, '-------------------'
            print*, 'Calculating nonlinear term decompostion only (mode 3)'
            if (ifmeaneddy) then
                print*, "Warning: Eddy/mean decomposition is not calculated in this mode. "
            endif
            print*, '-------------------'
            cmode = 3
            if (ifdecompdebug) then
               cmode = -3
            endif
        elseif (.not. ifcurl .and. .not. ifdecomp .and.       ifmeaneddy .and. yrnm_clm == "") then
            print*, '-------------------'
            print*, "Eddy/mean decomposition will be calculated based on given decomposed zeta equation files"
            print*, "    Skipping the main loop ..."
            print*, '-------------------'
        else
            print*, '-------------------'
            print*, "Do nothing. Exit"
            print*, '-------------------'
            stop
        endif

        MVALUE = ieee_value(MVALUE, IEEE_QUIET_NAN)
        print*, 'Missing value: ', MVALUE

        !-------------------------------------------
        ! * Relative domain boundaries *
        !-------------------------------------------
        xl = xl_reg - xl_ref + 1
        xr = xr_reg - xl_ref + 1
        yd = yd_reg - yd_ref + 1
        yu = yu_reg - yd_ref + 1

        nx = xr - xl + 1
        ny = yu - yd + 1

        !-------------------------------------------
        ! * Interpreting date for the main loop *
        !-------------------------------------------
        if (len(trim(yrnm_clm)) > 0) then
            nyr = 1
            allocate(yrlist(nyr))
            yrlist = -1
        elseif (count(yrlist_in/=0) /= 0 .and. all(yrlist_in >= 0)) then
            nyr = count(yrlist_in/=0)
            allocate(yrlist(nyr))
            yrlist = yrlist_in(1:nyr)
        elseif (yrst /= 0 .and. yrst /= 0) then
            nyr = yred - yrst + 1
            allocate(yrlist(nyr))
            yrlist = (/ (iyr, iyr = yrst, yred) /)
        else
            print*, "Input missing: either year list or yrnm_clm is needed! Aborting"
            stop
        endif

        if (any(mnlist_in < 0)) then
            nmn = 1
            allocate(mnlist(nmn))
            mnlist = 0
        elseif (count(mnlist_in /= 0) /= 0 .and. all(mnlist_in >= 0)) then
            nmn = count(mnlist_in/=0)
            allocate(mnlist(nmn))
            mnlist = mnlist_in(1:nmn)
        else
            nmn = mned - mnst + 1
            allocate(mnlist(nmn))
            mnlist = (/ (imn, imn = mnst, mned) /)
        endif

        if (any(dalist_in < 0) .or. any(mnlist_in < 0)) then
            nda = 1
            allocate(dalist(nda))
            dalist = 0
        elseif (count(dalist_in /= 0) /= 0 .and. all(dalist_in >= 0)) then
            nda = count(dalist_in/=0)
            allocate(dalist(nda))
            dalist = dalist_in(1:nda)
        else
            nda = daed - dast + 1
            allocate(dalist(nda))
            dalist = (/ (ida, ida = dast, daed) /)
        endif

        idoy = 0
        do imn = 1, nmn, 1
           do ida = 1, nda, 1
               if (mnlist(imn) > 0 .and. dalist(ida) > eom(mnlist(imn))) exit
               idoy = idoy + 1
               if (dalist(ida) == 0 .and. mnlist(imn) /= 0) then
                   doylist_full(idoy) = -(1 + sum(eom(1:mnlist(imn)-1)))
               else
                   doylist_full(idoy) = dalist(ida) + sum(eom(1:mnlist(imn)-1))
              endif
          enddo
        enddo
        ndoy = idoy
        allocate(doylist(ndoy))
        doylist = doylist_full(1:ndoy)

        slyr  = 1
        sldoy = 1

        !-------------------------------------------
        ! * Modified loop and average yr doy list for the special case *
        !-------------------------------------------
        if (yrnm_clm == "" .and. cmode == 1) then
             write(yrnm_clm, '(I0.4, A, I0.4)') yrlist(1), '-', yrlist(nyr)
             call get_mmdd((/abs(doylist(1)), abs(doylist(ndoy))/), '', avnm_clm)

             slyr  = nyr
             sldoy = ndoy
        endif

        !-------------------------------------------
        ! * Interpreting mean/eddy decomposition sections (of a year) *
        !-------------------------------------------
        if (meanfreq=="m") then
            ! nda_sec = 0
            nsec = 12
            allocate(seclist(nsec, 2))
            do isec = 1, nsec
                seclist(isec, 1) = 1 + sum(eom(1:isec)) - eom(isec)
                seclist(isec, 2) = sum(eom(1:isec))
            enddo
        elseif (meanfreq=="a") then
            ! nda_sec = 365
            nsec = 1
            allocate(seclist(nsec, 2))
            seclist(1, 1) = 1
            seclist(1, 2) = 365
        elseif (count(seclist_in_st/=0) /= 0 .and. all(seclist_in_st >= 0) .and. &
                count(seclist_in_ed/=0) /= 0 .and. all(seclist_in_ed >= 0)) then
            nsec = min(count(seclist_in_st/=0), count(seclist_in_ed/=0))
            allocate(seclist(nsec, 2))
            seclist(:, 1) = seclist_in_st(1:nsec)
            seclist(:, 2) = seclist_in_ed(1:nsec)
        elseif (seced >= secst /= 0 .and. secst > 0) then
            nsec = 1
            allocate(seclist(nsec, 2))
            seclist(1, 1) = secst
            seclist(1, 2) = seced
        else
            nsec = ceiling(365. / nda_sec)
            allocate(seclist(nsec, 2))
            do isec = 1, nsec
                seclist(isec, 1) = (isec - 1) * nda_sec + 1
                seclist(isec, 2) = min(isec * nda_sec, 365)
            enddo
        endif

        if (len(trim(fn_out_sfx)) == 0) then
          write(fn_out_sfx, '(A, A)') '_', trim(subreg)
        endif

        if (len(trim(fn_out_dlm)) == 0) then
          write(fn_out_dlm, '(A)') trim(fn_in_dlm)
        endif

        !-------------------------------------------
        ! * Print basic info *
        !-------------------------------------------
        print*, '-------------------'
        print*, "Domain name :", subreg
        print*, "Domain size : ", "nx: ", nx, "ny: ", ny
        print*, '-------------------'
        if (yrnm_clm /= "") then
            print*, "Using or creating a climatology input file averaged over ", yrnm_clm
            print*, '-------------------'
        endif
        print*, "Year: " , yrlist
        ! print*, "Month: ", mnlist
        ! print*, "Day: "  , dalist
        print*, "Days of year: "  , doylist

        print*, "Mean eddy decompsoition : ", ifmeaneddy
        print*, "Sections : ", seclist(:, 1)
        print*, "Sections : ", seclist(:, 2)
    endsubroutine

    subroutine loadave_vars_sf(cmode, yrlist, doylist)
        use ncio, only : nc_read
        ! use popload, only : load_current_day, find_daily_file
        implicit none

        integer, intent(in) :: cmode, yrlist(:), doylist(:)
        character :: fn_in*300, yyyymmdd*20
        integer :: iyr, idoy
        integer :: nn, mm, dd
        real(kind=kd_r), dimension(nx, ny, nz, 1) :: WORK
        real(kind=kd_r), dimension(nx, ny, 1 , 1) :: WORK2

        nn = count(doylist/=0) * count(yrlist/=0)
        print*, '  '
        print*, 'Loading and averaing variables (velocity, mom. terms) from input files over ', nn, ' day(s)'

        do iyr = 1, count(yrlist/=0)
            do idoy = 1, count(doylist/=0)
                call get_yyyymmdd(yrlist(iyr), doylist(idoy), '', '', fn_in_dlm, yyyymmdd)

                write(fn_in, '(A, A, A, A, A)') &
                    trim(fn_in_dir), trim(fn_in_pfx), trim(yyyymmdd), trim(fn_in_sfx), '.nc'

                print*, fn_in

                call nc_read(fn_in, 'UVEL' , WORK, (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); uc = uc + WORK(:, :, :, 1) / nn
                call nc_read(fn_in, 'VVEL' , WORK, (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); vc = vc + WORK(:, :, :, 1) / nn
                if ( abs(cmode) == 0 .or. abs(cmode) == 1 .or. abs(cmode) == 2 ) then
                    call nc_read(fn_in, 'ADVU'   , WORK , (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); advx   = advx   + WORK(:, :, :, 1) / nn
                    call nc_read(fn_in, 'ADVV'   , WORK , (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); advy   = advy   + WORK(:, :, :, 1) / nn
                    call nc_read(fn_in, 'GRADX'  , WORK , (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); gradx  = gradx  + WORK(:, :, :, 1) / nn
                    call nc_read(fn_in, 'GRADY'  , WORK , (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); grady  = grady  + WORK(:, :, :, 1) / nn
                    call nc_read(fn_in, 'HDIFFU' , WORK , (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); hdiffx = hdiffx + WORK(:, :, :, 1) / nn
                    call nc_read(fn_in, 'HDIFFV' , WORK , (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); hdiffy = hdiffy + WORK(:, :, :, 1) / nn
                    call nc_read(fn_in, 'VDIFFU' , WORK , (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); vdiffx = vdiffx + WORK(:, :, :, 1) / nn
                    call nc_read(fn_in, 'VDIFFV' , WORK , (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); vdiffy = vdiffy + WORK(:, :, :, 1)
                    call nc_read(fn_in, 'SSH'    , WORK2, (/xl, yd, 1, 1/), (/nx, ny,  1, 1/)); ssh(:, :, 1) = ssh(:, :, 1) + WORK2(:, :, 1, 1) / nn
                endif
                if ( abs(cmode) == 1 .or. abs(cmode) == 2 .or. abs(cmode) == 3 ) then
                    call nc_read(fn_in, 'WVEL'   , WORK , (/xl, yd, 1, 1/), (/nx, ny, nz, 1/)); wc = wc + WORK(:, :, :, 1) / nn
                endif
            enddo
        enddo

        where(abs(uc) > 1e10) uc = 0.
        where(abs(vc) > 1e10) vc = 0.
        print*, 'uc: ', uc(xi_dp, yi_dp, zi_dpst:zi_dped)
        print*, 'vc: ', vc(xi_dp, yi_dp, zi_dpst:zi_dped)
        if ( abs(cmode) == 0 .or. abs(cmode) == 1 .or. abs(cmode) == 2 ) then
            where(abs(advx  ) > 1e10) advx   = 0.
            where(abs(advy  ) > 1e10) advy   = 0.
            where(abs(gradx ) > 1e10) gradx  = 0.
            where(abs(grady ) > 1e10) grady  = 0.
            where(abs(hdiffx) > 1e10) hdiffx = 0.
            where(abs(hdiffy) > 1e10) hdiffy = 0.
            where(abs(vdiffx) > 1e10) vdiffx = 0.
            where(abs(vdiffy) > 1e10) vdiffy = 0.
            where(abs(ssh   ) > 1e10) ssh    = 0.
            print*, 'Advx', advx(xi_dp, yi_dp, zi_dpst:zi_dped)
            print*, 'Advy', advy(xi_dp, yi_dp, zi_dpst:zi_dped)
            print*, 'Gradx', gradx(xi_dp, yi_dp, zi_dpst:zi_dped)
            print*, 'Grady', grady(xi_dp, yi_dp, zi_dpst:zi_dped)
            print*, 'Hdiffx', hdiffx(xi_dp, yi_dp, zi_dpst:zi_dped)
            print*, 'Hdiffy', hdiffy(xi_dp, yi_dp, zi_dpst:zi_dped)
            print*, 'Vdiffx', vdiffx(xi_dp, yi_dp, zi_dpst:zi_dped)
            print*, 'Vdiffy', vdiffy(xi_dp, yi_dp, zi_dpst:zi_dped)
            print*, 'SSH', ssh(xi_dp, yi_dp, 1)
        endif
        if ( abs(cmode) == 1 .or. abs(cmode) == 2 .or. abs(cmode) == 3 ) then
            where(abs(wc) > 1e10) wc = 0.
            print*, 'wc: ', wc(xi_dp, yi_dp, zi_dpst:zi_dped)
        endif
    endsubroutine

    subroutine load_mean_sf(fn_zetam)
        use ncio, only : nc_read
        ! use popload, only : load_current_day, find_daily_file
        implicit none
        character(len=*), intent(in) :: fn_zetam
        real(kind=kd_r), dimension(nx, ny, nz, 1) :: WORK

        print*, '  '
        print*, 'Loading mean field variables from input files'

        call nc_read(fn_zetam, 'advu'      , WORK); advu         = WORK(:, :, :, 1)
        call nc_read(fn_zetam, 'advv'      , WORK); advv         = WORK(:, :, :, 1)
        call nc_read(fn_zetam, 'advw'      , WORK); advw         = WORK(:, :, :, 1)
        call nc_read(fn_zetam, 'advVx'     , WORK); advVx        = WORK(:, :, :, 1)
        call nc_read(fn_zetam, 'advVy'     , WORK); advVy        = WORK(:, :, :, 1)
        call nc_read(fn_zetam, 'advVz'     , WORK); advVz        = WORK(:, :, :, 1)
        call nc_read(fn_zetam, 'curlmet'   , WORK); curlmet      = WORK(:, :, :, 1)
        call nc_read(fn_zetam, 'errnldcmp' , WORK); err_nldecomp = WORK(:, :, :, 1)

        print*, 'advu', advu(xi_dp, yi_dp, zi_dpst:zi_dped)
        print*, 'advv', advv(xi_dp, yi_dp, zi_dpst:zi_dped)
        print*, 'advw', advw(xi_dp, yi_dp, zi_dpst:zi_dped)
        print*, 'advVx', advVx(xi_dp, yi_dp, zi_dpst:zi_dped)
        print*, 'advVy', advVy(xi_dp, yi_dp, zi_dpst:zi_dped)
        print*, 'advVz', advVz(xi_dp, yi_dp, zi_dpst:zi_dped)
        print*, 'err_nldecomp', err_nldecomp(xi_dp, yi_dp, zi_dpst:zi_dped)
        print*, 'Curlmet: ', curlmet(xi_dp, yi_dp, zi_dpst:zi_dped)
    endsubroutine

    subroutine loadave_zeta_sf(yrlist, doylist)
        use ncio, only : nc_read
        implicit none
        integer, intent(in) :: yrlist(:), doylist(:)
        character :: fn_zeta*300, yyyymmdd*20
        integer :: iyr, idoy
        integer :: nn, mm, dd
        real(kind=kd_r), dimension(nx, ny, nz, 1) :: WORK
        real(kind=kd_r), dimension(nx, ny, 1 , 1) :: WORK2

        nn = count(doylist/=0) * count(yrlist/=0)
        print*, '  '
        print*, 'Loading and averaing variables (velocity, mom. terms) from input files over ', nn, ' day(s)'

        do iyr = 1, count(yrlist/=0)
            do idoy = 1, count(doylist/=0)
                call get_yyyymmdd(yrlist(iyr), doylist(idoy), '', '', fn_out_dlm, yyyymmdd)

                write(fn_zeta, '(A, A, A, A, A)') &
                    trim(fn_out_dir), trim(fn_out_pfx), trim(yyyymmdd), trim(fn_out_sfx), '.nc'
                print*, fn_zeta

                call nc_read(fn_zeta, 'curlnonl' , WORK); curlnonl      = curlnonl      + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'betav'    , WORK); betav         = betav         + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'stretchp' , WORK); stretchp      = stretchp      + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'errcor'   , WORK); err_cor       = err_cor       + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'curlpgrad', WORK); curlpgrad     = curlpgrad     + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'curlhdiff', WORK); curlhdiff     = curlhdiff     + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'curlvdiff', WORK); curlvdiff     = curlvdiff     + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'res'      , WORK); res           = res           + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'errnlsub' , WORK); err_nlsub     = err_nlsub     + WORK(:, :, :, 1) / nn

                call nc_read(fn_zeta, 'advu'     , WORK); advuT         = advuT         + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'advv'     , WORK); advvT         = advvT         + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'advw'     , WORK); advwT         = advwT         + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'advVx'    , WORK); advVxT        = advVxT        + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'advVy'    , WORK); advVyT        = advVyT        + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'advVz'    , WORK); advVzT        = advVzT        + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'curlmet'  , WORK); curlmetT      = curlmetT      + WORK(:, :, :, 1) / nn
                call nc_read(fn_zeta, 'errnldcmp', WORK); err_nldecompT = err_nldecompT + WORK(:, :, :, 1) / nn
            enddo
        enddo
    endsubroutine

    subroutine get_mmdd(sec, meanfreq, mmdd)
        implicit none
        integer, intent(in) :: sec(2)
        character(len = *), intent(in) :: meanfreq
        character(len = *), intent(inout) :: mmdd
        integer :: mm, dd

        if (len(trim(meanfreq)) == 0) then
            write(mmdd, '(A, I0.3, A, I0.3)') 'd', sec(1), '-', sec(2)
        elseif (meanfreq == "m") then
            call doy2date(sec(1), mm, dd)
            write(mmdd, '(I0.2)') mm
        elseif (meanfreq == "a") then
            write(mmdd, '(A)') 'ann'
        endif
    endsubroutine

    ! get yyyymmdd from year and doy. The string is overriden by yrnm_clm and avnm_clm
    subroutine get_yyyymmdd(iyr, idoy, yrnm_clm, avnm_clm, dlm, yyyymmdd)
        implicit none
        integer, intent(in) :: iyr, idoy
        character(len=*), intent(inout) :: yyyymmdd
        character(len=*), intent(in) :: yrnm_clm, avnm_clm, dlm
        character :: yyyy*9, mmdd*9
        integer :: mm, dd

        if (len(trim(yrnm_clm)) > 0) then
            write(yyyy, '(A)') trim(yrnm_clm)
        else
            write(yyyy, '(I0.4)') iyr
        endif

        if (len(trim(avnm_clm)) > 0) then
            write(mmdd, '(A, A)') trim(dlm), trim(avnm_clm)
        else
            if     (idoy == 0) then
                write(mmdd, '(A)') ''
            elseif (idoy <  0) then
                call doy2date(-idoy, mm, dd)
                write(mmdd, '(A, I0.2)') trim(dlm), mm
            else
                call doy2date( idoy, mm, dd)
                write(mmdd, '(A, I0.2, A, I0.2)') trim(dlm), mm, trim(dlm), dd
            endif
        endif
        write(yyyymmdd, '(A, A)') trim(yyyy), trim(mmdd)
    endsubroutine

    subroutine doy2date(idoy, mm, dd)
        implicit none
        integer, intent(in) :: idoy
        integer, intent(out) :: mm, dd
        integer :: imn, cumdays

        do imn = 1, 12
            cumdays = sum(eom(1:imn))
            if (idoy <= cumdays) exit
        enddo

        mm = imn
        dd = idoy - sum(eom(1:mm)) + eom(mm)
    endsubroutine

    subroutine output_me(fn_me)
        use netcdf
        implicit none
        character(len = *), intent(in) :: fn_me
        integer :: ncid, stat_create, stat_defdim, stat_defvar, stat_putatt, stat_inqvar, &
                   stat_getvar, stat_putvar, stat_io
        integer :: dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time
        integer :: varid_out_lon, varid_out_lat, varid_out_dep, varid_out_time

        print*, '  '
        print*, 'Creating output file ', trim(fn_me)
        print*, '  Start netcdf define ...'

        stat_create = nf90_create(trim(fn_me), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid)

        ! Dimension
        stat_defdim = nf90_def_dim(ncid, "nlon", nx, dimid_out_lon)
        stat_defdim = nf90_def_dim(ncid, "nlat", ny, dimid_out_lat)
        stat_defdim = nf90_def_dim(ncid, "z_t" , nz, dimid_out_dep)
        stat_defdim = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dimid_out_time)

        ! Coordinates
        stat_defvar = nf90_def_var(ncid, "TLONG", NF90_FLOAT, &
           (/dimid_out_lon, dimid_out_lat/), varid_out_lon)
        stat_defvar = nf90_def_var(ncid, "TLAT",  NF90_FLOAT, &
           (/dimid_out_lon, dimid_out_lat/), varid_out_lat)
        stat_defvar = nf90_def_var(ncid, "z_t" ,  NF90_FLOAT, dimid_out_dep , varid_out_dep )
        stat_defvar = nf90_def_var(ncid, "time" , NF90_FLOAT, dimid_out_time, varid_out_time)

        ! Variables
        stat_defvar = nf90_def_var(ncid, "curlnonl" , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlnonl )
        stat_putatt = nf90_put_att(ncid, varid_out_curlnonl , "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_curlnonl , "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_curlnonl , "long_name", "Curl of nonlinear term (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_curlnonl , "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "curlpgrad", nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlpgrad)
        stat_putatt = nf90_put_att(ncid, varid_out_curlpgrad, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_curlpgrad, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_curlpgrad, "long_name", "Curl of pressure gradient term (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_curlpgrad, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "res"      , nc_xtype, &
            (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_res)
        stat_putatt = nf90_put_att(ncid, varid_out_res      , "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_res      , "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_res      , "long_name", "Residual (lhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_res      , "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "curlhdiff", nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlhdiff)
        stat_putatt = nf90_put_att(ncid, varid_out_curlhdiff, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_curlhdiff, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_curlhdiff, "long_name", "Curl of horizontal diffusion (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_curlhdiff, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "curlvdiff", nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlvdiff)
        stat_putatt = nf90_put_att(ncid, varid_out_curlvdiff, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_curlvdiff, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_curlvdiff, "long_name", "Curl of vertical diffusion (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_curlvdiff, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "betav"    , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_betav    )
        stat_putatt = nf90_put_att(ncid, varid_out_betav    , "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_betav    , "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_betav    , "long_name", "Advection of planetary vorticity term (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_betav    , "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "stretchp",  nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_stretchp  )
        stat_putatt = nf90_put_att(ncid, varid_out_stretchp , "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_stretchp , "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_stretchp , "long_name", "Planetary vorticity stretching term (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_stretchp , "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "errcor"  ,  nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errcor    )
        stat_putatt = nf90_put_att(ncid, varid_out_errcor   , "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_errcor   , "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_errcor   , "long_name", "Error from decomposing curl(-fv, fu) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_errcor   , "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "errnlsub" , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnlsub)
        stat_putatt = nf90_put_att(ncid, varid_out_errnlsub, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_errnlsub, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_errnlsub, "long_name", "Error from nonlinear term due to calculating offline (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_errnlsub, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advu_m"     , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advu)
        stat_putatt = nf90_put_att(ncid, varid_out_advu, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advu, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advu, "long_name", "Mean advection of relative vorticity by zonal velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advu, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advv_m"     , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advv)
        stat_putatt = nf90_put_att(ncid, varid_out_advv, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advv, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advv, "long_name", "Mean advection of relative vorticity by meridional velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advv, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advw_m"     , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advw)
        stat_putatt = nf90_put_att(ncid, varid_out_advw, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advw, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advw, "long_name", "Mean advection of relative vorticity by vertical velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advw, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advVx_m" , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVx)
        stat_putatt = nf90_put_att(ncid, varid_out_advVx, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advVx, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advVx, "long_name", "Mean twisting of zonal voricity by vertical velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advVx, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advVy_m"     , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVy)
        stat_putatt = nf90_put_att(ncid, varid_out_advVy, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advVy, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advVy, "long_name", "Mean twisting of meridional voricity by vertical velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advVy, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advVz_m" , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVz)
        stat_putatt = nf90_put_att(ncid, varid_out_advVz, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advVz, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advVz, "long_name", "Mean twisting of vertical voricity by vertical velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advVz, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "curlmet_m"  , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlmet)
        stat_putatt = nf90_put_att(ncid, varid_out_curlmet, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_curlmet, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_curlmet, "long_name", "Mean curl of metric term (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_curlmet, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "errnldcmp_m", nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnldcmp)
        stat_putatt = nf90_put_att(ncid, varid_out_errnldcmp, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_errnldcmp, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_errnldcmp, "long_name", "Mean error from nonlinear term due to decomposition (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_errnldcmp, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advu_e"     , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advue)
        stat_putatt = nf90_put_att(ncid, varid_out_advue, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advue, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advue, "long_name", "Eddy advection of relative vorticity by zonal velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advue, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advv_e"     , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advve)
        stat_putatt = nf90_put_att(ncid, varid_out_advve, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advve, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advve, "long_name", "Eddy advection of relative vorticity by meridional velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advve, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advw_e"     , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advwe)
        stat_putatt = nf90_put_att(ncid, varid_out_advwe, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advwe, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advwe, "long_name", "Eddy advection of relative vorticity by vertical velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advwe, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advVx_e" , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVxe)
        stat_putatt = nf90_put_att(ncid, varid_out_advVxe, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advVxe, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advVxe, "long_name", "Eddy twisting of zonal voricity by vertical velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advVxe, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advVy_e"     , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVye)
        stat_putatt = nf90_put_att(ncid, varid_out_advVye, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advVye, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advVye, "long_name", "Eddy twisting of meridional voricity by vertical velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advVye, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "advVz_e" , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVze)
        stat_putatt = nf90_put_att(ncid, varid_out_advVze, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_advVze, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_advVze, "long_name", "Eddy twisting of vertical voricity by vertical velocity (flux form) (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_advVze, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "curlmet_e"  , nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlmete)
        stat_putatt = nf90_put_att(ncid, varid_out_curlmete, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_curlmete, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_curlmete, "long_name", "Eddy curl of metric term (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_curlmete, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "errnldcmp_e", nc_xtype, &
          (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnldcmpe)
        stat_putatt = nf90_put_att(ncid, varid_out_errnldcmpe, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varid_out_errnldcmpe, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varid_out_errnldcmpe, "long_name", "Eddy error from nonlinear term due to decomposition (rhs)")
        stat_putatt = nf90_put_att(ncid, varid_out_errnldcmpe, "missing_value", MVALUE)

        stat_create = nf90_enddef(ncid)
        print*, "    Finished netcdf define!", stat_create

        ! Writing cooordinates
        stat_putvar = nf90_put_var(ncid, varid_out_lat,  tlat , &
           start = (/1, 1/), count = (/nx, ny/))
        stat_putvar = nf90_put_var(ncid, varid_out_lon,  tlong, &
           start = (/1, 1/), count = (/nx, ny/))
        stat_putvar = nf90_put_var(ncid, varid_out_dep,  z_t)

        print*, '  '
        print*, '  Start writing zeta eddy mean file ...'
        stat_putvar = nf90_put_var(ncid, varid_out_curlnonl  , curlnonl      , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_betav     , betav         , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_stretchp  , stretchp      , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_errcor    , err_cor       , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_curlpgrad , curlpgrad     , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_curlhdiff , curlhdiff     , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_curlvdiff , curlvdiff     , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_res       , res           , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_errnlsub  , err_nlsub     , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))

        stat_putvar = nf90_put_var(ncid, varid_out_advu      , advu          , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_advv      , advv          , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_advw      , advw          , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_advVx     , advVx         , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_advVy     , advVy         , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_advVz     , advVz         , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_curlmet   , curlmet       , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_errnldcmp , err_nldecomp  , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))

        stat_putvar = nf90_put_var(ncid, varid_out_advue     , advuT  - advu , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_advve     , advvT  - advv , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_advwe     , advwT  - advw , &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_advVxe    , advVxT - advVx, &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_advVye    , advVyT - advVy, &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_advVze    , advVzT - advVz, &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_curlmete  , curlmetT - curlmet, &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        stat_putvar = nf90_put_var(ncid, varid_out_errnldcmpe, err_nldecompT - err_nldecomp, &
            start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))

        stat_io = nf90_close(ncid)
        print*, '  Finished writing zeta mean eddy equation file ', stat_io
    endsubroutine

    ! subroutine load_vars(iyr, imn, idy)
    !   use popload, only : load_current_day, find_daily_file
    !   implicit none
    !   integer, intent(in) :: iyr, imn, idy
    !   character(len = 300) :: fn_in
    !   real(kind=kd_r), dimension(nx, ny, nz, 1) :: WORK
    !   real(kind=kd_r), dimension(nx, ny, 1 , 1) :: WORK2
    !
    !   print*, '  '
    !   print*, 'Loading variables (velocity, mom. terms) from input files'
    !
    !
    !   call nc_read(fn_in, 'UVEL' , WORK)
    !   uc = WORK(:, :, :, 1)
    !   call nc_read(fn_in, 'VVEL' , WORK)
    !   vc = WORK(:, :, :, 1)
    !
    !   where(abs(uc    ) > 1e10) uc     = 0.
    !   where(abs(vc    ) > 1e10) vc     = 0.
    !
    !   if ( abs(cmode) == 0 .or. abs(cmode) == 1 .or. abs(cmode) == 2 ) then
    !       call nc_read(fn_in, 'ADVU'   , WORK); advx   = WORK(:, :, :, 1)
    !       call nc_read(fn_in, 'ADVV'   , WORK); advy   = WORK(:, :, :, 1)
    !
    !       call nc_read(fn_in, 'GRADX'  , WORK); gradx  = WORK(:, :, :, 1)
    !       call nc_read(fn_in, 'GRADY'  , WORK); grady  = WORK(:, :, :, 1)
    !
    !       call nc_read(fn_in, 'HDIFFU' , WORK); hdiffx = WORK(:, :, :, 1)
    !       call nc_read(fn_in, 'HDIFFV' , WORK); hdiffy = WORK(:, :, :, 1)
    !
    !       call nc_read(fn_in, 'VDIFFU' , WORK); vdiffx = WORK(:, :, :, 1)
    !       call nc_read(fn_in, 'VDIFFV' , WORK); vdiffy = WORK(:, :, :, 1)
    !
    !       call nc_read(fn_in, 'SSH', WORK2); ssh(:, :, 1) = WORK2(:, :, 1, 1)
    !
    !       where(abs(advx  ) > 1e10) advx   = 0.
    !       where(abs(advy  ) > 1e10) advx   = 0.
    !       where(abs(gradx ) > 1e10) gradx  = 0.
    !       where(abs(grady ) > 1e10) grady  = 0.
    !       where(abs(hdiffx) > 1e10) hdiffx = 0.
    !       where(abs(hdiffy) > 1e10) hdiffy = 0.
    !       where(abs(vdiffx) > 1e10) vdiffx = 0.
    !       where(abs(vdiffy) > 1e10) vdiffy = 0.
    !       where(abs(ssh   ) > 1e10) ssh    = 0.
    !   endif
    !
    !   if ( abs(cmode) == 1 .or. abs(cmode) == 2 .or. abs(cmode) == 3 ) then
    !       call nc_read(fn_in, 'WVEL' , WORK); wc = WORK(:, :, :, 1)
    !       where(abs(wc) > 1e10) wc = 0.
    !   endif
    !
    !   call find_daily_file(fn_in_dir, fn_in_pfx, iyr, imn, idy, fn_in)
    !   print*, 'Reading from file : ', trim(fn_in)
    !
    !   call load_current_day(iyr, imn, idy, fn_in, 'UVEL'  , WORK)
    !   uc = WORK(:, :, :, 1)
    !   call load_current_day(iyr, imn, idy, fn_in, 'VVEL'  , WORK)
    !   vc = WORK(:, :, :, 1)
    !   where(abs(uc    ) > 1e10) uc     = 0.
    !   where(abs(vc    ) > 1e10) vc     = 0.
    !
    !   if ( abs(cmode) == 0 .or. abs(cmode) == 1 .or. abs(cmode) == 2 ) then
    !       call load_current_day(iyr, imn, idy, fn_in, 'ADVU'  , WORK)
    !       advx = WORK(:, :, :, 1)
    !       call load_current_day(iyr, imn, idy, fn_in, 'ADVV'  , WORK)
    !       advy = WORK(:, :, :, 1)
    !
    !       call load_current_day(iyr, imn, idy, fn_in, 'GRADX' , WORK)
    !       gradx = WORK(:, :, :, 1)
    !       call load_current_day(iyr, imn, idy, fn_in, 'GRADY' , WORK)
    !       grady = WORK(:, :, :, 1)
    !
    !       call load_current_day(iyr, imn, idy, fn_in, 'HDIFFU', WORK)
    !       hdiffx = WORK(:, :, :, 1)
    !       call load_current_day(iyr, imn, idy, fn_in, 'HDIFFV', WORK)
    !       hdiffy = WORK(:, :, :, 1)
    !
    !       call load_current_day(iyr, imn, idy, fn_in, 'VDIFFU', WORK)
    !       vdiffx = WORK(:, :, :, 1)
    !       call load_current_day(iyr, imn, idy, fn_in, 'VDIFFV', WORK)
    !       vdiffy = WORK(:, :, :, 1)
    !
    !       call load_current_day(iyr, imn, idy, fn_in, 'SSH'   , WORK2)
    !       ssh(:, :, 1) = WORK2(:, :, 1, 1)
    !
    !       where(abs(advx  ) > 1e10) advx   = 0.
    !       where(abs(advy  ) > 1e10) advx   = 0.
    !       where(abs(gradx ) > 1e10) gradx  = 0.
    !       where(abs(grady ) > 1e10) grady  = 0.
    !       where(abs(hdiffx) > 1e10) hdiffx = 0.
    !       where(abs(hdiffy) > 1e10) hdiffy = 0.
    !       where(abs(vdiffx) > 1e10) vdiffx = 0.
    !       where(abs(vdiffy) > 1e10) vdiffy = 0.
    !       where(abs(ssh   ) > 1e10) ssh    = 0.
    !   endif
    !
    !   if ( abs(cmode) == 1 .or. abs(cmode) == 2 .or. abs(cmode) == 3 ) then
    !       call load_current_day(iyr, imn, idy, fn_in, 'WVEL'  , WORK)
    !       wc = WORK(:, :, :, 1)
    !       where(abs(wc) > 1e10) wc = 0.
    !   endif
    !
    !   print*, 'uc: ', uc(xi_dp, yi_dp, zi_dpst:zi_dped)
    !   print*, 'vc: ', vc(xi_dp, yi_dp, zi_dpst:zi_dped)
    !   print*, 'Advx', advx(xi_dp, yi_dp, zi_dpst:zi_dped)
    !   print*, 'Advy', advy(xi_dp, yi_dp, zi_dpst:zi_dped)
    !   print*, 'Gradx', gradx(xi_dp, yi_dp, zi_dpst:zi_dped)
    !   print*, 'Grady', grady(xi_dp, yi_dp, zi_dpst:zi_dped)
    !   print*, 'Hdiffx', hdiffx(xi_dp, yi_dp, zi_dpst:zi_dped)
    !   print*, 'Hdiffy', hdiffy(xi_dp, yi_dp, zi_dpst:zi_dped)
    !   print*, 'Vdiffx', vdiffx(xi_dp, yi_dp, zi_dpst:zi_dped)
    !   print*, 'Vdiffy', vdiffy(xi_dp, yi_dp, zi_dpst:zi_dped)
    !   print*, 'SSH', ssh(xi_dp, yi_dp, 1)
    ! endsubroutine

        ! subroutine create_outputfiles(fn_zeta, ncid_zeta, fn_error, ncid_error)

    subroutine create_outputfiles(cmode, fn_zeta, fn_error, ncid_zeta, ncid_error)
        use netcdf
        implicit none
        integer, intent(in) :: cmode
        character(len = *), intent(in) :: fn_zeta, fn_error
        integer, intent(inout) :: ncid_zeta, ncid_error
        integer :: stat_create, stat_defdim, stat_defvar, stat_putatt, stat_inqvar, &
                   stat_getvar, stat_putvar, stat_io
        integer :: dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time
        integer :: varid_out_lon, varid_out_lat, varid_out_dep, varid_out_time

        print*, '  '
        print*, 'Creating output file ', trim(fn_zeta)
        print*, '  Start netcdf define ...'

        stat_create = nf90_create(trim(fn_zeta), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid_zeta)

        ! Dimension
        stat_defdim = nf90_def_dim(ncid_zeta, "nlon", nx, dimid_out_lon)
        stat_defdim = nf90_def_dim(ncid_zeta, "nlat", ny, dimid_out_lat)
        stat_defdim = nf90_def_dim(ncid_zeta, "z_t" , nz, dimid_out_dep)
        stat_defdim = nf90_def_dim(ncid_zeta, "time", NF90_UNLIMITED, dimid_out_time)

        ! Coordinates
        stat_defvar = nf90_def_var(ncid_zeta, "TLONG", NF90_FLOAT, &
           (/dimid_out_lon, dimid_out_lat/), varid_out_lon)
        stat_defvar = nf90_def_var(ncid_zeta, "TLAT",  NF90_FLOAT, &
           (/dimid_out_lon, dimid_out_lat/), varid_out_lat)
        stat_defvar = nf90_def_var(ncid_zeta, "z_t" ,  NF90_FLOAT, dimid_out_dep , varid_out_dep )
        stat_defvar = nf90_def_var(ncid_zeta, "time" , NF90_FLOAT, dimid_out_time, varid_out_time)

        ! Variables
        if (abs(cmode) == 0 .or. abs(cmode) == 1 .or. abs(cmode) == 2) then
            stat_defvar = nf90_def_var(ncid_zeta, "curlnonl" , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlnonl )
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonl , "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonl , "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonl , "long_name", "Curl of nonlinear term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonl , "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "curlpgrad", nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlpgrad)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlpgrad, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlpgrad, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlpgrad, "long_name", "Curl of pressure gradient term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlpgrad, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "res"      , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_res)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_res      , "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_res      , "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_res      , "long_name", "Residual (lhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_res      , "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "curlhdiff", nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlhdiff)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlhdiff, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlhdiff, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlhdiff, "long_name", "Curl of horizontal diffusion (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlhdiff, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "curlvdiff", nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlvdiff)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlvdiff, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlvdiff, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlvdiff, "long_name", "Curl of vertical diffusion (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlvdiff, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "betav"    , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_betav    )
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_betav    , "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_betav    , "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_betav    , "long_name", "Advection of planetary vorticity term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_betav    , "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "stretchp",  nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_stretchp  )
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_stretchp , "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_stretchp , "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_stretchp , "long_name", "Planetary vorticity stretching term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_stretchp , "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "errcor"  ,  nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errcor    )
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errcor   , "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errcor   , "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errcor   , "long_name", "Error from decomposing curl(-fv, fu) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errcor   , "missing_value", MVALUE)
        endif

        if (abs(cmode) == 1) then
            stat_defvar = nf90_def_var(ncid_zeta, "curlnonl_m"  ,  nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlnonlc)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonlc, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonlc, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonlc, "long_name", "Mean curl of nonlinear term offline (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlnonlc, "missing_value", MVALUE)
        endif

        if (abs(cmode == 1) .or. (abs(cmode == 2) .and. yrnm_clm /= "")) then
            stat_defvar = nf90_def_var(ncid_zeta, "curlnonl_e" , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnlsub)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "long_name", "Eddy curl of nonlinear term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "missing_value", MVALUE)
        endif

        if (abs(cmode == 2) .and. yrnm_clm == "") then
            stat_defvar = nf90_def_var(ncid_zeta, "errnlsub" , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnlsub)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "long_name", "Error from nonlinear term due to calculating offline (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnlsub, "missing_value", MVALUE)
        endif

        if (abs(cmode) == 2 .or. abs(cmode) == 3) then
            stat_defvar = nf90_def_var(ncid_zeta, "advu"     , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advu)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advu, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advu, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advu, "long_name", "Advection of relative vorticity by zonal velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advu, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "advv"     , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advv)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advv, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advv, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advv, "long_name", "Advection of relative vorticity by meridional velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advv, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "advw"     , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advw)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advw, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advw, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advw, "long_name", "Advection of relative vorticity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advw, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "advVx" , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVx)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVx, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVx, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVx, "long_name", "Twisting of zonal voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVx, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "advVy"     , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVy)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVy, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVy, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVy, "long_name", "Twisting of meridional voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVy, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "advVz"     , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_advVz)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVz, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVz, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVz, "long_name", "Twisting of vertical voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_advVz, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "curlmet"  , nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_curlmet)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlmet, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlmet, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlmet, "long_name", "Curl of metric term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_curlmet, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "errnldcmp", nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnldcmp)
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnldcmp, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnldcmp, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnldcmp, "long_name", "Error from nonlinear term due to decomposition (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varid_out_errnldcmp, "missing_value", MVALUE)
        endif

        stat_create = nf90_enddef(ncid_zeta)
        print*, "    Finished netcdf define!", stat_create

        ! Writing cooordinates
        stat_putvar = nf90_put_var(ncid_zeta, varid_out_lat,  tlat , &
           start = (/1, 1/), count = (/nx, ny/))
        stat_putvar = nf90_put_var(ncid_zeta, varid_out_lon,  tlong, &
           start = (/1, 1/), count = (/nx, ny/))
        stat_putvar = nf90_put_var(ncid_zeta, varid_out_dep,  z_t)

        if (cmode < 0) then
            print*, '  '
            print*, 'Creating output error file ', trim(fn_error)

            stat_create = nf90_create(trim(fn_error), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid_error)

            stat_defdim = nf90_def_dim(ncid_error, "nlon", nx, dimid_out_lon)
            stat_defdim = nf90_def_dim(ncid_error, "nlat", ny, dimid_out_lat)
            stat_defdim = nf90_def_dim(ncid_error, "z_t" , nz, dimid_out_dep)
            stat_defdim = nf90_def_dim(ncid_error, "time", NF90_UNLIMITED, dimid_out_time)

            stat_defvar = nf90_def_var(ncid_error, "TLONG", NF90_FLOAT, &
               (/dimid_out_lon, dimid_out_lat/), varid_out_lon)
            stat_defvar = nf90_def_var(ncid_error, "TLAT",  NF90_FLOAT, &
               (/dimid_out_lon, dimid_out_lat/), varid_out_lat)
            stat_defvar = nf90_def_var(ncid_error, "z_t" ,  NF90_FLOAT, dimid_out_dep , varid_out_dep )
            stat_defvar = nf90_def_var(ncid_error, "time" , NF90_FLOAT, dimid_out_time, varid_out_time)

            stat_defvar = nf90_def_var(ncid_error, "rr_rev",   nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnl_rev)
            stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_rev, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_rev, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_rev, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_error, "rr_cha",   nc_xtype, &
              (/dimid_out_lon, dimid_out_lat, dimid_out_dep, dimid_out_time/), varid_out_errnl_cha)
            stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_cha, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_cha, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_error, varid_out_errnl_cha, "missing_value", MVALUE)

            stat_create = nf90_enddef(ncid_error)
            print*, "    Finished netcdf define!", stat_create

            stat_putvar = nf90_put_var(ncid_error, varid_out_lat,  tlat , &
               start = (/1, 1/), count = (/nx, ny/))
            stat_putvar = nf90_put_var(ncid_error, varid_out_lon,  tlong, &
               start = (/1, 1/), count = (/nx, ny/))
            stat_putvar = nf90_put_var(ncid_error, varid_out_dep,  z_t)
        else
            ncid_error = 0
        endif
    endsubroutine

    subroutine write_outputfiles(cmode, ncid_zeta, ncid_error)
        use netcdf
        implicit none
        integer, intent(in) :: cmode, ncid_zeta, ncid_error
        integer :: stat_putvar

        print*, '  '
        print*, '  Start writing zeta equation file ...'

        if (abs(cmode) == 0 .or. abs(cmode) == 1 .or. abs(cmode) == 2) then
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlnonl , curlnonl , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_betav    , betav    , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_stretchp , stretchp , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_errcor   , err_cor  , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlpgrad, curlpgrad, &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlhdiff, curlhdiff, &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlvdiff, curlvdiff, &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_res      , res      , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        endif

        if (abs(cmode) == 1) then
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlnonlc , curladv + curlmet , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        endif

        if (abs(cmode == 1) .or. abs(cmode == 2)) then
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_errnlsub , err_nlsub   , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        endif

        if (abs(cmode) == 2 .or. abs(cmode) == 3) then
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_advu     , advu        , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_advv     , advv        , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_advw     , advw        , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_advVx    , advVx       , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_advVy    , advVy       , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_advVz    , advVz       , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_curlmet  , curlmet     , &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varid_out_errnldcmp, err_nldecomp, &
                start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
        endif
        print*, '  Finished writing zeta equation file '

        if (cmode < 0) then
            print*, '  '
            print*, '  Start writing error file ...'

            stat_putvar = nf90_put_var(ncid_error, varid_out_errnl_rev, rr_rev, &
               start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            stat_putvar = nf90_put_var(ncid_error, varid_out_errnl_cha, rr_cha, &
               start = (/1, 1, 1, 1/), count = (/nx, ny, nz, 1/))
            print*, '  Finished writing error file '
        endif
    endsubroutine

    subroutine close_outputfiles(cmode, ncid_zeta, ncid_error)
        use netcdf
        implicit none
        integer, intent(in) :: cmode, ncid_zeta, ncid_error
        integer :: stat_io

        print*, '  '
        print*, 'Closing output files ...'
        stat_io = nf90_close(ncid_zeta)
        if (cmode < 0) then
          stat_io = nf90_close(ncid_error)
        endif
    endsubroutine

    subroutine init_meaneddy_vars()
        deallocate(curladv, curlmet)

        call init_zetavars_output(1)
        allocate(advuT(nx, ny, nz), advvT(nx, ny, nz), advwT(nx, ny, nz), &
            advVxT(nx, ny, nz), advVyT(nx, ny, nz), advVzT(nx, ny, nz), &
            curlmetT(nx, ny, nz), err_nldecompT(nx, ny, nz))
        advuT = 0. ; advvT = 0.; advwT = 0.
        advVxT = 0.; advVyT = 0.; advVzT = 0.
        curlmetT = 0.; err_nldecompT = 0.
    endsubroutine
endmodule
