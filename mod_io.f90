module io
    use params
    use zeta, only : uc, vc, wc, ssh, advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy, &
                     ueu, uev, vnu, vnv, wtu, wtv, &
                     tlat, tlong, z_t, &
                     curlnonl, betav, stretchp, err_cor, curlpgrad, curlhdiff, curlvdiff, res, &
                     curladv, curlmet, err_nlsub, advu, advv, advw, advVx, advVy, advVz, err_nldecomp, &
                     curladvf
    implicit none
    private
    public :: load_params, get_yyyymmdd, loadave_mom_sf, output_sf, load_mean_sf, &
              loadave_vor_sf, init_zetavars_mean, release_zetavars_mean, output_me, check_meanfile

    type filename
        ! Filename = dir + pfx + YYYY + dlm + MM + dlm + DD + sfx
        ! YYYY + dlm + MM + dlm + DD can be overriden by wildcards.
        character(len = 300) :: fn         ! filename
        character(len = 300) :: dir        ! directory
        character(len = 100) :: pfx        ! prefix
        character(len = 20 ) :: sfx        ! suffix
        character(len = 1  ) :: dlm        ! delimiter in dates
    endtype

    type timelist
        ! Year part of the input climatology file
        character(len=9) :: yrnm_clm = ""
        ! Month+Day part of the input climatology file.
        character(len=9) :: avnm_clm = ""
        ! list of years to loop
        integer, allocatable :: yrlist(:)
        ! list of doy to loop
        integer, allocatable :: doylist(:)
        ! Wildcard for the year part of mean file if ifmeanclm = T. It is used to override the default "yrst-yred".
        character(len=9) :: menm_clm = ""
        ! list of sections for mean eddy decomposition
        integer, allocatable :: seclist(:,:)
        ! list of section names
        character(len=9), dimension(:), allocatable :: meannm
        integer :: ndoy, nyr, nsec
    endtype

    ! netCDF handles
    type nc_varid
        integer :: curlnonl , betav    , stretchp , errcor    , &
                   curlpgrad, curlhdiff, curlvdiff, res       , &
                   curlnonlc, errnlsub ,                        &
                   advu     , advv     , advw     , errnldcmp , &
                   advVx    , advVy    , advVz    , curlmet   , &
                   advum    , advvm    , advwm    , errnldcmpm, &
                   advVxm   , advVym   , advVzm   , curlmetm
    endtype

    !!--------------------------------------------------------------------------
    ! * Calculation modes *
    logical, public :: ifcurl = .False., ifdecomp = .False.
    logical, public :: ifmeaneddy = .False., ifmeanclm = .FALSE., ifdecomposed = .True.

    !!--------------------------------------------------------------------------
    type(filename), public :: fn_mom, fn_vor, fn_vorm, fn_vore
    type(timelist), public :: T

    !!--------------------------------------------------------------------------
    ! * Mean field variables (eddy field is the difference between these and the total) *
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
        advu_m, advv_m, advw_m, advVx_m, advVy_m, advVz_m, curlmet_m, err_nldecomp_m, curlnonl_m

    !!--------------------------------------------------------------------------
    ! * Formats for screen display logs *
    character(len = 100) :: fmts_vel, fmtm_vel, fmts_mom, fmtm_mom, fmts_flx, fmtm_flx, fmts_vor, fmtm_vor

    contains
    subroutine load_params(func, func_m, func_me)
        use, intrinsic :: IEEE_ARITHMETIC
        implicit none
        ! Function mode (for mod_zeta) for the main loop, mean and eddy
        character(len=*), intent(inout) :: func, func_m, func_me
        integer :: nml_error
        integer :: iyr, imn, nmn, ida, nda, idoy, isec

        !!--------------------------------------------------------------------------
        ! * Date for the input files *
        ! Date format
        ! Two options in input.nml for yr, mn and da (listed by prioity)
        ! (a): list
        ! (b): starting and ending (st, ed)
        ! For year, additonal option is the name for climatology (yrnm_clm). If none is provide, then fatal error.
        ! For month and days, if neither list or st/ed is provided, then assume 12 months and 31* days
        ! If any negative number appears in mnlist, dalist, then assume average at this level.
        !    For example, mnlist_in = -1 indicates input file is annually averaged (therefore, no mnlist AND dalist needed).
        ! Starting and ending year
        integer :: yrst = 0, yred = 0
        ! Starting and ending month
        integer :: mnst = 1, mned = 12
        ! Starting and ending day
        integer :: dast = 1, daed = 31
        ! list of yr
        integer :: yrlist_in(60) = 0
        ! list of month
        integer :: mnlist_in(12) = 0 ! default is zero, = 12 months
        integer, allocatable :: mnlist(:)
        ! list of days
        integer :: dalist_in(31) = 0 ! default is zero, = 31 days
        integer, allocatable :: dalist(:)
        ! list of days in year
        integer :: doylist_full(365) = 0

        ! Wildcard for the Year part of the input climatology file.
        character(len=9) :: yrnm_clm = ""
        ! Wildcard for the Month+Day part of the input climatology file.
        character(len=9) :: avnm_clm = ""

        !!--------------------------------------------------------------------------
        ! * Mean field frequency *
        ! Date format
        ! Four options are available (listed by decreasing priority)
        ! (a) meanfreq: period name ("a" for annual, "m" for monthly)
        ! (b) seclist_in_st, seclist_in_ed: list of starting and ending date
        ! (c) secst, seced: starting and ending day (of a year) for a single section
        ! (d) nda_sec: Fixed number of days in each section
        character(len = 5) :: meanfreq = ""
        integer, dimension(365) :: seclist_in_st = 0, seclist_in_ed = 0
        integer :: secst = 0, seced = 0
        integer :: nda_sec = 365
        ! Wildcard for the year part of mean file if ifmeanclm = T. It is used to override the default "yrst-yred".
        character(len = 9) :: menm_clm = ""

        !!--------------------------------------------------------------------------
        ! * Input file info *
        ! Input file directory (suppose all in the same directory)
        character(len = 300) :: fn_mom_dir = "/lustre/atlas1/cli115/proj-shared/jritchie/SIO/MODEL_DATA/yelpatch60/IOC/DATA/DAT_IO_all/nday1/"
        ! Input file prefix (before "YYYY")
        character(len = 100) :: fn_mom_pfx = "ia_top_tx0.1_v2_yel_patc_1948_intel.pop.h.nday1."
        ! Input file suffix (after "DD", or subregion)
        character(len = 20 ) :: fn_mom_sfx = ".IO"
        ! Input file date delimiter (between YYYY, MM, DD)
        character(len = 1  ) :: fn_mom_dlm = '-'

        !!----------------------------------------------------------------------------
        ! * Output file info *
        ! Vorticity output file
        character(len = 300) :: fn_vor_dir = "/lustre/atlas/proj-shared/cli115/hewang/data/"
        character(len = 100) :: fn_vor_pfx = "vort_bud_"
        character(len = 20 ) :: fn_vor_sfx = ""     ! Default is subregion
        character(len = 1  ) :: fn_vor_dlm = ""     ! Default is the same as input

        ! Mean/Eddy output files
        character(len = 300) :: fn_vorm_dir = ""    ! Default is the same as fn_vor
        character(len = 100) :: fn_vorm_pfx = ""    ! Default is the same as fn_vor
        character(len = 20 ) :: fn_vorm_sfx = "_m"  ! Default is fn_vor_sfx + _m
        character(len = 1  ) :: fn_vorm_dlm = ""    ! Default is the same as input

        character(len = 300) :: fn_vore_dir = ""    ! Default is the same as fn_vor
        character(len = 100) :: fn_vore_pfx = ""    ! Default is the same as fn_vor
        character(len = 20 ) :: fn_vore_sfx = "_me" ! Default is fn_vor_sfx + _me
        character(len = 1  ) :: fn_vore_dlm = ""    ! Default is the same as input

        namelist /calcmode/ ifcurl, ifdecomp
        namelist /memode/ ifmeaneddy, ifmeanclm, ifdecomposed, menm_clm, &
                          meanfreq, nda_sec, seclist_in_st, seclist_in_ed, secst, seced
        namelist /time/ yrst, yred, mnst, mned, dast, daed, &
                        yrnm_clm, avnm_clm, yrlist_in, mnlist_in, dalist_in
        namelist /bnds/ B
        namelist /grid_files/ fngrid
        namelist /input_files/ fn_mom_dir, fn_mom_pfx, fn_mom_sfx, fn_mom_dlm
        namelist /output_files/ fn_vor_dir , fn_vor_pfx , fn_vor_sfx , fn_vor_dlm , &
                                fn_vorm_dir, fn_vorm_pfx, fn_vorm_sfx, fn_vorm_dlm, &
                                fn_vore_dir, fn_vore_pfx, fn_vore_sfx, fn_vore_dlm

        write(*, '(A)') '-----------------------------------------------------'
        open(101, file=fn_namelist, status="old", iostat=nml_error)
        read(101, nml=calcmode, iostat=nml_error)
          write(*, '(A30, I)') "input <calcmode>: ", nml_error
        read(101, nml=memode, iostat=nml_error)
          write(*, '(A30, I)') "input <memode>: ", nml_error
        read(101, nml=time, iostat=nml_error)
          write(*, '(A30, I)') "input <time>: ", nml_error
        read(101, nml=bnds, iostat=nml_error)
          write(*, '(A30, I)') "input <bnds>: ", nml_error
        read(101, nml=grid_files, iostat=nml_error)
          write(*, '(A30, I)') "input <grid>: ", nml_error
        read(101, nml=input_files, iostat=nml_error)
          write(*, '(A30, I)') "input <input>: ", nml_error
        read(101, nml=output_files, iostat=nml_error)
          write(*, '(A30, I)') "input <output>: ", nml_error
        close(101)
        write(*, *)

        !-------------------------------------------
        ! * Interpreting calculation mode *
        !-------------------------------------------
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') "Calculation options from user:"
        write(*, '(A40, L)'), "Curl of momentum equation? ", ifcurl
        write(*, '(A40, L)'), "Decomposition of the nonlinear term? ", ifdecomp
        write(*, '(A40, L)'), "Mean/Eddy decomposition? ", ifmeaneddy
        if (ifmeaneddy) then
            write(*, '(A40, L)'), "  with full decomposition for nonlinear term? ", ifdecomposed
            write(*, '(A40, L)'), "  use climatology as mean? ", ifmeanclm
        endif
        write(*, *)

        if (.not. ifcurl .and. .not. ifdecomp .and. .not. ifmeaneddy) then
            write(*, '(A)') "Do nothing. Exit"
            stop
        endif

        !-------------------------------------------
        ! * Translating "if" options into func, func_m, func_me *
        !-------------------------------------------
        write(*, '(A)') '-----------------------------------------------------'
        if (yrnm_clm /= "") then
            write(*, '(4A)') "yrnm_clm = ", trim(yrnm_clm), ". Input files are climatology averaged over ", trim(yrnm_clm)
            write(*, *)
            func_m = ""
            func_me = ""
            if (ifmeaneddy .and. .not. ifmeanclm) then
                write(*, '(A)') "  Warning: input files are climatology and 'ifmeanclm' is by definiation turned on for mean/eddy decomposition."
                write(*, *)
                ifmeanclm = .True.
            endif
            if     (      ifcurl .and. .not. ifdecomp .and. .not. ifmeaneddy) then
                write(*, '(A)') 'Curl of momentum equation for each input file'
                func = "c"
            elseif (      ifcurl .and. .not. ifdecomp .and.       ifmeaneddy) then
                write(*, '(A)') 'Curl of momentum equation and mean/eddy decomposition for the nonlinear term'
                write(*, '(A)') '  ERR_nlsub == curlnonl_eddy'
                func = "came"
            elseif (      ifcurl .and.       ifdecomp                       ) then
                write(*, '(A)') 'Curl of momentum equation and full decomposition for the nonlinear term'
                write(*, '(A)') '  ERR_nlsub == curlnonl_eddy and only mean terms have decompostion.'
                func = "cdme"
                if (.not. ifmeaneddy) then
                    write(*,'(A)') "  Note: 'ifmeaneddy' was set to False. A decomposition of nonlinear terms of the climatology &
                             is essentially equivalent to a mean/eddy decomposition."
                endif
            elseif (.not. ifcurl .and. .not. ifdecomp) then
                write(*, '(A)') "Both 'ifcurl' and 'ifdecomp' are set to False. Cannot do mean/eddy decomposition for climatology! Exit."
                stop
            endif
        else
            if     (      ifcurl .and. .not. ifdecomp .and. .not. ifmeaneddy) then
                write(*, '(A)') 'Curl of momentum equation for each input file'
                func = "c"
                func_m = ""
                func_me = ""
            elseif (      ifcurl .and. .not. ifdecomp .and.       ifmeaneddy) then
                write(*, '(A)') 'Curl of momentum equation and mean/eddy decomposition for the nonlinear term'
                func = "c"
                func_m = "am"
                func_me = func
            elseif (      ifcurl .and.       ifdecomp                       ) then
                write(*, '(A)') 'Curl of momentum equation and decomposition of the nonlinear term'
                func = "cdme"
                if (ifmeaneddy) then
                   write(*, '(A)'), '  and Mean/eddy decomposition'
                   func_m = "dm"
                   func_me = func
                else
                   func_m = ""
                   func_me = ""
                endif
            elseif (.not. ifcurl .and. .not. ifdecomp .and.       ifmeaneddy) then
                write(*, '(A)') "Mean/eddy decomposition will be calculated based on given decomposed zeta equation files"
                ! write(*, '(A)') "  Skipping the main loop ..."
                if (ifdecomposed) then
                    write(*, '(A)') "Zeta equation files include decomposed nonlinear terms."
                    func = ""
                    func_m = "dm"
                    func_me = "cdme"
                else
                    write(*, '(A)') "Zeta equation files does not include decomposed nonlinear terms."
                    func = ""
                    func_m = "am"
                    func_me = "c"
                endif
            endif
        endif
        if (.not. ifcurl .and. ifdecomp ) then
            write(*, '(A)'), 'Decomposition of the nonlinear term only.'
            func = "dm"
            func_m = ""
            func_me = ""
            if (ifmeaneddy) then
                write(*,'(A)') "  Warning: Cannot do mean/eddy decomposition without curl of momentum equation."
            endif
        endif
        write(*, *)
        write(*, '(A, A)') '  func    = ', trim(func)
        write(*, '(A, A)') '  func_m  = ', trim(func_m)
        write(*, '(A, A)') '  func_me = ', trim(func_me)
        write(*, *)

        !-------------------------------------------
        ! * Set missing value *
        !-------------------------------------------
        MVALUE = ieee_value(MVALUE, IEEE_QUIET_NAN)
        !print*, 'Missing value: ', MVALUE

        !-------------------------------------------
        ! * Relative domain boundaries *
        !-------------------------------------------
        B%xl = B%xl_reg - B%xl_ref + 1
        B%xr = B%xr_reg - B%xl_ref + 1
        B%yd = B%yd_reg - B%yd_ref + 1
        B%yu = B%yu_reg - B%yd_ref + 1

        B%nx = B%xr - B%xl + 1
        B%ny = B%yu - B%yd + 1
        B%nz = nzgl

        !-------------------------------------------
        ! * Interpreting date for the input files *
        !-------------------------------------------
        if (len(trim(yrnm_clm)) > 0) then
            T%nyr = 1
            allocate(T%yrlist(T%nyr))
            T%yrlist = -1
        elseif (count(yrlist_in/=0) /= 0 .and. all(yrlist_in >= 0)) then
            T%nyr = count(yrlist_in/=0)
            allocate(T%yrlist(T%nyr))
            T%yrlist = yrlist_in(1:T%nyr)
        elseif (yrst /= 0 .and. yrst /= 0) then
            T%nyr = yred - yrst + 1
            allocate(T%yrlist(T%nyr))
            T%yrlist = (/ (iyr, iyr = yrst, yred) /)
        else
            write(*, '(A)') "ERROR: No year information (yrlist or yrnm_clm) provided! Exit."
            stop
        endif

        if (any(mnlist_in < 0) .or. len(trim(avnm_clm)) > 0) then  !annual mean or arbitary averaged period
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

        if (any(dalist_in < 0) .or. any(mnlist_in < 0) .or. len(trim(avnm_clm)) > 0) then  !monthly or annual mean or arbitary averaged period
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

        !-------------------------------------------
        !    Translate mnlist and dalist to doylist (day of year)
        !-------------------------------------------
        ! Annual mean: doy_list = (/ 0 /)
        ! Monthly mean: doy_list = (/ -first day of the month /)
        ! Daily mean: doy_list = doy
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
        T%ndoy = idoy
        allocate(T%doylist(T%ndoy))
        T%doylist = doylist_full(1:T%ndoy)

        !-------------------------------------------
        !    Translate doylist to seclist
        !-------------------------------------------
        ! The main loop is through doylist, which is further devided into "sec" (sections).
        ! When mean/eddy decomposition is turned on, secs are based on the frequency/period。
        ! A sanity check will be needed to compare doylist and sections.
        if (trim(func_me)=="") then
            T%nsec = T%ndoy
            allocate(T%seclist(T%nsec, 2))
            T%seclist(:, 1) = (/ (T%doylist(ida), ida = 1, T%ndoy) /)
            T%seclist(:, 2) = T%seclist(:, 1)
        else
            if (meanfreq=="m") then
                T%nsec = size(mnlist)
                allocate(T%seclist(T%nsec, 2))
                allocate(T%meannm(T%nsec))
                do isec = 1, T%nsec
                    T%seclist(isec, 1) = 1 + sum(eom(1:mnlist(isec)-1))
                    T%seclist(isec, 2) = sum(eom(1:mnlist(isec)))
                    write(T%meannm(isec), '(I0.2)') mnlist(isec)
                enddo
            elseif (meanfreq=="a") then
                T%nsec = 1
                allocate(T%seclist(T%nsec, 2))
                allocate(T%meannm(T%nsec))
                T%seclist(1, 1) = 1
                T%seclist(1, 2) = 365
                write(T%meannm(1), '(A)') 'ann'
            elseif (count(seclist_in_st/=0) /= 0 .and. all(seclist_in_st >= 0) .and. &
                    count(seclist_in_ed/=0) /= 0 .and. all(seclist_in_ed >= 0)) then
                T%nsec = min(count(seclist_in_st/=0), count(seclist_in_ed/=0))
                allocate(T%seclist(T%nsec, 2))
                allocate(T%meannm(T%nsec))
                do isec = 1, T%nsec
                    T%seclist(isec, 1) = seclist_in_st(isec)
                    T%seclist(isec, 2) = seclist_in_ed(isec)
                    write(T%meannm(isec), '(A, I0.3, A, I0.3)') 'd', T%seclist(isec, 1), '-', T%seclist(isec, 2)
                enddo
            elseif (seced >= secst /= 0 .and. secst > 0) then
                T%nsec = 1
                allocate(T%seclist(T%nsec, 2))
                allocate(T%meannm(T%nsec))
                T%seclist(1, 1) = secst
                T%seclist(1, 2) = seced
                write(T%meannm(1), '(A, I0.3, A, I0.3)') 'd', secst, '-', seced
            else !
            ! If fixed number of day is specifed, min(doyslist) and max(doylist) are
            !   used as lower and upper bounds.
                T%nsec = (T%doylist(T%ndoy)-T%doylist(1))/nda_sec + 1
                allocate(T%seclist(T%nsec, 2))
                allocate(T%meannm(T%nsec))
                do isec = 1, T%nsec
                    T%seclist(isec, 1) = (isec - 1) * nda_sec + T%doylist(1)
                    T%seclist(isec, 2) = min(isec * nda_sec, T%doylist(T%ndoy))
                    write(T%meannm(isec), '(A, I0.3, A, I0.3)') &
                          'd', T%seclist(isec, 1), '-', T%seclist(isec, 2)
                enddo
            endif
        endif

        T%yrnm_clm = yrnm_clm
        T%avnm_clm = avnm_clm
        T%menm_clm = menm_clm

        !-------------------------------------------
        !  * Constructing filenames
        !-------------------------------------------
        !
        if (len(trim(fn_vor_sfx)) == 0) then
            write(fn_vor_sfx, '(A, A)') '_', trim(B%subreg)
        endif

        ! if (ifdecomp) then
        !     write(fn_vor_sfx, '(A, A)') trim(fn_vor_sfx), '_decomp'
        !     if (ifmeaneddy) then
        !         write(fn_vorm_sfx, '(A, A)') trim(fn_vorm_sfx), '_decomp'
        !         write(fn_vore_sfx, '(A, A)') trim(fn_vore_sfx), '_decomp'
        !     endif
        ! endif

        if (len(trim(fn_vor_dlm)) == 0) then
            write(fn_vor_dlm, '(A)') trim(fn_mom_dlm)
        endif
        if (len(trim(fn_vorm_dlm)) == 0) then
            write(fn_vorm_dlm, '(A)') trim(fn_mom_dlm)
        endif
        if (len(trim(fn_vore_dlm)) == 0) then
            write(fn_vore_dlm, '(A)') trim(fn_mom_dlm)
        endif
        if (len(trim(fn_vorm_dir)) == 0) then
            write(fn_vorm_dir, '(A)') trim(fn_vor_dir)
        endif
        if (len(trim(fn_vore_dir)) == 0) then
            write(fn_vore_dir, '(A)') trim(fn_vor_dir)
        endif
        if (len(trim(fn_vorm_pfx)) == 0) then
            write(fn_vorm_pfx, '(A)') trim(fn_vor_pfx)
        endif
        if (len(trim(fn_vore_pfx)) == 0) then
            write(fn_vore_pfx, '(A)') trim(fn_vor_pfx)
        endif

        fn_mom%dir = fn_mom_dir
        fn_mom%pfx = fn_mom_pfx
        fn_mom%sfx = fn_mom_sfx
        fn_mom%dlm = fn_mom_dlm

        fn_vor%dir = fn_vor_dir
        fn_vor%pfx = fn_vor_pfx
        fn_vor%sfx = fn_vor_sfx
        fn_vor%dlm = fn_vor_dlm

        fn_vorm%dir = fn_vorm_dir
        fn_vorm%pfx = fn_vorm_pfx
        write(fn_vorm%sfx, '(A, A)') trim(fn_vor_sfx), trim(fn_vorm_sfx)
        fn_vorm%dlm = fn_vorm_dlm

        fn_vore%dir = fn_vore_dir
        fn_vore%pfx = fn_vore_pfx
        write(fn_vore%sfx, '(A, A)') trim(fn_vor_sfx), trim(fn_vore_sfx)
        fn_vore%dlm = fn_vore_dlm

        !-------------------------------------------
        ! * Display format *
        !-------------------------------------------
        write(fmts_mom, '(A, A, A)' ) '(A20, ', trim(fmt_exp), ')'
        write(fmtm_mom, '(A, I2, A, A)') '(A20, ', B%zi_dped - B%zi_dpst + 1, trim(fmt_exp), ')'
        write(fmts_flx, '(A, A, A)' ) '(A20, ', trim(fmt_exp), ')'
        write(fmtm_flx, '(A, I2, A, A)') '(A20, ', B%zi_dped - B%zi_dpst + 1, trim(fmt_exp), ')'
        write(fmts_vel, '(A, A, A)' ) '(A20, ', trim(fmt_flt), ')'
        write(fmtm_vel, '(A, I2, A, A)') '(A20, ', B%zi_dped - B%zi_dpst + 1, trim(fmt_flt), ')'
        write(fmts_vor, '(A, A, A)' ) '(A20, ', trim(fmt_exp), ')'
        write(fmtm_vor, '(A, I2, A, A)') '(A20, ', B%zi_dped - B%zi_dpst + 1, trim(fmt_exp), ')'

        !-------------------------------------------
        ! * Print basic info *
        !-------------------------------------------
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A, A)') "Domain name: ", trim(B%subreg)
        write(*, '(A, I3, A, I3)') "xl: ", B%xl, " yd: ", B%yd
        write(*, '(A, I3, A, I3)') "Domain size (Nx x Ny): ", B%nx, " x ", B%ny
        write(*, *)

        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)', advance="no") "List of years: "
        do iyr = 1, T%nyr
            write(*, '(I4, A)', advance="no") T%yrlist(iyr), ' '
        enddo
        write(*, *)
        write(*, '(A)') "List of sections: "
        write(*, '(A)', advance="no") "  Start: "
        do ida = 1, T%nsec
            write(*, '(I3, A)', advance="no") T%seclist(ida, 1), ' '
        enddo
        write(*, *)
        write(*, '(A)', advance="no") "  End: "
        do ida = 1, T%nsec
            write(*, '(I3, A)', advance="no") T%seclist(ida, 2), ' '
        enddo
        write(*, *)
    endsubroutine

    ! Load and average momentum terms from single timestep files
    subroutine loadave_mom_sf(func, yrlist, doylist, yrnm_clm, avnm_clm, fn)
        use ncio, only : nc_read
        ! use popload, only : load_current_day, find_daily_file
        implicit none
        character(len=*), intent(in) :: func
        integer, intent(in) :: yrlist(:), doylist(:)
        character(len=*), intent(in) :: yrnm_clm, avnm_clm
        type(filename), intent(in) :: fn
        character :: fn_mom*300, yyyymmdd*20
        integer :: iyr, idoy
        integer :: nn, mm, dd
        real(kind=kd_r), dimension(B%nx, B%ny, B%nz, 1) :: WORK
        real(kind=kd_r), dimension(B%nx, B%ny, 1 , 1) :: WORK2

        ! nn = count(doylist/=0) * count(yrlist/=0)
        nn = size(doylist) * size(yrlist)

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A, I3, A)') 'Loading and averaing variables (velocity, mom. terms) from input files over ', nn, ' day(s)'

        do iyr = 1, size(yrlist)
            do idoy = 1, size(doylist)
                call get_yyyymmdd(yrlist(iyr), doylist(idoy),  yrnm_clm, avnm_clm, fn%dlm, yyyymmdd)
                write(fn_mom, '(A, A, A, A, A)') &
                    trim(fn%dir), trim(fn%pfx), trim(yyyymmdd), trim(fn%sfx), '.nc'

                write(*, '(A, A)') '  Load from file: ', trim(fn_mom)

                if (index(func, "a") /= 0 .or. index(func, "m") /= 0 .or. &
                    index(func, "d") /= 0 .or. index(func, "c") /= 0) then
                    call nc_read(fn_mom, 'UVEL' , WORK, (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); uc = uc + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'VVEL' , WORK, (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); vc = vc + WORK(:, :, :, 1) / nn
                endif
                if (index(func, "c") /= 0) then
                    call nc_read(fn_mom, 'ADVU'   , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); advx   = advx   + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'ADVV'   , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); advy   = advy   + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'GRADX'  , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); gradx  = gradx  + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'GRADY'  , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); grady  = grady  + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'HDIFFU' , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); hdiffx = hdiffx + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'HDIFFV' , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); hdiffy = hdiffy + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'VDIFFU' , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); vdiffx = vdiffx + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'VDIFFV' , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); vdiffy = vdiffy + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'SSH'    , WORK2, (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny,    1, 1/)); ssh(:, :, 1) = ssh(:, :, 1) + WORK2(:, :, 1, 1) / nn
                endif
                if (index(func, "a") /= 0 .or. index(func, "d") /= 0) then
                    call nc_read(fn_mom, 'WVEL'   , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); wc = wc + WORK(:, :, :, 1) / nn
                endif
                if (index(func, "f") /= 0) then
                    call nc_read(fn_mom, 'UEU'   , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); ueu = ueu + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'UEV'   , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); uev = uev + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'VNU'   , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); vnu = vnu + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'VNV'   , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); vnv = vnv + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'WTU'   , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); wtu = wtu + WORK(:, :, :, 1) / nn
                    call nc_read(fn_mom, 'WTV'   , WORK , (/B%xl, B%yd, 1, 1/), (/B%nx, B%ny, B%nz, 1/)); wtv = wtv + WORK(:, :, :, 1) / nn
                endif
            enddo
        enddo

        if (index(func, "a") /= 0 .or. index(func, "m") /= 0 .or. &
            index(func, "d") /= 0 .or. index(func, "c") /= 0) then
            where(abs(uc) > 1e10) uc = 0.
            where(abs(vc) > 1e10) vc = 0.
            write(*, fmtm_vel) 'UVEL: ', uc(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vel) 'VVEL: ', vc(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
        if (index(func, "c") /= 0) then
            where(abs(advx  ) > 1e10) advx   = 0.
            where(abs(advy  ) > 1e10) advy   = 0.
            where(abs(gradx ) > 1e10) gradx  = 0.
            where(abs(grady ) > 1e10) grady  = 0.
            where(abs(hdiffx) > 1e10) hdiffx = 0.
            where(abs(hdiffy) > 1e10) hdiffy = 0.
            where(abs(vdiffx) > 1e10) vdiffx = 0.
            where(abs(vdiffy) > 1e10) vdiffy = 0.
            where(abs(ssh   ) > 1e10) ssh    = 0.
            write(*, fmtm_mom) 'ADVX: ', advx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'ADVY: ', advy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'GRADX: ', gradx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'GRADY: ', grady(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'HDIFFU: ', hdiffx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'HDIFFV: ', hdiffy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'VDIFFU: ', vdiffx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'VDIFFV: ', vdiffy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmts_mom) 'SSH: ', ssh(B%xi_dp, B%yi_dp, 1)
        endif
        if (index(func, "a") /= 0 .or. index(func, "d") /= 0) then
            where(abs(wc) > 1e10) wc = 0.
            write(*, fmtm_vel) 'WVEL: ', wc(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
        if (index(func, "f") /= 0) then
            where(abs(ueu) > 1e10) ueu = 0.
            where(abs(uev) > 1e10) uev = 0.
            where(abs(vnu) > 1e10) vnu = 0.
            where(abs(vnv) > 1e10) vnv = 0.
            where(abs(wtu) > 1e10) wtu = 0.
            where(abs(wtv) > 1e10) wtv = 0.
            write(*, fmtm_mom) 'UEU: ', ueu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'UEV: ', uev(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'VNU: ', vnu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'VNV: ', vnv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'WTU: ', wtu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_mom) 'WTV: ', wtv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
    endsubroutine

    subroutine init_zetavars_mean(func)
        character(len=*), intent(in) :: func
        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Initializing mean field variables'

        if (index(func, "d") /= 0) then
            allocate(advu_m(B%nx, B%ny, B%nz), advv_m(B%nx, B%ny, B%nz), advw_m(B%nx, B%ny, B%nz), &
                     advVx_m(B%nx, B%ny, B%nz), advVy_m(B%nx, B%ny, B%nz), advVz_m(B%nx, B%ny, B%nz), &
                     curlmet_m(B%nx, B%ny, B%nz), err_nldecomp_m(B%nx, B%ny, B%nz))
            advu_m = 0. ; advv_m = 0.; advw_m = 0.
            advVx_m = 0.; advVy_m = 0.; advVz_m = 0.
            curlmet_m = 0.; err_nldecomp_m = 0.
        endif
        if (index(func, "f") /= 0) then
            allocate(curlnonl_m(B%nx, B%ny, B%nz))
            curlnonl_m = 0.
        endif
    endsubroutine

    subroutine release_zetavars_mean(func)
        character(len=*), intent(in) :: func
        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Release mean field variables'

        if (index(func, "d") /= 0) then
            deallocate(advu_m, advv_m, advw_m, advVx_m, advVy_m, advVz_m, curlmet_m, err_nldecomp_m)
        endif

        if (index(func, "a") /= 0) then
            deallocate(curlnonl_m)
        endif
    endsubroutine

    subroutine load_mean_sf(func, fn_zetam)
        use ncio, only : nc_read
        ! use popload, only : load_current_day, find_daily_file
        implicit none
        character(len=*), intent(in) :: func, fn_zetam
        real(kind=kd_r), dimension(B%nx, B%ny, B%nz, 1) :: WORK

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A, A)') 'Loading mean field variables from file ', trim(fn_zetam)

        if (index(func, "d") /= 0) then
            call nc_read(fn_zetam, 'advu'      , WORK); advu_m         = WORK(:, :, :, 1)
            call nc_read(fn_zetam, 'advv'      , WORK); advv_m         = WORK(:, :, :, 1)
            call nc_read(fn_zetam, 'advw'      , WORK); advw_m         = WORK(:, :, :, 1)
            call nc_read(fn_zetam, 'advVx'     , WORK); advVx_m        = WORK(:, :, :, 1)
            call nc_read(fn_zetam, 'advVy'     , WORK); advVy_m        = WORK(:, :, :, 1)
            call nc_read(fn_zetam, 'advVz'     , WORK); advVz_m        = WORK(:, :, :, 1)
            call nc_read(fn_zetam, 'curlmet'   , WORK); curlmet_m      = WORK(:, :, :, 1)
            call nc_read(fn_zetam, 'errnldcmp' , WORK); err_nldecomp_m = WORK(:, :, :, 1)

            write(*, fmtm_vor) 'advu: ', advu_m(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'advv: ', advv_m(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'advw: ', advw_m(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'advVx: ', advVx_m(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'advVy: ', advVy_m(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'advVz: ', advVz_m(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'err_nldecomp: ', err_nldecomp_m(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'curlmet: ', curlmet_m(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif

        if (index(func, "a") /= 0) then
            call nc_read(fn_zetam, 'curlnonl_m', WORK); curlnonl_m = WORK(:, :, :, 1)
            write(*, fmtm_vor) 'curlnonl_m: ', curlnonl_m(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
    endsubroutine

    ! Load and average vorticity terms from single files
    subroutine loadave_vor_sf(func, yrlist, doylist, yrnm_clm, avnm_clm, fn)
        use ncio, only : nc_read
        implicit none
        character(len=*), intent(in) :: func
        integer, intent(in) :: yrlist(:), doylist(:)
        character(len=*), intent(in) :: yrnm_clm, avnm_clm
        type(filename), intent(in) :: fn
        character :: fn_vor*300, yyyymmdd*20
        integer :: iyr, idoy
        integer :: nn, mm, dd
        real(kind=kd_r), dimension(B%nx, B%ny, B%nz, 1) :: WORK
        real(kind=kd_r), dimension(B%nx, B%ny, 1 , 1) :: WORK2

        nn = size(doylist) * size(yrlist)

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A, I3, A)') 'Loading and averaing variables (vorticity) from input files over ', nn, ' day(s)'

        do iyr = 1, size(yrlist)
            do idoy = 1, size(doylist)
                call get_yyyymmdd(yrlist(iyr), doylist(idoy), yrnm_clm, avnm_clm, fn%dlm, yyyymmdd)

                write(fn_vor, '(A, A, A, A, A)') &
                    trim(fn%dir), trim(fn%pfx), trim(yyyymmdd), trim(fn%sfx), '.nc'

                write(*, '(A, A)') '  Load from file: ', trim(fn_vor)

                if (index(func, "c") /= 0) then
                    call nc_read(fn_vor, 'curlnonl' , WORK); curlnonl      = curlnonl      + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'betav'    , WORK); betav         = betav         + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'stretchp' , WORK); stretchp      = stretchp      + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'errcor'   , WORK); err_cor       = err_cor       + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'curlpgrad', WORK); curlpgrad     = curlpgrad     + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'curlhdiff', WORK); curlhdiff     = curlhdiff     + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'curlvdiff', WORK); curlvdiff     = curlvdiff     + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'res'      , WORK); res           = res           + WORK(:, :, :, 1) / nn
                endif

                ! if (abs(cmode) == 1) then
                !     ! call nc_read(fn_vor, 'curlnonl_m' , WORK); curlnonlc     = curlnonlc     + WORK(:, :, :, 1) / nn
                ! endif
                !
                ! if (abs(cmode) == 1 .or. (abs(cmode) == 2 .and. yrnm_clm /= "")) then
                !       ! call nc_read(fn_vor, 'curlnonl_m' , WORK); curlnonlc     = curlnonlc     + WORK(:, :, :, 1) / nn
                ! endif
                !
                ! if (abs(cmode) == 2 .and. yrnm_clm == "") then
                !     call nc_read(fn_vor, 'errnlsub' , WORK); err_nlsub     = err_nlsub     + WORK(:, :, :, 1) / nn
                ! endif

                if (index(func, "d") /= 0) then
                    call nc_read(fn_vor, 'advu'     , WORK); advu         = advu         + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'advv'     , WORK); advv         = advv         + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'advw'     , WORK); advw         = advw         + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'advVx'    , WORK); advVx        = advVx        + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'advVy'    , WORK); advVy        = advVy        + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'advVz'    , WORK); advVz        = advVz        + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'curlmet'  , WORK); curlmet      = curlmet      + WORK(:, :, :, 1) / nn
                    call nc_read(fn_vor, 'errnldcmp', WORK); err_nldecomp = err_nldecomp + WORK(:, :, :, 1) / nn
                endif
            enddo
        enddo

        if (index(func, "c") /= 0) then
            write(*, fmtm_vor) 'curlnonl: ' , curlnonl (B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'betav: '    , betav    (B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'stretchp: ' , stretchp (B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'errcor: '   , err_cor  (B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'curlpgrad: ', curlpgrad(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'curlhdiff: ', curlhdiff(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'curlvdiff: ', curlvdiff(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'res: '      , res      (B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif

        if (index(func, "d") /= 0) then
            write(*, fmtm_vor) 'advu: ', advu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'advv: ', advv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'advw: ', advw(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'advVx: ', advVx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'advVy: ', advVy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'advVz: ', advVz(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'curlmet: ', curlmet(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) 'errnldcmp: ', err_nldecomp(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
    endsubroutine

    ! Get yyyymmdd from year and doy.
    ! iyr and idoy is overriden by yrnm_clm and avnm_clm.
    ! if idoy < 0, the string returns only the month part
    ! if idoy = 0, the string returns only the year part
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

    ! Converting month and day to day of the year
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

    subroutine output_sf(func, fn_zeta)
        implicit none
        character(len=*), intent(in) :: func, fn_zeta
        integer :: ncid_zeta
        type(nc_varid) :: varids

        call create_output(func, fn_zeta, ncid_zeta, varids)
        call write_output(func, ncid_zeta, varids)
        call close_output(func, ncid_zeta)
    endsubroutine

    subroutine output_me(func_m, fn_me)
        use netcdf
        implicit none
        character(len = *), intent(in) :: func_m, fn_me
        integer :: ncid, stat_create, stat_defdim, stat_defvar, stat_putatt, stat_inqvar, &
                   stat_getvar, stat_putvar, stat_io
        integer :: dimid_lon, dimid_lat, dimid_dep, dimid_time
        integer :: varid_lon, varid_lat, varid_dep, varid_time
        type(nc_varid) :: varids

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A, A)') 'Creating output file: ', trim(fn_me)
        write(*, '(A)') '  Start netcdf define ...'

        stat_create = nf90_create(trim(fn_me), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid)

        ! Dimension
        stat_defdim = nf90_def_dim(ncid, "nlon", B%nx, dimid_lon)
        stat_defdim = nf90_def_dim(ncid, "nlat", B%ny, dimid_lat)
        stat_defdim = nf90_def_dim(ncid, "z_t" , B%nz, dimid_dep)
        stat_defdim = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dimid_time)

        ! Coordinates
        stat_defvar = nf90_def_var(ncid, "TLONG", NF90_FLOAT, &
           (/dimid_lon, dimid_lat/), varid_lon)
        stat_defvar = nf90_def_var(ncid, "TLAT",  NF90_FLOAT, &
           (/dimid_lon, dimid_lat/), varid_lat)
        stat_defvar = nf90_def_var(ncid, "z_t" ,  NF90_FLOAT, dimid_dep , varid_dep )
        stat_defvar = nf90_def_var(ncid, "time" , NF90_FLOAT, dimid_time, varid_time)

        ! Variables
        stat_defvar = nf90_def_var(ncid, "curlnonl" , nc_xtype, &
          (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlnonl )
        stat_putatt = nf90_put_att(ncid, varids%curlnonl , "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varids%curlnonl , "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varids%curlnonl , "long_name", "Curl of nonlinear term (rhs)")
        stat_putatt = nf90_put_att(ncid, varids%curlnonl , "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "curlpgrad", nc_xtype, &
          (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlpgrad)
        stat_putatt = nf90_put_att(ncid, varids%curlpgrad, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varids%curlpgrad, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varids%curlpgrad, "long_name", "Curl of pressure gradient term (rhs)")
        stat_putatt = nf90_put_att(ncid, varids%curlpgrad, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "res"      , nc_xtype, &
            (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%res)
        stat_putatt = nf90_put_att(ncid, varids%res      , "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varids%res      , "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varids%res      , "long_name", "Residual (lhs)")
        stat_putatt = nf90_put_att(ncid, varids%res      , "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "curlhdiff", nc_xtype, &
          (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlhdiff)
        stat_putatt = nf90_put_att(ncid, varids%curlhdiff, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varids%curlhdiff, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varids%curlhdiff, "long_name", "Curl of horizontal diffusion (rhs)")
        stat_putatt = nf90_put_att(ncid, varids%curlhdiff, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "curlvdiff", nc_xtype, &
          (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlvdiff)
        stat_putatt = nf90_put_att(ncid, varids%curlvdiff, "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varids%curlvdiff, "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varids%curlvdiff, "long_name", "Curl of vertical diffusion (rhs)")
        stat_putatt = nf90_put_att(ncid, varids%curlvdiff, "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "betav"    , nc_xtype, &
          (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%betav    )
        stat_putatt = nf90_put_att(ncid, varids%betav    , "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varids%betav    , "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varids%betav    , "long_name", "Advection of planetary vorticity term (rhs)")
        stat_putatt = nf90_put_att(ncid, varids%betav    , "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "stretchp",  nc_xtype, &
          (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%stretchp  )
        stat_putatt = nf90_put_att(ncid, varids%stretchp , "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varids%stretchp , "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varids%stretchp , "long_name", "Planetary vorticity stretching term (rhs)")
        stat_putatt = nf90_put_att(ncid, varids%stretchp , "missing_value", MVALUE)

        stat_defvar = nf90_def_var(ncid, "errcor"  ,  nc_xtype, &
          (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%errcor    )
        stat_putatt = nf90_put_att(ncid, varids%errcor   , "Units", "1/s^2")
        stat_putatt = nf90_put_att(ncid, varids%errcor   , "coordinates", "TLONG TLAT z_t time")
        stat_putatt = nf90_put_att(ncid, varids%errcor   , "long_name", "Error from decomposing curl(-fv, fu) (rhs)")
        stat_putatt = nf90_put_att(ncid, varids%errcor   , "missing_value", MVALUE)

        if (index(func_m, "d") /=0) then
            stat_defvar = nf90_def_var(ncid, "errnlsub" , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%errnlsub)
            stat_putatt = nf90_put_att(ncid, varids%errnlsub, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%errnlsub, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%errnlsub, "long_name", "Error from nonlinear term due to calculating offline (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%errnlsub, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advu_m"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advum)
            stat_putatt = nf90_put_att(ncid, varids%advum, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advum, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advum, "long_name", "Mean advection of relative vorticity by zonal velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advum, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advv_m"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advvm)
            stat_putatt = nf90_put_att(ncid, varids%advvm, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advvm, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advvm, "long_name", "Mean advection of relative vorticity by meridional velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advvm, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advw_m"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advwm)
            stat_putatt = nf90_put_att(ncid, varids%advwm, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advwm, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advwm, "long_name", "Mean advection of relative vorticity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advwm, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advVx_m" , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advVxm)
            stat_putatt = nf90_put_att(ncid, varids%advVxm, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advVxm, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advVxm, "long_name", "Mean twisting of zonal voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advVxm, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advVy_m"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advVym)
            stat_putatt = nf90_put_att(ncid, varids%advVym, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advVym, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advVym, "long_name", "Mean twisting of meridional voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advVym, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advVz_m" , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advVzm)
            stat_putatt = nf90_put_att(ncid, varids%advVzm, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advVzm, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advVzm, "long_name", "Mean twisting of vertical voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advVzm, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "curlmet_m"  , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlmetm)
            stat_putatt = nf90_put_att(ncid, varids%curlmetm, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%curlmetm, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%curlmetm, "long_name", "Mean curl of metric term (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%curlmetm, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "errnldcmp_m", nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%errnldcmpm)
            stat_putatt = nf90_put_att(ncid, varids%errnldcmpm, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%errnldcmpm, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%errnldcmpm, "long_name", "Mean error from nonlinear term due to decomposition (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%errnldcmpm, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advu_e"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advu)
            stat_putatt = nf90_put_att(ncid, varids%advu, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advu, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advu, "long_name", "Eddy advection of relative vorticity by zonal velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advu, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advv_e"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advv)
            stat_putatt = nf90_put_att(ncid, varids%advv, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advv, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advv, "long_name", "Eddy advection of relative vorticity by meridional velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advv, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advw_e"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advw)
            stat_putatt = nf90_put_att(ncid, varids%advw, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advw, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advw, "long_name", "Eddy advection of relative vorticity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advw, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advVx_e" , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advVx)
            stat_putatt = nf90_put_att(ncid, varids%advVx, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advVx, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advVx, "long_name", "Eddy twisting of zonal voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advVx, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advVy_e"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advVy)
            stat_putatt = nf90_put_att(ncid, varids%advVy, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advVy, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advVy, "long_name", "Eddy twisting of meridional voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advVy, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "advVz_e" , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advVz)
            stat_putatt = nf90_put_att(ncid, varids%advVz, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%advVz, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%advVz, "long_name", "Eddy twisting of vertical voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%advVz, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "curlmet_e"  , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlmet)
            stat_putatt = nf90_put_att(ncid, varids%curlmet, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%curlmet, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%curlmet, "long_name", "Eddy curl of metric term (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%curlmet, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "errnldcmp_e", nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%errnldcmp)
            stat_putatt = nf90_put_att(ncid, varids%errnldcmp, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%errnldcmp, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%errnldcmp, "long_name", "Eddy error from nonlinear term due to decomposition (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%errnldcmp, "missing_value", MVALUE)
        endif

        if (index(func_m, "a") /=0) then
            stat_defvar = nf90_def_var(ncid, "curlnonl_m"  ,  nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlnonlc)
            stat_putatt = nf90_put_att(ncid, varids%curlnonlc, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%curlnonlc, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%curlnonlc, "long_name", "Mean curl of nonlinear term offline (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%curlnonlc, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid, "curlnonl_e" , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%errnlsub)
            stat_putatt = nf90_put_att(ncid, varids%errnlsub, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid, varids%errnlsub, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid, varids%errnlsub, "long_name", "Eddy curl of nonlinear term (rhs)")
            stat_putatt = nf90_put_att(ncid, varids%errnlsub, "missing_value", MVALUE)
        endif

        stat_create = nf90_enddef(ncid)
        write(*, '(A, I1)') "  Finished netcdf define!", stat_create

        ! Writing cooordinates
        stat_putvar = nf90_put_var(ncid, varid_lat,  tlat , &
           start = (/1, 1/), count = (/B%nx, B%ny/))
        stat_putvar = nf90_put_var(ncid, varid_lon,  tlong, &
           start = (/1, 1/), count = (/B%nx, B%ny/))
        stat_putvar = nf90_put_var(ncid, varid_dep,  z_t)

        write(*, *)
        write(*, '(A)') '  Start writing zeta mean/eddy file'
        stat_putvar = nf90_put_var(ncid, varids%curlnonl  , curlnonl      , &
            start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        stat_putvar = nf90_put_var(ncid, varids%betav     , betav         , &
            start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        stat_putvar = nf90_put_var(ncid, varids%stretchp  , stretchp      , &
            start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        stat_putvar = nf90_put_var(ncid, varids%errcor    , err_cor       , &
            start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        stat_putvar = nf90_put_var(ncid, varids%curlpgrad , curlpgrad     , &
            start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        stat_putvar = nf90_put_var(ncid, varids%curlhdiff , curlhdiff     , &
            start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        stat_putvar = nf90_put_var(ncid, varids%curlvdiff , curlvdiff     , &
            start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        stat_putvar = nf90_put_var(ncid, varids%res       , res           , &
            start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))

      if (index(func_m, "d") /= 0) then
            stat_putvar = nf90_put_var(ncid, varids%errnlsub  , err_nlsub     , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))

            stat_putvar = nf90_put_var(ncid, varids%advum      , advu_m          , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%advvm      , advv_m          , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%advwm      , advw_m          , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%advVxm     , advVx_m         , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%advVym     , advVy_m         , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%advVzm     , advVz_m         , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%curlmetm   , curlmet_m       , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%errnldcmpm , err_nldecomp_m  , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))

            stat_putvar = nf90_put_var(ncid, varids%advu     , advu  - advu_m , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%advv     , advv  - advv_m , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%advw     , advw  - advw_m, &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%advVx    , advVx - advVx_m, &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%advVy    , advVy - advVy_m, &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%advVz    , advVz - advVz_m, &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%curlmet  , curlmet - curlmet_m, &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%errnldcmp, err_nldecomp - err_nldecomp_m, &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        endif

        if (index(func_m, "a") /=0) then
            stat_putvar = nf90_put_var(ncid, varids%errnlsub  , curlnonl - curlnonl_m     , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid, varids%curlnonlc , curlnonl_m     , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        endif

        stat_io = nf90_close(ncid)
        write(*, '(A, I1)') '  Finished writing zeta mean/eddy file ', stat_io
    endsubroutine

    subroutine create_output(func, fn_zeta, ncid_zeta, varids)
        use netcdf
        implicit none
        character(len=*), intent(in) :: func, fn_zeta
        integer, intent(inout) :: ncid_zeta
        type(nc_varid), intent(inout) :: varids
        integer :: stat_create, stat_defdim, stat_defvar, stat_putatt, stat_inqvar, &
                   stat_getvar, stat_putvar, stat_io
        integer :: dimid_lon, dimid_lat, dimid_dep, dimid_time
        integer :: varid_lon, varid_lat, varid_dep, varid_time

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A, A)') 'Creating output file: ', trim(fn_zeta)
        write(*, '(A)') '  Start netcdf define ...'

        stat_create = nf90_create(trim(fn_zeta), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid_zeta)

        ! Dimension
        stat_defdim = nf90_def_dim(ncid_zeta, "nlon", B%nx, dimid_lon)
        stat_defdim = nf90_def_dim(ncid_zeta, "nlat", B%ny, dimid_lat)
        stat_defdim = nf90_def_dim(ncid_zeta, "z_t" , B%nz, dimid_dep)
        stat_defdim = nf90_def_dim(ncid_zeta, "time", NF90_UNLIMITED, dimid_time)

        ! Coordinates
        stat_defvar = nf90_def_var(ncid_zeta, "TLONG", NF90_FLOAT, &
           (/dimid_lon, dimid_lat/), varid_lon)
        stat_defvar = nf90_def_var(ncid_zeta, "TLAT",  NF90_FLOAT, &
           (/dimid_lon, dimid_lat/), varid_lat)
        stat_defvar = nf90_def_var(ncid_zeta, "z_t" ,  NF90_FLOAT, dimid_dep , varid_dep )
        stat_defvar = nf90_def_var(ncid_zeta, "time" , NF90_FLOAT, dimid_time, varid_time)

        ! Variables
        if (index(func, "c") /= 0) then
            stat_defvar = nf90_def_var(ncid_zeta, "curlnonl" , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlnonl )
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlnonl , "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlnonl , "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlnonl , "long_name", "Curl of nonlinear term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlnonl , "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "curlpgrad", nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlpgrad)
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlpgrad, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlpgrad, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlpgrad, "long_name", "Curl of pressure gradient term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlpgrad, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "res"      , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%res)
            stat_putatt = nf90_put_att(ncid_zeta, varids%res      , "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%res      , "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%res      , "long_name", "Residual (lhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%res      , "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "curlhdiff", nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlhdiff)
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlhdiff, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlhdiff, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlhdiff, "long_name", "Curl of horizontal diffusion (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlhdiff, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "curlvdiff", nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlvdiff)
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlvdiff, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlvdiff, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlvdiff, "long_name", "Curl of vertical diffusion (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlvdiff, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "betav"    , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%betav    )
            stat_putatt = nf90_put_att(ncid_zeta, varids%betav    , "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%betav    , "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%betav    , "long_name", "Advection of planetary vorticity term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%betav    , "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "stretchp",  nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%stretchp  )
            stat_putatt = nf90_put_att(ncid_zeta, varids%stretchp , "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%stretchp , "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%stretchp , "long_name", "Planetary vorticity stretching term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%stretchp , "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "errcor"  ,  nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%errcor    )
            stat_putatt = nf90_put_att(ncid_zeta, varids%errcor   , "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errcor   , "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errcor   , "long_name", "Error from decomposing curl(-fv, fu) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errcor   , "missing_value", MVALUE)
        endif

        if (index(func, "a") /= 0) then
            stat_defvar = nf90_def_var(ncid_zeta, "curlnonl_m"  ,  nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlnonlc)
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlnonlc, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlnonlc, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlnonlc, "long_name", "Mean curl of nonlinear term offline (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlnonlc, "missing_value", MVALUE)
        endif

        if (index(func, "e") /= 0 .and. T%yrnm_clm /= "") then
            stat_defvar = nf90_def_var(ncid_zeta, "curlnonl_e" , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%errnlsub)
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnlsub, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnlsub, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnlsub, "long_name", "Eddy curl of nonlinear term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnlsub, "missing_value", MVALUE)
        endif

        if (index(func, "e") /= 0 .and. T%yrnm_clm == "") then
            stat_defvar = nf90_def_var(ncid_zeta, "errnlsub" , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%errnlsub)
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnlsub, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnlsub, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnlsub, "long_name", "Error from nonlinear term due to calculating offline (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnlsub, "missing_value", MVALUE)
        endif

        if (index(func, "d") /= 0) then
            stat_defvar = nf90_def_var(ncid_zeta, "advu"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advu)
            stat_putatt = nf90_put_att(ncid_zeta, varids%advu, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advu, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advu, "long_name", "Advection of relative vorticity by zonal velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advu, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "advv"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advv)
            stat_putatt = nf90_put_att(ncid_zeta, varids%advv, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advv, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advv, "long_name", "Advection of relative vorticity by meridional velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advv, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "advw"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advw)
            stat_putatt = nf90_put_att(ncid_zeta, varids%advw, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advw, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advw, "long_name", "Advection of relative vorticity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advw, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "advVx" , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advVx)
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVx, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVx, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVx, "long_name", "Twisting of zonal voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVx, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "advVy"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advVy)
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVy, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVy, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVy, "long_name", "Twisting of meridional voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVy, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "advVz"     , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%advVz)
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVz, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVz, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVz, "long_name", "Twisting of vertical voricity by vertical velocity (flux form) (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%advVz, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "curlmet"  , nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%curlmet)
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlmet, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlmet, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlmet, "long_name", "Curl of metric term (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%curlmet, "missing_value", MVALUE)

            stat_defvar = nf90_def_var(ncid_zeta, "errnldcmp", nc_xtype, &
              (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), varids%errnldcmp)
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnldcmp, "Units", "1/s^2")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnldcmp, "coordinates", "TLONG TLAT z_t time")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnldcmp, "long_name", "Error from nonlinear term due to decomposition (rhs)")
            stat_putatt = nf90_put_att(ncid_zeta, varids%errnldcmp, "missing_value", MVALUE)
        endif

        stat_create = nf90_enddef(ncid_zeta)
        write(*, '(A, I1)') "  Finished netcdf define: ", stat_create

        ! Writing cooordinates
        stat_putvar = nf90_put_var(ncid_zeta, varid_lat,  tlat , &
           start = (/1, 1/), count = (/B%nx, B%ny/))
        stat_putvar = nf90_put_var(ncid_zeta, varid_lon,  tlong, &
           start = (/1, 1/), count = (/B%nx, B%ny/))
        stat_putvar = nf90_put_var(ncid_zeta, varid_dep,  z_t)
    endsubroutine

    subroutine write_output(func, ncid_zeta, varids)
        use netcdf
        implicit none
        character(len=*), intent(in) :: func
        integer, intent(in) :: ncid_zeta
        type(nc_varid), intent(in) :: varids
        integer :: stat_putvar

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Start writing zeta equation file ...'

        if (index(func, "c") /=0) then
            stat_putvar = nf90_put_var(ncid_zeta, varids%curlnonl , curlnonl , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%betav    , betav    , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%stretchp , stretchp , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%errcor   , err_cor  , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%curlpgrad, curlpgrad, &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%curlhdiff, curlhdiff, &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%curlvdiff, curlvdiff, &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%res      , res      , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        endif

        if (index(func, "am") /=0) then
            stat_putvar = nf90_put_var(ncid_zeta, varids%curlnonlc , curladv + curlmet , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        endif

        if (index(func, "e") /=0) then
            stat_putvar = nf90_put_var(ncid_zeta, varids%errnlsub , err_nlsub   , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        endif

        if (index(func, "d") /=0) then
            stat_putvar = nf90_put_var(ncid_zeta, varids%advu     , advu        , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%advv     , advv        , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%advw     , advw        , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%advVx    , advVx       , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%advVy    , advVy       , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%advVz    , advVz       , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%curlmet  , curlmet     , &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
            stat_putvar = nf90_put_var(ncid_zeta, varids%errnldcmp, err_nldecomp, &
                start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
        endif
        write(*, '(A)') '  Finished writing'
    endsubroutine

    subroutine close_output(func, ncid_zeta)
        use netcdf
        implicit none
        character(len=*), intent(in) :: func
        integer, intent(in) :: ncid_zeta
        integer :: stat_io

        write(*, *)
        write(*, '(A)') 'Closing output files'
        stat_io = nf90_close(ncid_zeta)
    endsubroutine

    ! Check if files for mean exists and contains required variables depending on the calculation mode
    ! Returns a logical variable
    subroutine check_meanfile(func_m, fn, filestat)
        use netcdf
        implicit none
        character(len=*), intent(in) :: func_m, fn
        logical, intent(inout) :: filestat
        integer :: iostat, ncid, varid
        integer, dimension(:), allocatable :: varstat

        filestat = .False.
        iostat = nf90_open(trim(fn), NF90_NOWRITE, ncid)

        if (iostat == nf90_noerr) then
            if (index(func_m, 'a') /= 0) then
                allocate(varstat(1))
                varstat(1) = nf90_inq_varid(ncid = ncid, name = 'curlnonl_m', varid = varid)
            endif
            if (index(func_m, 'd') /= 0) then
                allocate(varstat(9))
                varstat(1) = nf90_inq_varid(ncid = ncid, name = 'advu', varid = varid)
                varstat(2) = nf90_inq_varid(ncid = ncid, name = 'advv', varid = varid)
                varstat(3) = nf90_inq_varid(ncid = ncid, name = 'advw', varid = varid)
                varstat(4) = nf90_inq_varid(ncid = ncid, name = 'advu', varid = varid)
                varstat(5) = nf90_inq_varid(ncid = ncid, name = 'advVx', varid = varid)
                varstat(6) = nf90_inq_varid(ncid = ncid, name = 'advVy', varid = varid)
                varstat(7) = nf90_inq_varid(ncid = ncid, name = 'advVz', varid = varid)
                varstat(8) = nf90_inq_varid(ncid = ncid, name = 'curlmet', varid = varid)
                varstat(9) = nf90_inq_varid(ncid = ncid, name = 'errnldcmp', varid = varid)
            endif
            iostat = nf90_close(ncid)

            if (allocated(varstat) .and. all(varstat == nf90_noerr)) then
                filestat = .True.
            endif
        endif
    endsubroutine
endmodule
