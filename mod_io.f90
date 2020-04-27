module io
    use params
    use zeta, only : uc, vc, wc, ssh, advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy, &
                     ueu, uev, vnu, vnv, wtu, wtv, &
                     tlat, tlong, z_t, vgrp, zetavar, vargrp
                    !  curlnonl, betav, stretchp, err_cor, curlpgrad, curlhdiff, curlvdiff, res, &
                    !  curladv, curlmet, curladvu, curladvv, curladvw, err_nlsub, advu, advv, advw, advVx, advVy, advVz, err_nldecomp, &
                    !  curladvf, & 
                    !  advu_x, advv_x, advw_x, advVx_x, advVy_x, advVz_x, errnl_ux, errnl_vx, errnl_wx, &
                    !  advu_y, advv_y, advw_y, advVx_y, advVy_y, advVz_y, errnl_uy, errnl_vy, errnl_wy
             
    implicit none
    private
    public :: load_params, get_yyyymmdd, loadave_mom_sf, output_sf,  &
              loadave_vor_sf, output_me, check_meanfile

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

    !!--------------------------------------------------------------------------
    ! * Calculation modes *
    logical, public :: ifcurl = .False., ifdecomp = .False., flxtwi = .True.
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

        namelist /calcmode/ ifcurl, ifdecomp, flxtwi
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
        if (ifdecomp) then 
            write(*, '(A40, L)'), "  Using flux form twisting term? ", flxtwi
        endif
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
                write(*, '(A)') "  WARNING: input files are climatology and 'ifmeanclm' is by definiation turned on for mean/eddy decomposition."
                write(*, *)
                ifmeanclm = .True.
            endif
            if     (      ifcurl .and. .not. ifdecomp) then
                if (ifmeaneddy) then 
                    write(*, '(A)') 'Curl of momentum equation and mean/eddy decomposition for the nonlinear term'
                    write(*, '(A)') '  ERR_nlsub == curlnonl_eddy'
                    func = "came"
                else
                    write(*, '(A)') 'Curl of momentum equation for each input file'
                    func = "c"
                endif
            elseif (      ifcurl .and.       ifdecomp) then
                write(*, '(A)') 'Curl of momentum equation and full decomposition for the nonlinear term'
                write(*, '(A)') '  ERR_nlsub == curlnonl_eddy and only mean terms have decompostion.'
                func = "cdme"
                if (.not. ifmeaneddy) then
                    write(*,'(A)') "  Note: 'ifmeaneddy' was set to False. A decomposition of nonlinear terms of the climatology &
                             is essentially equivalent to a mean/eddy decomposition."
                endif
            elseif (.not. ifcurl .and. .not. ifdecomp) then 
                write(*, '(A)') "Both 'ifcurl' and 'ifdecomp' are set to False. Nothing to do."
                stop
            ! See below for (.not. ifcurl .and. ifdecomp)
            endif
        else
            if     (      ifcurl .and. .not. ifdecomp) then
                func = "c"
                if (ifmeaneddy) then 
                    write(*, '(A)') 'Curl of momentum equation and mean/eddy decomposition for the nonlinear term'
                    func_m = "am"
                    func_me = func
                else 
                    write(*, '(A)') 'Curl of momentum equation for each input file'
                    func_m = ""
                    func_me = ""
                endif
            elseif (      ifcurl .and.       ifdecomp) then
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
            elseif (.not. ifcurl .and. .not. ifdecomp) then
                if (ifmeaneddy) then
                    write(*, '(A)') "Mean/eddy decomposition will be calculated based on given decomposed zeta equation files"
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
                else 
                    write(*, '(A)') " 'ifcurl', 'ifdecomp' and 'ifmeaneddy' are set to False. Nothing to do."
                    stop
                endif
            endif
        endif
        if (.not. ifcurl .and. ifdecomp ) then
            write(*, '(A)'), 'Decomposition of the nonlinear term only.'
            func = "dm"
            func_m = ""
            func_me = ""
            if (ifmeaneddy) then
                write(*,'(A)') "  WARNING: Cannot do mean/eddy decomposition without curl of momentum equation."
            endif
        endif

        ! Inserting "#" for the nonflux twisting option
        if (index(func, "d") /=0 .and. .not. flxtwi) then
            func = func(1:index(func, "d"))//"#"//func(index(func, "d")+1:9)
        endif
        if (index(func_m, "d") /=0 .and. .not. flxtwi) then
            func_m = func_m(1:index(func_m, "d"))//"#"//func_m(index(func_m, "d")+1:9)
        endif        
        if (index(func_me, "d") /=0 .and. .not. flxtwi) then
            func_me = func_me(1:index(func_me, "d"))//"#"//func_me(index(func_me, "d")+1:9)
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
        ! Monthly mean: doy_list = (/ -number of the month /)
        ! Daily mean: doy_list = doy
        idoy = 0
        do imn = 1, nmn, 1
           do ida = 1, nda, 1
               if (mnlist(imn) > 0 .and. dalist(ida) > eom(mnlist(imn))) exit
               idoy = idoy + 1
               if (dalist(ida) == 0 .and. mnlist(imn) /= 0) then
                   doylist_full(idoy) = -mnlist(imn)
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
        ! When mean/eddy decomposition is turned on, secs are based on its frequency/period.
        ! For daily input:
        !   seclist overrides doylist with "m", "a", secst/seced and seclist_in_st/ed options.
        !   It uses min(doyslist) and max(doylist) as lower and upper bounds with nda_sec option
        ! For monthly input:
        !   Only "a" and "m" options are accepted. Other options are invalid.
        !   seclist overrides doylist with "a" option, using the full 12 month. 
        ! For annual input:
        !   Only "a" option is accepted. 
        if (trim(func_me)=="") then
            T%nsec = T%ndoy
            allocate(T%seclist(T%nsec, 2))
            T%seclist(:, 1) = (/ (T%doylist(ida), ida = 1, T%ndoy) /)
            T%seclist(:, 2) = T%seclist(:, 1)
        else
            if (all(T%doylist > 0)) then
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
                elseif (seced >= secst .and. secst > 0) then
                    T%nsec = 1
                    allocate(T%seclist(T%nsec, 2))
                    allocate(T%meannm(T%nsec))
                    T%seclist(1, 1) = secst
                    T%seclist(1, 2) = seced
                    write(T%meannm(1), '(A, I0.3, A, I0.3)') 'd', secst, '-', seced
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
                else
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
            elseif (all(T%doylist == 0)) then
                if (meanfreq=="a") then 
                    T%nsec = 1
                    allocate(T%seclist(T%nsec, 2))
                    allocate(T%meannm(T%nsec))
                    T%seclist(1, 1) = 1
                    T%seclist(1, 2) = 365             
                    write(T%meannm(1), '(A)') 'ann' 
                else
                    write(*, '(A)') "ERROR: Mean/Eddy frequency conflicts annual input."
                    STOP                    
                endif
            else
                if (meanfreq=="a") then 
                    T%nsec = 1
                    allocate(T%seclist(T%nsec, 2))
                    allocate(T%meannm(T%nsec))
                    T%seclist(1, 1) = -1
                    T%seclist(1, 2) = -12             
                    write(T%meannm(1), '(A)') 'ann' 
                elseif (meanfreq=="m") then 
                    T%nsec = T%ndoy
                    allocate(T%seclist(T%nsec, 2))
                    allocate(T%meannm(T%nsec))
                    T%seclist(:, 1) = (/ (T%doylist(ida), ida = 1, T%ndoy) /)
                    T%seclist(:, 2) = T%seclist(:, 1)
                    do isec = 1, T%nsec
                        write(T%meannm(isec), '(I0.2)') -T%seclist(isec, 1)
                    enddo
                else
                    write(*, '(A)') "ERROR: Mean/Eddy frequency conflicts monthly input."
                    STOP                    
                endif
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

        if (len(trim(fn_vor_dlm))  == 0) fn_vor_dlm  = trim(fn_mom_dlm)

        if (len(trim(fn_vorm_dlm)) == 0) fn_vorm_dlm = trim(fn_mom_dlm)
        if (len(trim(fn_vore_dlm)) == 0) fn_vore_dlm = trim(fn_mom_dlm)

        if (len(trim(fn_vorm_dir)) == 0) fn_vorm_dir = trim(fn_vor_dir)
        if (len(trim(fn_vore_dir)) == 0) fn_vore_dir = trim(fn_vor_dir)

        if (len(trim(fn_vorm_pfx)) == 0) fn_vorm_pfx = trim(fn_vor_pfx)
        if (len(trim(fn_vore_pfx)) == 0) fn_vore_pfx = trim(fn_vor_pfx)

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
        fn_vorm%sfx = trim(fn_vor_sfx) // trim(fn_vorm_sfx)
        fn_vorm%dlm = fn_vorm_dlm

        fn_vore%dir = fn_vore_dir
        fn_vore%pfx = fn_vore_pfx
        fn_vore%sfx = trim(fn_vor_sfx) // trim(fn_vore_sfx)
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
        write(*, '(A9)', advance="no") "  Start: "
        do isec = 1, T%nsec
            write(*, '(I3, A)', advance="no") T%seclist(isec, 1), ' '
        enddo
        write(*, *)
        write(*, '(A9)', advance="no") "  End: "
        do isec = 1, T%nsec
            write(*, '(I3, A)', advance="no") T%seclist(isec, 2), ' '
        enddo
        write(*, *)

        if (ifmeaneddy) then 
            write(*, *)
            write(*, '(A)') "  Mean/eddy section names: "
            do isec = 1, T%nsec
                write(*, '(A, A)', advance="no") trim(T%meannm(isec)), ' '
            enddo        
            write(*, *)
        endif 
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
                fn_mom = trim(fn%dir) // trim(fn%pfx) // trim(yyyymmdd) // trim(fn%sfx) // '.nc'

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

    ! Load and average vorticity terms from single files
    subroutine loadave_vor_sf(func, yrlist, doylist, yrnm_clm, avnm_clm, fn)
        use ncio, only : nc_read
        implicit none
        character(len=*), intent(in) :: func
        integer, intent(in) :: yrlist(:), doylist(:)
        character(len=*), intent(in) :: yrnm_clm, avnm_clm
        type(filename), intent(in) :: fn
        character :: fn_vor*300, yyyymmdd*20
        integer :: iyr, idoy, ig, iv
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
                fn_vor = trim(fn%dir) // trim(fn%pfx) // trim(yyyymmdd) // trim(fn%sfx) // '.nc'
                write(*, '(A, A)') '  Load from file: ', trim(fn_vor)
    
                do ig = 1, size(vgrp)
                    if (vgrp(ig)%key) then 
                        do iv = 1, size(vgrp(ig)%vlist)
                            call nc_read(fn_vor, trim(vgrp(ig)%vlist(iv)%name), WORK); 
                            vgrp(ig)%vlist(iv)%value = vgrp(ig)%vlist(iv)%value + WORK(:, :, :, 1) / nn
                        enddo
                    endif
                enddo
            enddo
        enddo
    
        do ig = 1, size(vgrp)
            if (vgrp(ig)%key) then 
                do iv = 1, size(vgrp(ig)%vlist)
                    write(*, fmtm_vor) trim(vgrp(ig)%vlist(iv)%name)// ': ' , &
                        vgrp(ig)%vlist(iv)%value(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
                enddo
            endif
        enddo
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
                write(mmdd, '(A, I0.2)') trim(dlm), -idoy
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

    ! Outputing single timestep file (Wrapper of the create_output, write_output and close_output)
    subroutine output_sf(fn_zeta)
        implicit none
        character(len=*), intent(in) :: fn_zeta
        integer :: ncid_zeta

        call create_output(fn_zeta, ncid_zeta)
        call write_output(ncid_zeta)
        call close_output(ncid_zeta)
    endsubroutine

    subroutine create_output(fn_zeta, ncid_zeta)
        use netcdf
        implicit none
        character(len=*), intent(in) :: fn_zeta
        integer, intent(inout) :: ncid_zeta
        integer :: dimid_lon, dimid_lat, dimid_dep, dimid_time
        integer :: varid_lon, varid_lat, varid_dep, varid_time
        integer :: stat, ig, iv
    
        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
    
        stat = nf90_create(trim(fn_zeta), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid_zeta)
        write(*, '(A, A, A, I2)') 'Creating output file: ', trim(fn_zeta), ': ', stat
    
        ! Dimension
        stat = nf90_def_dim(ncid_zeta, "nlon", B%nx, dimid_lon)
        stat = nf90_def_dim(ncid_zeta, "nlat", B%ny, dimid_lat)
        stat = nf90_def_dim(ncid_zeta, "z_t" , B%nz, dimid_dep)
        stat = nf90_def_dim(ncid_zeta, "time", NF90_UNLIMITED, dimid_time)
    
        ! Coordinates
        stat = nf90_def_var(ncid_zeta, "TLONG", NF90_FLOAT, (/dimid_lon, dimid_lat/), varid_lon)
        stat = nf90_def_var(ncid_zeta, "TLAT",  NF90_FLOAT, (/dimid_lon, dimid_lat/), varid_lat)
        stat = nf90_def_var(ncid_zeta, "z_t" ,  NF90_FLOAT, dimid_dep , varid_dep )
        stat = nf90_def_var(ncid_zeta, "time" , NF90_FLOAT, dimid_time, varid_time)
    
        ! For climatological input
        if (trim(T%yrnm_clm) /= "") then
            ! curladv
            vgrp(2)%vlist(1)%long_name = '(mean) ' // trim(vgrp(2)%vlist(1)%long_name)
            ! curlmet
            vgrp(4)%vlist(1)%long_name = '(mean) ' // trim(vgrp(4)%vlist(1)%long_name)
            ! errsub
            vgrp(7)%vlist(1)%long_name = '(eddy) curl of nonlinear (rhs)'
            ! decomposed adv
            do iv = 1, size(vgrp(5)%vlist)
                vgrp(5)%vlist(iv)%long_name = '(mean) ' // trim(vgrp(5)%vlist(iv)%long_name)
            enddo
        endif
    
        do ig = 1, size(vgrp)
            if ( vgrp(ig)%key ) then 
                do iv = 1, size(vgrp(ig)%vlist)
                    stat = nf90_def_var(ncid_zeta, vgrp(ig)%vlist(iv)%name, nc_xtype, &
                      (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), vgrp(ig)%vlist(iv)%varid)
                    stat = nf90_put_att(ncid_zeta, vgrp(ig)%vlist(iv)%varid, "Units", trim(vgrp(ig)%vlist(iv)%Units))
                    stat = nf90_put_att(ncid_zeta, vgrp(ig)%vlist(iv)%varid, "coordinates", trim(vgrp(ig)%vlist(iv)%coordinates))
                    stat = nf90_put_att(ncid_zeta, vgrp(ig)%vlist(iv)%varid, "long_name", trim(vgrp(ig)%vlist(iv)%long_name))
                    stat = nf90_put_att(ncid_zeta, vgrp(ig)%vlist(iv)%varid, "missing_value", MVALUE)
                enddo
            endif
        enddo
        stat = nf90_enddef(ncid_zeta)
        write(*, '(A, I2)') "  Finished netcdf define: ", stat
    
        ! Writing cooordinates
        stat = nf90_put_var(ncid_zeta, varid_lat, tlat , start = (/1, 1/), count = (/B%nx, B%ny/))
        stat = nf90_put_var(ncid_zeta, varid_lon, tlong, start = (/1, 1/), count = (/B%nx, B%ny/))
        stat = nf90_put_var(ncid_zeta, varid_dep,  z_t)
    endsubroutine

    subroutine write_output(ncid_zeta)
        use netcdf
        implicit none
        integer, intent(in) :: ncid_zeta
        integer :: stat, ig, iv
    
        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Writing zeta equation file'
    
        do ig = 1, size(vgrp)
            if ( vgrp(ig)%key ) then 
                do iv = 1, size(vgrp(ig)%vlist)
                    stat = nf90_put_var(ncid_zeta, vgrp(ig)%vlist(iv)%varid, vgrp(ig)%vlist(iv)%value, &
                        start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
                    write(*, '(A2, A9, A, I2)') '  ', trim(vgrp(ig)%vlist(iv)%name), ': ', stat
                enddo
            endif
        enddo
    endsubroutine
    
    subroutine close_output(ncid_zeta)
        use netcdf
        implicit none
        integer, intent(in) :: ncid_zeta
        integer :: stat

        stat = nf90_close(ncid_zeta)
        write(*, '(A, I2)') 'Closing output file: ', stat
    endsubroutine

    ! Calculating and outputing mean/eddy file
    subroutine output_me(fn_m, fn_me)
        use netcdf
        use ncio, only : nc_read
        implicit none
        character(len = *), intent(in) :: fn_m, fn_me
        type(zetavar), dimension(:), target :: vl_aM(1), vl_mM(1), vl_dM(7)
        type(vargrp), dimension(3) :: vgrp_m
        integer, dimension(:), parameter :: idx_regular(2) = (/1, 7/), idx_me(3) = (/2, 4, 5/)
        integer :: dimid_lon, dimid_lat, dimid_dep, dimid_time
        integer :: varid_lon, varid_lat, varid_dep, varid_time
        integer :: ncid, stat, ig, iv, idx
        character(len = 10) :: vnm
        real(kind=kd_r), dimension(B%nx, B%ny, B%nz, 1) :: WORK
    
        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
    
        ! Create a vgrp_m for the mean components. vgrp is used for the eddy part
        !   For simplicity, debug groups are not included
        vl_aM(1)%name = trim(vgrp(2)%vlist(1)%name) // '_m'
        vl_aM(1)%long_name = '(mean) ' // trim(vgrp(2)%vlist(1)%long_name)
    
        vl_mM(1)%name = trim(vgrp(4)%vlist(1)%name) // '_m'
        vl_mM(1)%long_name = '(mean) ' // trim(vgrp(4)%vlist(1)%long_name)
    
        do iv = 1, size(vl_dM)
            vl_dM(iv)%name = trim(vgrp(5)%vlist(iv)%name) // '_m'
            vl_dM(iv)%long_name = '(mean) ' // trim(vgrp(5)%vlist(iv)%long_name)
        enddo
    
        vgrp_m(1)%name = "aM"; vgrp_m(1)%vlist => vl_aM
        vgrp_m(2)%name = "mM"; vgrp_m(2)%vlist => vl_mM
        vgrp_m(3)%name = "dM"; vgrp_m(3)%vlist => vl_dM
        do idx = 1, size(idx_me)
            ig = idx_me(idx)
            if ( vgrp(ig)%key ) vgrp_m(idx)%key = .True.
        enddo

        ! Modify nonlinear term names
        ! curladv
        vgrp(2)%vlist(1)%name = trim(vgrp(2)%vlist(1)%name) // '_e'
        vgrp(2)%vlist(1)%long_name = '(eddy) ' // trim(vgrp(2)%vlist(1)%long_name)
        ! curlmet
        vgrp(4)%vlist(1)%name = trim(vgrp(4)%vlist(1)%name) // '_e'
        vgrp(4)%vlist(1)%long_name = '(eddy) ' // trim(vgrp(4)%vlist(1)%long_name)
        ! errsub
        vgrp(7)%vlist(1)%long_name = '(eddy) ' // trim(vgrp(7)%vlist(1)%long_name)
        ! decomposed adv
        do iv = 1, size(vgrp(5)%vlist)
            vgrp(5)%vlist(iv)%name = trim(vgrp(5)%vlist(iv)%name) // '_e'
            vgrp(5)%vlist(iv)%long_name = '(eddy) ' // trim(vgrp(5)%vlist(iv)%long_name)
        enddo
    
        stat = nf90_create(trim(fn_me), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid)
        write(*, '(A, A, A, I2)') 'Creating output file: ', trim(fn_me), ': ', stat
        write(*, '(A, A)') '  Using mean field file: ', trim(fn_m)
    
        ! Dimension
        stat = nf90_def_dim(ncid, "nlon", B%nx, dimid_lon)
        stat = nf90_def_dim(ncid, "nlat", B%ny, dimid_lat)
        stat = nf90_def_dim(ncid, "z_t" , B%nz, dimid_dep)
        stat = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dimid_time)
    
        ! Coordinates
        stat = nf90_def_var(ncid, "TLONG", NF90_FLOAT, &
           (/dimid_lon, dimid_lat/), varid_lon)
        stat = nf90_def_var(ncid, "TLAT",  NF90_FLOAT, &
           (/dimid_lon, dimid_lat/), varid_lat)
        stat = nf90_def_var(ncid, "z_t" ,  NF90_FLOAT, dimid_dep , varid_dep )
        stat = nf90_def_var(ncid, "time" , NF90_FLOAT, dimid_time, varid_time)
       
        do ig = 1, size(vgrp)
            if ( vgrp(ig)%key ) then 
                do iv = 1, size(vgrp(ig)%vlist)
                    stat = nf90_def_var(ncid, vgrp(ig)%vlist(iv)%name, nc_xtype, &
                      (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), vgrp(ig)%vlist(iv)%varid)
                    stat = nf90_put_att(ncid, vgrp(ig)%vlist(iv)%varid, "Units", trim(vgrp(ig)%vlist(iv)%Units))
                    stat = nf90_put_att(ncid, vgrp(ig)%vlist(iv)%varid, "coordinates", trim(vgrp(ig)%vlist(iv)%coordinates))
                    stat = nf90_put_att(ncid, vgrp(ig)%vlist(iv)%varid, "long_name", vgrp(ig)%vlist(iv)%long_name)
                    stat = nf90_put_att(ncid, vgrp(ig)%vlist(iv)%varid, "missing_value", MVALUE)
                enddo
            endif
        enddo
    
        do ig = 1, size(vgrp_m)
            if ( vgrp_m(ig)%key ) then 
                do iv = 1, size(vgrp_m(ig)%vlist)
                    stat = nf90_def_var(ncid, vgrp_m(ig)%vlist(iv)%name, nc_xtype, &
                      (/dimid_lon, dimid_lat, dimid_dep, dimid_time/), vgrp_m(ig)%vlist(iv)%varid)
                    stat = nf90_put_att(ncid, vgrp_m(ig)%vlist(iv)%varid, "Units", trim(vgrp_m(ig)%vlist(iv)%Units))
                    stat = nf90_put_att(ncid, vgrp_m(ig)%vlist(iv)%varid, "coordinates", trim(vgrp_m(ig)%vlist(iv)%coordinates))
                    stat = nf90_put_att(ncid, vgrp_m(ig)%vlist(iv)%varid, "long_name", vgrp_m(ig)%vlist(iv)%long_name)
                    stat = nf90_put_att(ncid, vgrp_m(ig)%vlist(iv)%varid, "missing_value", MVALUE)
                enddo
            endif 
        enddo
    
        stat = nf90_enddef(ncid)
        write(*, '(A, I2)') "  Finished netcdf define: ", stat
    
        ! Writing cooordinates
        stat = nf90_put_var(ncid, varid_lat, tlat , start = (/1, 1/), count = (/B%nx, B%ny/))
        stat = nf90_put_var(ncid, varid_lon, tlong, start = (/1, 1/), count = (/B%nx, B%ny/))
        stat = nf90_put_var(ncid, varid_dep, z_t)
    
        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Writing zeta mean/eddy file'
    
        ! Terms other than nonlinear components
        do idx = 1, size(idx_regular)
            ig = idx_regular(idx)
            if ( vgrp(ig)%key ) then 
                do iv = 1, size(vgrp(ig)%vlist)
                    stat = nf90_put_var(ncid, vgrp(ig)%vlist(iv)%varid, vgrp(ig)%vlist(iv)%value, &
                        start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
                    write(*, '(A2, A9, A, I2)') '  ', trim(vgrp(ig)%vlist(iv)%name), ': ', stat
                enddo
            endif 
        enddo
    
        ! mean and eddy nonlinear
        do idx = 1, size(idx_me)
            ig = idx_me(idx)
            if ( vgrp(ig)%key ) then 
                do iv = 1, size(vgrp(ig)%vlist)
                    vnm = trim(vgrp(ig)%vlist(iv)%name)
                    call nc_read(fn_m, vnm(1:len(trim(vnm))-2), WORK);
                    stat = nf90_put_var(ncid, vgrp_m(idx)%vlist(iv)%varid, WORK(:,:,:,1), &
                        start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
                    write(*, '(A2, A9, A, I2)') '  ', trim(vgrp_m(idx)%vlist(iv)%name), ': ', stat
                    stat = nf90_put_var(ncid, vgrp(ig)%vlist(iv)%varid, vgrp(ig)%vlist(iv)%value - WORK(:,:,:,1), &
                        start = (/1, 1, 1, 1/), count = (/B%nx, B%ny, B%nz, 1/))
                    write(*, '(A2, A9, A, I2)') '  ', trim(vgrp(ig)%vlist(iv)%name), ': ', stat
                enddo
            endif
        enddo
    
        stat = nf90_close(ncid)
        write(*, '(A, I2)') 'Closing output file: ', stat
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
                allocate(varstat(2))
                varstat(1) = nf90_inq_varid(ncid = ncid, name = 'curladv', varid = varid)
                varstat(2) = nf90_inq_varid(ncid = ncid, name = 'curlmet', varid = varid)
            endif
            if (index(func_m, 'd') /= 0) then
                allocate(varstat(9))
                varstat(1) = nf90_inq_varid(ncid = ncid, name = 'advu', varid = varid)
                varstat(2) = nf90_inq_varid(ncid = ncid, name = 'advv', varid = varid)
                varstat(3) = nf90_inq_varid(ncid = ncid, name = 'advw', varid = varid)
                varstat(4) = nf90_inq_varid(ncid = ncid, name = 'twix', varid = varid)
                varstat(5) = nf90_inq_varid(ncid = ncid, name = 'twiy', varid = varid)
                varstat(6) = nf90_inq_varid(ncid = ncid, name = 'twiz', varid = varid)
                varstat(7) = nf90_inq_varid(ncid = ncid, name = 'erradv', varid = varid)
                varstat(8) = nf90_inq_varid(ncid = ncid, name = 'curlmet', varid = varid)
            endif
            iostat = nf90_close(ncid)

            if (allocated(varstat) .and. all(varstat == nf90_noerr)) then
                filestat = .True.
            endif
        endif
    endsubroutine
endmodule
