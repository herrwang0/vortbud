module zeta
    use params, only : kd_r, MVALUE, B, fngrid, fmt_exp, fmt_flt
    implicit none
    private

    ! public :: load_const, zeta_equation, decomp_curladv, init_zetavars, release_zetavars, create_outputfiles
    ! public :: load_const, calc_zeta, init_zetavars, release_zetavars
    public :: load_const, calc_zeta, init_zetavars_input, init_zetavars_output, release_zetavars_input, release_zetavars_output

    ! Input variables
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: uc, vc, wc
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
        ssh, advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
        ueu, uev, vnu, vnv, wtu, wtv

    ! Output variables
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
        curlnonl, betav, stretchp, err_cor, curlpgrad, curlhdiff, curlvdiff, res
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
        curladv, curlmet, err_nlsub, advu, advv, advw, advVx, advVy, advVz, err_nldecomp, curladvf
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: curladvu, curladvv, curladvw
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
        advu_x, advv_x, advw_x, advVx_x, advVy_x, advVz_x, errnl_ux, errnl_vx, errnl_wx, &
        advu_y, advv_y, advw_y, advVx_y, advVy_y, advVz_y, errnl_uy, errnl_vy, errnl_wy

    ! Constants and grid info
    real(kind=kd_r), allocatable, public :: tlat(:,:), tlong(:,:), z_t(:)
    real(kind=kd_r), dimension(:, :)   , allocatable :: ulat, ulong, dxu, dyu, tarea, uarea, huw , hus
    real(kind=kd_r), dimension(:)      , allocatable :: z_w, dz
    real(kind=kd_r), dimension(:, :, :), allocatable :: dzu, dzt
    real(kind=kd_r), dimension(:, :)   , allocatable :: fcor, fcort
    real(kind=kd_r) :: grav
    ! real(kind=kd_r), dimension(:, :)   , allocatable :: dxue, dyue, tareae, &
    !                                                     dxun, dyun, tarean
    logical, dimension(:, :, :), allocatable :: umask, tmask

    real(kind=kd_r), dimension(:,:,:), allocatable :: ue, vn, wt, ume, vme, umn, vmn, umt, vmt
    ! real(kind=kd_r), dimension(:,:,:), allocatable, public :: rrxx, rrxy, rryx, rryy, rrzx, rrzy
    ! real(kind=kd_r), dimension(:,:,:), allocatable, public :: rcuv, rcuu, rcvv, rcvu, rcwv, rcwu
    ! real(kind=kd_r), dimension(:,:,:), allocatable, public :: rsx1, rsx2, rsy1, rsy2, rsz1, rsz2
    ! real(kind=kd_r), dimension(:,:,:), allocatable, public :: vDdivDx, uDdivDy

    ! print format
    character(len = 100) :: fmts_vel, fmtm_vel, fmts_vor, fmtm_vor

    contains
    !!--------------------------------------------------------------------------
    ! * Input functions
    ! null (empty): do nothing.
    ! curl ("c"): curl of zeta equations terms, including two (+1) tersm from Coriolis
    !                (nonlinear, curl of pgrad, hdiff, vdiff, residual,
    !                 betav from Coriolis, stretching from Coriolis, error with Coriolis decomposition)
    ! offline adv ("a"): offline calculation of nonlinear advection term (curladv)
    !   ("a-"): Debug mode, outputting curladvu, curladvv, curladvw instead of curladv
    ! offline met ("m"): offline calculation of nonlinear metric term (curlmet)
    ! decomp adv ("d"): decomposition of nonlinear advection term.
    !   ("d"): flux form twisting term (default)
    !   ("d#"): traditional twisting term
    ! online adv from flux ("f"): "online" calculation of nonlinear advection term through momentum fluxes
    ! + err_nlsub ("e"): difference between curlnonl online and offline (should not be used as a standalone function)
    !!--------------------------------------------------------------------------

    subroutine calc_zeta(func)
        implicit none
        character(len=*), intent(in), optional :: func
        character(len=10) :: func_c

        write(fmts_vor, '(A, A, A)' ) '(A20, ', trim(fmt_exp), ')'
        write(fmtm_vor, '(A, I2, A, A)') '(A20, ', B%zi_dped - B%zi_dpst + 1, trim(fmt_exp), ')'
        write(fmts_vel, '(A, A, A)' ) '(A20, ', trim(fmt_flt), ')'
        write(fmtm_vel, '(A, I2, A, A)') '(A20, ', B%zi_dped - B%zi_dpst + 1, trim(fmt_flt), ')'

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)'), 'Start calculating voriticity equation'

        if (present(func) .and. trim(func) /= "") then
            func_c = func
        else
            call warning_msg("noinput")
            return
        endif
        write(*, '(2A)') "  Function code: ", trim(func_c)

        if (index(func_c, "c") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            write(*, '(A)') '  Start calculating zcurl of momentum equations'
            call zeta_equation()
        endif
        if (index(func_c, "m") /= 0) then
            call calc_curlmet()
        endif
        if (index(func_c, "a") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            if (index(func_c, "a-") /= 0) then
                write(*, '(A)') '  Calculating curl of advection term (offline, w/ u,v,w commponents)'
                call calc_curladv(.True.)
            else
                write(*, '(A)') '  Calculating curl of advection term (offline)'
                call calc_curladv(.False.)
            endif
        endif
        if (index(func_c, "d") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            if (index(func_c, "d#") /= 0) then
                if (index(func_c, "d#-") /= 0) then
                    write(*, '(A)') '  Decomposing curl of advection term (w/ nonflux twisting term, w/ subcomponents)'
                    call decomp_curladv(.False., .True.)
                else
                    write(*, '(A)') '  Decomposing curl of advection term (w/ nonflux twisting term)'
                    call decomp_curladv(.False., .False.)
                endif
            else
                if (index(func_c, "d-") /= 0) then
                    write(*, '(A)') '  Decomposing curl of advection term (w/ flux twisting term, w/ subcomponents)'
                    call decomp_curladv(.True., .True.)
                else
                    write(*, '(A)') '  Decomposing curl of advection term (w/ flux twisting term)'
                    call decomp_curladv(.True., .False.)
                endif
            endif
        endif
        if (index(func_c, "e") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            write(*, '(A)') '  Calculating error term from offline calculation of nonlinear terms'
            call calc_errnlsub()
        endif
        if (index(func_c, "f") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            write(*, '(A)') '  Calculating curl of advection terms from momentum fluxes'
            call calc_curladv_flx()
        endif

        if (index(func_c, "a") /= 0 .or. index(func_c, "d") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            write(*, '(A)') '  Verifying offline calculations'
            call verify_nonl()
        endif
    endsubroutine

    subroutine load_const()
        use ncio, only : nc_read
        implicit none
        real(kind=kd_r) :: omega
        real(kind=kd_r) :: htn(B%nx, B%ny), hte(B%nx, B%ny), WORK(B%nx, B%ny, B%nz, 1)
        integer :: ix, iy

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(2A)') "Loading grid info from ", trim(fngrid%grid)
        write(*, '(2A)') "  and ", trim(fngrid%cons)
        write(*, '(2A)') "  and ", trim(fngrid%dz)


        allocate(tlat(B%nx, B%ny), tlong(B%nx, B%ny), z_t(B%nz), z_w(B%nz), dz(B%nz))
        allocate(ulat(B%nx, B%ny), ulong(B%nx, B%ny), dxu(B%nx, B%ny), dyu(B%nx, B%ny), tarea(B%nx, B%ny))
        allocate(dzt(B%nx, B%ny, B%nz), dzu(B%nx, B%ny, B%nz))
        allocate(tmask(B%nx, B%ny, B%nz), umask(B%nx, B%ny, B%nz))
        tmask = .false.
        umask = .false.

        call nc_read(fngrid%grid, 'TLAT' , WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        tlat  = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'TLONG', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        tlong = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'ULAT' , WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        ulat  = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'ULONG', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        ulong = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'TAREA', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        tarea = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'DXU'  , WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        dxu   = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'DYU'  , WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        dyu   = WORK(:, :, 1, 1)

        call nc_read(fngrid%grid, 'z_t'  , WORK, (/1/), (/B%nz/))
        z_t   = WORK(1:B%nz, 1, 1, 1)
        call nc_read(fngrid%grid, 'z_w'  , WORK, (/1/), (/B%nz/))
        z_w   = WORK(1:B%nz, 1, 1, 1)

        dz(1:B%nz-1) = z_w(2:) - z_w(1:B%nz-1)
        dz(B%nz) = dz(B%nz-1)

        call nc_read(fngrid%cons, 'omega', WORK, (/1/), (/1/))
        omega = WORK(1, 1, 1, 1)

        call nc_read(fngrid%cons, 'grav' , WORK, (/1/), (/1/))
        grav = WORK(1, 1, 1, 1)

        call nc_read(fngrid%dz, 'DZT', WORK, (/B%xl_reg, B%yd_reg, 1/), (/B%nx, B%ny, B%nz/))
        dzt = WORK(:, :, :, 1)
        call nc_read(fngrid%dz, 'DZU', WORK, (/B%xl_reg, B%yd_reg, 1/), (/B%nx, B%ny, B%nz/))
        dzu = WORK(:, :, :, 1)

        call nc_read(fngrid%dz, 'TMASK', WORK, (/B%xl_reg, B%yd_reg, 1/), (/B%nx, B%ny, B%nz/))
        where (abs(WORK(:, :, :, 1)) > 1e-10) tmask = .true.
        call nc_read(fngrid%dz, 'UMASK', WORK, (/B%xl_reg, B%yd_reg, 1/), (/B%nx, B%ny, B%nz/))
        where (abs(WORK(:, :, :, 1)) > 1e-10) umask = .true.

        allocate(fcor(B%nx, B%ny), fcort(B%nx, B%ny))
        do iy = 1, B%ny
            do ix = 1, B%nx
                fcor (ix, iy) = 2 * omega * sind(ulat(ix, iy))
                fcort(ix, iy) = 2 * omega * sind(tlat(ix, iy))
            enddo
        enddo

        ! do iy = 1, B%ny
        !     do ix = 1, B%nx
        !        where (abs(dzu(ix, iy, :)) < 1e-10) dzu(ix, iy, :) = dz
        !     enddo
        ! enddo

        write(*, *)
        write(*, *) "Calculating weight functions for decomposing nonlinear term"
        write(*, '(2A)') "  Loading grid info from ", trim(fngrid%grid)
        write(*, '(2A)') "  and ", trim(fngrid%dz)

        allocate(uarea(B%nx, B%ny), huw(B%nx, B%ny), hus(B%nx, B%ny))
        ! allocate(dxue(B%nx, B%ny), dyue(B%nx, B%ny), tareae(B%nx, B%ny), &
        !          dxun(B%nx, B%ny), dyun(B%nx, B%ny), tarean(B%nx, B%ny))

        ! call nc_read(fngrid%grid, 'HTN', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        ! htn = WORK(:, :, 1, 1)
        ! call nc_read(fngrid%grid, 'HTE', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        ! hte = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'HUW', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        huw = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'HUS', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        hus = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'UAREA', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        uarea = WORK(:, :, 1, 1)

        ! dyue(1:B%nx-1, :) = huw(2:B%nx, :)
        ! dxue(1:B%nx-1, :) = htn(2:B%nx, :)
        ! tareae = hus * hte

        ! dyun(:, 1:B%ny-1) = hte(:, 2:B%ny)
        ! dxun(:, 1:B%ny-1) = hus(:, 2:B%ny)
        ! tarean = htn * huw
    endsubroutine

    subroutine calc_velw()
        use popfun, only : t2u
        implicit none
        integer :: iz

        allocate(ue (B%nx, B%ny, B%nz), vn (B%nx, B%ny, B%nz), wt (B%nx, B%ny, B%nz), &
                 ume(B%nx, B%ny, B%nz), vme(B%nx, B%ny, B%nz), &
                 umn(B%nx, B%ny, B%nz), vmn(B%nx, B%ny, B%nz), &
                 umt(B%nx, B%ny, B%nz), vmt(B%nx, B%ny, B%nz))
        ue  = MVALUE; vn  = MVALUE; wt  = MVALUE
        ume = MVALUE; vme = MVALUE; umn = MVALUE; vmn = MVALUE; umt = MVALUE; vmt = MVALUE

        do iz = 1, B%nz
            ue(1:B%nx-1, 2:B%ny-1, iz) = &
              (0.25  * (uc(1:B%nx-1, 2:B%ny-1, iz) * dyu(1:B%nx-1, 2:B%ny-1) &
                                                   * dzu(1:B%nx-1, 2:B%ny-1, iz) + &
                        uc(2:B%nx  , 2:B%ny-1, iz) * dyu(2:B%nx  , 2:B%ny-1) &
                                                   * dzu(2:B%nx  , 2:B%ny-1, iz)) + &
               0.125 * (uc(2:B%nx  , 1:B%ny-2, iz) * dyu(2:B%nx  , 1:B%ny-2) &
                                                   * dzu(2:B%nx  , 1:B%ny-2, iz) + &
                        uc(1:B%nx-1, 1:B%ny-2, iz) * dyu(1:B%nx-1, 1:B%ny-2) &
                                                   * dzu(1:B%nx-1, 1:B%ny-2, iz) + &
                        uc(2:B%nx  , 3:B%ny  , iz) * dyu(2:B%nx  , 3:B%ny  ) &
                                                   * dzu(2:B%nx  , 3:B%ny  , iz) + &
                        uc(1:B%nx-1, 3:B%ny  , iz) * dyu(1:B%nx-1, 3:B%ny  ) &
                                                   * dzu(1:B%nx-1, 3:B%ny  , iz)))

            vn(2:B%nx-1, 1:B%ny-1, iz) = &
              (0.25  * (vc(2:B%nx-1, 1:B%ny-1, iz) * dxu(2:B%nx-1, 1:B%ny-1) &
                                                   * dzu(2:B%nx-1, 1:B%ny-1, iz) + &
                        vc(2:B%nx-1, 2:B%ny  , iz) * dxu(2:B%nx-1, 2:B%ny  ) &
                                                   * dzu(2:B%nx-1, 2:B%ny  , iz)) + &
               0.125 * (vc(1:B%nx-2, 2:B%ny  , iz) * dxu(1:B%nx-2, 2:B%ny  ) &
                                                   * dzu(1:B%nx-2, 2:B%ny  , iz) + &
                        vc(1:B%nx-2, 1:B%ny-1, iz) * dxu(1:B%nx-2, 1:B%ny-1) &
                                                   * dzu(1:B%nx-2, 1:B%ny-1, iz) + &
                        vc(3:B%nx  , 2:B%ny  , iz) * dxu(3:B%nx  , 2:B%ny  ) &
                                                   * dzu(3:B%nx  , 2:B%ny  , iz) + &
                        vc(3:B%nx  , 1:B%ny-1, iz) * dxu(3:B%nx  , 1:B%ny-1) &
                                                   * dzu(3:B%nx  , 1:B%ny-1, iz)))
            if (iz == 1) then
                wt(:, :, iz) = t2u(wc(:, :, 1), tarea, uarea)
            else
                wt(2:B%nx-1, 2:B%ny-1, iz) = wt(2:B%nx-1, 2:B%ny-1, iz-1) + &
                  ((vn(2:B%nx-1, 2:B%ny-1, iz-1) - vn(2:B%nx-1, 1:B%ny-2, iz-1)) + &
                   (ue(2:B%nx-1, 2:B%ny-1, iz-1) - ue(1:B%nx-2, 2:B%ny-1, iz-1))) &
                  / uarea(2:B%nx-1, 2:B%ny-1)
            endif
        enddo

        ume(1:B%nx-1, :, :) = (uc(1:B%nx-1, :, :) + uc(2:B%nx, :, :)) / 2.
        vme(1:B%nx-1, :, :) = (vc(1:B%nx-1, :, :) + vc(2:B%nx, :, :)) / 2.

        umn(:, 1:B%ny-1, :) = (uc(:, 1:B%ny-1, :) + uc(:, 2:B%ny, :)) / 2.
        vmn(:, 1:B%ny-1, :) = (vc(:, 1:B%ny-1, :) + vc(:, 2:B%ny, :)) / 2.

        umt(:, :, 1) = uc(:, :, 1)
        umt(:, :, 2:B%nz) = (uc(:, :, 1:B%nz-1) + uc(:, :, 2:B%nz)) / 2.
        vmt(:, :, 1) = vc(:, :, 1)
        vmt(:, :, 2:B%nz) = (vc(:, :, 1:B%nz-1) + vc(:, :, 2:B%nz)) / 2.

        ! write(*, fmtm_vel) 'ue: ', ue(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vel) 'vn: ', vn(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vel) 'wt: ', wt(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        write(*, fmtm_vel) 'ume: ', ume(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vel) 'vme: ', vme(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vel) 'umn: ', umn(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vel) 'vmn: ', vmn(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vel) 'umt: ', umt(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vel) 'vmt: ', vmt(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
    endsubroutine

    subroutine zeta_equation()
        use popfun, only : zcurl, u2t
        use derives, only : dd_xsw_chain, dd_ysw_chain
        implicit none
        real(kind=kd_r), dimension(:, :, :), allocatable :: curlcor, corx, cory
        real(kind=kd_r), dimension(B%nx, B%ny) :: pgradsfx, pgradsfy
        real(kind=kd_r), dimension(B%nx, B%ny) :: stretchpx, betavx, stretchpy, betavy, &
                                                  errcorx, errcory, ONES
        integer :: iz

        !!----------------------------------------------------------------------
        ! Coriolis terms in momentum equqation
        write(*, *)
        write(*, '(A)'), '    Calculating Coriolis terms in momentum equations'
        ! write(*, *) "f", fcor(B%xi_dp, B%yi_dp), fcort(B%xi_dp, B%yi_dp)

        allocate(corx(B%nx, B%ny, B%nz), cory(B%nx, B%ny, B%nz))
        do iz = 1, B%nz
            corx(:, :, iz) = - fcor * vc(:, :, iz)
            cory(:, :, iz) =   fcor * uc(:, :, iz)
        enddo
        write(*, fmtm_vor) 'Corx: ', corx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'Cory: ', cory(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        !!----------------------------------------------------------------------
        ! Barotropic and baroclinic pressure gradient
        write(*, *)
        write(*, '(A)'), '    Calculating and combining barotropic and baroclinic pressure gradient'

        ! Note that the last row and last column of the pressure terms (at u cells)
        !   are set to be zero, as no SSH is available there. Therefore the values
        !   for the values of curlpgrad in the last row and last colunm are invalid.
        pgradsfx = 0.
        pgradsfy = 0.
        pgradsfx(1:B%nx-1, 1:B%ny-1) = ((ssh(2:B%nx,   1:B%ny-1, 1) + ssh(2:B%nx,   2:B%ny  , 1))/2 - &
                                        (ssh(1:B%nx-1, 1:B%ny-1, 1) + ssh(1:B%nx-1, 2:B%ny  , 1))/2) / dxu(1:B%nx-1, 1:B%ny-1) * (-grav)
        pgradsfy(1:B%nx-1, 1:B%ny-1) = ((ssh(1:B%nx-1, 2:B%ny  , 1) + ssh(2:B%nx,   2:B%ny  , 1))/2 - &
                                        (ssh(1:B%nx-1, 1:B%ny-1, 1) + ssh(2:B%nx,   1:B%ny-1, 1))/2) / dyu(1:B%nx-1, 1:B%ny-1) * (-grav)

        do iz = 1, B%nz
            gradx(:, :, iz) = pgradsfx - gradx(:, :, iz)
            grady(:, :, iz) = pgradsfy - grady(:, :, iz)
        enddo

        where(umask) gradx = 0.
        where(umask) grady = 0.

        write(*, fmtm_vor) 'pgradx: ', gradx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'pgrady: ', grady(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        !!----------------------------------------------------------------------
        ! Curl of nonlinear, cor, pgrad, h/v diffusion and res terms
        write(*, *)
        write(*, '(A)') '    Calculating curl of nonlinear, Coriolis, pgrad, horizontal/vertical diffusion and residual terms'
        allocate(curlcor(B%nx, B%ny, B%nz))
        curlcor = 0.
        do iz = 1, B%nz
            curlnonl (:, :, iz) = zcurl(-advx (:, :, iz), -advy (:, :, iz), &
                                        dxu*dzu(:, :, iz), dyu*dzu(:, :, iz), tarea*dzt(:, :, iz))
            curlcor  (:, :, iz) = zcurl(-corx (:, :, iz), -cory (:, :, iz), &
                                        dxu*dzu(:, :, iz), dyu*dzu(:, :, iz), tarea*dzt(:, :, iz))
            curlpgrad(:, :, iz) = zcurl(gradx (:, :, iz), grady (:, :, iz), &
                                        dxu*dzu(:, :, iz), dyu*dzu(:, :, iz), tarea*dzt(:, :, iz))
            curlhdiff(:, :, iz) = zcurl(hdiffx(:, :, iz), hdiffy(:, :, iz), &
                                        dxu*dzu(:, :, iz), dyu*dzu(:, :, iz), tarea*dzt(:, :, iz))
            curlvdiff(:, :, iz) = zcurl(vdiffx(:, :, iz), vdiffy(:, :, iz), &
                                        dxu*dzu(:, :, iz), dyu*dzu(:, :, iz), tarea*dzt(:, :, iz))
        enddo
        deallocate(corx, cory)
        res = curlpgrad + curlhdiff + curlvdiff + curlnonl + curlcor

        where(tmask)
            curlnonl  = MVALUE
            curlcor   = MVALUE
            curlpgrad = MVALUE
            curlhdiff = MVALUE
            curlvdiff = MVALUE
            res       = MVALUE
        endwhere

        write(*, fmtm_vor) 'curlnonl: ' , curlnonl (B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'curlcor: '  , curlcor  (B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'curlpgrad: ', curlpgrad(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'curlhdiff: ', curlhdiff(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'curlvdiff: ', curlvdiff(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'res: '      , res      (B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        !!----------------------------------------------------------------------
        ! Decompositing of curl of Coriolis term
        write(*, *)
        write(*, '(A)') '    Decomposing curl of Coriolis term'

        ONES = 1.
        do iz = 1, B%nz
            call dd_xsw_chain(fcor, fcort, &
                              uc(:,:,iz) * dyu * dzu(:,:,iz), u2t(uc(:,:,iz), dyu*dzu(:,:,iz), ONES), &
                              tarea * dzt(:,:,iz), &
                              stretchpx, betavx, errcorx)
            call dd_ysw_chain(fcor, fcort, &
                              vc(:,:,iz) * dxu * dzu(:,:,iz), u2t(vc(:,:,iz), dxu*dzu(:,:,iz), ONES), &
                              tarea * dzt(:,:,iz), &
                              stretchpy, betavy, errcory)

            betav   (:, :, iz) = -(betavx + betavy)
            stretchp(:, :, iz) = -(stretchpx + stretchpy)
            err_cor (: ,:, iz) = -(errcorx + errcory)
        enddo

        where(tmask)
            betav    = MVALUE
            stretchp = MVALUE
            err_cor  = MVALUE
        endwhere

        !!----------------------------------------------------------------------
        ! Verifying decomposition and clean up
        write(*, *)
        write(*, '(A)') '    Verify decomposition of curl(-fu, fv) ...'
        do iz = B%zi_dpst, B%zi_dped
            write(*, '(A, I02)') '      iz = ',  iz
            write(*, fmts_vor) 'Curlcor: '    , curlcor  (B%xi_dp, B%yi_dp, iz)
            write(*, fmts_vor) 'bv + -fdwdz: ', betav(B%xi_dp, B%yi_dp, iz) + stretchp(B%xi_dp, B%yi_dp, iz)
            write(*, fmts_vor) 'bv: '         , betav(B%xi_dp, B%yi_dp, iz)
            write(*, fmts_vor) '-fdwdz: '     , stretchp(B%xi_dp, B%yi_dp, iz)
            write(*, fmts_vor) 'errcor: '     , err_cor(B%xi_dp, B%yi_dp, iz)
            write(*, fmts_vor) 'Diff: '       , curlcor (B%xi_dp, B%yi_dp, iz) - &
                                                betav   (B%xi_dp, B%yi_dp, iz) - &
                                                stretchp(B%xi_dp, B%yi_dp, iz) - &
                                                err_cor (B%xi_dp, B%yi_dp, iz)
            write(*, *)
        enddo
        deallocate(curlcor)
    endsubroutine

    subroutine calc_errnlsub()
        if (allocated(curladv)) then
            err_nlsub = curlnonl - curladv - curlmet
        elseif (allocated(advu)) then
            err_nlsub = curlnonl - advu - advv - advw - advVx - advVy - advVz - err_nldecomp - curlmet
        endif
    endsubroutine

    subroutine calc_curlmet()
        use popfun, only : zcurl
        implicit none
        real(kind=kd_r), dimension(B%nx, B%ny) :: kxu, kyu, metx, mety
        integer :: iz

        write(*, *)
        write(*, '(A)') '  ---------------------------------------------------'
        write(*, '(A)') '  Calculating curl of metric terms (offline)'
        kxu = MVALUE
        kyu = MVALUE

        kxu(2:B%nx, :) = (huw(2:B%nx, :) - huw(1:B%nx-1, :)) / uarea(2:B%nx, :)
        kyu(:, 2:B%ny) = (hus(:, 2:B%ny) - hus(:, 1:B%ny-1)) / uarea(:, 2:B%ny)

        do iz = 1, B%nz
            metx = uc(:, :, iz) * vc(:, :, iz) * kyu - vc(:, :, iz) * vc(:, :, iz) * kxu
            mety = uc(:, :, iz) * vc(:, :, iz) * kxu - uc(:, :, iz) * uc(:, :, iz) * kyu

            curlmet(:, :, iz) = zcurl(-metx, -mety, dxu*dzu(:, :, iz), dyu*dzu(:, :, iz), tarea*dzt(:, :, iz))
            curlmet(1, :, iz) = MVALUE
            curlmet(:, 1, iz) = MVALUE
        enddo
        where(tmask) curlmet = MVALUE
        ! write(*, fmtm_vor) 'curlmet', curlmet(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
    endsubroutine

    subroutine calc_curladv(debug)
        use popfun, only : zcurl
        implicit none
        logical, intent(in) :: debug
        ! real(kind=kd_r), dimension(B%nx, B%ny) :: advx, advy
        real(kind=kd_r), dimension(B%nx, B%ny) :: advxu, advxv, advxw, advyu, advyv, advyw
        ! real(kind=kd_r), dimension(B%nx, B%ny) :: WORKx, WORKy, WORKz
        integer :: iz

        ! Calculating derived velocity at walls
        write(*, *)
        write(*, '(A)') '    -------------------------------------------------'
        write(*, '(A)') '    Calculating velocity at the walls'
        call calc_velw()

        do iz = 1, B%nz
            advxu = 0.; advxv = 0.; advxw = 0. 
            advyu = 0.; advyv = 0.; advyw = 0. 
            advxu(2:B%nx, 2:B%ny) = (ue(2:B%nx  , 2:B%ny  , iz) * ume (2:B%nx  , 2:B%ny  , iz) - &
                                     ue(1:B%nx-1, 2:B%ny  , iz) * ume (1:B%nx-1, 2:B%ny  , iz))  &
                                    / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)
            advxv(2:B%nx, 2:B%ny) = (vn(2:B%nx  , 2:B%ny  , iz) * umn (2:B%nx  , 2:B%ny  , iz) - &
                                     vn(2:B%nx  , 1:B%ny-1, iz) * umn (2:B%nx  , 1:B%ny-1, iz))  &
                                    / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)
            advxw(2:B%nx, 2:B%ny) = wt(2:B%nx, 2:B%ny, iz) * umt(2:B%nx, 2:B%ny, iz) / dzu(2:B%nx, 2:B%ny, iz)

            advyu(2:B%nx, 2:B%ny) = (ue(2:B%nx  , 2:B%ny  , iz) * vme (2:B%nx  , 2:B%ny  , iz) - &
                                     ue(1:B%nx-1, 2:B%ny  , iz) * vme (1:B%nx-1, 2:B%ny  , iz))  &
                                    / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)
            advyv(2:B%nx, 2:B%ny) = (vn(2:B%nx  , 2:B%ny  , iz) * vmn (2:B%nx  , 2:B%ny  , iz) - &
                                     vn(2:B%nx  , 1:B%ny-1, iz) * vmn (2:B%nx  , 1:B%ny-1, iz))  &
                                    / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)
            advyw(2:B%nx, 2:B%ny) = wt(2:B%nx, 2:B%ny, iz) * vmt(2:B%nx, 2:B%ny, iz) / dzu(2:B%nx, 2:B%ny, iz)

            if (iz < B%nz) then
                advxw(2:B%nx, 2:B%ny) = advxw(2:B%nx, 2:B%ny) &
                    - wt(2:B%nx, 2:B%ny, iz+1) * umt(2:B%nx, 2:B%ny, iz+1) / dzu(2:B%nx, 2:B%ny, iz)
                advyw(2:B%nx, 2:B%ny) = advyw(2:B%nx, 2:B%ny) &
                    - wt(2:B%nx, 2:B%ny, iz+1) * vmt(2:B%nx, 2:B%ny, iz+1) / dzu(2:B%nx, 2:B%ny, iz)
            endif

            where(umask(:,:,iz)) ! For cases where dzu = 0.
                advxu = 0.; advxv = 0.; advxw = 0. 
                advyu = 0.; advyv = 0.; advyw = 0. 
            endwhere

            curladv(:, :, iz) = zcurl(-advxu-advxv-advxw, -advyu-advyu-advyv, &
            dxu*dzu(:,:,iz), dyu*dzu(:,:,iz), tarea*dzt(:,:,iz))
            if (debug) then 
                curladvu(:, :, iz) = zcurl(-advxu, -advyu, dxu*dzu(:,:,iz), dyu*dzu(:,:,iz), tarea*dzt(:,:,iz))
                curladvv(:, :, iz) = zcurl(-advxv, -advyv, dxu*dzu(:,:,iz), dyu*dzu(:,:,iz), tarea*dzt(:,:,iz))
                curladvw(:, :, iz) = zcurl(-advxw, -advyw, dxu*dzu(:,:,iz), dyu*dzu(:,:,iz), tarea*dzt(:,:,iz))
            endif
        enddo

        where(tmask) curladv = MVALUE
        curladv(1:2, :, :) = MVALUE  ! advxu/v, advyu/v, advwu/v at ix = 1 and iy = 1 is NaN, 
        curladv(:, 1:2, :) = MVALUE  !   which are needed for curl at ix = 2 and iy = 2
        curladv(B%nx, :, :) = MVALUE ! ue, vn at ix = nx and iy = ny is NaN
        curladv(:, B%ny, :) = MVALUE
        if (debug) then
            where(tmask)
               curladvu = MVALUE; curladvv = MVALUE; curladvw = MVALUE
            endwhere
            curladvu(1:2, :, :) = MVALUE; curladvv(1:2, :, :) = MVALUE; curladvw(1:2, :, :) = MVALUE
            curladvu(:, 1:2, :) = MVALUE; curladvv(:, 1:2, :) = MVALUE; curladvw(:, 1:2, :) = MVALUE
            curladvu(B%nx, :, :) = MVALUE; curladvv(B%nx, :, :) = MVALUE; curladvw(B%nx, :, :) = MVALUE
            curladvu(:, B%ny, :) = MVALUE; curladvv(:, B%ny, :) = MVALUE; curladvw(:, B%ny, :) = MVALUE
        endif
        deallocate(ue, vn, wt, ume, vme, umn, vmn, umt, vmt)
        ! write(*, fmtm_vor), 'curladv', curladv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
    endsubroutine

    subroutine calc_curladv_flx()
        use popfun, only : zcurl
        implicit none
        real(kind=kd_r), dimension(B%nx, B%ny) :: advx, advy
        real(kind=kd_r), dimension(B%nx, B%ny) :: ueu_c, uev_c, vnu_c, vnv_c, wtu_c, wtv_c, wtu_c2, wtv_c2
        integer :: iz

        do iz = 1, B%nz
            advx = 0.
            advy = 0.

            ueu_c = ueu(:, :, iz) * dzu(:, :, iz) * uarea
            uev_c = uev(:, :, iz) * dzu(:, :, iz) * uarea
            vnu_c = vnu(:, :, iz) * dzu(:, :, iz) * uarea
            vnv_c = vnv(:, :, iz) * dzu(:, :, iz) * uarea
            wtu_c = wtu(:, :, iz) * dz(iz)
            wtv_c = wtv(:, :, iz) * dz(iz)

            advx(2:B%nx, 2:B%ny) = (ueu_c(2:B%nx  , 2:B%ny  ) - ueu_c(1:B%nx-1, 2:B%ny  ) + &
                                    vnu_c(2:B%nx  , 2:B%ny  ) - vnu_c(2:B%nx  , 1:B%ny-1))  &
                                    / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)       &
                                    + wtu_c(2:B%nx, 2:B%ny) / dzu(2:B%nx, 2:B%ny, iz)

            advy(2:B%nx, 2:B%ny) = (uev_c(2:B%nx  , 2:B%ny  ) - uev_c(1:B%nx-1, 2:B%ny  ) + &
                                    vnv_c(2:B%nx  , 2:B%ny  ) - vnv_c(2:B%nx  , 1:B%ny-1))  &
                                    / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)       &
                                   + wtv_c(2:B%nx, 2:B%ny) / dzu(2:B%nx, 2:B%ny, iz)

            if (iz < B%nz) then
                wtu_c2 = wtu(:, :, iz+1) * dz(iz+1)
                wtv_c2 = wtv(:, :, iz+1) * dz(iz+1)
                advx(2:B%nx, 2:B%ny) = advx(2:B%nx, 2:B%ny) &
                                       - wtu_c2(2:B%nx, 2:B%ny) / dzu(2:B%nx, 2:B%ny, iz)
                advy(2:B%nx, 2:B%ny) = advy(2:B%nx, 2:B%ny) &
                                       - wtv_c2(2:B%nx, 2:B%ny) / dzu(2:B%nx, 2:B%ny, iz)
            endif

            advx(1, :) = MVALUE
            advx(:, 1) = MVALUE
            advy(1, :) = MVALUE
            advy(:, 1) = MVALUE

            where(umask(:,:,iz)) ! For cases where ucell = 0.
                advx = 0.
                advy = 0.
            endwhere

            curladvf(:, :, iz) = zcurl(-advx, -advy, dxu*dzu(:,:,iz), dyu*dzu(:,:,iz), tarea*dzt(:,:,iz))
            curladvf(1, :, iz) = MVALUE
            curladvf(:, 1, iz) = MVALUE
        enddo
    endsubroutine

    subroutine decomp_curladv(twif, debug)
        use derives
        use popfun, only : u2t
        implicit none
        logical, intent(in) :: twif, debug
        real(kind=kd_r), dimension(B%nx, B%ny) :: ONES, WORK, u1, u2, test1, test2, test3
        real(kind=kd_r), dimension(B%nx, B%ny) :: F_uudxdy, F_vvdxdy, F_wwdxdy, F_wvdxdz, F_wudydz, wm
        real(kind=kd_r), dimension(B%nx, B%ny, 2) :: u10du2, u20du1, u10du2_zx, u20du1_zx, u10du2_zy, u20du1_zy
        real(kind=kd_r), dimension(B%nx, B%ny, 2) :: F_wdzx, F_wdzy
        integer :: iz

        write(*, *)
        write(*, '(A)') '    -------------------------------------------------'
        write(*, '(A)') '    Calculating velocity at the walls'
        call calc_velw()

        !!----------------------------------------------------------------------
        ! Decomposition
        write(*, *)
        write(*, '(A)') '    Calculating decomposition'
        ONES = 1.

        do iz = 1, B%nz
          !! Chain rule: advection terms
          ! advu
            ! d [d(uv)/dx] / dx = d (udv/dx) / dx + d (vdu/dx) / dx
            u1 = ue(:,:,iz)
            u2 = vme(:,:,iz)/dxu
            where (umask(:,:,iz))
                u1 = 0.; u2 = 0.
            endwhere
            call dd_xw_chain(u1, mean_xw(u1), u2, mean_xw(u2), ONES, u10du2(:,:,2), u20du1(:,:,2))
            u1 = shift_xe(ue(:,:,iz))
            u2 = shift_xe(vme(:,:,iz))/dxu
            where (umask(:,:,iz))
                u1 = 0.; u2 = 0.
            endwhere
            call dd_xw_chain(u1, mean_xw(u1), u2, mean_xw(u2), ONES, u10du2(:,:,1), u20du1(:,:,1))

            WORK = mean_ys(u10du2(:,:,2) - u10du2(:,:,1)) / tarea / dzt(:,:,iz)
            advu(:, :, iz) = advu(:, :, iz) + WORK
            if (debug) advu_x(:, :, iz) = WORK

            if (twif) then
                WORK = mean_ys(u20du1(:,:,2) - u20du1(:,:,1)) / tarea / dzt(:,:,iz)
                advVx(:, :, iz) = advVx(:, :, iz) + WORK
                if (debug) advVx_x(:, :, iz) = WORK
            else
                WORK = mean_ys(mean_xw(dd_xw(ue(:,:,iz), ONES)) * mean_xw(dd_xw(vme(:,:,iz), dxu))) / tarea / dzt(:,:,iz)
                advVx(:, :, iz) = advVx(:, :, iz) + WORK
                if (debug) advVx_x(:, :, iz) = WORK
                WORK = mean_ys(mean_xw(mean_xw(vme(:,:,iz))/dxu) * dd_xw(dd_xw(ue(:,:,iz), ONES), ONES)) / tarea / dzt(:,:,iz)
                err_nldecomp(:,:,iz) = err_nldecomp(:,:,iz) + WORK
                if (debug) errnl_ux(:, :, iz) = WORK                  
            endif

            ! d [d(uu)/dy] / dx = d (udu/dy) / dx + d (udu/dy) / dx
            u1 = ue(:,:,iz)
            u2 = ume(:,:,iz)/dyu
            where (umask(:,:,iz))
                u1 = 0.; u2 = 0.
            endwhere            
            call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2(:,:,2), u20du1(:,:,2))
            u1 = shift_xe(ue(:,:,iz))
            u2 = shift_xe(ume(:,:,iz))/dyu
            where (umask(:,:,iz))
                u1 = 0.; u2 = 0.
            endwhere
            call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2(:,:,1), u20du1(:,:,1))
            WORK = mean_xw(u10du2(:,:,2) - u10du2(:,:,1)) / tarea / dzt(:,:,iz)
            advu(:, :, iz) = advu(:, :, iz) - WORK
            if (debug) advu_y(:, :, iz) = -WORK

            if (twif) then
                WORK = mean_xw(u20du1(:,:,2) - u20du1(:,:,1)) / tarea / dzt(:,:,iz)
                advVx(:, :, iz) = advVx(:, :, iz) - WORK
                if (debug) advVx_y(:, :, iz) = -WORK
            else
                WORK = mean_xw(mean_xw(dd_ys(ue(:,:,iz), ONES)) * mean_ys(dd_xw(ume(:,:,iz), dyu))) / tarea / dzt(:,:,iz)
                advVx(:, :, iz) = advVx(:, :, iz) - WORK
                if (debug) advVx_y(:, :, iz) = -WORK
                WORK = mean_xw(mean_ys(mean_xw(ume(:,:,iz))/dyu) * dd_xw(dd_ys(ue(:,:,iz), ONES), ONES)) / tarea / dzt(:,:,iz)
                err_nldecomp(:,:,iz) = err_nldecomp(:,:,iz) - WORK
                if (debug) errnl_uy(:, :, iz) = -WORK
            endif

          ! advv
            ! d [d(vv)/dx] / dy = d (vdvdx) / dy + d (vdvdx) / dy
            u1 = vn(:,:,iz)
            u2 = vmn(:,:,iz)/dxu
            where (umask(:,:,iz))
                u1 = 0.; u2 = 0.
            endwhere
            call dd_xw_chain(u1, mean_xw(u1), u2, mean_xw(u2), ONES, u10du2(:,:,2), u20du1(:,:,2))
            u1 = shift_yn(vn(:,:,iz))
            u2 = shift_yn(vmn(:,:,iz))/dxu
            where (umask(:,:,iz))
                u1 = 0.; u2 = 0.
            endwhere
            call dd_xw_chain(u1, mean_xw(u1), u2, mean_xw(u2), ONES, u10du2(:,:,1), u20du1(:,:,1))
            WORK = mean_ys(u10du2(:,:,2) - u10du2(:,:,1)) / tarea / dzt(:,:,iz)
            advv(:, :, iz) = advv(:, :, iz) + WORK
            if (debug) advv_x(:, :, iz) = WORK

            if (twif) then
                WORK = mean_ys(u20du1(:,:,2) - u20du1(:,:,1)) / tarea / dzt(:,:,iz)
                advVy(:, :, iz) = advVy(:, :, iz) + WORK
                if (debug) advVy_x(:, :, iz) = WORK
            else
                WORK = mean_ys(mean_ys(dd_xw(vn(:,:,iz), ONES)) * mean_xw(dd_ys(vmn(:,:,iz), dxu))) / tarea / dzt(:,:,iz)
                advVy(:, :, iz) = advVy(:, :, iz) + WORK
                if (debug) advVy_x(:, :, iz) = WORK
                WORK = mean_ys(mean_xw(mean_ys(vmn(:,:,iz))/dxu) * dd_ys(dd_xw(vn(:,:,iz), ONES), ONES)) / tarea / dzt(:,:,iz)
                err_nldecomp(:,:,iz) = err_nldecomp(:,:,iz) + WORK
                if (debug) errnl_vx(:, :, iz) = WORK
            endif

            ! d [d(vu)/dy] / dy = d (vdu/dy) / dy + d (udv/dy) / dy
            u1 = vn(:,:,iz)
            u2 = umn(:,:,iz)/dyu
            call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2(:,:,2), u20du1(:,:,2))
            u1 = shift_yn(vn(:,:,iz))
            u2 = shift_yn(umn(:,:,iz))/dyu
            call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2(:,:,1), u20du1(:,:,1))
            WORK = mean_xw(u10du2(:,:,2) - u10du2(:,:,1)) / tarea / dzt(:,:,iz)
            advv(:, :, iz) = advv(:, :, iz) - WORK
            if (debug) advv_y(:, :, iz) = -WORK

            if (twif) then
                WORK = mean_xw(u20du1(:,:,2) - u20du1(:,:,1)) / tarea / dzt(:,:,iz)
                advVy(:, :, iz) = advVy(:, :, iz) - WORK
                if (debug) advVy_y(:, :, iz) = -WORK
            else
                WORK = mean_xw(mean_ys(dd_ys(vn(:,:,iz), ONES)) * mean_ys(dd_ys(umn(:,:,iz), dyu))) / tarea / dzt(:,:,iz)
                advVy(:, :, iz) = advVy(:, :, iz) - WORK
                if (debug) advVy_y(:, :, iz) = -WORK
                WORK = mean_xw(mean_ys(mean_ys(umn(:,:,iz))/dyu) * dd_ys(dd_ys(vn(:,:,iz), ONES), ONES)) / tarea / dzt(:,:,iz)
                err_nldecomp(:,:,iz) = err_nldecomp(:,:,iz) - WORK
                if (debug) errnl_vy(:, :, iz) = -WORK
            endif

          ! advw
            ! d(wv) / dx  &  d(wu) / dy
            if (iz == 1) then
                u1 = wt(:,:,iz)
                u2 = vmt(:,:,iz)*dyu
                where (umask(:,:,iz)) u1 = 0.
                where (umask(:,:,iz)) u2 = 0.
                call dd_xw_chain(u1, mean_xw(u1), u2, mean_xw(u2), ONES, u10du2_zx(:,:,2), u20du1_zx(:,:,2))
                u2 = umt(:,:,iz)*dxu
                where (umask(:,:,iz)) u2 = 0.
                call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2_zy(:,:,2), u20du1_zy(:,:,2))
            endif

            if (iz == B%nz) then
                u10du2_zx(:,:,1) = 0.
                u20du1_zx(:,:,1) = 0.
                u10du2_zy(:,:,1) = 0.
                u20du1_zy(:,:,1) = 0.
            else
                u1 = wt(:,:,iz+1)
                u2 = vmt(:,:,iz+1)*dyu
                where (umask(:,:,iz)) u1 = 0.
                where (umask(:,:,iz)) u2 = 0.
                call dd_xw_chain(u1, mean_xw(u1), u2, mean_xw(u2), ONES, u10du2_zx(:,:,1), u20du1_zx(:,:,1))
                u2 = umt(:,:,iz+1)*dxu
                where (umask(:,:,iz)) u2 = 0.
                call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2_zy(:,:,1), u20du1_zy(:,:,1))
            endif

            WORK = mean_ys(u10du2_zx(:,:,2) - u10du2_zx(:,:,1)) / tarea / dzt(:, :, iz)
            advw(:, :, iz) = advw(:, :, iz) + WORK
            if (debug) advw_x(:, :, iz) = WORK

            WORK = mean_xw(u10du2_zy(:,:,2) - u10du2_zy(:,:,1)) / tarea / dzt(:, :, iz)
            advw(:, :, iz) = advw(:, :, iz) - WORK
            if (debug) advw_y(:, :, iz) = -WORK

            if (twif) then
                WORK = mean_ys(u20du1_zx(:,:,2) - u20du1_zx(:,:,1)) / tarea / dzt(:, :, iz)
                advVz(:, :, iz) = advVz(:, :, iz) + WORK
                if (debug) advVz_x(:, :, iz) = WORK

                WORK = mean_xw(u20du1_zy(:,:,2) - u20du1_zy(:,:,1)) / tarea / dzt(:, :, iz)
                advVz(:, :, iz) = advVz(:, :, iz) - WORK
                if (debug) advVz_y(:, :, iz) = -WORK
            else
                if (iz == B%nz) then
                    WORK = mean_ys(mean_xw((vmt(:,:,iz) - 0.) * dyu) * dd_xw((wt(:,:,iz) + 0.)/2, ONES)) / tarea / dzt(:,:,iz) 
                    advVz(:, :, iz) = advVz(:, :, iz) + WORK
                    if (debug) advVz_x(:, :, iz) = WORK

                    WORK = mean_xw(mean_ys((umt(:,:,iz) - 0.) * dxu) * dd_ys((wt(:,:,iz) + 0.)/2, ONES)) / tarea / dzt(:,:,iz) 
                    advVz(:, :, iz) = advVz(:, :, iz) - WORK
                    if (debug) advVz_y(:, :, iz) = -WORK

                    WORK = mean_ys(dd_xw(wt(:,:,iz) - 0., ONES) * mean_xw((vmt(:,:,iz) + 0.)/2 * dyu)) / tarea / dzt(:,:,iz)
                    err_nldecomp(:,:,iz) = err_nldecomp(:,:,iz) + WORK
                    if (debug) errnl_wx(:, :, iz) = WORK

                    WORK = mean_xw(dd_ys(wt(:,:,iz) - 0., ONES) * mean_ys((umt(:,:,iz) + 0.)/2 * dxu)) / tarea / dzt(:,:,iz)
                    err_nldecomp(:,:,iz) = err_nldecomp(:,:,iz) - WORK
                    if (debug) errnl_wy(:, :, iz) = -WORK
               else
                    WORK = mean_ys(mean_xw((vmt(:,:,iz) - vmt(:,:,iz+1)) * dyu) * dd_xw((wt(:,:,iz) + wt(:,:,iz+1))/2, ONES)) & 
                           / tarea / dzt(:,:,iz) 
                    advVz(:, :, iz) = advVz(:, :, iz) + WORK
                    if (debug) advVz_x(:, :, iz) = WORK

                    WORK = mean_xw(mean_ys((umt(:,:,iz) - umt(:,:,iz+1)) * dxu) * dd_ys((wt(:,:,iz) + wt(:,:,iz+1))/2, ONES)) & 
                           / tarea / dzt(:,:,iz) 
                    advVz(:, :, iz) = advVz(:, :, iz) - WORK
                    if (debug) advVz_y(:, :, iz) = -WORK

                    WORK = mean_ys(dd_xw(wt(:,:,iz) - wt(:,:,iz+1), ONES) * mean_xw((vmt(:,:,iz) + vmt(:,:,iz+1))/2 * dyu)) &
                           / tarea / dzt(:,:,iz)
                    err_nldecomp(:,:,iz) = err_nldecomp(:,:,iz) + WORK
                    if (debug) errnl_wx(:, :, iz) = WORK

                    WORK = mean_xw(dd_ys(wt(:,:,iz) - wt(:,:,iz+1), ONES) * mean_ys((umt(:,:,iz) + umt(:,:,iz+1))/2 * dxu)) &
                           / tarea / dzt(:,:,iz)
                    err_nldecomp(:,:,iz) = err_nldecomp(:,:,iz) - WORK
                    if (debug) errnl_wy(:, :, iz) = -WORK
                endif
            endif
           
            u10du2_zx(:,:,2) = u10du2_zx(:,:,1)
            u10du2_zy(:,:,2) = u10du2_zy(:,:,1)
            u20du1_zx(:,:,2) = u20du1_zx(:,:,1)
            u20du1_zy(:,:,2) = u20du1_zy(:,:,1)

            ! Add five terms to compose flux of 3D voricity twisting
            if (twif) then 
                F_uudxdy = dd_ys(dd_xw(0.5 * uc(:, :, iz) * uc(:, :, iz) * dzu(:, :, iz), ONES), &
                            tarea * dzt(:,:,iz))
                F_vvdxdy = dd_ys(dd_xw(0.5 * vc(:, :, iz) * vc(:, :, iz) * dzu(:, :, iz), ONES), &
                            tarea * dzt(:,:,iz))

                if (iz == B%nz) then
                    wm = (wt(:, :, iz) + 0.) / 2
                else
                    wm = (wt(:, :, iz) + wt(:, :, iz+1)) / 2
                endif
                F_wwdxdy = dd_ys(dd_xw(0.5 * wm * wm * dzu(:, :, iz), ONES), tarea * dzt(:,:,iz))

                if (iz == 1) then
                    F_wdzx(:,:,2) = dd_xsw(wt(:, :, iz) * vmt(:, :, iz) * dyu, tarea)
                    F_wdzy(:,:,2) = dd_ysw(wt(:, :, iz) * umt(:, :, iz) * dxu, tarea)
                endif

                if (iz == B%nz) then
                    F_wdzx(:,:,1) = 0.
                    F_wdzy(:,:,1) = 0.
                else
                    F_wdzx(:,:,1) = dd_xsw(wt(:, :, iz+1) * vmt(:, :, iz+1) * dyu, ONES)
                    F_wdzy(:,:,1) = dd_ysw(wt(:, :, iz+1) * umt(:, :, iz+1) * dxu, ONES)
                endif
                F_wvdxdz = (F_wdzx(:,:,2) - F_wdzx(:,:,1)) / tarea / dzt(:, :, iz)
                F_wudydz = (F_wdzy(:,:,2) - F_wdzy(:,:,1)) / tarea / dzt(:, :, iz)

                advVx(:, :, iz) = advVx(:, :, iz) + F_wvdxdz + F_vvdxdy - F_wwdxdy + F_uudxdy
                advVy(:, :, iz) = advVy(:, :, iz) - F_wudydz - F_vvdxdy + F_wwdxdy - F_uudxdy
                advVz(:, :, iz) = advVz(:, :, iz) + F_wudydz - F_wvdxdz

                ! stepping down in z direction
                F_wdzx(:,:,2) = F_wdzx(:,:,1)
                F_wdzy(:,:,2) = F_wdzy(:,:,1)
            endif
        enddo

        ! RHS
        advu = -advu; advv = -advv; advw = -advw;
        advVx = -advVx; advVy = -advVy; advVz = -advVz

        where(tmask)
            advu  = MVALUE; advv  = MVALUE; advw  = MVALUE
            advVx = MVALUE; advVy = MVALUE; advVz = MVALUE
        endwhere
        advu (1:2, :, :) = MVALUE; advv (1:2, :, :) = MVALUE; advw (1:2, :, :) = MVALUE
        advVx(1:2, :, :) = MVALUE; advVy(1:2, :, :) = MVALUE; advVz(1:2, :, :) = MVALUE
        advu (:, 1:2, :) = MVALUE; advv (:, 1:2, :) = MVALUE; advw (:, 1:2, :) = MVALUE
        advVx(:, 1:2, :) = MVALUE; advVy(:, 1:2, :) = MVALUE; advVz(:, 1:2, :) = MVALUE
        advu (B%nx, :, :) = MVALUE; advv (B%nx, :, :) = MVALUE; advw (B%nx, :, :) = MVALUE
        advVx(B%nx, :, :) = MVALUE; advVy(B%nx, :, :) = MVALUE; advVz(B%nx, :, :) = MVALUE
        advu (:, B%ny, :) = MVALUE; advv (:, B%ny, :) = MVALUE; advw (:, B%ny, :) = MVALUE
        advVx(:, B%ny, :) = MVALUE; advVy(:, B%ny, :) = MVALUE; advVz(:, B%ny, :) = MVALUE

        if (twif) then
            where(tmask) err_nldecomp = MVALUE
            err_nldecomp(1:2, :, :) = MVALUE; err_nldecomp(B%nx, :, :) = MVALUE
            err_nldecomp(:, 1:2, :) = MVALUE; err_nldecomp(:, B%nx, :) = MVALUE
        endif 

        write(*, fmtm_vor) 'advu: ', advu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        if (debug) then
            write(*, fmtm_vor) '  advu_x: ', advu_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  advu_y: ', advu_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
        write(*, fmtm_vor) 'advv: ', advv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        if (debug) then
            write(*, fmtm_vor) '  advv_x: ', advv_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  advv_y: ', advv_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
        write(*, fmtm_vor) 'advw: ', advw(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        if (debug) then
            write(*, fmtm_vor) '  advw_x: ', advw_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  advw_y: ', advw_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif

        write(*, fmtm_vor) 'advVx: ', advVx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        if (debug) then
            write(*, fmtm_vor) '  advVx_x : ', advVx_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  errnl_ux: ', errnl_ux(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  advVx_y : ', advVx_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  errnl_uy: ', errnl_uy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
        write(*, fmtm_vor) 'advVy: ', advVy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        if (debug) then
            write(*, fmtm_vor) '  advVy_x : ', advVy_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  errnl_vx: ', errnl_vx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  advVy_y : ', advVy_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  errnl_vy: ', errnl_vy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
        write(*, fmtm_vor) 'advVz: ', advVz(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        if (debug) then
            write(*, fmtm_vor) '  advVz_x : ', advVz_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  errnl_wx: ', errnl_wx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  advVz_y : ', advVz_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  errnl_wy: ', errnl_wy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif

        ! write(*, fmtm_vor) 'advX: ', advu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped) + advVx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'advY: ', advv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped) + advVy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        ! write(*, fmtm_vor) 'vDdivDx', vDdivDx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'uDdivDy', uDdivDy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsx1', rsx1(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsx2', rsx2(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsy1', rsy1(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsy2', rsy2(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsz1', rsz1(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsz2', rsz2(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        deallocate(ue, vn, wt, ume, vme, umn, vmn, umt, vmt)
    endsubroutine

    subroutine verify_nonl()
        implicit none
        integer :: iz
        real(kind=kd_r) :: total

        do iz = B%zi_dpst, B%zi_dped
            write(*, *)
            if (allocated(advu)) then
                total = advu(B%xi_dp, B%yi_dp, iz) + advv(B%xi_dp, B%yi_dp, iz) + advw(B%xi_dp, B%yi_dp, iz) + &
                        advVx(B%xi_dp, B%yi_dp, iz) + advVy(B%xi_dp, B%yi_dp, iz) + advVz(B%xi_dp, B%yi_dp, iz)
                write(*, '(A, I02)') '      iz = ',  iz
                write(*, fmts_vor) 'All components: ',  total
                write(*, fmts_vor) 'All errors: ', err_nldecomp(B%xi_dp, B%yi_dp, iz)
                if (allocated(curladv)) then
                    write(*, fmts_vor) 'Diff: ', curladv(B%xi_dp, B%yi_dp, iz) - total - err_nldecomp(B%xi_dp, B%yi_dp, iz)
                endif
            endif

            if (allocated(curladv)) then
                write(*, fmts_vor) 'Curladv: ', curladv(B%xi_dp, B%yi_dp, iz)
            endif

            if (allocated(curlmet)) then
                write(*, fmts_vor) 'Curlmet: ', curlmet(B%xi_dp, B%yi_dp, iz)
            endif

            if (allocated(err_nlsub)) then
                write(*, fmts_vor) 'err_nlsub: ', err_nlsub(B%xi_dp, B%yi_dp, iz)
            endif

            if (allocated(curladvf)) then
                write(*, fmts_vor) 'Curladv (f): ', curladvf(B%xi_dp, B%yi_dp, iz)
            endif
        enddo
    endsubroutine

    subroutine init_zetavars_input(func)
        implicit none
        character(len=*), intent(in), optional :: func
        character(len=10) :: func_c

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Initializing input variables'
        if (present(func) .and. trim(func) /= "") then
            func_c = func
        else
            call warning_msg("noinput")
            return
        endif

        if (index(func_c, "a") /= 0 .or. index(func_c, "m") /= 0 .or. &
            index(func_c, "d") /= 0 .or. index(func_c, "c") /= 0) then
            allocate(uc(B%nx, B%ny, B%nz), vc(B%nx, B%ny, B%nz))
            uc = 0.; vc = 0.
        endif
        if (index(func_c, "c") /= 0) then
            allocate(advx  (B%nx, B%ny, B%nz), advy  (B%nx, B%ny, B%nz), &
                     gradx (B%nx, B%ny, B%nz), grady (B%nx, B%ny, B%nz), &
                     hdiffx(B%nx, B%ny, B%nz), hdiffy(B%nx, B%ny, B%nz), &
                     vdiffx(B%nx, B%ny, B%nz), vdiffy(B%nx, B%ny, B%nz), &
                     ssh   (B%nx, B%ny, 1 ))
            advx   = 0.; advy   = 0.; gradx  = 0.; grady  = 0.
            hdiffx = 0.; hdiffy = 0.; vdiffx = 0.; vdiffy = 0.
            ssh = 0.
        endif
        if (index(func_c, "a") /= 0 .or. index(func_c, "d") /= 0) then
            allocate(wc(B%nx, B%ny, B%nz))
            wc = 0.
        endif
        if (index(func_c, "f") /= 0) then
            allocate(ueu(B%nx, B%ny, B%nz), uev(B%nx, B%ny, B%nz), &
                     vnu(B%nx, B%ny, B%nz), vnv(B%nx, B%ny, B%nz), &
                     wtu(B%nx, B%ny, B%nz), wtv(B%nx, B%ny, B%nz))
            ueu = 0.; uev = 0.; vnu = 0.; vnv = 0.; wtu = 0.; wtv = 0.
        endif
    endsubroutine

    subroutine init_zetavars_output(func)
        implicit none
        character(len=*), intent(in), optional :: func
        character(len=10) :: func_c

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Initializing output variables'
        if (present(func) .and. trim(func) /= "") then
            func_c = func
        else
            call warning_msg("noinput")
            return
        endif

        if (index(func_c, "c") /= 0) then
            allocate(curlnonl (B%nx, B%ny, B%nz), curlpgrad(B%nx, B%ny, B%nz), res(B%nx, B%ny, B%nz), &
                     curlhdiff(B%nx, B%ny, B%nz), curlvdiff(B%nx, B%ny, B%nz), &
                     betav(B%nx, B%ny, B%nz), stretchp(B%nx, B%ny, B%nz), err_cor(B%nx, B%ny, B%nz))
            curlnonl = 0.; curlpgrad = 0.; curlhdiff = 0.; curlvdiff = 0.; res = 0.
            betav    = 0.; stretchp  = 0.; err_cor  = 0.
        endif
        if (index(func_c, "e") /= 0) then
            allocate(err_nlsub(B%nx, B%ny, B%nz))
            err_nlsub = 0.
        endif
        if (index(func_c, "m") /= 0) then
            allocate(curlmet(B%nx, B%ny, B%nz))
            curlmet = 0.
        endif
        if (index(func_c, "a") /= 0) then
            allocate(curladv(B%nx, B%ny, B%nz))
            curladv = 0.
            if (index(func_c, "a-") /= 0) then
                allocate(curladvu(B%nx, B%ny, B%nz), &
                         curladvv(B%nx, B%ny, B%nz), &
                         curladvw(B%nx, B%ny, B%nz))
                curladvu = 0.; curladvv = 0.; curladvw = 0.
            endif
        endif
        if (index(func_c, "f") /= 0) then
            allocate(curladvf(B%nx, B%ny, B%nz))
            curladvf = 0.
        endif
        if (index(func_c, "d") /= 0) then
            allocate(advu (B%nx, B%ny, B%nz), advv (B%nx, B%ny, B%nz), advw (B%nx, B%ny, B%nz), &
                     advVx(B%nx, B%ny, B%nz), advVy(B%nx, B%ny, B%nz), advVz(B%nx, B%ny, B%nz), err_nldecomp(B%nx, B%ny, B%nz))
            advu  = 0.; advv  = 0.; advw  = 0.
            advVx = 0.; advVy = 0.; advVz = 0.; err_nldecomp = 0.
            if (index(func_c, "d-") /= 0) then 
                allocate(advu_x  (B%nx, B%ny, B%nz), advv_x  (B%nx, B%ny, B%nz), advw_x  (B%nx, B%ny, B%nz), &
                         advVx_x (B%nx, B%ny, B%nz), advVy_x (B%nx, B%ny, B%nz), advVz_x (B%nx, B%ny, B%nz), &
                         errnl_ux(B%nx, B%ny, B%nz), errnl_vx(B%nx, B%ny, B%nz), errnl_wx(B%nx, B%ny, B%nz), &
                         advu_y  (B%nx, B%ny, B%nz), advv_y  (B%nx, B%ny, B%nz), advw_y  (B%nx, B%ny, B%nz), &
                         advVx_y (B%nx, B%ny, B%nz), advVy_y (B%nx, B%ny, B%nz), advVz_y (B%nx, B%ny, B%nz), &
                         errnl_uy(B%nx, B%ny, B%nz), errnl_vy(B%nx, B%ny, B%nz), errnl_wy(B%nx, B%ny, B%nz))
                advu_x = 0.; advv_x = 0; advw_x = 0.; advVx_x = 0.; advVy_x = 0.; advVz_x = 0.
                errnl_ux = 0.; errnl_vx = 0.; errnl_wx = 0.
                advu_y = 0.; advv_y = 0; advw_y = 0.; advVx_y = 0.; advVy_y = 0.; advVz_y = 0.;
                errnl_uy = 0.; errnl_vy = 0.; errnl_wy = 0.  
            endif                   
        endif
    endsubroutine

    subroutine release_zetavars_input(func)
        implicit none
        character(len=*), intent(in), optional :: func
        character(len=10) :: func_c

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Releasing input variables'
        if (present(func) .and. trim(func) /= "") then
            func_c = func
        else
            call warning_msg("noinput")
            return
        endif

        if (index(func_c, "a") /= 0 .or. index(func_c, "m") /= 0 .or. &
            index(func_c, "d") /= 0 .or. index(func_c, "c") /= 0) then
            deallocate(uc, vc)
        endif
        if (index(func_c, "c") /= 0) then
            deallocate(advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy, ssh)
        endif
        if (index(func_c, "a") /= 0 .or. index(func_c, "d") /= 0) then
            deallocate(wc)
        endif
        if (index(func_c, "f") /= 0) then
            deallocate(ueu, uev, vnu, vnv, wtu, wtv)
        endif
    endsubroutine

    subroutine release_zetavars_output(func)
        implicit none
        character(len=*), intent(in), optional :: func
        character(len=10) :: func_c

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Releasing output variables'
        if (present(func) .and. trim(func) /= "") then
            func_c = func
        else
            call warning_msg("noinput")
            return
        endif

        if (index(func_c, "c") /= 0) then
            deallocate(curlnonl, curlpgrad, res, curlhdiff, curlvdiff, betav, stretchp, err_cor)
        endif
        if (index(func_c, "e") /= 0) then
            deallocate(err_nlsub)
        endif
        if (index(func_c, "m") /= 0) then
            deallocate(curlmet)
        endif
        if (index(func_c, "a") /= 0) then
            deallocate(curladv)
            if (index(func_c, "a-") /= 0) then
                deallocate(curladvu, curladvv, curladvw)
            endif
        endif
        if (index(func_c, "f") /= 0) then
            deallocate(curladvf)
        endif
        if (index(func_c, "d") /= 0) then
            deallocate(advu, advv, advw, advVx, advVy, advVz, err_nldecomp)
            if (index(func_c, "d-") /= 0) then 
                deallocate(advu_x, advv_x, advw_x, advVx_x, advVy_x, advVz_x, errnl_ux, errnl_vx, errnl_wx, &
                           advu_y, advv_y, advw_y, advVx_y, advVy_y, advVz_y, errnl_uy, errnl_vy, errnl_wy)
           endif
        endif
    endsubroutine

    subroutine warning_msg(incident)
        character(len=*), intent(in) :: incident
        select case (trim(incident))
            case("noinput")
                write(*, '(A)') "Neither cmode or func is given. Exit."
        endselect
    endsubroutine
endmodule
