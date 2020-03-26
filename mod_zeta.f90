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

    ! Constants and grid info
    real(kind=kd_r), allocatable, public :: tlat(:,:), tlong(:,:), z_t(:)
    real(kind=kd_r), dimension(:, :)   , allocatable :: ulat, ulong, dxu, dyu, tarea, uarea, huw , hus
    real(kind=kd_r), dimension(:)      , allocatable :: z_w, dz
    real(kind=kd_r), dimension(:, :, :), allocatable :: dzu, dzt
    real(kind=kd_r), dimension(:, :)   , allocatable :: fcor, fcort
    real(kind=kd_r) :: grav
    real(kind=kd_r), dimension(:, :)   , allocatable :: dxue, dyue, tareae, &
                                                        dxun, dyun, tarean
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
            write(*, '(A)') '  Calculating curl of advection term (offline)'
            call calc_curladv()
        endif
        if (index(func_c, "d") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            if (index(func_c, "d#") /= 0) then
                write(*, '(A)') '  Decomposing curl of advection term (w/ nonflux twisting term)'
                call decomp_curladv(.False.)
            else
                write(*, '(A)') '  Decomposing curl of advection term (w/ flux twisting term)'
                call decomp_curladv(.True.)
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
        allocate(dxue(B%nx, B%ny), dyue(B%nx, B%ny), tareae(B%nx, B%ny), &
                 dxun(B%nx, B%ny), dyun(B%nx, B%ny), tarean(B%nx, B%ny))

        call nc_read(fngrid%grid, 'HTN', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        htn = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'HTE', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        hte = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'HUW', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        huw = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'HUS', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        hus = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'UAREA', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        uarea = WORK(:, :, 1, 1)

        dyue(1:B%nx-1, :) = huw(2:B%nx, :)
        dxue(1:B%nx-1, :) = htn(2:B%nx, :)
        tareae = hus * hte

        dyun(:, 1:B%ny-1) = hte(:, 2:B%ny)
        dxun(:, 1:B%ny-1) = hus(:, 2:B%ny)
        tarean = htn * huw
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

    subroutine calc_curladv()
        use popfun, only : zcurl
        implicit none
        real(kind=kd_r), dimension(B%nx, B%ny) :: advx, advy
        real(kind=kd_r), dimension(B%nx, B%ny) :: WORKx, WORKy, WORKz
        integer :: iz

        ! Calculating derived velocity at walls
        write(*, *)
        write(*, '(A)') '    -------------------------------------------------'
        write(*, '(A)') '    Calculating velocity at the walls'
        call calc_velw()

        do iz = 1, B%nz
            advx = 0.
            advy = 0.
            advx(2:B%nx, 2:B%ny) = (ue(2:B%nx  , 2:B%ny  , iz) * ume (2:B%nx  , 2:B%ny  , iz) - &
                                    ue(1:B%nx-1, 2:B%ny  , iz) * ume (1:B%nx-1, 2:B%ny  , iz) + &
                                    vn(2:B%nx  , 2:B%ny  , iz) * umn (2:B%nx  , 2:B%ny  , iz) - &
                                    vn(2:B%nx  , 1:B%ny-1, iz) * umn (2:B%nx  , 1:B%ny-1, iz))  &
                                   / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)            &
                                   + wt(2:B%nx, 2:B%ny, iz) * umt(2:B%nx, 2:B%ny, iz) / dzu(2:B%nx, 2:B%ny, iz)

            advy(2:B%nx, 2:B%ny) = (ue(2:B%nx  , 2:B%ny  , iz) * vme (2:B%nx  , 2:B%ny  , iz) - &
                                    ue(1:B%nx-1, 2:B%ny  , iz) * vme (1:B%nx-1, 2:B%ny  , iz) + &
                                    vn(2:B%nx  , 2:B%ny  , iz) * vmn (2:B%nx  , 2:B%ny  , iz) - &
                                    vn(2:B%nx  , 1:B%ny-1, iz) * vmn (2:B%nx  , 1:B%ny-1, iz))  &
                                  / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)            &
                                   + wt(2:B%nx, 2:B%ny, iz) * vmt(2:B%nx, 2:B%ny, iz) / dzu(2:B%nx, 2:B%ny, iz)

            if (iz < B%nz) then
                advx(2:B%nx, 2:B%ny) = advx(2:B%nx, 2:B%ny) &
                                       - wt(2:B%nx, 2:B%ny, iz+1) * umt(2:B%nx, 2:B%ny, iz+1) / dzu(2:B%nx, 2:B%ny, iz)
                advy(2:B%nx, 2:B%ny) = advy(2:B%nx, 2:B%ny) &
                                       - wt(2:B%nx, 2:B%ny, iz+1) * vmt(2:B%nx, 2:B%ny, iz+1) / dzu(2:B%nx, 2:B%ny, iz)
            endif

            advx(1, :) = MVALUE
            advx(:, 1) = MVALUE
            advy(1, :) = MVALUE
            advy(:, 1) = MVALUE

            where(umask(:,:,iz)) ! For cases where dzu = 0.
                advx = 0.
                advy = 0.
            endwhere

            curladv(:, :, iz) = zcurl(-advx, -advy, dxu*dzu(:,:,iz), dyu*dzu(:,:,iz), tarea*dzt(:,:,iz))
            curladv(1, :, iz) = MVALUE
            curladv(1, :, iz) = MVALUE

            ! if (iz >=  B%zi_dpst .and. iz <= B%zi_dped) then
            !      write(*, '(A, I02)') '      iz = ',  iz
            !      WORKx(2:B%nx, 2:B%ny) = (ue(2:B%nx  , 2:B%ny  , iz) * ume (2:B%nx  , 2:B%ny  , iz) - &
            !               ue(1:B%nx-1, 2:B%ny  , iz) * ume (1:B%nx-1, 2:B%ny  , iz))  &
            !                        / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)
            !
            !      WORKy(2:B%nx, 2:B%ny) = (ue(2:B%nx  , 2:B%ny  , iz) * vme (2:B%nx  , 2:B%ny  , iz) - &
            !               ue(1:B%nx-1, 2:B%ny  , iz) * vme (1:B%nx-1, 2:B%ny  , iz))  &
            !                        / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)
            !      WORKz = zcurl(-WORKx, -WORKy, dxu*dzu(:,:,iz), dyu*dzu(:,:,iz), tarea*dzt(:,:,iz))
            !      write(*, fmts_vor) 'x: ', WORKz(B%xi_dp, B%yi_dp)
            !
            !      WORKx(2:B%nx, 2:B%ny) = (vn(2:B%nx  , 2:B%ny  , iz) * umn (2:B%nx  , 2:B%ny  , iz) - &
            !               vn(2:B%nx  , 1:B%ny-1, iz) * umn (2:B%nx  , 1:B%ny-1, iz))  &
            !                        / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)
            !
            !      WORKy(2:B%nx, 2:B%ny) = (vn(2:B%nx  , 2:B%ny  , iz) * vmn (2:B%nx  , 2:B%ny  , iz) - &
            !               vn(2:B%nx  , 1:B%ny-1, iz) * vmn (2:B%nx  , 1:B%ny-1, iz))  &
            !                        / dzu(2:B%nx, 2:B%ny, iz) / uarea(2:B%nx, 2:B%ny)
            !      WORKz = zcurl(-WORKx, -WORKy, dxu*dzu(:,:,iz), dyu*dzu(:,:,iz), tarea*dzt(:,:,iz))
            !      write(*, fmts_vor) 'y: ', WORKz(B%xi_dp, B%yi_dp)
            !
            ! endif
        enddo
        where(tmask) curladv = MVALUE

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
            curladvf(1, :, iz) = MVALUE
        enddo
    endsubroutine

    subroutine decomp_curladv(twif)
        use derives
        use popfun, only : u2t
        implicit none
        logical, intent(in) :: twif
        real(kind=kd_r), dimension(B%nx, B%ny, 2) :: u10du2, u20du1
        real(kind=kd_r), dimension(B%nx, B%ny) :: u1, u2, u10, u20, test1, test2, test3
        real(kind=kd_r), dimension(B%nx, B%ny) :: ONES, var_F, var_B, &
              u10du2_F, u20du1_F, u10du2_B, u20du1_B, &
              u10du2_xt, u10du2_yt, u10du2_xb, u10du2_yb, u20du1_xt, u20du1_yt, u20du1_xb, u20du1_yb
        real(kind=kd_r), dimension(B%nx, B%ny) :: WORK, uudxdy, vvdxdy, wwdxdy, wvdxdz, wudydz, wm
        real(kind=kd_r), dimension(B%nx, B%ny, 2) :: WORK4zx, WORK4zy
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
            ! d(uv) / dx
            u1 = ue(:,:,iz)
            u2 = vme(:,:,iz)/dxu
            call dd_xw_chain(u1, mean_xw(u1), u2, mean_xw(u2), ONES, u10du2(:,:,2), u20du1(:,:,2))
            u1 = shift_xe(ue(:,:,iz))
            u2 = shift_xe(vme(:,:,iz))/dxu
            call dd_xw_chain(u1, mean_xw(u1), u2, mean_xw(u2), ONES, u10du2(:,:,1), u20du1(:,:,1))
            advu(:, :, iz) = advu(:, :, iz) + mean_ys(u10du2(:,:,2) - u10du2(:,:,1)) / tarea / dzt(:,:,iz)

            ! if (twif) then
            !     advVx(:, :, iz) = advVx(:, :, iz) + mean_ys(u20du1(:,:,2) - u20du1(:,:,1)) / tarea / dzt(:,:,iz)
            ! else
            !     u1 = mean_xw(vme(:,:,iz))/dxu
            !     u2 = dd_xw(ue(:,:,iz), ONES)
            !     call dd_xw_chain(u1, mean_xw(u1), u2, mean_xw(u2), ONES, u10du2(:,:,1), u20du1(:,:,1))
            !     advVx(:, :, iz) = advVx(:, :, iz) + mean_ys(u20du1(:,:,1))
            !     err_nldecomp(:,:,iz) = err_nldecomp(:,:,iz) + mean_ys(u10du2(:,:,1))
            ! endif

            test1 = mean_ys(u20du1(:,:,2) - u20du1(:,:,1)) / tarea / dzt(:,:,iz)
            test2 = mean_ys(mean_xw(mean_xw(vme(:,:,iz))/dxu) * dd_xw(dd_xw(ue(:,:,iz), ONES), ONES))) / tarea / dzt(:,:,iz)
            test3 = mean_ys(mean_xw(dd_xw(ue(:,:,iz), ONES)) * mean_xw(dd_xw(vme(:,:,iz), dxu))) / tarea / dzt(:,:,iz)

            print*, test1(B%xi_dp, B%yi_dp)
            print*, test2(B%xi_dp, B%yi_dp)
            print*, test3(B%xi_dp, B%yi_dp)
            print*, test1(B%xi_dp, B%yi_dp) - test2(B%xi_dp, B%yi_dp) - test3(B%xi_dp, B%yi_dp)
            ! advu(:, :, iz) = advu(:, :, iz) + 
            ! mean_sw(ue(:,:,iz)) * dd_xw(vme(:,:,iz)/dxu, ONES) - mean_sw(shift_xe(ue(:,:,iz))) * dd_xw(shift_xe(vme(:,:,iz))/dxu, ONES)

            ! d(uu) / dy
            u1 = ue(:,:,iz)
            u2 = ume(:,:,iz)/dyu
            call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2(:,:,2), u20du1(:,:,2))
            u1 = shift_xe(ue(:,:,iz))
            u2 = shift_xe(ume(:,:,iz))/dyu
            call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2(:,:,1), u20du1(:,:,1))
            advu(:, :, iz) = advu(:, :, iz) - mean_xw(u10du2(:,:,2) - u10du2(:,:,1)) / tarea / dzt(:,:,iz)

            ! if (twif) then
            !     advVx(:, :, iz) = advVx(:, :, iz) - mean_xw(u20du1(:,:,2) - u20du1(:,:,1)) / tarea / dzt(:,:,iz)
            ! else
            !     u1 = mean_ys(ume(:,:,iz))/dxu

          ! advv
            ! d(vv) / dx
            var_F = vmn(:,:,iz)/dxu
            var_B = shift_yn(vmn(:,:,iz))/dxu
            ! var_B(:,1) = MVALUE
            call dd_xw_chain(vn(:,:,iz), mean_xw(vn(:,:,iz)), var_F, mean_xw(var_F), &
                             ONES, u10du2_F, u20du1_F)
            call dd_xw_chain(shift_yn(vn(:,:,iz)), mean_xw(shift_yn(vn(:,:,iz))), var_B, mean_xw(var_B), &
                             ONES, u10du2_B, u20du1_B)

            advv  (:, :, iz) = advv  (:, :, iz) + mean_ys(u10du2_F - u10du2_B) / tarea / dzt(:,:,iz)
            advVy (:, :, iz) = advVy (:, :, iz) + mean_ys(u20du1_F - u20du1_B) / tarea / dzt(:,:,iz)

            ! d(vu) / dy
            var_F = umn(:,:,iz)/dyu
            var_B = shift_yn(umn(:,:,iz))/dyu
            ! var_B(:,1) = MVALUE
            call dd_ys_chain(vn(:,:,iz), mean_ys(vn(:,:,iz)), var_F, mean_ys(var_F), &
                             ONES, u10du2_F, u20du1_F)
            call dd_ys_chain(shift_yn(vn(:,:,iz)), mean_ys(shift_yn(vn(:,:,iz))), var_B, mean_ys(var_B), &
                             ONES, u10du2_B, u20du1_B)


            advv  (:, :, iz) = advv  (:, :, iz) - mean_xw(u10du2_F - u10du2_B) / tarea / dzt(:,:,iz)
            advVy (:, :, iz) = advVy (:, :, iz) - mean_xw(u20du1_F - u20du1_B) / tarea / dzt(:,:,iz)

          ! advw
            ! d(wv) / dx  &  d(wu) / dy
            if (iz == 1) then
                call dd_xw_chain(wt(:,:,iz), mean_xw(wt(:,:,iz)), vmt(:,:,iz)*dyu, mean_xw(vmt(:,:,iz)*dyu), &
                                 ONES, u10du2_xt, u20du1_xt)
                call dd_ys_chain(wt(:,:,iz), mean_ys(wt(:,:,iz)), umt(:,:,iz)*dxu, mean_ys(umt(:,:,iz)*dxu), &
                                 ONES, u10du2_yt, u20du1_yt)
            endif

            if (iz < B%nz) then
                call dd_xw_chain(wt(:,:,iz+1), mean_xw(wt(:,:,iz+1)), vmt(:,:,iz+1)*dyu, mean_xw(vmt(:,:,iz+1)*dyu), &
                                 ONES, u10du2_xb, u20du1_xb)
                call dd_ys_chain(wt(:,:,iz+1), mean_ys(wt(:,:,iz+1)), umt(:,:,iz+1)*dxu, mean_ys(umt(:,:,iz+1)*dxu), &
                                 ONES, u10du2_yb, u20du1_yb)
            else
                u10du2_xb = 0.
                u20du1_xb = 0.
                u10du2_yb = 0.
                u20du1_yb = 0.
            endif

            advw (:, :, iz) = (mean_ys(u10du2_xt - u10du2_xb) - mean_xw(u10du2_yt - u10du2_yb)) / tarea / dzt(:, :, iz)
            advVz(:, :, iz) = (mean_ys(u20du1_xt - u20du1_xb) - mean_xw(u20du1_yt - u20du1_yb)) / tarea / dzt(:, :, iz)

            u10du2_xt = u10du2_xb
            u10du2_yt = u10du2_yb
            u20du1_xt = u20du1_xb
            u20du1_yt = u20du1_yb

          ! Add terms to compose flux of 3D voricity twisting
            ! call ddx_w(0.5 * uc(:, :, iz) * uc(:, :, iz) * dzu(:, :, iz), ONES, WORK)
            ! call ddy_s(WORK, tarea * dzt(:,:,iz), uudxdy)

            ! call ddx_w(0.5 * vc(:, :, iz) * vc(:, :, iz) * dzu(:, :, iz), ONES, WORK)
            ! call ddy_s(WORK, tarea * dzt(:,:,iz), vvdxdy)
            
            uudxdy = dd_ys(dd_xw(0.5 * uc(:, :, iz) * uc(:, :, iz) * dzu(:, :, iz), ONES), &
                           tarea * dzt(:,:,iz))
            vvdxdy = dd_ys(dd_xw(0.5 * vc(:, :, iz) * vc(:, :, iz) * dzu(:, :, iz), ONES), &
                           tarea * dzt(:,:,iz))

            if (iz == B%nz) then
                wm = (wt(:, :, iz) + 0.) / 2
            else
                wm = (wt(:, :, iz) + wt(:, :, iz+1)) / 2
            endif
            ! call ddx_w(0.5 * wm * wm * dzu(:, :, iz), ONES, WORK)
            ! call ddy_s(WORK, tarea * dzt(:,:,iz), wwdxdy)

            wwdxdy = dd_ys(dd_xw(0.5 * wm * wm * dzu(:, :, iz), ONES), tarea * dzt(:,:,iz))

            if (iz == 1) then
                ! call ddx(wt(:, :, iz) * vmt(:, :, iz), dyu, tarea, WORK4zx(:,:,1))
                ! call ddy(wt(:, :, iz) * umt(:, :, iz), dxu, tarea, WORK4zy(:,:,1))
                WORK4zx(:,:,1) = dd_xsw(wt(:, :, iz) * vmt(:, :, iz) * dyu, tarea)
                WORK4zy(:,:,1) = dd_ysw(wt(:, :, iz) * umt(:, :, iz) * dxu, tarea)
            endif

            if (iz == B%nz) then
                WORK4zx(:,:,2) = 0.
                WORK4zy(:,:,2) = 0.
            else
                ! call ddx(wt(:, :, iz+1) * vmt(:, :, iz+1), dyu, tarea, WORK4zx(:,:,2))
                ! call ddy(wt(:, :, iz+1) * umt(:, :, iz+1), dxu, tarea, WORK4zy(:,:,2))
                WORK4zx(:,:,2) = dd_xsw(wt(:, :, iz+1) * vmt(:, :, iz+1) * dyu, tarea)
                WORK4zy(:,:,2) = dd_ysw(wt(:, :, iz+1) * umt(:, :, iz+1) * dxu, tarea)
            endif
            wvdxdz = (WORK4zx(:,:,1) - WORK4zx(:,:,2)) / dzt(:, :, iz)
            wudydz = (WORK4zy(:,:,1) - WORK4zy(:,:,2)) / dzt(:, :, iz)

            advVx(:, :, iz) = advVx(:, :, iz) + wvdxdz + vvdxdy - wwdxdy + uudxdy
            advVy(:, :, iz) = advVy(:, :, iz) - wudydz - vvdxdy + wwdxdy - uudxdy
            advVz(:, :, iz) = advVz(:, :, iz) - wvdxdz + wudydz

            ! stepping down in z direction
            WORK4zx(:,:,1) = WORK4zx(:,:,2)
            WORK4zy(:,:,1) = WORK4zy(:,:,2)
        enddo

        ! RHS
        advu = -advu; advv = -advv; advw = -advw;
        advVx = -advVx; advVy = -advVy; advVz = -advVz

        where(tmask)
            advu  = MVALUE
            advv  = MVALUE
            advw  = MVALUE
            advVx = MVALUE
            advVy = MVALUE
            advVz = MVALUE
        endwhere

        write(*, fmtm_vor) 'advu: ', advu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'advv: ', advv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'advw: ', advw(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcuv', rcuv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcuu', rcuu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcvv', rcvv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcvu', rcvu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcwv', rcwv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcwu', rcwu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        write(*, fmtm_vor) 'advVx: ', advVx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'advVy: ', advVy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'advVz: ', advVz(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        write(*, fmtm_vor) 'advX: ', advu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped) + advVx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'advY: ', advv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped) + advVy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

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
        endif
        if (index(func_c, "f") /= 0) then
            deallocate(curladvf)
        endif
        if (index(func_c, "d") /= 0) then
            deallocate(advu, advv, advw, advVx, advVy, advVz, err_nldecomp)
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
