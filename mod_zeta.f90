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
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
        rr_rev, rr_cha

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
    ! * Basic functions
    ! curl ("c"): curl of zeta equations terms, including two (+1) tersm from Coriolis
    !                (nonlinear, curl of pgrad, hdiff, vdiff, residual,
    !                 betav from Coriolis, stretching from Coriolis, error with Coriolis decomposition)
    ! offline adv ("a"): offline calculation of nonlinear advection term (curladv)
    ! offline met ("m"): offline calculation of nonlinear metric term (curlmet)
    ! decomp adv ("d"): decomposition of nonlinear advection term
    ! online adv from flux ("f"): "online" calculation of nonlinear advection term through momentum fluxes
    ! + err_nlsub ("e"): difference between curlnonl online and offline (should not be used as a standalone function)

    ! * Calculation mode
    ! cmode = 0: "c"
    ! cmode = 1: "c" + "a" + "m" + "e"
    ! cmode = 2: "c" + "d" + "m" + "e"
    ! cmode = 3: "d" + "m"
    ! cmode = 4: "a" + "m"
    ! cmode = 5: "f"
    ! cmode = -1: null mode
    !!--------------------------------------------------------------------------
    subroutine calc_zeta(cmode, func)
        implicit none
        integer, intent(in), optional :: cmode
        character, intent(in), optional :: func
        character(len=10) :: func_c

        write(fmts_vor, '(A, A, A)' ) '(A20, ', trim(fmt_exp), ')'
        write(fmtm_vor, '(A, I2, A, A)') '(A20, ', B%zi_dped - B%zi_dpst + 1, trim(fmt_exp), ')'
        write(fmts_vel, '(A, A, A)' ) '(A20, ', trim(fmt_flt), ')'
        write(fmtm_vel, '(A, I2, A, A)') '(A20, ', B%zi_dped - B%zi_dpst + 1, trim(fmt_flt), ')'

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)'), 'Start calculating voriticity equation'

        if (present(cmode)) then
            if (cmode == -1) then
                call warning_msg("null")
                return
            endif
            call decode_cmode(cmode, func_c)
        elseif (present(func)) then
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
            write(*, '(A)') '  Decomposing curl of advection term'
            call decomp_curladv()
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

        write(*, *)
        write(*, '(A)') '  ---------------------------------------------------'
        write(*, '(A)') '  Verifying offline calculations'
        call verify_nonl()
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
        use derives, only : ddx_chain, ddy_chain
        implicit none
        real(kind=kd_r), dimension(:, :, :), allocatable :: curlcor, corx, cory
        real(kind=kd_r), dimension(B%nx, B%ny) :: pgradsfx, pgradsfy
        real(kind=kd_r), dimension(B%nx, B%ny) :: stretchpx, betavx, stretchpy, betavy, &
                                              WORKrr, WORKrrcx, WORKrrcy, u0, v0, ONES
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

        pgradsfx = 0.
        pgradsfy = 0.
        pgradsfx(1:B%nx-1, 1:B%ny-1) = -grav * ((ssh(2:B%nx,   1:B%ny-1, 1) + ssh(2:B%nx,   2:B%ny  , 1))/2 - &
                                            (ssh(1:B%nx-1, 1:B%ny-1, 1) + ssh(1:B%nx-1, 2:B%ny  , 1))/2) / dxu(1:B%nx-1, 1:B%ny-1)
        pgradsfy(1:B%nx-1, 1:B%ny-1) = -grav * ((ssh(1:B%nx-1, 2:B%ny  , 1) + ssh(2:B%nx,   2:B%ny  , 1))/2 - &
                                            (ssh(1:B%nx-1, 1:B%ny-1, 1) + ssh(2:B%nx,   1:B%ny-1, 1))/2) / dyu(1:B%nx-1, 1:B%ny-1)

        do iz = 1, B%nz
            gradx(:, :, iz) = pgradsfx - gradx(:, :, iz)
            grady(:, :, iz) = pgradsfy - grady(:, :, iz)
        enddo

        where(umask) gradx = 0.
        where(umask) grady = 0.

        write(*, fmtm_vor) 'pgradx', gradx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'pgrady', grady(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        !!----------------------------------------------------------------------
        ! Curl of nonlinear, cor, pgrad, h/v diffusion and res terms
        write(*, *)
        write(*, '(A)') '    Calculating curl of nonlinear, Coriolis, pgrad, horizontal/vertical diffusion and residual terms'
        allocate(curlcor(B%nx, B%ny, B%nz))
        curlcor = 0.
        do iz = 1, B%nz
            curlnonl (:, :, iz) = zcurl(-advx (:, :, iz), -advy (:, :, iz), dxu, dyu, tarea)
            curlcor  (:, :, iz) = zcurl(-corx (:, :, iz), -cory (:, :, iz), dxu, dyu, tarea)
            curlpgrad(:, :, iz) = zcurl(gradx (:, :, iz), grady (:, :, iz), dxu, dyu, tarea)
            curlhdiff(:, :, iz) = zcurl(hdiffx(:, :, iz), hdiffy(:, :, iz), dxu, dyu, tarea)
            curlvdiff(:, :, iz) = zcurl(vdiffx(:, :, iz), vdiffy(:, :, iz), dxu, dyu, tarea)
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
            u0 = u2t(uc(:,:,iz), dyu, ONES)
            v0 = u2t(vc(:,:,iz), dxu, ONES)

            call ddx_chain(fcor, fcort, uc(:,:,iz) * dyu, u0, ONES, tarea, &
                   stretchpx, betavx, WORKrrcx, WORKrr)
            call ddy_chain(fcor, fcort, vc(:,:,iz) * dxu, v0, ONES, tarea, &
                   stretchpy, betavy, WORKrrcy, WORKrr)

            betav   (:, :, iz) = -(betavx + betavy)
            stretchp(:, :, iz) = -(stretchpx + stretchpy)
            err_cor (: ,:, iz) = -(WORKrrcx + WORKrrcy)
        enddo

        where(tmask)
            betav    = MVALUE
            stretchp = MVALUE
            err_cor  = MVALUE
        endwhere

        !!----------------------------------------------------------------------
        ! Verifying decomposition and clean up
        write(*, *)
        write(*, '(A)') '    Clean up and verifying decomposition of curl(-fu, fv) ...'
        do iz = B%zi_dpst, B%zi_dped
            write(*, fmts_vor) 'Curlcor: ', curlcor  (B%xi_dp, B%yi_dp, iz)
            write(*, fmts_vor) 'bv + -fdwdz: ', betav(B%xi_dp, B%yi_dp, iz) + stretchp(B%xi_dp, B%yi_dp, iz)
            write(*, fmts_vor) 'errcor: ' , err_cor(B%xi_dp, B%yi_dp, iz)
            write(*, fmts_vor) 'Diff: '   , curlcor (B%xi_dp, B%yi_dp, iz) - &
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

            curlmet(:, :, iz) = zcurl(-metx, -mety, dxu, dyu, tarea)
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

            curladv(:, :, iz) = zcurl(-advx, -advy, dxu, dyu, tarea)
            curladv(1, :, iz) = MVALUE
            curladv(1, :, iz) = MVALUE
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

            curladvf(:, :, iz) = zcurl(-advx, -advy, dxu, dyu, tarea)
            curladvf(1, :, iz) = MVALUE
            curladvf(1, :, iz) = MVALUE
        enddo
    endsubroutine

    subroutine decomp_curladv()
        use derives
        use popfun, only : u2t
        implicit none
        real(kind=kd_r) :: WORK(B%nx, B%ny), WORKB(B%nx, B%ny), WORKTx(B%nx, B%ny), WORKTy(B%nx, B%ny), &
                           advthencur(B%nx, B%ny), curthenadv(B%nx, B%ny)
        real(kind=kd_r), dimension(B%nx, B%ny) :: ue0, vn0, wt0, ume0, vme0, umn0, vmn0, wm, &
                              WORK1, WORK2, WORK3, WORK1a, WORK2a, WORK3a, DUMMY, ONES
        real(kind=kd_r), dimension(B%nx, B%ny) :: DueDx_vme0, DueDy_ume0, DvnDx_vmn0, DvnDy_umn0
        real(kind=kd_r), dimension(B%nx, B%ny, 2) :: vmt0, umt0, &
                WORK1zx, WORK1zy, WORK2zx, WORK2zy, WORK3zx, WORK3zy, WORK4zx, WORK4zy
        real(kind=kd_r), dimension(B%nx, B%ny) :: uudxdy, vvdxdy, wwdxdy, wvdxdz, wudydz
        real(kind=kd_r), dimension(B%nx, B%ny, B%nz) :: calc_total
        integer :: iz

        write(*, *)
        write(*, '(A)') '    -------------------------------------------------'
        write(*, '(A)') '    Calculating velocity at the walls'
        call calc_velw()

        !!----------------------------------------------------------------------
        ! Calculating error from reversing operators
        !    *Efficiency of this function can be improved by direct calculation. See doc for derivation.
        write(*, *)
        write(*, '(A)') '    Calculating error term from reversing operators'
        do iz = 1, B%nz
            ! d(duv/dx) / dx
            call ddx_w(ue(:,:,iz) * vme(:,:,iz), uarea * dzu(:,:,iz), WORK)
            where(umask(:,:,iz)) WORK = 0.
            call ddx(WORK, dyu, tarea, advthencur)

            call ddx(ue(:,:,iz) * vme(:,:,iz), dyue, tareae, WORK)
            call ddx_w(WORK, tarea * dzt(:,:,iz), curthenadv)

            ! rrxx(:, :, iz) = advthencur - curthenadv
            rr_rev(:, :, iz) = rr_rev(:, :, iz) + (advthencur - curthenadv)


            ! d(duu/dx) / dy
            call ddx_w(ue(:,:,iz) * ume(:,:,iz), uarea * dzu(:,:,iz), WORK)
            where(umask(:,:,iz)) WORK = 0.
            call ddy(WORK, dxu, tarea, advthencur)

            call ddy(ue(:,:,iz) * ume(:,:,iz), dxue, tareae, WORK)
            call ddx_w(WORK, tarea * dzt(:,:,iz), curthenadv)

            ! rrxy(:, :, iz) = advthencur - curthenadv
            rr_rev(:, :, iz) = rr_rev(:, :, iz) - (advthencur - curthenadv)


            ! d(dvv/dy) / dx
            call ddy_s(vn(:,:,iz) * vmn(:,:,iz), uarea * dzu(:,:,iz), WORK)
            where(umask(:,:,iz)) WORK = 0.
            call ddx(WORK, dyu, tarea, advthencur)

            call ddx(vn(:,:,iz) * vmn(:,:,iz), dyun, tarean, WORK)
            call ddy_s(WORK, tarea * dzt(:,:,iz), curthenadv)

            ! rryx(:, :, iz) = advthencur - curthenadv
            rr_rev(:, :, iz) = rr_rev(:, :, iz) + (advthencur - curthenadv)


            ! d(dvu/dy) / dy
            call ddy_s(vn(:,:,iz) * umn(:,:,iz), uarea * dzu(:,:,iz), WORK)
            where(umask(:,:,iz)) WORK = 0.
            call ddy(WORK, dxu, tarea, advthencur)

            call ddy(vn(:,:,iz) * umn(:,:,iz), dxun, tarean, WORK)
            call ddy_s(WORK, tarea * dzt(:,:,iz), curthenadv)

            ! rryy(:, :, iz) = advthencur - curthenadv
            rr_rev(:, :, iz) = rr_rev(:, :, iz) - (advthencur - curthenadv)


            ! d(dwv/dz) / dx
            if (iz == B%nz) then
                WORK = (wt(:,:,iz) * vmt(:,:,iz)) / dzu(:,:,iz)
            else
                WORK = (wt(:,:,iz) * vmt(:,:,iz) - wt(:,:,iz+1) * vmt(:,:,iz+1)) / dzu(:,:,iz)
            endif
            where(umask(:,:,iz)) WORK = 0.
            call ddx(WORK, dyu, tarea, advthencur)

            if (iz == 1) then
                call ddx(wt(:,:,iz) * vmt(:,:,iz), dyu, tarea, WORKTx)
            endif
            if (iz == B%nz) then
                WORKB = 0.
            else
                call ddx(wt(:,:,iz+1) * vmt(:,:,iz+1), dyu, tarea, WORKB)
            endif
            curthenadv = (WORKTx - WORKB) / dzt(:,:,iz)
            WORKTx = WORKB

            ! rrzx(:, :, iz) = advthencur - curthenadv
            rr_rev(:, :, iz) = rr_rev(:, :, iz) + (advthencur - curthenadv)

            ! d(dwu/dz) / dy
            if (iz == B%nz) then
                WORK = (wt(:,:,iz) * umt(:,:,iz)) / dzu(:,:,iz)
            else
                WORK = (wt(:,:,iz) * umt(:,:,iz) - wt(:,:,iz+1) * umt(:,:,iz+1)) / dzu(:,:,iz)
            endif
            where(umask(:,:,iz)) WORK = 0.
            call ddy(WORK, dxu, tarea, advthencur)

            if (iz == 1) then
                call ddy(wt(:,:,iz) * umt(:,:,iz), dxu, tarea, WORKTy)
            endif
            if (iz == B%nz) then
                WORKB = 0.
            else
                call ddy(wt(:,:,iz+1) * umt(:,:,iz+1), dxu, tarea, WORKB)
            endif
            curthenadv = (WORKTy - WORKB) / dzt(:,:,iz)
            WORKTy = WORKB

            ! rrzy(:, :, iz) = advthencur - curthenadv
            rr_rev(:, :, iz) = rr_rev(:, :, iz) - (advthencur - curthenadv)
        enddo

        where(abs(dzt) < 1e-10) rr_rev = MVALUE
        rr_rev = -rr_rev ! RHS
        write(*, fmtm_vor) 'rr_rev: ', rr_rev(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rrxx', rrxx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rrxy', rrxy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rryx', rryx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rryy', rryy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rrzx', rrzx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rrzy', rrzy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        !!----------------------------------------------------------------------
        ! Decomposition
        write(*, *)
        write(*, '(A)') '    Calculating decomposition'
        ONES = 1.

        do iz = 1, B%nz
          !! Chain rule: advection terms
          ! advu
            ue0 = u2t(ue(:,:,iz), ONES, ONES)
            vme0 = u2t(vme(:,:,iz), dyue, ONES)
            ume0 = u2t(ume(:,:,iz), dxue, ONES)

            ! d(uv) / dx
            call ddx_chain(ue(:,:,iz), ue0, vme(:,:,iz) * dyue, vme0, &
                   ONES, tareae, WORK1, WORK2, WORK3, DUMMY)
            call ddx_w(WORK1, tarea * dzt(:,:,iz), WORK1a)
            call ddx_w(WORK2, tarea * dzt(:,:,iz), WORK2a)
            call ddx_w(WORK3, tarea * dzt(:,:,iz), WORK3a)

            advu  (:, :, iz) = advu  (:, :, iz) + WORK1a
            advVx (:, :, iz) = advVx (:, :, iz) + WORK2a
            rr_cha(:, :, iz) = rr_cha(:, :, iz) + WORK3a

            ! d(uu) / dy
            call ddy_chain(ue(:,:,iz), ue0, ume(:,:,iz) * dxue, ume0, &
                   ONES, tareae, WORK1, WORK2, WORK3, DUMMY)
            call ddx_w(WORK1, tarea * dzt(:,:,iz), WORK1a)
            call ddx_w(WORK2, tarea * dzt(:,:,iz), WORK2a)
            call ddx_w(WORK3, tarea * dzt(:,:,iz), WORK3a)

            advu  (:, :, iz) = advu  (:, :, iz) - WORK1a
            advVx (:, :, iz) = advVx (:, :, iz) - WORK2a
            rr_cha(:, :, iz) = rr_cha(:, :, iz) - WORK3a


          ! advv
            vn0 = u2t(vn(:,:,iz), ONES, ONES)
            vmn0 = u2t(vmn(:,:,iz), dyun, ONES)
            umn0 = u2t(umn(:,:,iz), dxun, ONES)

            ! d(vv) / dx
            call ddx_chain(vn(:,:,iz), vn0, vmn(:,:,iz) * dyun, vmn0, &
                   ONES, tarean, WORK1, WORK2, WORK3, DUMMY)
            call ddy_s(WORK1, tarea * dzt(:,:,iz), WORK1a)
            call ddy_s(WORK2, tarea * dzt(:,:,iz), WORK2a)
            call ddy_s(WORK3, tarea * dzt(:,:,iz), WORK3a)

            advv  (:, :, iz) = advv  (:, :, iz) + WORK1a
            advVy (:, :, iz) = advVy (:, :, iz) + WORK2a
            rr_cha(:, :, iz) = rr_cha(:, :, iz) + WORK3a

            ! d(vu) / dy
            call ddy_chain(vn(:,:,iz), vn0, umn(:,:,iz) * dxun, umn0, &
                   ONES, tarean, WORK1, WORK2, WORK3, DUMMY)
            call ddy_s(WORK1, tarea * dzt(:,:,iz), WORK1a)
            call ddy_s(WORK2, tarea * dzt(:,:,iz), WORK2a)
            call ddy_s(WORK3, tarea * dzt(:,:,iz), WORK3a)

            advv  (:, :, iz) = advv  (:, :, iz) - WORK1a
            advVy (:, :, iz) = advVy (:, :, iz) - WORK2a
            rr_cha(:, :, iz) = rr_cha(:, :, iz) - WORK3a


          ! advw
            ! d(wv) / dx  &  d(wu) / dy
            if (iz == 1) then
                wt0 = u2t(wt(:,:,iz), uarea, tarea)
                vmt0(:,:,1) = u2t(vmt(:,:,iz), dyu, ONES)
                umt0(:,:,1) = u2t(umt(:,:,iz), dxu, ONES)

                call ddx_chain(wt(:,:,iz), wt0, vmt(:,:,iz) * dyu, vmt0(:,:,1), &
                       ONES, tarea, WORK1zx(:,:,1), WORK2zx(:,:,1), WORK3zx(:,:,1), DUMMY)
                call ddy_chain(wt(:,:,iz), wt0, umt(:,:,iz) * dxu, umt0(:,:,1), &
                       ONES, tarea, WORK1zy(:,:,1), WORK2zy(:,:,1), WORK3zy(:,:,1), DUMMY)
            endif

            if (iz == B%nz) then
                WORK1zx(:,:,2) = 0.
                WORK2zx(:,:,2) = 0.
                WORK3zx(:,:,2) = 0.
                WORK1zy(:,:,2) = 0.
                WORK2zy(:,:,2) = 0.
                WORK3zy(:,:,2) = 0.
            else
                wt0 = u2t(wt(:,:,iz+1), uarea, tarea)
                vmt0(:,:,2) = u2t(vmt(:,:,iz+1), dyu, ONES)
                umt0(:,:,2) = u2t(umt(:,:,iz+1), dxu, ONES)

                call ddx_chain(wt(:,:,iz+1), wt0, vmt(:,:,iz+1) * dyu, vmt0(:,:,2), &
                       ONES, tarea, WORK1zx(:,:,2), WORK2zx(:,:,2), WORK3zx(:,:,2), DUMMY)
                call ddy_chain(wt(:,:,iz+1), wt0, umt(:,:,iz+1) * dxu, umt0(:,:,2), &
                       ONES, tarea, WORK1zy(:,:,2), WORK2zy(:,:,2), WORK3zy(:,:,2), DUMMY)
            endif

            advw (:, :, iz) = ((WORK1zx(:,:,1) - WORK1zx(:,:,2)) - (WORK1zy(:,:,1) - WORK1zy(:,:,2))) / dzt(:, :, iz)
            advVz(:, :, iz) = ((WORK2zx(:,:,1) - WORK2zx(:,:,2)) - (WORK2zy(:,:,1) - WORK2zy(:,:,2))) / dzt(:, :, iz)
            rr_cha(:, :, iz) = rr_cha(:, :, iz) + &
              ((WORK3zx(:,:,1) - WORK3zx(:,:,2)) - (WORK3zy(:,:,1) - WORK3zy(:,:,2))) / dzt(:, :, iz)


          ! Add terms to compose flux of 3D voricity twisting
            call ddx_w(uc(:, :, iz) * uc(:, :, iz) * dzu(:, :, iz), ONES, WORK1)
            call ddy_s(WORK1, tarea * dzt(:,:,iz), uudxdy)

            call ddx_w(vc(:, :, iz) * vc(:, :, iz) * dzu(:, :, iz), ONES, WORK1)
            call ddy_s(WORK1, tarea * dzt(:,:,iz), vvdxdy)

            if (iz == B%nz) then
                wm = (wt(:, :, iz) + 0.) / 2
            else
                wm = (wt(:, :, iz) + wt(:, :, iz+1)) / 2
            endif
            call ddx_w(wm * wm * dzu(:, :, iz), ONES, WORK1)
            call ddy_s(WORK1, tarea * dzt(:,:,iz), wwdxdy)

            if (iz == 1) then
                call ddx(wt(:, :, iz) * vmt(:, :, iz), dyu, tarea, WORK4zx(:,:,1))
                call ddy(wt(:, :, iz) * umt(:, :, iz), dxu, tarea, WORK4zy(:,:,1))
            endif

            if (iz == B%nz) then
                WORK4zx(:,:,2) = 0.
                WORK4zy(:,:,2) = 0.
            else
                call ddx(wt(:, :, iz+1) * vmt(:, :, iz+1), dyu, tarea, WORK4zx(:,:,2))
                call ddy(wt(:, :, iz+1) * umt(:, :, iz+1), dxu, tarea, WORK4zy(:,:,2))
            endif
            wvdxdz = (WORK4zx(:,:,1) - WORK4zx(:,:,2)) / dzt(:, :, iz)
            wudydz = (WORK4zy(:,:,1) - WORK4zy(:,:,2)) / dzt(:, :, iz)

            advVx(:, :, iz) = advVx(:, :, iz) + wvdxdz + vvdxdy - wwdxdy + uudxdy
            advVy(:, :, iz) = advVy(:, :, iz) - wudydz - vvdxdy + wwdxdy - uudxdy
            advVz(:, :, iz) = advVz(:, :, iz) - wvdxdz + wudydz

            ! stepping down in z direction
            WORK1zx(:,:,1) = WORK1zx(:,:,2)
            WORK1zy(:,:,1) = WORK1zy(:,:,2)
            WORK2zx(:,:,1) = WORK2zx(:,:,2)
            WORK2zy(:,:,1) = WORK2zy(:,:,2)
            WORK3zx(:,:,1) = WORK3zx(:,:,2)
            WORK3zy(:,:,1) = WORK3zy(:,:,2)

            WORK4zx(:,:,1) = WORK4zx(:,:,2)
            WORK4zy(:,:,1) = WORK4zy(:,:,2)
        enddo

        ! RHS
        advu = -advu; advv = -advv; advw = -advw; rr_cha = -rr_cha
        advVx = -advVx; advVy = -advVy; advVz = -advVz

        where(tmask)
            advu  = MVALUE
            advv  = MVALUE
            advw  = MVALUE
            rr_cha = MVALUE
            advVx = MVALUE
            advVy = MVALUE
            advVz = MVALUE
        endwhere

        write(*, fmtm_vor) 'advu: ', advu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'advv: ', advv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'advw: ', advw(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'rr_cha: ', rr_cha(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcuv', rcuv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcuu', rcuu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcvv', rcvv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcvu', rcvu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcwv', rcwv(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rcwu', rcwu(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        write(*, fmtm_vor) 'advVx: ', advVx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'advVy: ', advVy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        write(*, fmtm_vor) 'advVz: ', advVz(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        ! write(*, fmtm_vor) 'vDdivDx', vDdivDx(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'uDdivDy', uDdivDy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsx1', rsx1(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsx2', rsx2(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsy1', rsy1(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsy2', rsy2(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsz1', rsz1(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        ! write(*, fmtm_vor) 'rsz2', rsz2(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)

        deallocate(ue, vn, wt, ume, vme, umn, vmn, umt, vmt)
        err_nldecomp = rr_rev + rr_cha
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

    subroutine init_zetavars_input(cmode, func)
        implicit none
        integer, intent(in), optional :: cmode
        character(len=*), intent(in), optional :: func
        character(len=10) :: func_c

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Initializing input variables'
        if (present(cmode)) then
            if (cmode == -1) then
                call warning_msg("null")
                return
            else
                call decode_cmode(cmode, func_c)
                if (func_c == "") then
                    call warning_msg("invalid")
                    return
                endif
            endif
        elseif (present(func)) then
            if (func == "") then
                call warning_msg("null")
                return
            else
                func_c = func
            endif
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

    subroutine init_zetavars_output(cmode, func)
        implicit none
        integer, intent(in), optional :: cmode
        character(len=*), intent(in), optional :: func
        character(len=10) :: func_c

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Initializing output variables'
        if (present(cmode)) then
            if (cmode == -1) then
                call warning_msg("null")
                return
            else
                call decode_cmode(cmode, func_c)
                if (func_c == "") then
                    call warning_msg("invalid")
                    return
                endif
            endif
        elseif (present(func)) then
            if (func == "") then
                call warning_msg("null")
                return
            else
                func_c = func
            endif
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
            allocate(rr_rev(B%nx, B%ny, B%nz), rr_cha(B%nx, B%ny, B%nz))
            rr_rev = 0.; rr_cha = 0.
        endif
    endsubroutine

    subroutine release_zetavars_input(cmode, func)
        implicit none
        integer, intent(in), optional :: cmode
        character(len=*), intent(in), optional :: func
        character(len=10) :: func_c

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Releasing input variables'
        if (present(cmode)) then
            if (cmode == -1) then
                call warning_msg("null")
                return
            else
                call decode_cmode(cmode, func_c)
                if (func_c == "") then
                    call warning_msg("invalid")
                    return
                endif
            endif
        elseif (present(func)) then
            if (func == "") then
                call warning_msg("null")
                return
            else
                func_c = func
            endif
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

    subroutine release_zetavars_output(cmode, func)
        implicit none
        integer, intent(in), optional :: cmode
        character(len=*), intent(in), optional :: func
        character(len=10) :: func_c

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Initializing output variables'
        if (present(cmode)) then
            if (cmode == -1) then
                call warning_msg("null")
                return
            else
                call decode_cmode(cmode, func_c)
                if (func_c == "") then
                    call warning_msg("invalid")
                    return
                endif
            endif
        elseif (present(func)) then
            if (func == "") then
                call warning_msg("null")
                return
            else
                func_c = func
            endif
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
            deallocate(advu, advv, advw, advVx, advVy, advVz, err_nldecomp, rr_rev, rr_cha)
        endif
    endsubroutine

    subroutine decode_cmode(cmode, func)
        implicit none
        integer, intent(in) :: cmode
        character(len=*), intent(inout) :: func
        select case (cmode)
            case (0)
                func = "c"
            case (1)
                func = "came"
            case (2, -2)
                func = "cdme"
            case (3, -3)
                func = "dm"
            case (4)
                func = "am"
            case (5)
                func = "f"
            case default
                func = ""
        endselect
    endsubroutine

    subroutine warning_msg(incident)
        character(len=*), intent(in) :: incident
        select case (trim(incident))
            case("null")
                write(*, '(A)') "cmode = -1 or func = "". Do nothing."
            case("invalid")
                write(*, '(A)') "Invalid cmode or func. Exit."
            case("noinput")
                write(*, '(A)') "Neither cmode or func is given. Exit."
        endselect
    endsubroutine
endmodule
