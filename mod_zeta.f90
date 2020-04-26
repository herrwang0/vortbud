module zeta
    use params, only : kd_r, MVALUE, B, fngrid, fmt_exp, fmt_flt
    implicit none
    private
    public :: zetavar, vargrp, load_const, calc_zeta, &
        init_zetavars_input, init_zetavars_output, release_zetavars_input, release_zetavars_output

    ! Output variable type w/ meta for netCDF
    type zetavar
        character(len = 10) :: name
        real(kind=kd_r), dimension(:,:,:), allocatable :: value
        real(kind=kd_r) :: missing_value
        character(len = 300) :: long_name
        character(len = 10)  :: Units = "1/s^2"
        character(len = 300) :: coordinates = "TLONG TLAT z_t time" 
        integer :: varid 
    endtype zetavar

    ! Wrapper for variables groups
    type vargrp
        character(len = 10) :: name
        type(zetavar), dimension(:), pointer :: vlist
        logical :: key = .False.
    endtype vargrp

    ! Input variables
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: uc, vc, wc
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
        ssh, advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy
    real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
        ueu, uev, vnu, vnv, wtu, wtv

    ! Output variables
    ! Variables are gathered in groups corresponding to each individual calc mode
    ! vgrp is a wrapper over the groups, which includes the group's name and if it is activated 
    type(zetavar), dimension(:), target :: &
        vl_c(8), vl_a(1), vl_aD(3), vl_m(1), vl_d(7), vl_dD(18), vl_e(1), vl_f(1)
    ! Hash table: 1 = c, 2 = a, 3 = a-, 4 = m, 5 = d, 6 = d-, 7 = e, 8 = f
    type(vargrp), dimension(8), target, public :: vgrp

    real(kind=kd_r), dimension(:, :, :), pointer :: &
        curlnonl, betav, stretchp, err_cor, curlpgrad, curlhdiff, curlvdiff, res
    real(kind=kd_r), dimension(:, :, :), pointer :: &
        curladv, curlmet, err_nlsub, advu, advv, advw, twix, twiy, twiz, err_adv, curladvf
    real(kind=kd_r), dimension(:, :, :), pointer :: curladvu, curladvv, curladvw
    real(kind=kd_r), dimension(:, :, :), pointer :: &
        advu_x, advv_x, advw_x, twix_x, twiy_x, twiz_x, erax_x, eray_x, eraz_x, &
        advu_y, advv_y, advw_y, twix_y, twiy_y, twiz_y, erax_y, eray_y, eraz_y

    ! Constants and grid info
    real(kind=kd_r), allocatable, public :: tlat(:,:), tlong(:,:), z_t(:)
    real(kind=kd_r), dimension(:, :)   , allocatable :: ulat, ulong, dxu, dyu, tarea, uarea, huw , hus
    real(kind=kd_r), dimension(:)      , allocatable :: z_w, dz
    real(kind=kd_r), dimension(:, :, :), allocatable :: dzu, dzt
    real(kind=kd_r), dimension(:, :)   , allocatable :: fcor, fcort
    real(kind=kd_r) :: grav
    logical, dimension(:, :, :), allocatable :: umask, tmask
    real(kind=kd_r), dimension(:,:,:), allocatable :: ue, vn, wt, ume, vme, umn, vmn, umt, vmt

    ! print format
    character(len = 100) :: fmts_vel, fmtm_vel, fmts_vor, fmtm_vor

    contains
    !!--------------------------------------------------------------------------
    ! * Input functions
    ! curl ("c"): curl of zeta equations terms, including three terms from Coriolis
    !   (nonlinear, curl of pgrad, hdiff, vdiff, residual, betav from Coriolis,
    !    stretching from Coriolis, error with Coriolis decomposition)
    ! offline adv ("a"): offline calculation of nonlinear advection term (curladv)
    !   ("a-"): Debug mode, outputting additionally curladvu, curladvv, curladvw
    ! offline met ("m"): offline calculation of nonlinear metric term (curlmet)
    ! decomp adv ("d"): decomposition of nonlinear advection term.
    !   ("d"): flux form twisting term (default)
    !   ("d#"): traditional twisting term
    !   ("d-"/"d#-"): Debug mode, outputting the original x and y components of each term
    ! online adv from flux ("f"): "online" calculation of the advection term from momentum fluxes
    ! + err_nlsub ("e"): difference between curlnonl online and offline 
    !   (should not be used as a standalone function)
    !!--------------------------------------------------------------------------

    subroutine calc_zeta(func)
        implicit none
        character(len=*), intent(in) :: func
        integer :: ig

        write(fmts_vor, '(A, A, A)' ) '(A20, ', trim(fmt_exp), ')'
        write(fmtm_vor, '(A, I2, A, A)') '(A20, ', B%zi_dped - B%zi_dpst + 1, trim(fmt_exp), ')'
        write(fmts_vel, '(A, A, A)' ) '(A20, ', trim(fmt_flt), ')'
        write(fmtm_vel, '(A, I2, A, A)') '(A20, ', B%zi_dped - B%zi_dpst + 1, trim(fmt_flt), ')'

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)'), 'Calculating voriticity equation'
        write(*, '(2A)') "  Function code: ", trim(func)

        res       => vl_c(1)%value
        curlnonl  => vl_c(2)%value
        curlpgrad => vl_c(3)%value
        curlhdiff => vl_c(4)%value
        curlvdiff => vl_c(5)%value
        betav     => vl_c(6)%value
        stretchp  => vl_c(7)%value
        err_cor   => vl_c(8)%value
    
        curladv   => vl_a(1)%value
        curladvu  => vl_aD(1)%value; curladvv  => vl_aD(2)%value; curladvw  => vl_aD(3)%value
    
        curlmet   => vl_m(1)%value
    
        advu      => vl_d(1)%value;  advv      => vl_d(2)%value;  advw      => vl_d(3)%value
        twix      => vl_d(4)%value;  twiy      => vl_d(5)%value;  twiz      => vl_d(6)%value
        err_adv   => vl_d(7)%value
    
        advu_x    => vl_dD( 1)%value; advu_y    => vl_dD( 2)%value
        advv_x    => vl_dD( 3)%value; advv_y    => vl_dD( 4)%value
        advw_x    => vl_dD( 5)%value; advw_y    => vl_dD( 6)%value
        twix_x    => vl_dD( 7)%value; twix_y    => vl_dD( 8)%value
        twiy_x    => vl_dD( 9)%value; twiy_y    => vl_dD(10)%value
        twiz_x    => vl_dD(11)%value; twiz_y    => vl_dD(12)%value
        erax_x    => vl_dD(13)%value; erax_y    => vl_dD(14)%value
        eray_x    => vl_dD(15)%value; eray_y    => vl_dD(16)%value
        eraz_x    => vl_dD(17)%value; eraz_y    => vl_dD(18)%value
    
        err_nlsub => vl_e(1)%value
        curladvf  => vl_f(1)%value
    
        if (index(func, "c") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            write(*, '(A)') '  Start calculating zcurl of momentum equations'
            call zeta_equation()
        endif
        if (index(func, "m") /= 0) then
            call calc_curlmet()
        endif
        if (index(func, "a") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            if (index(func, "a-") /= 0) then
                write(*, '(A)') '  Calculating curl of advection term (offline, w/ u,v,w commponents)'
                call calc_curladv(.True.)
            else
                write(*, '(A)') '  Calculating curl of advection term (offline)'
                call calc_curladv(.False.)
            endif
        endif
        if (index(func, "d") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'

            if   (index(func, "d#-") /= 0) then
                write(*, '(A)') 'ERROR: "d#-" found in calc mode, replace it with "d-#"'
                stop
            endif
            if     (index(func, "d-#") /= 0) then
                write(*, '(A)') '  Decomposing curl of advection term (w/ nonflux twisting term, w/ subcomponents)'
                call decomp_curladv(.False., .True.)
            elseif (index(func, "d-") /= 0 .and. index(func, "d-#") == 0) then
                write(*, '(A)') '  Decomposing curl of advection term (w/ flux twisting term, w/ subcomponents)'
                call decomp_curladv(.True., .True.)
            elseif (index(func, "d#") /= 0) then
                write(*, '(A)') '  Decomposing curl of advection term (w/ nonflux twisting term)'
                call decomp_curladv(.False., .False.)
            else
                write(*, '(A)') '  Decomposing curl of advection term (w/ flux twisting term)'
                call decomp_curladv(.True., .False.)
            endif
        endif
        if (index(func, "e") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            write(*, '(A)') '  Calculating error term from offline calculation of nonlinear terms'
            call calc_errnlsub()
        endif
        if (index(func, "f") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            write(*, '(A)') '  Calculating curl of advection terms from momentum fluxes'
            call calc_curladv_flx()
        endif

        if (index(func, "a") /= 0 .or. index(func, "d") /= 0) then
            write(*, *)
            write(*, '(A)') '  ---------------------------------------------------'
            write(*, '(A)') '  Verifying offline calculations'
            call verify_nonl()
        endif

        nullify(curlnonl, betav, stretchp, err_cor, curlpgrad, curlhdiff, curlvdiff, res, &
                curladv, curlmet, err_nlsub, advu, advv, advw, twix, twiy, twiz, err_adv, curladvf, &
                curladvu, curladvv, curladvw, &
                advu_x, advv_x, advw_x, twix_x, twiy_x, twiz_x, erax_x, eray_x, eraz_x, &
                advu_y, advv_y, advw_y, twix_y, twiy_y, twiz_y, erax_y, eray_y, eraz_y)
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

        allocate(uarea(B%nx, B%ny), huw(B%nx, B%ny), hus(B%nx, B%ny))
        call nc_read(fngrid%grid, 'HUW', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        huw = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'HUS', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        hus = WORK(:, :, 1, 1)
        call nc_read(fngrid%grid, 'UAREA', WORK, (/B%xl_reg, B%yd_reg/), (/B%nx, B%ny/))
        uarea = WORK(:, :, 1, 1)
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
        curlnonl(1,:,:) = MVALUE; curlcor(1,:,:) = MVALUE; curlpgrad(1,:,:) = MVALUE; 
        curlhdiff(1,:,:) = MVALUE; curlvdiff(1,:,:) = MVALUE; res(1,:,:) = MVALUE;
        curlnonl(:,1,:) = MVALUE; curlcor(:,1,:) = MVALUE; curlpgrad(:,1,:) = MVALUE; 
        curlhdiff(:,1,:) = MVALUE; curlvdiff(:,1,:) = MVALUE; res(:,1,:) = MVALUE;
        curlpgrad(B%nx,:,:) = MVALUE; curlpgrad(:,B%ny,:) = MVALUE; 


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
        betav(1,:,:) = MVALUE; stretchp(1,:,:) = MVALUE; err_cor(1,:,:) = MVALUE;
        betav(:,1,:) = MVALUE; stretchp(:,1,:) = MVALUE; err_cor(:,1,:) = MVALUE; 

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
        if ( vgrp(2)%key ) then
            err_nlsub = curlnonl - curladv - curlmet
        elseif ( vgrp(5)%key ) then
            err_nlsub = curlnonl - advu - advv - advw - twix - twiy - twiz - err_adv - curlmet
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
        real(kind=kd_r), dimension(B%nx, B%ny) :: advxu, advxv, advxw, advyu, advyv, advyw
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

            curladv(:, :, iz) = zcurl(-advxu-advxv-advxw, -advyu-advyv-advyw, &
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
                twix(:, :, iz) = twix(:, :, iz) + WORK
                if (debug) twix_x(:, :, iz) = WORK
            else
                WORK = mean_ys(mean_xw(dd_xw(ue(:,:,iz), ONES)) * mean_xw(dd_xw(vme(:,:,iz), dxu))) / tarea / dzt(:,:,iz)
                twix(:, :, iz) = twix(:, :, iz) + WORK
                if (debug) twix_x(:, :, iz) = WORK
                WORK = mean_ys(mean_xw(mean_xw(vme(:,:,iz))/dxu) * dd_xw(dd_xw(ue(:,:,iz), ONES), ONES)) / tarea / dzt(:,:,iz)
                err_adv(:,:,iz) = err_adv(:,:,iz) + WORK
                if (debug) erax_x(:, :, iz) = WORK                  
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
                twix(:, :, iz) = twix(:, :, iz) - WORK
                if (debug) twix_y(:, :, iz) = -WORK
            else
                WORK = mean_xw(mean_xw(dd_ys(ue(:,:,iz), ONES)) * mean_ys(dd_xw(ume(:,:,iz), dyu))) / tarea / dzt(:,:,iz)
                twix(:, :, iz) = twix(:, :, iz) - WORK
                if (debug) twix_y(:, :, iz) = -WORK
                WORK = mean_xw(mean_ys(mean_xw(ume(:,:,iz))/dyu) * dd_xw(dd_ys(ue(:,:,iz), ONES), ONES)) / tarea / dzt(:,:,iz)
                err_adv(:,:,iz) = err_adv(:,:,iz) - WORK
                if (debug) erax_y(:, :, iz) = -WORK
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
                twiy(:, :, iz) = twiy(:, :, iz) + WORK
                if (debug) twiy_x(:, :, iz) = WORK
            else
                WORK = mean_ys(mean_ys(dd_xw(vn(:,:,iz), ONES)) * mean_xw(dd_ys(vmn(:,:,iz), dxu))) / tarea / dzt(:,:,iz)
                twiy(:, :, iz) = twiy(:, :, iz) + WORK
                if (debug) twiy_x(:, :, iz) = WORK
                WORK = mean_ys(mean_xw(mean_ys(vmn(:,:,iz))/dxu) * dd_ys(dd_xw(vn(:,:,iz), ONES), ONES)) / tarea / dzt(:,:,iz)
                err_adv(:,:,iz) = err_adv(:,:,iz) + WORK
                if (debug) eray_x(:, :, iz) = WORK
            endif

            ! d [d(vu)/dy] / dy = d (vdu/dy) / dy + d (udv/dy) / dy
            u1 = vn(:,:,iz)
            u2 = umn(:,:,iz)/dyu
            where (umask(:,:,iz))
                u1 = 0.; u2 = 0.
            endwhere
            call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2(:,:,2), u20du1(:,:,2))
            u1 = shift_yn(vn(:,:,iz))
            u2 = shift_yn(umn(:,:,iz))/dyu
            where (umask(:,:,iz))
                u1 = 0.; u2 = 0.
            endwhere
            call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2(:,:,1), u20du1(:,:,1))
            WORK = mean_xw(u10du2(:,:,2) - u10du2(:,:,1)) / tarea / dzt(:,:,iz)
            advv(:, :, iz) = advv(:, :, iz) - WORK
            if (debug) advv_y(:, :, iz) = -WORK

            if (twif) then
                WORK = mean_xw(u20du1(:,:,2) - u20du1(:,:,1)) / tarea / dzt(:,:,iz)
                twiy(:, :, iz) = twiy(:, :, iz) - WORK
                if (debug) twiy_y(:, :, iz) = -WORK
            else
                WORK = mean_xw(mean_ys(dd_ys(vn(:,:,iz), ONES)) * mean_ys(dd_ys(umn(:,:,iz), dyu))) / tarea / dzt(:,:,iz)
                twiy(:, :, iz) = twiy(:, :, iz) - WORK
                if (debug) twiy_y(:, :, iz) = -WORK
                WORK = mean_xw(mean_ys(mean_ys(umn(:,:,iz))/dyu) * dd_ys(dd_ys(vn(:,:,iz), ONES), ONES)) / tarea / dzt(:,:,iz)
                err_adv(:,:,iz) = err_adv(:,:,iz) - WORK
                if (debug) eray_y(:, :, iz) = -WORK
            endif

          ! advw  
            ! d(wv) / dx  &  d(wu) / dy
            u1 = wt(:,:,iz)
            u2 = vmt(:,:,iz)*dyu
            where (umask(:,:,iz)) u1 = 0.
            where (umask(:,:,iz)) u2 = 0.
            call dd_xw_chain(u1, mean_xw(u1), u2, mean_xw(u2), ONES, u10du2_zx(:,:,2), u20du1_zx(:,:,2))
            u2 = umt(:,:,iz)*dxu
            where (umask(:,:,iz)) u2 = 0.
            call dd_ys_chain(u1, mean_ys(u1), u2, mean_ys(u2), ONES, u10du2_zy(:,:,2), u20du1_zy(:,:,2))

            if (iz == B%nz) then
                u10du2_zx(:,:,1) = 0.
                u20du1_zx(:,:,1) = 0.
                u10du2_zy(:,:,1) = 0.
                u20du1_zy(:,:,1) = 0.
            else
                u1 = wt(:,:,iz+1)
                u2 = vmt(:,:,iz+1)*dyu
                where (umask(:,:,iz)) u1 = 0. ! We cannnot step down and reuse last layer's bottom interface as the masks differ
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
                twiz(:, :, iz) = twiz(:, :, iz) + WORK
                if (debug) twiz_x(:, :, iz) = WORK

                WORK = mean_xw(u20du1_zy(:,:,2) - u20du1_zy(:,:,1)) / tarea / dzt(:, :, iz)
                twiz(:, :, iz) = twiz(:, :, iz) - WORK
                if (debug) twiz_y(:, :, iz) = -WORK
            else
                if (iz == B%nz) then
                    WORK = mean_ys(mean_xw((vmt(:,:,iz) - 0.) * dyu) * dd_xw((wt(:,:,iz) + 0.)/2, ONES)) / tarea / dzt(:,:,iz) 
                    twiz(:, :, iz) = twiz(:, :, iz) + WORK
                    if (debug) twiz_x(:, :, iz) = WORK

                    WORK = mean_xw(mean_ys((umt(:,:,iz) - 0.) * dxu) * dd_ys((wt(:,:,iz) + 0.)/2, ONES)) / tarea / dzt(:,:,iz) 
                    twiz(:, :, iz) = twiz(:, :, iz) - WORK
                    if (debug) twiz_y(:, :, iz) = -WORK

                    WORK = mean_ys(dd_xw(wt(:,:,iz) - 0., ONES) * mean_xw((vmt(:,:,iz) + 0.)/2 * dyu)) / tarea / dzt(:,:,iz)
                    err_adv(:,:,iz) = err_adv(:,:,iz) + WORK
                    if (debug) eraz_x(:, :, iz) = WORK

                    WORK = mean_xw(dd_ys(wt(:,:,iz) - 0., ONES) * mean_ys((umt(:,:,iz) + 0.)/2 * dxu)) / tarea / dzt(:,:,iz)
                    err_adv(:,:,iz) = err_adv(:,:,iz) - WORK
                    if (debug) eraz_y(:, :, iz) = -WORK
               else
                    WORK = mean_ys(mean_xw((vmt(:,:,iz) - vmt(:,:,iz+1)) * dyu) * dd_xw((wt(:,:,iz) + wt(:,:,iz+1))/2, ONES)) & 
                           / tarea / dzt(:,:,iz) 
                    twiz(:, :, iz) = twiz(:, :, iz) + WORK
                    if (debug) twiz_x(:, :, iz) = WORK

                    WORK = mean_xw(mean_ys((umt(:,:,iz) - umt(:,:,iz+1)) * dxu) * dd_ys((wt(:,:,iz) + wt(:,:,iz+1))/2, ONES)) & 
                           / tarea / dzt(:,:,iz) 
                    twiz(:, :, iz) = twiz(:, :, iz) - WORK
                    if (debug) twiz_y(:, :, iz) = -WORK

                    WORK = mean_ys(dd_xw(wt(:,:,iz) - wt(:,:,iz+1), ONES) * mean_xw((vmt(:,:,iz) + vmt(:,:,iz+1))/2 * dyu)) &
                           / tarea / dzt(:,:,iz)
                    err_adv(:,:,iz) = err_adv(:,:,iz) + WORK
                    if (debug) eraz_x(:, :, iz) = WORK

                    WORK = mean_xw(dd_ys(wt(:,:,iz) - wt(:,:,iz+1), ONES) * mean_ys((umt(:,:,iz) + umt(:,:,iz+1))/2 * dxu)) &
                           / tarea / dzt(:,:,iz)
                    err_adv(:,:,iz) = err_adv(:,:,iz) - WORK
                    if (debug) eraz_y(:, :, iz) = -WORK
                endif
            endif

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

                twix(:, :, iz) = twix(:, :, iz) + F_wvdxdz + F_vvdxdy - F_wwdxdy + F_uudxdy
                twiy(:, :, iz) = twiy(:, :, iz) - F_wudydz - F_vvdxdy + F_wwdxdy - F_uudxdy
                twiz(:, :, iz) = twiz(:, :, iz) + F_wudydz - F_wvdxdz

                ! stepping down in z direction
                F_wdzx(:,:,2) = F_wdzx(:,:,1)
                F_wdzy(:,:,2) = F_wdzy(:,:,1)
            endif
        enddo

        ! RHS
        advu = -advu; advv = -advv; advw = -advw
        twix = -twix; twiy = -twiy; twiz = -twiz
        err_adv = -err_adv

        where(tmask)
            advu = MVALUE; advv = MVALUE; advw = MVALUE
            twix = MVALUE; twiy = MVALUE; twiz = MVALUE
        endwhere
        advu(1:2, : , :) = MVALUE; advv(1:2, :,  :) = MVALUE; advw(1:2, :,  :) = MVALUE
        twix(1:2, : , :) = MVALUE; twiy(1:2, :,  :) = MVALUE; twiz(1:2, :,  :) = MVALUE
        advu( :, 1:2, :) = MVALUE; advv( :, 1:2, :) = MVALUE; advw( :, 1:2, :) = MVALUE
        twix( :, 1:2, :) = MVALUE; twiy( :, 1:2, :) = MVALUE; twiz( :, 1:2, :) = MVALUE
        advu(B%nx, :, :) = MVALUE; advv(B%nx, :, :) = MVALUE; advw(B%nx, :, :) = MVALUE
        twix(B%nx, :, :) = MVALUE; twiy(B%nx, :, :) = MVALUE; twiz(B%nx, :, :) = MVALUE
        advu(:, B%ny, :) = MVALUE; advv(:, B%ny, :) = MVALUE; advw(:, B%ny, :) = MVALUE
        twix(:, B%ny, :) = MVALUE; twiy(:, B%ny, :) = MVALUE; twiz(:, B%ny, :) = MVALUE

        if (twif) then
            where(tmask) err_adv = MVALUE
            err_adv(1:2, :, :) = MVALUE; err_adv(B%nx, :, :) = MVALUE
            err_adv(:, 1:2, :) = MVALUE; err_adv(:, B%nx, :) = MVALUE
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

        write(*, fmtm_vor) 'twix: ', twix(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        if (debug) then
            write(*, fmtm_vor) '  twix_x: ', twix_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  erax_x: ', erax_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  twix_y: ', twix_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  erax_y: ', erax_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
        write(*, fmtm_vor) 'twiy: ', twiy(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        if (debug) then
            write(*, fmtm_vor) '  twiy_x: ', twiy_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  eray_x: ', eray_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  twiy_y: ', twiy_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  eray_y: ', eray_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
        write(*, fmtm_vor) 'twiz: ', twiz(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        if (debug) then
            write(*, fmtm_vor) '  twiz_x: ', twiz_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  eraz_x: ', eraz_x(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  twiz_y: ', twiz_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
            write(*, fmtm_vor) '  eraz_y: ', eraz_y(B%xi_dp, B%yi_dp, B%zi_dpst:B%zi_dped)
        endif
        deallocate(ue, vn, wt, ume, vme, umn, vmn, umt, vmt)
    endsubroutine

    subroutine verify_nonl()
        implicit none
        integer :: iz
        real(kind=kd_r) :: total

        do iz = B%zi_dpst, B%zi_dped
            write(*, *)
            if ( vgrp(5)%key ) then
                total = advu(B%xi_dp, B%yi_dp, iz) + advv(B%xi_dp, B%yi_dp, iz) + advw(B%xi_dp, B%yi_dp, iz) + &
                        twix(B%xi_dp, B%yi_dp, iz) + twiy(B%xi_dp, B%yi_dp, iz) + twiz(B%xi_dp, B%yi_dp, iz)
                write(*, '(A, I02)') '      iz = ',  iz
                write(*, fmts_vor) 'All components: ',  total
                write(*, fmts_vor) 'All errors: ', err_adv(B%xi_dp, B%yi_dp, iz)
                if ( vgrp(2)%key ) then
                    write(*, fmts_vor) 'Diff: ', curladv(B%xi_dp, B%yi_dp, iz) - total - err_adv(B%xi_dp, B%yi_dp, iz)
                endif
            endif

            if ( vgrp(2)%key ) then
                write(*, fmts_vor) 'Curladv: ', curladv(B%xi_dp, B%yi_dp, iz)
            endif

            if ( vgrp(4)%key ) then
                write(*, fmts_vor) 'Curlmet: ', curlmet(B%xi_dp, B%yi_dp, iz)
            endif

            if ( vgrp(7)%key ) then
                write(*, fmts_vor) 'err_nlsub: ', err_nlsub(B%xi_dp, B%yi_dp, iz)
            endif

            if ( vgrp(8)%key ) then
                write(*, fmts_vor) 'Curladv (f): ', curladvf(B%xi_dp, B%yi_dp, iz)
            endif
        enddo
    endsubroutine
    
    subroutine init_zetavars_input(func)
        implicit none
        character(len=*), intent(in) :: func

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Initializing input variables'

        if (index(func, "a") /= 0 .or. index(func, "m") /= 0 .or. &
            index(func, "d") /= 0 .or. index(func, "c") /= 0) then
            allocate(uc(B%nx, B%ny, B%nz), vc(B%nx, B%ny, B%nz))
            uc = 0.; vc = 0.
        endif
        if (index(func, "c") /= 0) then
            allocate(advx  (B%nx, B%ny, B%nz), advy  (B%nx, B%ny, B%nz), &
                     gradx (B%nx, B%ny, B%nz), grady (B%nx, B%ny, B%nz), &
                     hdiffx(B%nx, B%ny, B%nz), hdiffy(B%nx, B%ny, B%nz), &
                     vdiffx(B%nx, B%ny, B%nz), vdiffy(B%nx, B%ny, B%nz), &
                     ssh   (B%nx, B%ny, 1 ))
            advx   = 0.; advy   = 0.; gradx  = 0.; grady  = 0.
            hdiffx = 0.; hdiffy = 0.; vdiffx = 0.; vdiffy = 0.
            ssh = 0.
        endif
        if (index(func, "a") /= 0 .or. index(func, "d") /= 0) then
            allocate(wc(B%nx, B%ny, B%nz))
            wc = 0.
        endif
        if (index(func, "f") /= 0) then
            allocate(ueu(B%nx, B%ny, B%nz), uev(B%nx, B%ny, B%nz), &
                     vnu(B%nx, B%ny, B%nz), vnv(B%nx, B%ny, B%nz), &
                     wtu(B%nx, B%ny, B%nz), wtv(B%nx, B%ny, B%nz))
            ueu = 0.; uev = 0.; vnu = 0.; vnv = 0.; wtu = 0.; wtv = 0.
        endif
    endsubroutine

    subroutine init_zetavars_output(func)
        implicit none
        character(len=*), intent(in) :: func
        integer :: ig, iv 
    
        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Initializing output variables'
    
        ! Output variable list (vl)
        ! Group c: Simple curl of momentum terms
        vl_c(1)%name = "res"      ; vl_c(1)%long_name = "Residual (lhs)"
        vl_c(2)%name = "curlnonl" ; vl_c(2)%long_name = "Curl of nonlinear term (rhs)"
        vl_c(3)%name = "curlpgrad"; vl_c(3)%long_name = "Curl of pressure gradient (rhs)"
        vl_c(4)%name = "curlhdiff"; vl_c(4)%long_name = "Curl of horizontal diffusion (rhs)"
        vl_c(5)%name = "curlvdiff"; vl_c(5)%long_name = "Curl of vertical diffusion (rhs)"
        vl_c(6)%name = "betav"    ; vl_c(6)%long_name = "Planetary vorticity advection (rhs)"
        vl_c(7)%name = "stretchp" ; vl_c(7)%long_name = "Planetary vorticity stretching (rhs)"
        vl_c(8)%name = "errcor"   ; vl_c(8)%long_name = "Error from decomposing curl(-fv, fu) (rhs)"
        ! Group a: Offline curl of advection
        vl_a(1)%name = "curladv"  ; vl_a(1)%long_name = "Curl of advection (offline) (rhs)"
        ! Group aD: For group a debug
        vl_aD(1)%name = "curladvu"; vl_aD(2)%name = "curladvv"; vl_aD(3)%name = "curladvw"; 
        ! Group m: Offline curl of metric
        vl_m(1)%name = "curlmet"  ; vl_m(1)%long_name = "Curl of metric (offline) (rhs)"
        ! Group d: Offline decomposition of curl advection
        vl_d(1)%name = "advu"     ; vl_d(1)%long_name = "Relative vorticity advection zonal (rhs)"
        vl_d(2)%name = "advv"     ; vl_d(2)%long_name = "Relative vorticity advection meridional (rhs)"
        vl_d(3)%name = "advw"     ; vl_d(3)%long_name = "Relative vorticity advection vertical (rhs)"
        vl_d(4)%name = "twix"     ; vl_d(4)%long_name = "Relative vorticity tilting zonal (rhs)"
        vl_d(5)%name = "twiy"     ; vl_d(5)%long_name = "Relative vorticity tilting meridional (rhs)"
        vl_d(6)%name = "twiz"     ; vl_d(6)%long_name = "Relative vorticity stretching (rhs)"
        vl_d(7)%name = "erradv"   ; vl_d(7)%long_name = "Error from decomposing advection (rhs)"
        ! Group dD: For group d debug
        vl_dD( 1)%name = "advu_x"; vl_dD( 2)%name = "advu_y"
        vl_dD( 3)%name = "advv_x"; vl_dD( 4)%name = "advv_y"
        vl_dD( 5)%name = "advw_x"; vl_dD( 6)%name = "advw_y"
        vl_dD( 7)%name = "twix_x"; vl_dD( 8)%name = "twix_y"
        vl_dD( 9)%name = "twiy_x"; vl_dD(10)%name = "twiy_y"
        vl_dD(11)%name = "twiz_x"; vl_dD(12)%name = "twiz_y"
        vl_dD(13)%name = "erax_x"; vl_dD(14)%name = "erax_y"
        vl_dD(15)%name = "eray_x"; vl_dD(16)%name = "eray_y"
        vl_dD(17)%name = "eraz_x"; vl_dD(18)%name = "eraz_y"
        ! Group e: Difference between curlnonl and offline nonlinear
        vl_e(1)%name = "errsub"   ; vl_e(1)%long_name = "Error from offline calculated advection (rhs)"
        ! Group f: Curl of advection using output flues
        vl_f(1)%name = "curladvf" ; vl_f(1)%long_name = "Curl of advection using output fluxes (ueu, uev, vnu, vnv, wtu, wtv) (rhs)"
    
        ! Wrapper group to control on/off of the variable groups
        vgrp(1)%name = "c" ; vgrp(1)%vlist => vl_c
        vgrp(2)%name = "a" ; vgrp(2)%vlist => vl_a
        vgrp(3)%name = "a-"; vgrp(3)%vlist => vl_aD
        vgrp(4)%name = "m" ; vgrp(4)%vlist => vl_m
        vgrp(5)%name = "d" ; vgrp(5)%vlist => vl_d
        vgrp(6)%name = "d-"; vgrp(6)%vlist => vl_dD
        vgrp(7)%name = "e" ; vgrp(7)%vlist => vl_e
        vgrp(8)%name = "f" ; vgrp(8)%vlist => vl_f
    
        do ig = 1, size(vgrp)
            if (index(func, trim(vgrp(ig)%name)) /= 0) then 
                vgrp(ig)%key = .True.
                do iv = 1, size(vgrp(ig)%vlist)
                    allocate(vgrp(ig)%vlist(iv)%value(B%nx, B%ny, B%nz))
                    vgrp(ig)%vlist(iv)%value = 0.
                enddo
            endif
        enddo
    endsubroutine

    subroutine release_zetavars_input(func)
        implicit none
        character(len=*), intent(in) :: func

        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Releasing input variables'

        if (index(func, "a") /= 0 .or. index(func, "m") /= 0 .or. &
            index(func, "d") /= 0 .or. index(func, "c") /= 0) then
            deallocate(uc, vc)
        endif
        if (index(func, "c") /= 0) then
            deallocate(advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy, ssh)
        endif
        if (index(func, "a") /= 0 .or. index(func, "d") /= 0) then
            deallocate(wc)
        endif
        if (index(func, "f") /= 0) then
            deallocate(ueu, uev, vnu, vnv, wtu, wtv)
        endif
    endsubroutine

    subroutine release_zetavars_output()
        implicit none
        integer :: ig, iv
        
        write(*, *)
        write(*, '(A)') '-----------------------------------------------------'
        write(*, '(A)') 'Releasing output variables'

        do ig = 1, size(vgrp)
            if (vgrp(ig)%key) then 
                do iv = 1, size(vgrp(ig)%vlist)
                    deallocate(vgrp(ig)%vlist(iv)%value)
                enddo
            endif
        enddo
    endsubroutine
endmodule
