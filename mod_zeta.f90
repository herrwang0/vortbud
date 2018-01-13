module zeta
  use params, only : kd_r, nx, ny, nz, xi_dp, yi_dp, zi_dpst, zi_dped, MVALUE, &
                     xl_reg, yd_reg, fn_grid, fn_cons, fn_dz
  implicit none
  private
  public :: load_const, zeta_equation, decomp_curlnonl

  ! Constants and grid info
  real(kind=kd_r), dimension(:, :)   , allocatable, public :: ulat, tlat, ulong, tlong, dxu, dyu, tarea
  real(kind=kd_r), dimension(:)      , allocatable, public :: z_t
  real(kind=kd_r), dimension(:, :)   , allocatable, public :: fcor, fcort
  real(kind=kd_r), public :: grav

  real(kind=kd_r), dimension(:, :)   , allocatable, public :: dxue, dyue, tareae, uarea, &
                                                              dxun, dyun, tarean, huw , hus
  real(kind=kd_r), dimension(:, :, :), allocatable, public :: dzu, dzt
  logical, dimension(:, :, :), allocatable, public :: umask, tmask

  ! Work variables for zeta_equation and decomp_curlnonl
  real(kind=kd_r), dimension(: ,:, :), allocatable, public :: uc, vc, wc
  real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
      ssh, advx, advy, gradx, grady, hdiffx, hdiffy, vdiffx, vdiffy

  ! Output variables
  real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
      curlnonl, betav, stretchp, err_cor, curlpgrad, curlhdiff, curlvdiff, res
  real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
      advu, advv, advw, advVx, advVy, advVz, curlmet, err_nlsub, err_nldecomp
  real(kind=kd_r), dimension(:, :, :), allocatable, public :: &
      rr_rev, rr_cha

contains
subroutine load_const()
  use ncio, only : nc_read
  implicit none
  real(kind=kd_r) :: omega
  real(kind=kd_r) :: htn(nx, ny), hte(nx, ny), WORK(nx, ny, nz, 1)
  integer :: ix, iy

  print*, "  "
  print*, "Loading grid info from ", trim(fn_grid)
  print*, "  and ", trim(fn_cons)
  print*, "  and ", trim(fn_dz)

  allocate(tlat(nx, ny), tlong(nx, ny), z_t(nz))
  allocate(ulat(nx, ny), ulong(nx, ny), dxu(nx, ny), dyu(nx, ny), tarea(nx, ny))
  allocate(dzt(nx, ny, nz), dzu(nx, ny, nz))

  call nc_read(fn_grid, 'TLAT' , WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  tlat  = WORK(:, :, 1, 1)
  call nc_read(fn_grid, 'TLONG', WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  tlong = WORK(:, :, 1, 1)
  call nc_read(fn_grid, 'ULAT' , WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  ulat  = WORK(:, :, 1, 1)
  call nc_read(fn_grid, 'ULONG', WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  ulong = WORK(:, :, 1, 1)
  call nc_read(fn_grid, 'TAREA', WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  tarea = WORK(:, :, 1, 1)
  call nc_read(fn_grid, 'DXU'  , WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  dxu   = WORK(:, :, 1, 1)
  call nc_read(fn_grid, 'DYU'  , WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  dyu   = WORK(:, :, 1, 1)

  call nc_read(fn_grid, 'z_t'  , WORK, (/1/), (/nz/))
  z_t   = WORK(1:nz, 1, 1, 1)

  call nc_read(fn_cons, 'omega', WORK, (/1/), (/1/))
  omega = WORK(1, 1, 1, 1)

  call nc_read(fn_cons, 'grav' , WORK, (/1/), (/1/))
  grav = WORK(1, 1, 1, 1)

  call nc_read(fn_dz, 'DZT', WORK, (/xl_reg, yd_reg, 1/), (/nx, ny, nz/))
  dzt = WORK(:, :, :, 1)
  call nc_read(fn_dz, 'DZU', WORK, (/xl_reg, yd_reg, 1/), (/nx, ny, nz/))
  dzu = WORK(:, :, :, 1)

  allocate(fcor(nx, ny), fcort(nx, ny))
  do iy = 1, ny
    do ix = 1, nx
      fcor (ix, iy) = 2 * omega * sind(ulat(ix, iy))
      fcort(ix, iy) = 2 * omega * sind(tlat(ix, iy))
    enddo
  enddo


  allocate(tmask(nx, ny, nz), umask(nx, ny, nz))
  tmask = .false.
  umask = .false.
  where (abs(dzt) < 1e-10) tmask = .true.
  where (abs(dzu) < 1e-10) umask = .true.

  print*, "  "
  print*, "Calculating weight functions for decomposing nonlinear term"
  print*, "  Loading grid info from ", trim(fn_grid)
  print*, "  and ", trim(fn_dz)

  allocate(uarea(nx, ny), huw(nx, ny), hus(nx, ny))
  allocate(dxue(nx, ny), dyue(nx, ny), tareae(nx, ny), &
           dxun(nx, ny), dyun(nx, ny), tarean(nx, ny))

  call nc_read(fn_grid, 'HTN', WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  htn = WORK(:, :, 1, 1)
  call nc_read(fn_grid, 'HTE', WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  hte = WORK(:, :, 1, 1)
  call nc_read(fn_grid, 'HUW', WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  huw = WORK(:, :, 1, 1)
  call nc_read(fn_grid, 'HUS', WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  hus = WORK(:, :, 1, 1)
  call nc_read(fn_grid, 'UAREA', WORK, (/xl_reg, yd_reg/), (/nx, ny/))
  uarea = WORK(:, :, 1, 1)

  dyue(1:nx-1, :) = huw(2:nx, :)
  dxue(1:nx-1, :) = htn(2:nx, :)
  tareae = hus * hte

  dyun(:, 1:ny-1) = hte(:, 2:ny)
  dxun(:, 1:ny-1) = hus(:, 2:ny)
  tarean = htn * huw
endsubroutine

subroutine zeta_equation()
  implicit none
  real(kind=kd_r), dimension(:, :, :), allocatable :: curlcor, corx, cory

  print*, '  '
  print*, 'Start calculating zeta equation (zcurl of momentum equations): '

  ! calculating Coriolis terms with u and v
  call calc_cor()
  print*, 'Corx', corx(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'Cory', cory(xi_dp, yi_dp, zi_dpst:zi_dped)

  ! calculating and combining barotropic and baroclinic pressure gradient
  call calc_pgrad()

  ! calculating curl of nonlinear, cor, pgrad, h/v diffusion and res terms
  call calc_curlmom()
  print*, 'curlnonl' , curlnonl (xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'curlcor'  , curlcor  (xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'curlpgrad', curlpgrad(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'curlhdiff', curlhdiff(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'curlvdiff', curlvdiff(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'res'      , res      (xi_dp, yi_dp, zi_dpst:zi_dped)

  ! Decompositing curl of Coriolis term
  call decomp_curlcor()

  ! Verifying decomposition and clean up
  call verify()

  contains
  subroutine calc_cor()
    implicit none
    integer :: iz

    print*, '  '
    print*, '  Calculating Coriolis terms in momentum equations'
    ! print*, "f", fcor(xi_dp, yi_dp), fcort(xi_dp, yi_dp)

    allocate(corx(nx, ny, nz), cory(nx, ny, nz))
    do iz = 1, nz
      corx(:, :, iz) = - fcor * vc(:, :, iz)
      cory(:, :, iz) =   fcor * uc(:, :, iz)
    enddo
  endsubroutine

  subroutine calc_pgrad()
    implicit none
    integer :: iz
    real(kind=kd_r), dimension(nx, ny) :: pgradsfx, pgradsfy

    print*, '  '
    print*, '  Calculating and combining barotropic and baroclinic pressure gradient'

    pgradsfx = 0.
    pgradsfy = 0.
    pgradsfx(1:nx-1, 1:ny-1) = -grav * ((ssh(2:nx,   1:ny-1, 1) + ssh(2:nx,   2:ny  , 1))/2 - &
                                        (ssh(1:nx-1, 1:ny-1, 1) + ssh(1:nx-1, 2:ny  , 1))/2) / dxu(1:nx-1, 1:ny-1)
    pgradsfy(1:nx-1, 1:ny-1) = -grav * ((ssh(1:nx-1, 2:ny  , 1) + ssh(2:nx,   2:ny  , 1))/2 - &
                                        (ssh(1:nx-1, 1:ny-1, 1) + ssh(2:nx,   1:ny-1, 1))/2) / dyu(1:nx-1, 1:ny-1)

    do iz = 1, nz
      gradx(:, :, iz) = pgradsfx - gradx(:, :, iz)
      grady(:, :, iz) = pgradsfy - grady(:, :, iz)
    enddo

    where(umask) gradx = 0.
    where(umask) grady = 0.

    print*, 'pgradx', gradx(xi_dp, yi_dp, zi_dpst:zi_dped)
    print*, 'pgrady', grady(xi_dp, yi_dp, zi_dpst:zi_dped)
  endsubroutine

  subroutine calc_curlmom()
    use popfun, only : zcurl
    implicit none
    integer :: iz

    print*, '  '
    print*, '  Calculating curl of nonlinear, Coriolis, pgrad, horizontal/vertical diffusion and residual terms'
    allocate(curlcor(nx, ny, nz))
    curlcor = 0.
    do iz = 1, nz
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
  endsubroutine

  subroutine decomp_curlcor()
    use popfun, only : u2t
    use derives, only : ddx_chain, ddy_chain
    implicit none
    integer :: iz
    real(kind=kd_r), dimension(nx, ny) :: stretchpx, betavx, stretchpy, betavy, &
                                          WORKrr, WORKrrcx, WORKrrcy, u0, v0, ONES

    print*, '  '
    print*, '  Decomposing curl of Coriolis term'

    ONES = 1.
    do iz = 1, nz
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
  endsubroutine

  subroutine verify()
    implicit none
    integer :: iz

    print*, '  '
    print*, '  Clean up and verifying decomposition of curl(-fu, fv) ...'
    do iz = zi_dpst, zi_dped
      print*, 'Curlcor: '       , curlcor  (xi_dp, yi_dp, iz)
      print*, 'betav+stretchp: ', betav(xi_dp, yi_dp, iz) + stretchp(xi_dp, yi_dp, iz)
      print*, 'errcor: '        , err_cor(xi_dp, yi_dp, iz)
      print*, 'Diff: '          , curlcor (xi_dp, yi_dp, iz) - &
                                  betav   (xi_dp, yi_dp, iz) - &
                                  stretchp(xi_dp, yi_dp, iz) - &
                                  err_cor (xi_dp, yi_dp, iz)
      print*, ''
    enddo
    deallocate(curlcor)
  endsubroutine
endsubroutine

subroutine decomp_curlnonl()
  implicit none
  ! real(kind=kd_r), dimension(:,:,:), allocatable, public :: rrxx, rrxy, rryx, rryy, rrzx, rrzy
  ! real(kind=kd_r), dimension(:,:,:), allocatable, public :: rcuv, rcuu, rcvv, rcvu, rcwv, rcwu
  ! real(kind=kd_r), dimension(:,:,:), allocatable, public :: rsx1, rsx2, rsy1, rsy2, rsz1, rsz2
  ! real(kind=kd_r), dimension(:,:,:), allocatable, public :: vDdivDx, uDdivDy

  real(kind=kd_r), dimension(:,:,:), allocatable :: curladv
  real(kind=kd_r), dimension(:,:,:), allocatable :: ue, vn, wt, ume, vme, umn, vmn, umt, vmt

  print*, '  '
  print*, 'Start decomposing nonlinear term: '

  ! Calculating derived velocity at walls
  call calc_velw()
  ! print*, 'ue', ue(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'vn', vn(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'wt', wt(xi_dp, yi_dp, zi_dpst:zi_dped)
  !
  ! print*, 'ume', ume(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'vme', vme(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'umn', umn(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'vmn', vmn(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'umt', umt(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'vmt', vmt(xi_dp, yi_dp, zi_dpst:zi_dped)

  ! Calculating curladv offline
  call calc_curladv_offline()
  print*, 'curladv', curladv(xi_dp, yi_dp, zi_dpst:zi_dped)

  ! Calculating curlmet offline
  call calc_curlmet_offline()
  print*, 'curlmet', curlmet(xi_dp, yi_dp, zi_dpst:zi_dped)

  ! Calculating error from calculating nonlinear terms offline
  call calc_errnl_sub()
  if (allocated(err_nlsub)) then
    print*, 'err_nlsub', err_nlsub(xi_dp, yi_dp, zi_dpst:zi_dped)
  endif

  ! Calculating error from reversing operators
  call calc_errnl_rev()
  print*, 'rr_rev', rr_rev(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rrxx', rrxx(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rrxy', rrxy(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rryx', rryx(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rryy', rryy(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rrzx', rrzx(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rrzy', rrzy(xi_dp, yi_dp, zi_dpst:zi_dped)

  call calc_decomp()
  print*, 'advu', advu(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'advv', advv(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'advw', advw(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'rr_cha', rr_cha(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rcuv', rcuv(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rcuu', rcuu(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rcvv', rcvv(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rcvu', rcvu(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rcwv', rcwv(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rcwu', rcwu(xi_dp, yi_dp, zi_dpst:zi_dped)

  print*, 'advVx', advVx(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'advVy', advVy(xi_dp, yi_dp, zi_dpst:zi_dped)
  print*, 'advVz', advVz(xi_dp, yi_dp, zi_dpst:zi_dped)

  ! print*, 'vDdivDx', vDdivDx(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'uDdivDy', uDdivDy(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rsx1', rsx1(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rsx2', rsx2(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rsy1', rsy1(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rsy2', rsy2(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rsz1', rsz1(xi_dp, yi_dp, zi_dpst:zi_dped)
  ! print*, 'rsz2', rsz2(xi_dp, yi_dp, zi_dpst:zi_dped)

  ! Verifying decomposition and clean up
  call verify()

  contains
  subroutine calc_velw()
    use popfun, only : t2u
    implicit none
    integer :: iz

    print*, '  '
    print*, '  Calculating velocity at the walls'
    allocate(ue (nx, ny, nz), vn (nx, ny, nz), wt (nx, ny, nz), &
             ume(nx, ny, nz), vme(nx, ny, nz), &
             umn(nx, ny, nz), vmn(nx, ny, nz), &
             umt(nx, ny, nz), vmt(nx, ny, nz))
    ue  = MVALUE
    vn  = MVALUE
    wt  = MVALUE
    ume = MVALUE
    vme = MVALUE
    umn = MVALUE
    vmn = MVALUE
    umt = MVALUE
    vmt = MVALUE

    do iz = 1, nz
      ue(1:nx-1, 2:ny-1, iz) = &
          (0.25  * (uc(1:nx-1, 2:ny-1, iz) * dyu(1:nx-1, 2:ny-1) &
                                           * dzu(1:nx-1, 2:ny-1, iz) + &
                    uc(2:nx  , 2:ny-1, iz) * dyu(2:nx  , 2:ny-1) &
                                           * dzu(2:nx  , 2:ny-1, iz)) + &
           0.125 * (uc(2:nx  , 1:ny-2, iz) * dyu(2:nx  , 1:ny-2) &
                                           * dzu(2:nx  , 1:ny-2, iz) + &
                    uc(1:nx-1, 1:ny-2, iz) * dyu(1:nx-1, 1:ny-2) &
                                           * dzu(1:nx-1, 1:ny-2, iz) + &
                    uc(2:nx  , 3:ny  , iz) * dyu(2:nx  , 3:ny  ) &
                                           * dzu(2:nx  , 3:ny  , iz) + &
                    uc(1:nx-1, 3:ny  , iz) * dyu(1:nx-1, 3:ny  ) &
                                           * dzu(1:nx-1, 3:ny  , iz)))

      vn(2:nx-1, 1:ny-1, iz) = &
          (0.25  * (vc(2:nx-1, 1:ny-1, iz) * dxu(2:nx-1, 1:ny-1) &
                                           * dzu(2:nx-1, 1:ny-1, iz) + &
                    vc(2:nx-1, 2:ny  , iz) * dxu(2:nx-1, 2:ny  ) &
                                           * dzu(2:nx-1, 2:ny  , iz)) + &
           0.125 * (vc(1:nx-2, 2:ny  , iz) * dxu(1:nx-2, 2:ny  ) &
                                           * dzu(1:nx-2, 2:ny  , iz) + &
                    vc(1:nx-2, 1:ny-1, iz) * dxu(1:nx-2, 1:ny-1) &
                                           * dzu(1:nx-2, 1:ny-1, iz) + &
                    vc(3:nx  , 2:ny  , iz) * dxu(3:nx  , 2:ny  ) &
                                           * dzu(3:nx  , 2:ny  , iz) + &
                    vc(3:nx  , 1:ny-1, iz) * dxu(3:nx  , 1:ny-1) &
                                           * dzu(3:nx  , 1:ny-1, iz)))

      if (iz == 1) then
        wt(:, :, iz) = t2u(wc(:, :, 1), tarea, uarea)
      else
        wt(2:nx-1, 2:ny-1, iz) = wt(2:nx-1, 2:ny-1, iz-1) + &
          ((vn(2:nx-1, 2:ny-1, iz-1) - vn(2:nx-1, 1:ny-2, iz-1)) + &
           (ue(2:nx-1, 2:ny-1, iz-1) - ue(1:nx-2, 2:ny-1, iz-1))) &
          / uarea(2:nx-1, 2:ny-1)
      endif
    enddo

    ume(1:nx-1, :, :) = (uc(1:nx-1, :, :) + uc(2:nx, :, :)) / 2.
    vme(1:nx-1, :, :) = (vc(1:nx-1, :, :) + vc(2:nx, :, :)) / 2.

    umn(:, 1:ny-1, :) = (uc(:, 1:ny-1, :) + uc(:, 2:ny, :)) / 2.
    vmn(:, 1:ny-1, :) = (vc(:, 1:ny-1, :) + vc(:, 2:ny, :)) / 2.

    umt(:, :, 1) = uc(:, :, 1)
    umt(:, :, 2:nz) = (uc(:, :, 1:nz-1) + uc(:, :, 2:nz)) / 2.
    vmt(:, :, 1) = vc(:, :, 1)
    vmt(:, :, 2:nz) = (vc(:, :, 1:nz-1) + vc(:, :, 2:nz)) / 2.
  endsubroutine

  subroutine calc_curladv_offline()
    use popfun, only : zcurl
    implicit none
    real(kind=kd_r), dimension(nx, ny, nz) :: advx, advy
    integer :: iz

    print*, '  '
    print*, '  Calculating curl of advection terms (offline)'
    allocate(curladv(nx, ny, nz))
    curladv = 0.
    advx = 0.
    advy = 0.

    do iz = 1, nz
      advx(2:nx, 2:ny, iz) = (ue(2:nx  , 2:ny  , iz) * ume (2:nx  , 2:ny  , iz) - &
                              ue(1:nx-1, 2:ny  , iz) * ume (1:nx-1, 2:ny  , iz) + &
                              vn(2:nx  , 2:ny  , iz) * umn (2:nx  , 2:ny  , iz) - &
                              vn(2:nx  , 1:ny-1, iz) * umn (2:nx  , 1:ny-1, iz))  &
                             / dzu(2:nx, 2:ny, iz) / uarea(2:nx, 2:ny)            &
                             + wt(2:ny, 2:ny, iz) * umt(2:nx, 2:ny, iz) / dzu(2:nx, 2:ny, iz)

      advy(2:nx, 2:ny, iz) = (ue(2:nx  , 2:ny  , iz) * vme (2:nx  , 2:ny  , iz) - &
                              ue(1:nx-1, 2:ny  , iz) * vme (1:nx-1, 2:ny  , iz) + &
                              vn(2:nx  , 2:ny  , iz) * vmn (2:nx  , 2:ny  , iz) - &
                              vn(2:nx  , 1:ny-1, iz) * vmn (2:nx  , 1:ny-1, iz))  &
                             / dzu(2:nx, 2:ny, iz) / uarea(2:nx, 2:ny)            &
                             + wt(2:ny, 2:ny, iz) * vmt(2:nx, 2:ny, iz) / dzu(2:nx, 2:ny, iz)

      if (iz < nz) then
        advx(2:nx, 2:ny, iz) = advx(2:nx, 2:ny, iz) &
                        - wt(2:ny, 2:ny, iz+1) * umt(2:nx, 2:ny, iz+1) / dzu(2:nx, 2:ny, iz)
        advy(2:nx, 2:ny, iz) = advy(2:nx, 2:ny, iz) &
                        - wt(2:ny, 2:ny, iz+1) * vmt(2:nx, 2:ny, iz+1) / dzu(2:nx, 2:ny, iz)
      endif

      advx(1, :, iz) = MVALUE
      advx(:, 1, iz) = MVALUE
      advy(1, :, iz) = MVALUE
      advy(:, 1, iz) = MVALUE

      where(umask(:,:,iz)) ! For cases where dzu = 0.
        advx(:, :, iz) = 0.
        advy(:, :, iz) = 0.
      endwhere

      curladv(:, :, iz) = zcurl(-advx(:, :, iz), -advy(:, :, iz), dxu, dyu, tarea)
      curladv(1, :, iz) = MVALUE
      curladv(1, :, iz) = MVALUE
    enddo
    where(tmask) curladv = MVALUE
  endsubroutine

  subroutine calc_curlmet_offline()
    use popfun, only : zcurl
    implicit none
    real(kind=kd_r), dimension(nx, ny) :: kxu, kyu, metx, mety
    integer :: iz

    print*, '  '
    print*, '  Calculating curl of metric terms (offline)'
    kxu = MVALUE
    kyu = MVALUE

    kxu(2:nx, :) = (huw(2:nx, :) - huw(1:nx-1, :)) / uarea(2:nx, :)
    kyu(:, 2:ny) = (hus(:, 2:ny) - hus(:, 1:ny-1)) / uarea(:, 2:ny)

    do iz = 1, nz
      metx = uc(:, :, iz) * vc(:, :, iz) * kyu - vc(:, :, iz) * vc(:, :, iz) * kxu
      mety = uc(:, :, iz) * vc(:, :, iz) * kxu - uc(:, :, iz) * uc(:, :, iz) * kyu

      curlmet(:, :, iz) = zcurl(-metx, -mety, dxu, dyu, tarea)
      curlmet(1, :, iz) = MVALUE
      curlmet(:, 1, iz) = MVALUE
    enddo
    where(tmask) curlmet = MVALUE
  endsubroutine

  subroutine calc_errnl_sub()
    print*, '  '
    print*, '  Calculating error term from offline calculation'

    if (allocated(curlnonl) .and. allocated(err_nlsub)) then
      err_nlsub = curlnonl - curladv - curlmet
    endif
  endsubroutine

  !! Efficiency of this function can be improved by direct calculation. See doc for derivation.
  subroutine calc_errnl_rev()
    use derives
    implicit none
    integer :: iz
    real(kind=kd_r) :: WORK(nx, ny), WORKB(nx, ny), WORKTx(nx, ny), WORKTy(nx, ny), advthencur(nx, ny), curthenadv(nx, ny)
    ! allocate(rrxx(nx, ny, nz), rrxy(nx, ny, nz), rryx(nx, ny, nz), rryy(nx, ny, nz), &
    !          rrzx(nx, ny, nz), rrzy(nx, ny, nz))

    print*, '  '
    print*, '  Calculating error term from reversing operators'

    do iz = 1, nz
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
      if (iz == nz) then
        WORK = (wt(:,:,iz) * vmt(:,:,iz)) / dzu(:,:,iz)
      else
        WORK = (wt(:,:,iz) * vmt(:,:,iz) - wt(:,:,iz+1) * vmt(:,:,iz+1)) / dzu(:,:,iz)
      endif
      where(umask(:,:,iz)) WORK = 0.
      call ddx(WORK, dyu, tarea, advthencur)

      if (iz == 1) then
        call ddx(wt(:,:,iz) * vmt(:,:,iz), dyu, tarea, WORKTx)
      endif
      if (iz == nz) then
        WORKB = 0.
      else
        call ddx(wt(:,:,iz+1) * vmt(:,:,iz+1), dyu, tarea, WORKB)
      endif
      curthenadv = (WORKTx - WORKB) / dzt(:,:,iz)
      WORKTx = WORKB

      ! rrzx(:, :, iz) = advthencur - curthenadv
      rr_rev(:, :, iz) = rr_rev(:, :, iz) + (advthencur - curthenadv)

      ! d(dwu/dz) / dy
      if (iz == nz) then
        WORK = (wt(:,:,iz) * umt(:,:,iz)) / dzu(:,:,iz)
      else
        WORK = (wt(:,:,iz) * umt(:,:,iz) - wt(:,:,iz+1) * umt(:,:,iz+1)) / dzu(:,:,iz)
      endif
      where(umask(:,:,iz)) WORK = 0.
      call ddy(WORK, dxu, tarea, advthencur)

      if (iz == 1) then
        call ddy(wt(:,:,iz) * umt(:,:,iz), dxu, tarea, WORKTy)
      endif
      if (iz == nz) then
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
  endsubroutine

  subroutine calc_decomp()
    use derives
    use popfun, only : u2t
    implicit none
    integer :: iz
    real(kind=kd_r), dimension(nx, ny) :: ue0, vn0, wt0, ume0, vme0, umn0, vmn0, wm, &
                          WORK1, WORK2, WORK3, WORK1a, WORK2a, WORK3a, DUMMY, ONES
    real(kind=kd_r), dimension(nx, ny) :: DueDx_vme0, DueDy_ume0, DvnDx_vmn0, DvnDy_umn0
    real(kind=kd_r), dimension(nx, ny, 2) :: vmt0, umt0, &
            WORK1zx, WORK1zy, WORK2zx, WORK2zy, WORK3zx, WORK3zy, WORK4zx, WORK4zy
    real(kind=kd_r), dimension(nx, ny) :: uudxdy, vvdxdy, wwdxdy, wvdxdz, wudydz

    print*, '  '
    print*, '  Calculating decomposition ... '
    ONES = 1.

    do iz = 1, nz
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

      if (iz == nz) then
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

      if (iz == nz) then
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

      if (iz == nz) then
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
    advu = -advu
    advv = -advv
    advw = -advw
    rr_cha = -rr_cha

    advVx = -advVx
    advVy = -advVy
    advVz = -advVz

    where(tmask)
      advu  = MVALUE
      advv  = MVALUE
      advw  = MVALUE
      rr_cha = MVALUE
      advVx = MVALUE
      advVy = MVALUE
      advVz = MVALUE
    endwhere
  endsubroutine

  subroutine verify()
    implicit none
    real(kind=kd_r), dimension(nx, ny, nz) :: calc_total
    integer :: iz

    print*, '  '
    print*, '  Clean up and verifying sum of decomposition ...'
    err_nldecomp = rr_rev + rr_cha

    calc_total = advu + advv + advw + advVx + advVy + advVz
    do iz = zi_dpst, zi_dped
      print*, 'Curladv: ', curladv(xi_dp, yi_dp, iz)
      print*, 'All terms: ', calc_total(xi_dp, yi_dp, iz)
      print*, 'All errors: ', err_nldecomp(xi_dp, yi_dp, iz)
      print*, 'Diff: ', curladv(xi_dp, yi_dp, iz) - calc_total(xi_dp, yi_dp, iz) - err_nldecomp(xi_dp, yi_dp, iz)
      print*, ''
    enddo

    deallocate(curladv)
    deallocate(ue, vn, wt, ume, vme, umn, vmn, umt, vmt)
  endsubroutine
endsubroutine
endmodule
