Program main
  use params
  use control_eddy
  use zeta, only : load_const, zeta_equation, decomp_curlnonl

  implicit none
  integer :: iyr, isec, ida
  character(len = 15) :: datenm_output

!---------------------Initialization--------------------------------------------
  print*, 'Loading parameters ...'
  call load_params()
  print*, ' '

  print*, 'Loading grids and constants ...'
  call load_const()
  print*, ' '
!---------------------Main loop-------------------------------------------------
  print*, 'Starting looping ...'

  ! do isec = 1, nsec
  !   if individual
  !     do iyr = 1, nyr
  !       do ida = sec(isec,1), sec(isec,2)
  !         call loadave_vars(yrlist(iyr), ida, sec(isec,2)-sec(isec,1)+1)
  !       enddo
  !       call decomp_curlnonl()
  !       call write_outputfiles()
  !       call release_vars()
  !       call close_outputfiles()
  !     enddo
  !   else
  !     do iyr = 1, nyr
  !       do ida = sec(isec,1), sec(isec,2)
  !         call loadave_vars(yrlist(iyr), ida, sec(isec,2)-sec(isec,1)+1)
  !       enddo
  !     enddo
  !     call decomp_curlnonl()
  !     call write_outputfiles()
  !     call release_vars()
  !     call close_outputfiles()
  !   endif
  ! enddo

    do isec = 1, nsec
      if initidual
      do iyr = 1, nyr
        call init_zetavars()
        do ida = sec(isec,1), sec(isec,2)
          call loadave_vels(yrlist(iyr), ida, sec(isec,2)-sec(isec,1)+1)
        enddo

        call get_outputfn_mean()
        call create_outputfiles(fn, fn)
        call decomp_curlnonl()
        call write_outputfiles()
        call release_zetavars()
        call close_outputfiles()
      enddo

      else
      call init_zetavars()
      do iyr = 1, nyr
        do ida = sec(isec,1), sec(isec,2)
          call loadave_vels(yrlist(iyr), ida, (sec(isec,2)-sec(isec,1)+1) * nyr)
        enddo
      enddo

      call create_outputfiles(datenm_output)
      call decomp_curlnonl()
      call write_outputfiles()
      call release_vars()
      call close_outputfiles()
      endif
    enddo
