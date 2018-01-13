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
  do iyr = 1, nyr
    do isec = 1, nsec
      call get_yyyymmdd(yrlist(iyr), sec(isec,1), sec(isec,2), avenm, fn_out_dlm, datenm_output)
      call create_outputfiles(datenm_output)

      ! Initializing input/outut fields used by zeta module'
      call init_zetavars()
      ! Initializing outut fields'
      call init_outputvars()

      do ida = sec(isec,1), sec(isec,2)
        call loadave_vars(yrlist(iyr), ida, sec(isec,2)-sec(isec,1)+1)
      enddo

      call decomp_curlnonl()

      call write_outputfiles()
      call release_vars()
      call close_outputfiles()
      print*, "Finished ", yrlist(iyr), 'sec', isec, 'of ', nsec
    enddo
  enddo
endprogram
