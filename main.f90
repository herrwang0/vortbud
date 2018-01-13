Program main
  use params
  use control
  use zeta, only : load_const, zeta_equation, decomp_curlnonl

  implicit none
  integer :: iyr, imn, ida
  character(len = 15) :: datenm_input, datenm_output

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
    do imn = 1, nmn
      do ida = 1, nda
        if (mnlist(imn) > 0 .and. dalist(ida) > eom(mnlist(imn))) exit

        ! print*, 'Opening output files'
        call get_yyyymmdd(yrlist(iyr), mnlist(imn), dalist(ida), yrnm_clm, fn_out_dlm, datenm_output)
        call create_outputfiles(datenm_output)

        ! Initializing input/outut fields used by zeta module'
        call init_zetavars()

        ! print*, 'Loading variable fields from input files'
        call get_yyyymmdd(yrlist(iyr), mnlist(imn), dalist(ida), yrnm_clm, fn_in_dlm, datenm_input)
        call load_vars(datenm_input, yrlist(iyr), mnlist(imn), dalist(ida))

        call zeta_equation()
        if (cmode > 0) then
          call decomp_curlnonl()
        endif

        call write_outputfiles()
        call release_zetavars()
        call close_outputfiles()
        print*, "Finished ", datenm_output
      enddo
    enddo
  enddo
endprogram
