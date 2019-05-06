Vorticity budget for 0.1deg POP

Quick start:
Edit environment variables in makefile accordingly
Compile and get zeta.x
Generate an input file "input.nml" using the sample "sample_input.nml"
Run zeta.x

Algorithm:
The vertical vorticity equation now applies the full flux form, in that the stretching and tilting (twisting) terms by the relative vorticity are expressed as fluxes (as in Vallis (2006) p.167)
Additional terms added (and substracted) to merge the non-zero divergence into the twisting terms.
The only residuals from decomposing the nonlinear term come from reversing the operators and applying the chain rule.

Required files:
  Global grid and constants files:
  "ocn_static_grid.nc": 2-D fields ULAT, ULONG, UAREA, TLAT, TLONG, TAREA, DXU, DYU,
      z_t, z_w, HTN, HTE, HUS, HUW
  "ocn_static_dz": 3-D fields DZT, DZU, TMASK, UMASK
  "ocn_constant.nc": constants omega and grav

  Input files with momentum terms
  Depending on the calculation mode, variables needed include UVEL, VVEL, WVEL,
      SSH, ADVU, ADVV, GRADX, GRADY, HDIFFU, HDIFFV, VDIFFU, VDIFFV,
      UEU, UEV, VNU, VNV, WTU, WTV
