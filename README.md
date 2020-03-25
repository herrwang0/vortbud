# Vorticity budget for 0.1deg POP

## Quick start

1. Edit environment variables in makefile accordingly
2. Compile and get `zeta.x`
3. Generate an input file `input.nml` from the sample `sample_input.nml`
4. Run `zeta.x`

## Example

The testing example is using a small region in the Arabian Sea from Jan 1, 2009 in a 0.1deg POP/CICE simulation. The testing is perform on OLCF Rhea.  
Input: ./example/momentum_2009-01-01_as.nc  
Output: ./example/vorticity_2009-01-01_as.nc
In addition, it uses grid files stored at /ccs/home/hewang/pop0_1/ on Rhea
To test:  

1. copy ./example/input.nml to ./  
2. run ./zeta.x > log.out
3. compare sampled outputs in log.out with ./example/log.out

## Algorithm

<!-- The vertical vorticity equation now applies the full flux form, in that the stretching and tilting (twisting) terms by the relative vorticity are expressed as fluxes (as in Vallis (2006) p.167)
Additional terms added (and substracted) to merge the non-zero divergence into the twisting terms.
The only residuals from decomposing the nonlinear term come from reversing the operators and applying the chain rule. -->

## Required file

* Global grid and constants files
  * "ocn_static_grid.nc": 2-D fields `ULAT`, `ULONG`, `UAREA`, `TLAT`, `TLONG`, `TAREA`, `DXU`, `DYU`, `z_t`, `z_w`, `HTN`, `HTE`, `HUS`, `HUW`
  * "ocn_static_dz": 3-D fields `DZT`, `DZU`, `TMASK`, `UMASK`
  * "ocn_constant.nc": constants `omega` and `grav`

* Input files with momentum terms  
  Depending on the calculation mode, variables needed include 
  * `UVEL`, `VVEL`, `WVEL`
  * `SSH`, `ADVU`, `ADVV`, `GRADX`, `GRADY`, `HDIFFU`, `HDIFFV`, `VDIFFU`, `VDIFFV`,
  * `UEU`, `UEV`, `VNU`, `VNV`, `WTU`, `WTV`
