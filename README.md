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

## Input files (containing momentum fields) and namelist (input.nml) specifications:

* Input files are supposed to be netCDF files with time record length = 1. The fields can be of daily, monthly or annually frequency.

* It is recommended to use at least daily averaged fields for nonlinear decomposition. Technically, the decomposition can be done for any frequency, but large errors are expected with monthly and annual fields.

* `sec` controls the period of mean/eddy decomposition and overall output frequency. It can be set by a list of options (see comments in sample_input.nml for details), for an arbitrary average length (annually, monthly, fixed number of days etc) within a year.

* Mean/Eddy decomposition section frequency `meanfreq` can only be "a" for annual fields and "a" or "m" for monthly fields


### Climatology (`yrnm_clm` /= "")

If the daily, monthly and annual inpupt files are also averged across multiple years, the input files are considered to be "**climatology**".

* `yrnm_clm` is required. It is a *wildcard* option to replace the numerics for "year"  in the filename. it can be set to be anything, though it is suggested to follow the form of "YYYY-YYYY".
* When mean/eddy decomposition is turned on (`ifmeaneddy` = T), `ifmeanclm` is automatically set to T.
* When nonlinear decomposition is turned on (`ifdecomp` = T), it is calculated only for the mean fields. In this case, `ifmeaneddy` is equivalent of being turned on.  
    * In other words, when `ifdecomp` = F and `ifmeaneddy` = T, the outputs are curladv (mean), curlmet (mean), errsub (eddy); when `ifdecomp` = T, `ifmeaneddy` is ignored and the outputs are adv(u,v,w) (mean), twi(u,v,w) (mean), curlmet (mean), errsub (eddy)

### For **daily** averaged files

1. Curl of momentum terms (output at the same frequency as input)  
    `ifcurl` = T; `ifdecomp` = F;  
    `ifmeaneddy` = F  
    * func = c, func_m = NULL, func_me = NULL  
2. Curl + nonlinear decomposition (output at the same frequency as input)  
    `ifcurl` = T; `ifdecomp` = T; [`flxtwi` = T for flux twisting, and F for nonflux]  
    `ifmeaneddy` = F
    * func = cdme, func_m = NULL, func_me = NULL
3. Curl + mean/eddy decomposition (output at section frequency)  
    `ifcurl` = T; `ifdecomp` = F;  
    `ifmeaneddy` = F; `ifdecomposed` = F; [`ifmeanclm` = T for **clim** mean and F for **ann** mean]  
    * func = c, func_m = am, func_me = func
4. Curl + nonlinear decomposition + mean/eddy decomposition (output at section frequency)  
    `ifcurl` = T; `ifdecomp` = T; [`flxtwi` = T for flux twisting, and F for nonflux]  
    `ifmeaneddy` = T; `ifdecomposed` = T; [`ifmeanclm` = T for **clim** mean and F for **ann** mean]  
    *  func = cdme, func_m = dm, func_me = func

### For **monthly** averaged files

* `dalist_in` should be set to any single negative value, which sets doylist to the negative month number
* Nonlinear decomposition not recommened
* `meanfreq` = 'm' or 'a'

### For **annually** averaged files

* `dalist_in` and `mnlist_in` should be set to any single negative value, which sets doylist to 0
* Nonlinear decomposition not recommened
* `meanfreq` = 'a'

### Other Tips

* For mean/eddy decomposiion, the mean nonlinear fields are calculated first and saved, which means the files can be reused for re-calculations.

* For fully-decomposed nonlinear term, output frequency is also the same as the input. Mean/eddy decomposion then reads and averages these files.

* Resume calculation when using daily-averaged input with fully-decomposed nonlinear term:
    * Streamlines:
        1. Generate files containing mean nonlinear field
        2. Files containing curl of each term and decomposed nonlinear term of each day
        3. Load and average files from 2 and output mean (from 1) and eddy (difference averged 2 and 1)
        * If the calculation is interrupted during 2:  
            1. Turn off `ifmeaneddy` and reset lists of days, month and years and resume
            2. Turn off `ifcurl` and `ifdecomp` and turn on `ifmeaneddy`

## Algorithm
See the appendix of 10.1175/JPO-D-20-0223.1
<!-- The vertical vorticity equation now applies the full flux form, in that the stretching and tilting (twisting) terms by the relative vorticity are expressed as fluxes (as in Vallis (2006) p.167)
Additional terms added (and substracted) to merge the non-zero divergence into the twisting terms.
The only residuals from decomposing the nonlinear term come from reversing the operators and applying the chain rule. -->
