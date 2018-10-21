Vorticity decomposition for 0.1deg POP
make produces two executables
Run zeta.x to calculate vorticity decomposition for each input file
Run zeta_meaneddy.x to decompose mean eddy in the nonlinear term

v4.2
1. Curl operator is reversed to "zcurl_dz", in which dz is taken into consideration.
2. Add subroutine "find_daily_file" in mod_popload.f90. 

v4.1
Algorithm
The vertical vorticity equation now applies the full flux form, in that the stretching and tilting (twisting) terms by the relative vorticity are expressed as fluxes (as in ...).
Additional terms added (and substracted) to merge the non-zero divergence into the twisting terms.
The only residuals from decomposing the nonlinear term come from reversing the operators and applying the chain rule.
