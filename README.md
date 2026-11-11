# Mesoscale eddy analysis southwest of Taiwan

These six MATLAB files document how I analysed sea-level variability, Kuroshio intrusion, and the evolution of a tracked mesoscale eddy southwest of Taiwan.

| MATLAB file | What it computes |
|---|---|
| [KSI_analysis.m](ksi/KSI_analysis.m) | Computes the Kuroshio South China Sea Index from daily absolute dynamic topography and counts events by path class using the supplied event table. |
| [EOF_Taiwan_SW_2014_2023.m](eof/EOF_Taiwan_SW_2014_2023.m) | Decomposes area-weighted daily sea-level anomalies into spatial patterns and time coefficients for 2014–2023. |
| [RUN_eddy_fullpath_vorticity.m](eddy_diagnostics/RUN_eddy_fullpath_vorticity.m) | Maps planetary potential vorticity advection (APV) and the joint effect of baroclinicity and relief (JEBAR) at eight dates. Both diagnostics use upper-ocean data and full bathymetric depth. |
| [RUN_eddy_fullpath_conversion.m](eddy_diagnostics/RUN_eddy_fullpath_conversion.m) | Computes barotropic and baroclinic energy conversion in a 30-day window around each of the eight dates. |
| [Energy.m](energy_reservoirs/Energy.m) | Computes annual mean and eddy kinetic energy and includes the baseline potential-energy calculation. |
| [Energy_conversion.m](energy_reservoirs/Energy_conversion.m) | Computes the updated annual mean and eddy available potential energy with TEOS-10 stratification, alongside annual conversion diagnostics. |

## Running the scripts

The scripts were developed for MATLAB R2022b. To reproduce the figures, configure the local data and output paths and supply the relevant NetCDF fields, bathymetry, event table, and toolboxes (including GSW and M_Map). The two eddy-stage scripts also read a local track file named `typical_eddy_track.csv` for the eddy-centre markers. This repository contains the analysis code only.