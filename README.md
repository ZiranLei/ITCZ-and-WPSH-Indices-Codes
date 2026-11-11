# NCL scripts for calculating ITCZ latitude and WPSH indices

This repository contains the NCL scripts used to calculate the Intertropical Convergence Zone (ITCZ) latitude and the Western Pacific Subtropical High (WPSH) indices for the study:

*Isotopic signatures chart the East Asian rainband's 2,000-km trek driven by coupled tropical-subtropical circulations.*  

## Scripts

| File | Description | Input | Output |
|------|-------------|-------|--------|
| `calc_ITCZ_lat.ncl` | Mean latitudinal position of the ITCZ based on ERA5 500-hPa vertical velocity (ω | `Vertical_velocity.nc` | `ITCZ_lat_output.csv` |
| `calc_WPSH_NBI.ncl` | WPSH Northern Boundary Index: northernmost 5880-gpm crossing at 500 hPa | `Geopotential.nc` | `WPSH_North_Boundary_Index.csv` |
| `calc_WPSH_RI.ncl` | WPSH Ridge Line Index: latitude where 500-hPa zonal wind changes sign (u = 0) and geopotential ≥ 5880 gpm | `Geopotential.nc` , `U_component_wind.nc`  | `WPSH_Ridge_Line_Index.csv` |

## Requirements

- NCAR Command Language (NCL) **v6.6.2 or later**
- ERA5 500-hPa data (e.g., daily or monthly means)

## Usage

1. Edit the input/output paths at the top of each script (the `in_file` / `out_csv` lines).
2. Run directly:

```bash
ncl calc_ITCZ_lat.ncl
ncl calc_WPSH_NBI.ncl
ncl calc_WPSH_RI.ncl
