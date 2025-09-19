# Black Globe Temperature Model in Python

*This repository provides a Python version of the **globe-temperature (Tglobe)** component from James C. Liljegren’s WBGT model.*

> Note: Some default constants (globe diameter, convergence, and maximum iterations number) may be set to match the instrument used in my dataset. See the header in `bg_temp_model.py`.

### Sections included from Liljegren's WBGT model
- `solarposition(year, month, day, days_1900, latitude, longitude)`
- `calc_solar_parameters(year, month, day_fraction, lat, lon, solar_ghi)`
- Air/thermal property helper functions (viscosity, thermal conductivity, saturated vapor pressure, atmospheric emissivity, convective heat transfer coefficient)
- `Tglobe(Tair_K, rh_frac, Pair_mb, speed, solar, fdir, cza)`


*Not included in this model: natural/psychrometric wet-bulb, stability/height wind extrapolation (those are part of the full WBGT program, not this focused port).*
