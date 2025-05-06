# demotools

A set of tools for working with simulated SMBH populations.  The simulations used here come from [Pesce et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021ApJ...923..260P/abstract).

## Making source counts plots

The example script `plot_source_counts.py` shows how to produce SMBH source counts plots similar to those shown in Figure 4 of [Pesce et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021ApJ...923..260P/abstract).  Available observing frequencies are 86, 230, and 345 GHz, and the options for assumed black hole mass functions are [Shankar et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...690...20S/abstract), [Kormendy & Ho (2013)](https://ui.adsabs.harvard.edu/abs/2013ARA%26A..51..511K/abstract), and [Reines & Volonteri (2015)](https://ui.adsabs.harvard.edu/abs/2015ApJ...813...82R/abstract).

## Surveying SMBH mass/spin/shadow measurability thresholds

The example script `survey_measurement_thresholds.py` shows how to survey the space of SMBH shadow diameter and flux density to establish measurability thresholds for mass, spin, and shadow, following the approach described in [Pesce et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022Galax..10..109P/abstract).  The corresponding plots are made using `plot_measurement_thresholds.py`.
