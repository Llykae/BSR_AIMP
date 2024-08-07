set title "D-wave Phase Shift for electron scattering from He"
set xlabel "Energy (in eV)"
set ylabel "Phase-Shift (in radian)"

plot "./ScatteringResult/He/Scatteringtools_noDKH_noC.out" u ($1*27.21):4 w l lc rgb "blue" t'No DKH - Starting Point'
replot "./ScatteringResult/He/Scatteringtools_DKH_noC.out" u ($1*27.21):4 w l lc rgb "red" t'DKH - Starting Point'
replot "./ScatteringResult/He/Scatteringtools_noDKH_Corr.out" u ($1*27.21):4 w l dashtype 2 lc rgb "blue" t'No DKH - Corrected'
replot "./ScatteringResult/He/Scatteringtools_DKH_Corr.out" u ($1*27.21):4 w l dashtype 2 lc  rgb "red" t'DKH - Corrected'
replot "./data/electron/He/TCS/McEachran21.dat" u ($1*27.21):4 t'McEachran (2021)'
replot "./data/electron/He/TCS/Hudson.data" u 1:5 t'Hudson (1996)'
replot "./data/electron/He/TCS/Fon.data" u 1:5 t'Fon (1979)'
replot "./data/electron/He/TCS/Nesbet.data" u 1:5 t'Nesbet (1979)'
replot "./data/electron/He/TCS/Williams.data" u 1:5 t'Williams (1979)'




