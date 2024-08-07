set size 1.0,1.0
set title "S-wave Phase Shift for electron scattering from He"
set xlabel "Energy (in eV)"
set ylabel "Phase-Shift"
set xrange[0:20]

plot "./ScatteringResult/He/Scatteringtools_noDKH_noC.out" u ($1*27.21):($2+pi) w l lc rgb "blue" t'No DKH - Starting Point'
replot "./ScatteringResult/Scatteringtools.out" u ($1*27.21):($2+pi) w l lc rgb "magenta" t'ICI'
replot "./ScatteringResult/He/Scatteringtools_DKH_noC.out" u ($1*27.21):($2+pi) w l lc rgb "red" t'DKH - Starting Point'
replot "./ScatteringResult/He/Scatteringtools_noDKH_Corr.out" u ($1*27.21):($2+pi) w l dashtype 2 lc rgb "blue" t'No DKH - Corrected'
replot "./ScatteringResult/He/Scatteringtools_DKH_Corr.out" u ($1*27.21):($2+pi) w l dashtype 2 lc  rgb "red" t'DKH - Corrected'
replot "./data/electron/He/TCS/McEachran21.dat" u ($1*27.21):($2+pi) t'McEachran (2021)'
replot "./data/electron/He/TCS/Hudson.data" u 1:3 t'Hudson (1996)'
replot "./data/electron/He/TCS/Fon.data" u 1:3 t'Fon (1981)'
replot "./data/electron/He/TCS/Nesbet.data" u 1:3 t'Nesbet (1979)'
replot "./data/electron/He/TCS/Williams.data" u 1:3 t'Williams (1979)'
replot "./data/electron/He/TCS/AB.data" u 1:3 t'Andrick and Bitsch (1975)'
replot "./data/electron/He/TCS/Callaway.data" u 1:3 t'Callaway (1968)'
replot "./data/electron/He/TCS/Phu.data" u 1:3 t'Phu and Chang (1966)'

