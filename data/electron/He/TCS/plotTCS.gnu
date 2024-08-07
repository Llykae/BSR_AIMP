set title "TCS for electron scattering from He"
set xlabel "Energy (in eV)"
set ylabel "Total Cross Section"
set xrange[0:20]

plot "./ScatteringResult/He/Scatteringtools_noDKH_noC.out" u ($1*27.21):5 w l lc rgb "blue" t'No DKH - Starting Point'
replot "./ScatteringResult/He/Scatteringtools_DKH_noC.out" u ($1*27.21):5 w l lc rgb "red" t'DKH - Starting Point'
replot "./ScatteringResult/He/Scatteringtools_noDKH_Corr.out" u ($1*27.21):5 w l dashtype 2 lc rgb "blue" t'No DKH - Corrected'
replot "./ScatteringResult/He/Scatteringtools_DKH_Corr.out" u ($1*27.21):5 w l dashtype 2 lc  rgb "red" t'DKH - Corrected'
replot "./data/electron/He/TCS/Hudson.data" u 1:6 t'Hudson (1996)'
replot "./data/electron/He/TCS/Fon.data" u 1:6 t'Fon (1981)'
replot "./data/electron/He/TCS/Nesbet.data" u 1:6 t'Nesbet (1979)'
replot "./data/electron/He/TCS/Williams.data" u 1:6 t'Williams (1979)'
replot "./data/electron/He/TCS/AB.data" u 1:5 t'Andrick and Bitsch (1975)'
replot "./data/electron/He/TCS/Callaway.data" u 1:5 t'Callaway (1968)'
replot "./data/electron/He/TCS/Phu.data" u 1:5 t'Phu and Chang (1966)'

