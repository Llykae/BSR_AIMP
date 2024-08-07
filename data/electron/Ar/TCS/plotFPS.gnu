set title "F-wave Phase Shift for electron scattering from Ar"
set xlabel "Energy (in eV)"
set ylabel "Phase-Shift (in radian)"
set xrange[0:40]
set yrange[0:0.5]

plot "./ScatteringResult/Ar/Scatteringtools_noDKH_noC.out" u ($1*27.21):5 w l lc rgb "blue" t'No DKH - Starting Point'
replot "./ScatteringResult/Ar/Scatteringtools_DKH_noC.out" u ($1*27.21):5 w l lc rgb "red" t'DKH - Starting Point'
replot "./ScatteringResult/Ar/Scatteringtools_noDKH_Corr.out" u ($1*27.21):5 w l dashtype 2 lc rgb "blue" t'No DKH - Corrected'
replot "./ScatteringResult/Ar/Scatteringtools_DKH_Corr.out" u ($1*27.21):5 w l dashtype 2 lc  rgb "red" t'DKH - Corrected'
replot "./data/electron/Ar/TCS/McEachran97.dat" u ($1*$1*13.605):6 t'McEachran (1997)'
replot "./data/electron/Ar/TCS/Sienkiewicz87.dat" u 1:5 t'Sienkiewicz (1987)'
replot "./data/electron/Ar/TCS/Dasgupta85.dat" u ($1*$1*13.605):5 t'Dasgupta (1985)'
replot "./data/electron/Ar/TCS/McEachran80.dat-Fwave" u 2:3 t'McEachran (1980)'