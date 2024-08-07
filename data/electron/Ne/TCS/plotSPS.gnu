set size 1.0,1.0
set title "S-wave Phase Shift for electron scattering from Ne"
set xlabel "Energy (in eV)"
set ylabel "Phase-Shift"
set xrange[0:20]
set yrange[-1.2:0]

plot "./ScatteringResult/Ne/Scatteringtools_noDKH_noC.out" u ($1*27.21):2 w l lc rgb "blue" t'No DKH - Starting Point'
replot "./ScatteringResult/Scatteringtools.out" u ($1*27.21):2 w l lc rgb "magenta" t'ICI'
replot "./ScatteringResult/Ne/Scatteringtools_DKH_noC.out" u ($1*27.21):2 w l lc rgb "red" t'DKH - Starting Point'
replot "./ScatteringResult/Ne/Scatteringtools_noDKH_Corr.out" u ($1*27.21):2 w l dashtype 2 lc rgb "blue" t'No DKH - Corrected'
replot "./ScatteringResult/Ne/Scatteringtools_DKH_Corr.out" u ($1*27.21):2 w l dashtype 2 lc  rgb "red" t'DKH - Corrected'
replot "./ScatteringResult/Ne/Scatteringtools_newpola.out" u ($1*27.21):2 w l lc rgb "black" t'New pola - No Corrected'
replot "./ScatteringResult/Ne/Scatteringtools_newpola_Corr.out" u ($1*27.21):2 w l dashtype 2 lc rgb "black" t'New pola - Corrected'
replot "./data/electron/Ne/TCS/Saha_Th90.dat" u 2:3 t'Saha Simu (1990)'
replot "./data/electron/Ne/TCS/Saha_Exp90.dat" u 2:3  t'Saha Exp (1990)'
replot "./data/electron/Ne/TCS/Saha.dat" u ($1*$1*13.605):3 t'Saha (1989)'
replot "./data/electron/Ne/TCS/Gulley.dat" u ($1*27.21):2 t'Gulley (MERT5)'
