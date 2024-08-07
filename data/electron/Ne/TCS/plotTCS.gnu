set title "TCS for electron scattering from Ne"
set xlabel "Energy (in eV)"
set ylabel "Total Cross Section"

set logscale x

plot "./ScatteringResult/Ne/Scatteringtools_noDKH_noC.out" u ($1*27.21):5 w l lc rgb "blue" t'No DKH - Starting Point'
replot "./ScatteringResult/Ne/Scatteringtools_DKH_noC.out" u ($1*27.21):5 w l lc rgb "red" t'DKH - Starting Point'
replot "./ScatteringResult/Ne/Scatteringtools_noDKH_Corr.out" u ($1*27.21):5 w l dashtype 2 lc rgb "blue" t'No DKH - Corrected'
replot "./ScatteringResult/Ne/Scatteringtools_DKH_Corr.out" u ($1*27.21):5 w l dashtype 2 lc  rgb "red" t'DKH - Corrected'
replot "./ScatteringResult/Ne/Scatteringtools_newpola.out" u ($1*27.21):5 w l lc rgb "black" t'New pola - No Corrected'
replot "./ScatteringResult/Ne/Scatteringtools_newpola_Corr.out" u ($1*27.21):5 w l dashtype 2 lc rgb "black" t'New pola - Corrected'
replot "./data/electron/Ne/TCS/Saha.dat" u ($1*$1*13.605):8 t'Saha (1989)'
replot "./data/electron/Ne/TCS/Saha_Th90.dat" u 2:6 t'Saha Simu (1990)'
replot "./data/electron/Ne/TCS/Saha_Exp90.dat" u 2:6  t'Saha Exp (1990)'
replot "./data/electron/Ne/TCS/Gulley.dat" u ($1*27.21):5 t'Gulley (MERT5)'
