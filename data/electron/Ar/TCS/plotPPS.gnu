set title "P-wave Phase Shift for electron scattering from Ar"
set xlabel "Energy (in eV)"
set ylabel "Phase-Shift (in radian)"
set xrange[0:40]

plot "./ScatteringResult/Ar/Scatteringtools_noDKH_noC.out" u ($1*27.21):3 w l lc rgb "blue" t'No DKH - Starting Point'
replot "./ScatteringResult/Scatteringtools.out" u ($1*27.21):3 w l lc rgb "magenta" t'ICI'
replot "./ScatteringResult/Ar/Scatteringtools_DKH_noC.out" u ($1*27.21):3 w l lc rgb "red" t'DKH - Starting Point'
replot "./ScatteringResult/Ar/Scatteringtools_noDKH_Corr.out" u ($1*27.21):3 w l dashtype 2 lc rgb "blue" t'No DKH - Corrected'
replot "./ScatteringResult/Ar/Scatteringtools_DKH_Corr.out" u ($1*27.21):3 w l dashtype 2 lc  rgb "red" t'DKH - Corrected'
replot "./ScatteringResult/Ar/Scatteringtools_test.out" u ($1*27.21):3 w l dashtype 2 lc  rgb "black" t'Old - Corrected'
replot "./data/electron/Ar/TCS/McEachran97.dat" u ($1*$1*13.605):4 t'McEachran (1997)'
replot "./data/electron/Ar/TCS/Williams79.dat" u 1:4:5 with yerrorbars pointsize 0.5 t'Williams (1979)'
replot "./data/electron/Ar/TCS/Kurokawa_ps_extracted.out" u ($1*27.21):(-$3) t'Kurokawa - Extracted (2011)'
replot "./data/electron/Ar/TCS/Cheng2022.data" u ($1*$1*13.605):3 t'Cheng - Experimental (?) (2020)'
replot "./data/electron/Ar/TCS/ps_mixed.dat" u ($1*27.21):3 t'Erwan - Mixed Phase Shift'
replot "./data/electron/Ar/TCS/Ref-McEachran.dat" u ($2*27.21):4 t'Ref - McEachran'