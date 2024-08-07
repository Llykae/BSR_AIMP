set size 1.0,1.0
set title "S-wave Phase Shift for electron scattering from Ar"
set xlabel "Energy (in eV)"
set ylabel "Phase-Shift (in radian)"

#set logscale x
set xrange[0.0007:1]
set yrange[-0.25:0.06]

plot "./ScatteringResult/Ar/Scatteringtools_noDKH_noC.out" u ($1*27.21):2 w l lc rgb "blue" t'No DKH - Starting Point'
replot "./ScatteringResult/Scatteringtools.out" u ($1*27.21):2 w l lc rgb "magenta" t'correction 2'
replot "./ScatteringResult/Ar/Scatteringtools_ref.out" u ($1*27.21):2 w l lc rgb "cyan" t'ICI2'
replot "./ScatteringResult/Ar/Scatteringtools_DKH_noC.out" u ($1*27.21):2 w l lc rgb "red" t'DKH - Starting Point'
replot "./ScatteringResult/Ar/Scatteringtools_noDKH_Corr.out" u ($1*27.21):2 w l dashtype 2 lc rgb "blue" t'No DKH - Corrected'
replot "./ScatteringResult/Ar/Scatteringtools_DKH_Corr.out" u ($1*27.21):2 w l dashtype 2 lc  rgb "red" t'DKH - Corrected'
replot "./ScatteringResult/Ar/Scatteringtools_test.out" u ($1*27.21):2 w l dashtype 2 lc  rgb "black" t'Old - Corrected'
replot "./data/electron/Ar/TCS/Ref-McEachran.dat" u ($1*$1*13.605):3 t'McEachran (1997)'
replot "./data/electron/Ar/TCS/Williams79.dat" u 1:2:3 with yerrorbars pointsize 0.5 t'Williams (1979)' 
replot "./data/electron/Ar/TCS/Cheng2022.data" u ($1*$1*13.605):2 t'Cheng - Experimental (2020)'
replot "./data/electron/Ar/TCS/ps_mixed.dat" u ($1*27.21):2 t'Erwan - Mixed Phase Shift'
