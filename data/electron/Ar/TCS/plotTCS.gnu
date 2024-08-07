set title "TCS for electron scattering from Ar"
set xlabel "Energy (in eV)"
set ylabel "Total Cross Section (in a_0^2)"
set xrange[0:45]

set logscale x
set logscale y

plot "./ScatteringResult/Ar/Scatteringtools_noDKH_noC.out" u ($1*27.21):6 w l lc rgb "blue" t'No DKH - Starting Point'
replot "./ScatteringResult/Scatteringtools.out" u ($1*27.21):5 w l lc rgb "magenta" t'correction 2'
replot "./ScatteringResult/Ar/Scatteringtools_DKH_noC.out" u ($1*27.21):6 w l lc rgb "red" t'DKH - Starting Point'
replot "./ScatteringResult/Ar/Scatteringtools_noDKH_Corr.out" u ($1*27.21):6 w l dashtype 2 lc rgb "blue" t'No DKH - Corrected'
replot "./ScatteringResult/Ar/Scatteringtools_DKH_Corr.out" u ($1*27.21):6 w l dashtype 2 lc  rgb "red" t'DKH - Corrected'
replot "./ScatteringResult/Ar/Scatteringtools_ref.out" u ($1*27.21):5 w l lc rgb "cyan" t'ICI2'
replot "./data/electron/Ar/TCS/Williams79.dat" u 1:8:($9/(0.529*0.529)) with yerrorbars pointsize 0.5 t'Williams (1979)'
replot "./data/electron/Ar/TCS/Kurokawa2011.dat" u 1:($2/(0.529*0.529)) t'Kurokawa - Experimental (2011)'
replot "./data/electron/Ar/TCS/Cheng2022.data" u ($1*$1*13.605):5 t'Cheng - Fit Experimental (2020)'
replot "./data/electron/Ar/TCS/Sienkiewicz87.dat" u 1:6 t'Sienkiewicz (1987)'
replot "./data/electron/Ar/TCS/Dasgupta85.dat" u ($1*$1*13.605):6 t'Dasgupta (1985)'
replot "./data/electron/Ar/TCS/Ref-McEachran.dat" u ($2*27.21):6 t'Ref - McEachran'