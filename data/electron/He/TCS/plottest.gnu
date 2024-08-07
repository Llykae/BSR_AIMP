set title "D-Phase Shift for electron scattering from He"
set xlabel "Energy (in eV)"
set ylabel "Phase Shift"
set xrange[0:20]

plot "rmatrixdeph.out" u ($1*27.21):3 w l dashtype 0 t'Analytic Exchange potential (Erwan code)'
replot "rmSP.out" u ($1*27.21):3 w l t'Actual Exchange potential' 
