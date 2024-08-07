set title "P-wave Phase Shift for electron scattering from Ne"
set xlabel "Energy (in eV)"
set ylabel "Phase-Shift (in radian)"
set xrange[0:50]

plot "./output/TCS/rmatrixdeph.out" u ($1*27.21):3 w l t'Simu'
replot "./data/electron/Ne/Saha.dat" u ($1*$1*13.605):4 t'Saha (1989)'





