set title " Total Cross-Section for e^{-} - Ne"
set xlabel "Energy (in eV)"

set terminal jpeg size 1960,1200 font "{/:Bold} Helvetica,16"
set output "./output/TCS/Ne/Ne_multiplot.jpeg"

set xrange[0:30]

set multiplot layout 2,2

set title "S-wave Phase Shift for electron scattering from Ne"
set ylabel "Phase-Shift (in radian)"
plot "./output/TCS/rmatrixdeph.out" u ($1*27.21):2 w l t'Simulation with correction',\
"./output/TCS/rm_ne_wc.out" u ($1*27.21):2 with lines lt 0 lw 2 dt "-.-" lc 2 t'Simulation without correction',\
"./data/electron/Ne/TCS/Saha_Th90.dat" u 2:3 t'Saha Simu (1990)',\
"./data/electron/Ne/TCS/Saha_Exp90.dat" u 2:3  t'Saha Exp (1990)',\
"./data/electron/Ne/TCS/Saha.dat" u ($1*$1*13.605):3 t'Saha (1989)',\
"./data/electron/Ne/TCS/Dasgupta_EA.dat" u ($1*$1*13.605):3 t'Dasgupta EA (1984)',\
"./data/electron/Ne/TCS/Dasgupta_Ex.dat" u ($1*$1*13.605):3 t'Dasgupta Ex (1984)',\
"./data/electron/Ne/TCS/Dasgupta_PO.dat" u ($1*$1*13.605):3 t'Dasgupta PO (1984)',\
"./data/electron/Ne/TCS/Garbaty_066HFS.dat" u 1:2 t'Garbaty 2/3 HFS (1971)',\
"./data/electron/Ne/TCS/Garbaty_HFS.dat" u 2:3 t'Garbaty HFS (1971)',\
"./data/electron/Ne/TCS/Hoeper.dat" u 1:($2-2*pi) t'Hoeper (1967)'

set title "P-wave Phase Shift for electron scattering from Ne"
set ylabel "Phase-Shift (in radian)"
plot "./output/TCS/rmatrixdeph.out" u ($1*27.21):3 w l t'Simulation with correction',\
"./output/TCS/rm_ne_wc.out" u ($1*27.21):3 with lines lt 0 lw 2 dt "-.-" lc 2 t'Simulation without correction',\
"./data/electron/Ne/TCS/Saha_Exp90.dat" u 2:4 t'Saha (Experimental - 1990)',\
"./data/electron/Ne/TCS/Saha.dat" u ($1*$1*13.605):4 t'Saha (Simulation - 1989)'

set title "D-wave Phase Shift for electron scattering from Ne"
set ylabel "Phase-Shift (in radian)"
plot "./output/TCS/rmatrixdeph.out" u ($1*27.21):4 w l t'Simulation with correction',\
"./output/TCS/rm_ne_wc.out" u ($1*27.21):4 with lines lt 0 lw 2 dt "-.-" lc 2 t'Simulation without correction',\
"./data/electron/Ne/TCS/Saha_Exp90.dat" u 2:5 t'Saha (Experimental - 1990)',\
"./data/electron/Ne/TCS/Saha.dat" u ($1*$1*13.605):5 t'Saha (1989)'

set title "TCS for electron scattering from Ne"
set ylabel "TCS (in a²_0)"
plot "./output/TCS/rmatrixdeph.out" u ($1*27.21):6 w l t'Simulation with correction',\
"./output/TCS/rm_ne_wc.out" u ($1*27.21):6 with lines lt 0 lw 2 dt "-.-" lc 2 t'Simulation without correction',\
"./data/electron/Ne/TCS/Saha.dat" u ($1*$1*13.605):8 t'Saha (1989)'

unset multiplot
