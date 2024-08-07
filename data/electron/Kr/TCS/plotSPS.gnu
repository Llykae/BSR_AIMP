set size 1.0,1.0
set title "S-wave Phase Shift for electron scattering from Kr"
set xlabel "Energy (in eV)"
set ylabel "Phase-Shift"
set xrange[0:13]

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):2 w l t'Gervais - Present Work 2022',\
"./data/electron/Kr/TCS/Bell88.data" u 1:3 with lines lt 0 lw 2 dt "-.-" lc 2 t'Bell - Simulation (1988)',\
"./data/electron/Kr/TCS/McEachran87.data" u ($1*$1*13.605):2 with lines lt 0 lw 2 dt "-.-" lc 4 t'McEachran - Simulation (1987)',\
