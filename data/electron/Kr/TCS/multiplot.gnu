set title " Scattering Simulation for e^{-} - Kr"

set terminal jpeg size 1960,1200 font "{/:Bold} Helvetica,16"
set output "./output/TCS/Kr/Kr_multiplot_TCS.jpeg"
set key font ",10"

set multiplot layout 2,3

set title "S-wave Phase-Shift From 0 to 13 eV"
set xrange[0:13]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):2 w l t'Gervais - Present Work 2022',\
"./data/electron/Ar/TCS/Kurokawa_extracted.out" u ($1*27.21):(-$2) t'Kurokawa - Experimental Extraction (2011)',\
"./data/electron/Kr/TCS/Bell88.data" u 1:3 with lines lt 0 lw 2 dt "-.-" lc 2 t'Bell - Simulation (1988)',\
"./data/electron/Kr/TCS/McEachran87.data" u ($1*$1*13.605):2 with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1987)',\


set title "P-wave Phase-Shift From 0 to 20 eV"
set xrange[0:20]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):3 w l t'Gervais - Present Work 2022',\
"./data/electron/Ar/TCS/Kurokawa_extracted.out" u ($1*27.21):(-$3) t'Kurokawa - Experimental Extraction (2011)',\
"./data/electron/Kr/TCS/Bell88.data" u 1:4 with lines lt 0 lw 2 dt "-.-" lc 2 t'Bell - Simulation (1988)',\
"./data/electron/Kr/TCS/McEachran87.data" u ($1*$1*13.605):3 with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1987)',\


set title "D-wave Phase-Shift From 0 to 20 eV"
set xrange[0:20]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):4 w l t'Gervais - Present Work 2022',\
"./data/electron/Ar/TCS/Kurokawa_extracted.out" u ($1*27.21):4 t'Kurokawa - Experimental Extraction (2011)',\
"./data/electron/Kr/TCS/Bell88.data" u 1:5 with lines lt 0 lw 2 dt "-.-" lc 2 t'Bell - Simulation (1988)',\
"./data/electron/Kr/TCS/McEachran87.data" u ($1*$1*13.605):4 with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1987)',\


set title "F-wave Phase-Shift From 0 to 20 eV"
set xrange[0:20]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):5 w l t'Gervais - Present Work 2022',\
"./data/electron/Ar/TCS/Kurokawa_extracted.out" u ($1*27.21):5 t'Kurokawa - Experimental Extraction (2011)',\
"./data/electron/Kr/TCS/Bell88.data" u 1:6 with lines lt 0 lw 2 dt "-.-" lc 2 t'Bell - Simulation (1988)',\
"./data/electron/Kr/TCS/McEachran87.data" u ($1*$1*13.605):5 with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1987)',\


set title "Total Cross-Section below 3 eV"
set logscale x
set xrange[0.01:3]
set logscale y
set yrange[1:150]
set xlabel "Energy (eV)"
set ylabel "TCS (in a_0^2)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):6 w l t'Gervais - Present Work 2022',\
"./data/electron/Kr/TCS/Kurokawa2010.data" u 1:($2/(0.529*0.529)) with points pt 6 lc rgb "blue" t"Kurokawa - Experimental (2010)",\
"./data/electron/Kr/TCS/Bell88.data" u 1:($2*3.141592) with lines lt 0 lw 2 dt "-.-" lc 2 t'Bell - Simulation (1988)',\
"./data/electron/Kr/TCS/McEachran87.data" u ($1*$1*13.605):((4*3.141592/($1*$1))*(sin($2)*sin($2)+3*sin($3)*sin($3)+5*sin($4)*sin($4)+7*sin($5)*sin($5))) with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1987)',\


set title "Total Cross-Section from 0 to 30 eV"
set logscale x
set xrange[0.01:30]
set logscale y
set xlabel "Energy (eV)"
set ylabel "TCS (in a_0^2)"

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):6 w l t'Gervais - Present Work 2022',\
"./data/electron/Kr/TCS/Kurokawa2010.data" u 1:($2/(0.529*0.529)) with points pt 6 lc rgb "blue" t"Kurokawa - Experimental (2010)",\
"./data/electron/Kr/TCS/Bell88.data" u 1:($2*3.141592) with lines lt 0 lw 2 dt "-.-" lc 2 t'Bell - Simulation (1988)',\
"./data/electron/Kr/TCS/McEachran87.data" u ($1*$1*13.605):((4*3.141592/($1*$1))*(sin($2)*sin($2)+3*sin($3)*sin($3)+5*sin($4)*sin($4)+7*sin($5)*sin($5))) with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1987)',\

unset multiplot