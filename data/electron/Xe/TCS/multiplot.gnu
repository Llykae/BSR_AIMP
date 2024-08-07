set title " Scattering Simulation for e^{-} - Xe"

set terminal jpeg size 1960,1200 font "{/:Bold} Helvetica,16"
set output "./output/TCS/Xe/Xe_multiplot_TCS.jpeg"
set key font ",10"

set multiplot layout 2,3

set title "S-wave Phase-Shift From 0 to 8 eV"
set xrange[0:8]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):2 w l t'Gervais - Present Work 2022',\
"./data/electron/Xe/TCS/Kurokawa_extracted.out" u ($1*27.21):(-$2) t'Kurokawa - Experimental Extraction (2011)',\
"./data/electron/Xe/TCS/Yuan91.dat" u 1:2 with lines lt 0 lw 2 dt "-.-" lc 2 t'Yuan - Simulation (1991)',\
"./data/electron/Xe/TCS/McEachran84.dat" u ($1*$1*13.605):2 with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1984)',\


set title "P-wave Phase-Shift From 0 to 17 eV"
set xrange[0:17]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):3 w l t'Gervais - Present Work 2022',\
"./data/electron/Xe/TCS/Yuan91.dat" u 1:3 with lines lt 0 lw 2 dt "-.-" lc 2 t'Yuan - Simulation (1991)',\
"./data/electron/Xe/TCS/McEachran84.dat" u ($1*$1*13.605):3 with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1984)',\


set title "D-wave Phase-Shift From 0 to 20 eV"
set xrange[0:20]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):4 w l t'Gervais - Present Work 2022',\
"./data/electron/Xe/TCS/Yuan91.dat" u 1:4 with lines lt 0 lw 2 dt "-.-" lc 2 t'Yuan - Simulation (1991)',\
"./data/electron/Xe/TCS/McEachran84.dat" u ($1*$1*13.605):4 with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1984)',\


set title "F-wave Phase-Shift From 0 to 20 eV"
set xrange[0:20]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):5 w l t'Gervais - Present Work 2022',\
"./data/electron/Xe/TCS/Yuan91.dat" u 1:5 with lines lt 0 lw 2 dt "-.-" lc 2 t'Yuan - Simulation (1991)',\
"./data/electron/Xe/TCS/McEachran84.dat" u ($1*$1*13.605):5 with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1984)',\


set title "Total Cross-Section below 3 eV"
set logscale x
set xrange[0.1:3]
set logscale y
set yrange[4:100]
set xlabel "Energy (eV)"
set ylabel "TCS (in a_0^2)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):6 w l t'Gervais - Present Work 2022',\
"./data/electron/Xe/TCS/Kurokawa2011.dat" u 1:($2/(0.529*0.529)):(0.01*$3/(0.529*0.529)) with yerrorbars t"Kurokawa - Experimental (2011)",\
"./data/electron/Xe/TCS/Yuan91.dat" u 1:((4*3.141592/(sqrt($1/13.605)*sqrt($1/13.605)))*(sin($2)*sin($2)+3*sin($3)*sin($3)+5*sin($4)*sin($4)+7*sin($5)*sin($5))) with lines lt 0 lw 2 dt "-.-" lc 2 t'Yuan - Simulation (1991)',\
"./data/electron/Xe/TCS/McEachran84.dat" u ($1*$1*13.605):((4*3.141592/($1*$1))*(sin($2)*sin($2)+3*sin($3)*sin($3)+5*sin($4)*sin($4)+7*sin($5)*sin($5))) with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1984)',\


set title "Total Cross-Section from 0 to 70 eV"
set logscale x
set xrange[0.05:70]
set logscale y
set yrange[4:200]
set xlabel "Energy (eV)"
set ylabel "TCS (in a_0^2)"

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):6 w l t'Gervais - Present Work 2022',\
"./data/electron/Xe/TCS/Kurokawa2011.dat" u 1:($2/(0.529*0.529)):(0.01*$3/(0.529*0.529)) with yerrorbars t"Kurokawa - Experimental (2011)",\
"./data/electron/Xe/TCS/Yuan91.dat" u 1:((4*3.141592/(sqrt($1/13.605)*sqrt($1/13.605)))*(sin($2)*sin($2)+3*sin($3)*sin($3)+5*sin($4)*sin($4)+7*sin($5)*sin($5))) with lines lt 0 lw 2 dt "-.-" lc 2 t'Yuan - Simulation (1991)',\
"./data/electron/Xe/TCS/McEachran84.dat" u ($1*$1*13.605):((4*3.141592/($1*$1))*(sin($2)*sin($2)+3*sin($3)*sin($3)+5*sin($4)*sin($4)+7*sin($5)*sin($5))) with lines lt 0 lw 2 dt "---" lc 4 t'McEachran - Simulation (1984)',\

unset logscale
unset multiplot