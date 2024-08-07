set title " Differential-Cross at different energy for e^{-} - Ne"
set xlabel "Angle (degree)"
set ylabel "DCS (angström².sr^{-1})

set terminal jpeg size 1960,1200 font "{/:Bold} Helvetica,16"
set output "./output/DCS/Ne/Ne_multiplot2_DCS_6L.jpeg"

set multiplot layout 2,2

set title "Differential Cross-Section at 5 eV"
plot "./output/DCS/Ne/Without_Correction/DCS_5eV.out" u 1:($2*0.529*0.529) with lines lt 0 lw 2 dt " _ " lc 8  t'Simulation without Correction',\
"./data/electron/Ne/DCS/Khandker20.data" u 1:($2*0.529*0.529) with lines lt 0 lw 2 dt "... " lc 2 t'Khandker - Simulation (2020)',\
"./data/electron/Ne/DCS/Cho08_5eV.data" u 1:2:3 with yerrorbars lc rgb "red" t'Cho - Experimental (2008)',\
"./data/electron/Ne/DCS/Registrer84_5eV.data" u 1:2:($2*0.07) with yerrorbars t'Register - Experimental (1984)'

set title " Differential Cross-Section at 10 eV"
plot "./output/DCS/Ne/Without_Correction/DCS_10eV.out" u 1:($2*0.529*0.529) with lines lt 0 lw 2 dt "..." lc 8  t'Simulation without Correction',\
"./data/electron/Ne/DCS/Khandker20.data" u 1:($3*0.529*0.529) with lines lt 0 lw 2 dt "-.-" lc 2  t'Khandker - Simulation (2020)',\
"./data/electron/Ne/DCS/Cho08_10eV.data" u 1:2:3 with yerrorbars lc rgb "red" t'Cho - Experimental (2008)',\
"./data/electron/Ne/DCS/Registrer84_10eV.data" u 1:2:($2*0.07) with yerrorbars t'Register - Experimental (1984)'


set title " Differential Cross-Section at 15 eV"
plot "./output/DCS/Ne/Without_Correction/DCS_15eV.out" u 1:($2*0.529*0.529) with lines lt 0 lw 2 dt  "..." lc 8 t'Simulation without Correction',\
"./data/electron/Ne/DCS/Khandker20.data" u 1:($4*0.529*0.529) with lines lt 0 lw 2 dt "-.-" lc 2 t'Khandker - Simulation (2020)',\
"./data/electron/Ne/DCS/Registrer84_15eV.data" u 1:2:($2*0.07) with yerrorbars t'Register - Experimental (1984)'

set title " Differential Cross-Section at 20 eV"
plot "./output/DCS/Ne/Without_Correction/DCS_20eV.out" u 1:($2*0.529*0.529) with lines lt 0 lw 2 dt "..." lc 8 t'Simulation without Correction',\
"./data/electron/Ne/DCS/Khandker20.data" u 1:($5*0.529*0.529) with lines lt 0 lw 2 dt "-.-" lc 2 t'Khandker - Simulation (2020)',\
"./data/electron/Ne/DCS/Cho08_20eV.data" u 1:2:3 with yerrorbars lc rgb "red" t'Cho - Experimental (2008)',\
"./data/electron/Ne/DCS/Registrer84_20eV.data" u 1:2:($2*0.07) with yerrorbars t'Register - Experimental (1984)'

unset multiplot