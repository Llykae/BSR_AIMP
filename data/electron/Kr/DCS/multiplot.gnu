set title " Differential-Cross at different energy for e^{-} - Kr"
set xlabel "Angle (degree)"
set ylabel "DCS (a_0².sr^{-1})

set terminal jpeg size 1960,1200 font "{/:Bold} Helvetica,16"
set output "./output/DCS/Kr/Kr_multiplot_DCS_6L.jpeg"

set multiplot layout 2,2
set logscale y

set title "Differential Cross-Section at 5 eV"
plot "./output/DCS/Kr/DCS_5eV.out" u 1:2 w l t'Gervais - Present Work',\
"./data/electron/Kr/DCS/McEachran87.data" u 1:3 with lines lt 0 lw 2 dt "-.-" lc 2 t'McEachran - Simulation (1987)',\
"./data/electron/Kr/DCS/Linert2010.data" u 1:($2/(0.529*0.529)):($2*0.15/(0.529*0.529)) with yerrorbars t'Linert - Experimental (2010)',\
"./data/electron/Kr/DCS/Srivastava81.data" u 1:($3*0.529*0.529):($3*0.2*0.529*0.529) with yerrorbars lc rgb "red" t'Srivastava - Experimental (1981)',\

set title " Differential Cross-Section at 7.5 eV"
plot "./output/DCS/Kr/DCS_7.5eV.out" u 1:2 w l t'Gervais - Present Work',\
"./data/electron/Kr/DCS/McEachran87.data" u 1:4 with lines lt 0 lw 2 dt "-.-" lc 2 t'McEachran - Simulation (1987)',\
"./data/electron/Kr/DCS/Linert2010.data" u 1:($3/(0.529*0.529)):($3*0.15/(0.529*0.529)) with yerrorbar t'Linert - Experimental (2010)',\
"./data/electron/Kr/DCS/Srivastava81.data" u 1:($4*0.529*0.529):($4*0.2*0.529*0.529) with yerrorbars lc rgb "red" t'Srivastava - Experimental (1981)',\

set title " Differential Cross-Section at 10 eV"
plot "./output/DCS/Kr/DCS_10eV.out" u 1:2 w l t'Gervais - Present Work',\
"./data/electron/Kr/DCS/McEachran87.data" u 1:5 with lines lt 0 lw 2 dt "-.-" lc 2 t'McEachran - Simulation (1987)',\
"./data/electron/Kr/DCS/Linert2010.data" u 1:($4/(0.529*0.529)):($4*0.15/(0.529*0.529)) with yerrorbar t'Linert - Experimental (2010)',\
"./data/electron/Kr/DCS/Srivastava81.data" u 1:($5*0.529*0.529):($5*0.2*0.529*0.529) with yerrorbars lc rgb "red" t'Srivastava - Experimental (1981)',\

set title " Differential Cross-Section at 20 eV"
plot "./output/DCS/Kr/DCS_20eV.out" u 1:2 w l t'Gervais - Present Work',\
"./data/electron/Kr/DCS/McEachran87.data" u 1:7 with lines lt 0 lw 2 dt "-.-" lc 2 t'McEachran - Simulation (1987)',\
"./data/electron/Kr/DCS/Srivastava81.data" u 1:($7*0.529*0.529):($7*0.2*0.529*0.529) with yerrorbars lc rgb "red" t'Srivastava - Experimental (1981)',\

unset multiplot