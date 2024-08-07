set title " Differential-Cross at different energy for e^{-} - Ar"
set xlabel "Angle (degree)"
set ylabel "DCS (angström².sr^{-1})

set terminal jpeg size 1960,1200
#set terminal wxt size 2400,1000
set output "./output/DCS/Ar/Ar_multiplot_DCS_6L_logscale.jpeg"

set xrange[20:150]
set logscale y

set multiplot layout 2,3

set title " Differentiel Cross-Section at 3 eV"
plot "./output/DCS/Ar/DCS_3eV.out" u 1:($2*0.529*0.529) w l t'Simulation',\
"./data/electron/Ar/DCS/Nahar87.data" u 1:($2*0.529*0.529) with lines lt 0 lw 2 dt "-.-" lc 2  t'Nahar - Simulation (1987)',\
"./data/electron/Ar/DCS/McEachran.data" u 1:($2*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 3 t'McEachran',\
"./data/electron/Ar/DCS/Sienkiewicz.data" u 1:($2*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 4 t'Sienkiewicz',\
"./data/electron/Ar/DCS/Dasgupta.data" u 1:($2*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 5 t'Dasgupta'


set title " Differentiel Cross-Section at 5 eV"
plot "./output/DCS/Ar/DCS_5eV.out" u 1:($2*0.529*0.529) w l t'Simulation',\
"./data/electron/Ar/DCS/Nahar87.data" u 1:($3*0.529*0.529) with lines lt 0 lw 2 dt "-.-" lc 2 t'Nahar - Simulation (1987)',\
"./data/electron/Ar/DCS/McEachran.data" u 1:($3*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 3 t'McEachran',\
"./data/electron/Ar/DCS/Sienkiewicz.data" u 1:($3*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 4 t'Sienkiewicz',\
"./data/electron/Ar/DCS/Dasgupta.data" u 1:($3*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 5 t'Dasgupta'

set title " Differentiel Cross-Section at 10 eV"
plot "./output/DCS/Ar/DCS_10eV.out" u 1:($2*0.529*0.529) w l t'Simulation',\
"./data/electron/Ar/DCS/Panajotovic97.data" u 1:2:3 with yerrorbars t'Panajotovic - Experimental (1997)',\
"./data/electron/Ar/DCS/Nahar87.data" u 1:($4*0.529*0.529)  with lines lt 0 lw 2 dt "-.-" lc 2 t'Nahar - Simulation (1987)',\
"./data/electron/Ar/DCS/Srivastava81.data" u 1:4:5 with yerrorbars t'Srivastava - Experimental (1981)',\
"./data/electron/Ar/DCS/McEachran.data" u 1:($5*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 3 t'McEachran',\
"./data/electron/Ar/DCS/Sienkiewicz.data" u 1:($5*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 4 t'Sienkiewicz',\
"./data/electron/Ar/DCS/Dasgupta.data" u 1:($4*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 5 t'Dasgupta'

set title " Differentiel Cross-Section at 15 eV"
plot "./output/DCS/Ar/DCS_15eV.out" u 1:($2*0.529*0.529) w l t'Simulation',\
"./data/electron/Ar/DCS/Panajotovic97.data" u 1:4:5 with yerrorbars t'Panajotovic - Experimental (1997)',\
"./data/electron/Ar/DCS/Nahar87.data" u 1:($5*0.529*0.529) with lines lt 0 lw 2 dt "-.-" lc 2 t'Nahar - Simulation (1987)',\
"./data/electron/Ar/DCS/Srivastava81.data" u 1:6:7 with yerrorbars t'Srivastava - Experimental (1981)',\
"./data/electron/Ar/DCS/McEachran.data" u 1:($6*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 3 t'McEachran',\
"./data/electron/Ar/DCS/Sienkiewicz.data" u 1:($6*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 4 t'Sienkiewicz'

set title " Differentiel Cross-Section at 20 eV"
plot "./output/DCS/Ar/DCS_20eV.out" u 1:($2*0.529*0.529) w l t'Simulation',\
"./data/electron/Ar/DCS/Panajotovic97.data" u 1:6:7 with yerrorbars t'Panajotovic - Experimental (1997)',\
"./data/electron/Ar/DCS/Nahar87.data" u 1:($6*0.529*0.529) with lines lt 0 lw 2 dt "-.-" lc 2 t'Nahar - Simulation (1987)',\
"./data/electron/Ar/DCS/Srivastava81.data" u 1:8:9 with yerrorbars t'Srivastava - Experimental (1981)',\
"./data/electron/Ar/DCS/McEachran.data" u 1:($7*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 3 t'McEachran',\
"./data/electron/Ar/DCS/Sienkiewicz.data" u 1:($7*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 4 t'Sienkiewicz',\
"./data/electron/Ar/DCS/Dasgupta.data" u 1:($5*0.529*0.529) w l lt 0 lw 2 dt "-.-" lc 5 t'Dasgupta'

set title " Differentiel Cross-Section at 30 eV"
plot "./output/DCS/Ar/DCS_30eV.out" u 1:($2*0.529*0.529) w l t'Simulation',\
"./data/electron/Ar/DCS/Panajotovic97.data" u 1:8:9 with yerrorbars t'Panajotovic - Experimental (1997)',\
"./data/electron/Ar/DCS/Nahar87.data" u 1:($7*0.529*0.529) with lines lt 0 lw 2 dt "-.-" lc 2 t'Nahar - Simulation (1987)',\
"./data/electron/Ar/DCS/Srivastava81.data" u 1:10:11 with yerrorbars t'Srivastava - Experimental (1981)'

unset multiplot