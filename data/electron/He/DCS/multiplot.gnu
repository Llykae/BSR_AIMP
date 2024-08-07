set title " "
set xlabel "Angle (degree)"
set ylabel "DCS (angström/sr^{-1})

set terminal jpeg size 1960,1200
#set terminal wxt size 2400,1000
set output "./output/DCS/multiplot_DCS_6L.jpeg"

set multiplot layout 2,2

set title " Differentiel Cross-Section at 2 eV"
plot "./output/DCS/He/DCS_2eV.out" u 1:($2*0.529*0.529) w l t'Simulation',\
"./data/electron/He/DCS/Andrick75.data" u 1:2:3 with yerrorbars t'Andrick - Experimental (1975)'

set title " Differentiel Cross-Section at 5 eV"
plot "./output/DCS/He/DCS_5eV.out" u 1:($2*0.529*0.529) w l t'Simulation',\
"./data/electron/He/DCS/Brunger92.data" u 1:($4*0.1):($5*0.1) with yerrorbars t'Brunger - Experimental (1992)',\
"./data/electron/He/DCS/Andrick75.data" u 1:4:5 with yerrorbars t'Andrick - Experimental (1975)',\
"./data/electron/He/DCS/Abidzadeh05.data" u 1:($4*0.529*0.529) with lines lt 0  dt "-" t'Abidzadeh - Simulation (2005)'

set title " Differentiel Cross-Section at 12 eV"
plot "./output/DCS/He/DCS_12eV.out" u 1:($2*0.529*0.529) w l t'Simulation',\
"./data/electron/He/DCS/Andrick75.data" u 1:6:7 with yerrorbars t'Andrick - Experimental (1975)'

set title " Differentiel Cross-Section at 19 eV"
plot "./output/DCS/He/DCS_19eV.out" u 1:($2*0.529*0.529) w l t'Simulation',\
"./data/electron/He/DCS/Andrick75.data" u 1:8:9 with yerrorbars t'Andrick - Experimental (1975)'

#show output
unset multiplot