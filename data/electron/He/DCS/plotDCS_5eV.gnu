set title "Differential Cross-Section for electron scattering from He at 5eV"
set xlabel "Angle (degree)"
set ylabel "DCS (in angström/sr^{-1})"

plot "./output/DCS/He/DCS_5eV.out" u 1:($2*0.529*0.529) w l t'Simulation'
replot "./data/electron/He/DCS/Brunger92.data" u 1:($4*0.1):($5*0.1) with yerrorbars t'Brunger - Experimental (1992)'
replot "./data/electron/He/DCS/Andrick75.data" u 1:4:5 with yerrorbars t'Andrick - Experimental (1975)'