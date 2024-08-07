set title "Differential Cross-Section for electron scattering from He at 10eV"
set xlabel "Angle (degree)"
set ylabel "DCS (in angström/sr^{-1})"

plot "./output/DCS/He/DCS_10eV.out" u 1:($2*0.529*0.529) w l t'Simulation'
replot "./data/electron/He/DCS/Brunger92.data" u 1:($6*0.1):($7*0.1) with yerrorbars t'Brunger - Experimental (1992)'