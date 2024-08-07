set title "Differential Cross-Section for electron scattering from He at 18eV"
set xlabel "Angle (degree)"
set ylabel "DCS (in angström/sr^{-1})"

plot "./output/DCS/He/DCS_18eV.out" u 1:($2*0.529*0.529) w l t'Simulation'
replot "./data/electron/He/DCS/Brunger92.data" u 1:($8*0.1):($9*0.1) with yerrorbars t'Brunger - Experimental (1992)'