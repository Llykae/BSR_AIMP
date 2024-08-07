set title "Differential Cross-Section for electron scattering from He at 12eV"
set xlabel "Angle (degree)"
set ylabel "DCS (in angström/sr^{-1})"

plot "./output/DCS/He/DCS_12eV.out" u 1:($2*0.529*0.529) w l t'Simulation'
replot "./data/electron/He/DCS/Andrick75.data" u 1:6:7 with yerrorbars t'Andrick - Experimental (1975)'