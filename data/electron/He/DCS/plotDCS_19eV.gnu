set title "Differential Cross-Section for electron scattering from He at 19eV"
set xlabel "Angle (degree)"
set ylabel "DCS (in angström/sr^{-1})"

plot "./output/DCS/He/DCS_19eV.out" u 1:($2*0.529*0.529) w l t'Simulation'
replot "./data/electron/He/DCS/Andrick75.data" u 1:8:9 with yerrorbars t'Andrick - Experimental (1975)'