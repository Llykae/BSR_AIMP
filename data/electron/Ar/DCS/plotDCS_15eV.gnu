set title "Differential Cross-Section for electron scattering from Ar at 15eV"
set xlabel "Angle (degree)"
set ylabel "DCS (in angström².sr^{-1})"

plot "./output/DCS/Ar/DCS_15eV.out" u 1:($2*0.529*0.529) w l t'Simulation'
replot "./data/electron/Ar/DCS/Srivastava81.data" u 1:6:7 with yerrorbars t'Srivastava - Experimental (1981)'