set title "Differential Cross-Section for electron scattering from Ar at 20eV"
set xlabel "Angle (degree)"
set ylabel "DCS (in angström².sr^{-1})"

plot "./output/DCS/Ar/DCS_20eV.out" u 1:($2*0.529*0.529) w l t'Simulation'
replot "./data/electron/Ar/DCS/Srivastava81.data" u 1:8:9 with yerrorbars t'Srivastava - Experimental (1981)'