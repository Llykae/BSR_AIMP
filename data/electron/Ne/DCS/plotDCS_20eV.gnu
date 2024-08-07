set title "Differential Cross-Section for electron scattering from Ne at 20eV"
set xlabel "Angle (degree)"
set ylabel "DCS (in angström².sr^{-1})"

plot "./output/DCS/Ne/DCS_20eV.out" u 1:($2*0.529*0.529) w l t'Simulation'
replot "./data/electron/Ne/DCS/Khandker20.data" u 1:($5*0.529*0.529) t'Khandker - Simulation (2020)'
replot "./data/electron/Ne/DCS/Cho08_20eV.data" u 1:2:3 with yerrorbars t'Cho - Experimental (2008)'
replot "./data/electron/Ne/DCS/Registrer84_20eV.data" u 1:2:($2*0.07) with yerrorbars t'Register - Experimental (1984)'