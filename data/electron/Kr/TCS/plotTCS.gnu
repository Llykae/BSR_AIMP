set title "TCS for electron scattering from Kr"
set xlabel "Energy (in eV)"
set ylabel "Total Cross Section"
set xrange[0:50]

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):6 w l t'Gervais - Present Work 2022'
replot "./data/electron/Kr/TCS/Kurokawa2010.data" u 1:($2/(0.529*0.529)) t"Kurokawa (Experimental 2010)"
replot "./data/electron/Kr/TCS/McEachran87.data" u ($1*$1*13.605):((4*3.141592/($1*$1))*(sin($2)*sin($2)+3*sin($3)*sin($3)+5*sin($4)*sin($4)+7*sin($5)*sin($5))) t'McEachran (SIMU 1987)'
replot "./data/electron/Kr/TCS/Dabaneh80.data" u 1:($2/(0.529*0.529)) t"Dabaneh (Experimental 1980)"
