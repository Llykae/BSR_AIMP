set title "Partial TCS ratio for different l-wave"
set xlabel "Energy (eV)"
set ylabel "Ratio"

set xrange[0:20]
set yrange[0:1.05]

plot "./output/TCS/extraction/extraction_PS.out" u ($1*27.21):2 w l t"S-wave"
replot "./output/TCS/extraction/extraction_PS.out" u ($1*27.21):3 w l t"P-wave"
replot "./output/TCS/extraction/extraction_PS.out" u ($1*27.21):4 w l t"D-wave"
replot "./output/TCS/extraction/extraction_PS.out" u ($1*27.21):5 w l t"F-wave"
replot "./output/TCS/extraction/extraction_PS.out" u ($1*27.21):($2+$3+$4+$5) w l t"Norme"