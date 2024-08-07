set title " Scattering Simulation for e^{-} - Ar"

set terminal jpeg size 1960,1200 font "{/:Bold} Helvetica,16"
set output "./output/TCS/Ar/Ar_multiplot_TCS.jpeg"
set key font ",10"

set multiplot layout 2,3

set title "S-wave Phase-Shift From 0 to 16 eV"
set xrange[0:16]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):2 w l t'Present Work',\
"./data/electron/Ar/TCS/McEachran97.dat" u ($1*$1*13.605):3 t'McEachran - Simulation (1997)',\
"./data/electron/Ar/TCS/Sienkiewicz87.dat" u 1:2 t'Sienkiewicz - Simulation (1987)',\
"./data/electron/Ar/TCS/Williams79.dat" u 1:2:3 with yerrorbars pointsize 0.5 t'Williams - Experimental (1979)' ,\
"./data/electron/Ar/TCS/ps_mixed.dat" u ($1*27.21):2 t'Best Fit',\
"./data/electron/Ar/TCS/Ref-McEachran.dat" u ($2*27.21):3 t'Ref - McEachran'

set title "P-wave Phase-Shift From 0 to 20 eV"
set xrange[0:20]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):3 w l t'Present Work',\
"./data/electron/Ar/TCS/McEachran97.dat" u ($1*$1*13.605):4 t'McEachran - Simulation (1997)',\
"./data/electron/Ar/TCS/Sienkiewicz87.dat" u 1:3 t'Sienkiewicz - Simulation (1987)',\
"./data/electron/Ar/TCS/Williams79.dat" u 1:4:5 with yerrorbars pointsize 0.5 t'Williams - Experimental (1979)',\
"./data/electron/Ar/TCS/ps_mixed.dat" u ($1*27.21):3 t'Best Fit',\
"./data/electron/Ar/TCS/Ref-McEachran.dat" u ($2*27.21):4 t'Ref - McEachran'


set title "D-wave Phase-Shift From 0 to 16 eV"
set xrange[0:16]
set yrange[-0.5:2.0]
set xlabel "Energy (eV)"
set ylabel "Phase-Shift (in radian)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):4 w l t'Present Work',\
"./data/electron/Ar/TCS/McEachran97.dat" u ($1*$1*13.605):5 t'McEachran - Simulation (1997)',\
"./data/electron/Ar/TCS/Sienkiewicz87.dat" u 1:4 t'Sienkiewicz - Simulation (1987)',\
"./data/electron/Ar/TCS/Williams79.dat" u 1:6:7 with yerrorbars pointsize 0.5 t'Williams - Experimental (1979)',\
"./data/electron/Ar/TCS/ps_mixed.dat" u ($1*27.21):4 t'Best Fit',\
"./data/electron/Ar/TCS/Ref-McEachran.dat" u ($2*27.21):5 t'Ref - McEachran'

set title "Total Cross-Section below 3 eV"
set logscale x
set xrange[0.01:3]
set logscale y
set yrange[0.5:20]
set xlabel "Energy (eV)"
set ylabel "TCS (in a_0^2)" 

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):7 w l t'Present Work',\
"./data/electron/Ar/TCS/Kurokawa2011.data" u 1:($2/(0.529*0.529)):(0.01*$3/(0.529*0.529)) with yerrorbars  t'Kurokawa - Experimental (2011)',\
"./data/electron/Ar/TCS/McEachran97.dat" u ($1*$1*13.605):($9/(0.529*0.529)) t'McEachran - Simulation (1997)',\
"./data/electron/Ar/TCS/Sienkiewicz87.dat" u 1:6 t'Sienkiewicz - Simulation (1987)',\
"./data/electron/Ar/TCS/ps_mixed.dat" u ($1*27.21):5 t'Best Fit',\
"./data/electron/Ar/TCS/Ref-McEachran.dat" u ($2*27.21):6 t'Ref - McEachran'


set title "Total Cross-Section over 3 ev"
set logscale x
set xrange[3:70]
set logscale y
set yrange[20:100]
set xlabel "Energy (eV)"
set ylabel "TCS (in a_0^2)"

plot "./output/TCS/Scatteringtools.out" u ($1*27.21):7 w l t'Present Work',\
"./data/electron/Ar/TCS/Kurokawa2011.data" u 1:($2/(0.529*0.529)):(0.01*$3/(0.529*0.529)) with yerrorbars  t'Kurokawa - Experimental (2011)',\
"./data/electron/Ar/TCS/McEachran97.dat" u ($1*$1*13.605):($9/(0.529*0.529)) t'McEachran - Simulation (1997)',\
"./data/electron/Ar/TCS/Sienkiewicz87.dat" u 1:6 t'Sienkiewicz - Simulation (1987)',\
"./data/electron/Ar/TCS/ps_mixed.dat" u ($1*27.21):6 t'Best Fit',\
"./data/electron/Ar/TCS/Ref-McEachran.dat" u ($2*27.21):6 t'Ref - McEachran'


unset logscale
unset multiplot
