set title "Eigenfunction"
set xlabel "y"
set ylabel "phi"
plot "fort.11" u 1:2 w l t "Re(phi)", "fort.11" u 1:3 w l t "Im(phi)"
#replot "fort.12" u 1:2 w l, "fort.12" u 1:3 w l
#replot "fort.13" u 1:2 w l, "fort.13" u 1:3 w l
#replot "fort.14" u 1:2 w l, "fort.14" u 1:3 w l
