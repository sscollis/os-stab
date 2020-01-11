set title "Sweep in Frequncy"
set xlabel "{/Symbol w}"
set ylabel "Re({/Symbol a})"
set y2label "Im({/Symbol a})"
set ytics nomirror
set y2tics auto
plot "fort.17" u 1:2 w l axis x1y1 t "Re({/Symbol a})"
replot "fort.17" u 1:3 w l axis x1y2 t "Im({/Symbol a})"
