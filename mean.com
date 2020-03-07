set xlabel "y"
set ylabel "Profile"
set title "Mean flow profile"
plot "fort.10" u 1:2 w l t "u", "fort.10" u 1:3 w l t "u''"
