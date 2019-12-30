plot "fort.14" u 1:2 w l t "Re(d3phi/dx3)", "fort.14" u 1:3 w l t "Im(d3phi/dx3)"
