set title "Streamwise Velocity Spatial Eigenfunction"
set xlabel "y"
set ylabel "u(y)"
plot for [i=2:9] "fort.15" u 1:i w l t sprintf("u(t_%d)",i-1)
