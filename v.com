set title "Wall-normal Velocity Spatial Eigenfunction"
set xlabel "y"
set ylabel "v(y)"
plot for [i=2:9] "fort.16" u 1:i w l t sprintf("v(t_%d)",i-1)
