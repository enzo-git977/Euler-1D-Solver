set terminal pngcairo enhanced font 'arial,10' size 1200, 800
set output 'data_final/Solution_Exact_et_approche.png'
 
t_final =   0.0120
filename_final = sprintf('data_final/results_at_t=%.3fs.dat', t_final)
 
set multiplot layout 3,2 title "Exact vs Approximate Riemann solver"
 
set xlabel 'x (m)'
set ylabel 'Density '
set title 'Density Comparison'
plot 'data_final/Exact_Sod.dat' using 1:2 with lines lw 2 title 'Exact', \
     filename_final using 1:2 with lines lw 2 title sprintf('t=%.3fs', t_final)
set xlabel 'x (m)'
set ylabel 'Mach number (-)'
set title 'Mach Comparison'
plot 'data_final/Exact_Sod.dat' using 1:3 with lines lw 2 title 'Exact', \
     filename_final using 1:3 with lines lw 2 title sprintf('t=%.3fs', t_final)
set xlabel 'x (m)'
set ylabel 'Pressure '
set title 'Pressure Comparison'
plot 'data_final/Exact_Sod.dat' using 1:4 with lines lw 2 title 'Exact', \
     filename_final using 1:4 with lines lw 2 title sprintf('t=%.3fs', t_final)
set xlabel 'x (m)'
set ylabel 'Temperature (-)'
set title 'T adimensionne Comparison'
plot 'data_final/Exact_Sod.dat' using 1:5 with lines lw 2 title 'Exact', \
     filename_final using 1:5 with lines lw 2 title sprintf('t=%.3fs', t_final)
set xlabel 'x (m)'
set ylabel 'velocity '
set title 'velocity Comparison'
plot 'data_final/Exact_Sod.dat' using 1:6 with lines lw 2 title 'Exact', \
     filename_final using 1:6 with lines lw 2 title sprintf('t=%.3fs', t_final)
set xlabel 'x (m)'
set ylabel 's '
set title 'Entropy Comparison'
plot 'data_final/Exact_Sod.dat' using 1:7 with lines lw 2 title 'Exact', \
     filename_final using 1:7 with lines lw 2 title sprintf('t=%.3fs', t_final)
unset multiplot
