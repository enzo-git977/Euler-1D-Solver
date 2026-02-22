set terminal pngcairo enhanced font 'arial,10' size 1200, 800
set output 'data_final/Solution_Exact_et_apprcohe.png'

# Set time step to Tf
t_final = 0.20

# Generate the filename for the final time
filename_final = sprintf("data_final/results_at_t=%.2fs.dat", t_final)

set multiplot layout 3,2 title "Exact vs Approximate Riemann solver" 

# Plot Density (rho)
set xlabel 'x (m)'
set ylabel 'Density '
set title 'Density Comparison'
plot 'data_final/Exact_Sod.dat' using 1:2 with lines lw 2 title 'Exact', \
     filename_final using 1:2 with lines lw 2 title sprintf('t=%.2fs', t_final)

# Plot Mach number (Ma)
set xlabel 'x (m)'
set ylabel 'Mach number (-)'
set title 'Mach Comparison'
plot 'data_final/Exact_Sod.dat' using 1:3 with lines lw 2 title 'Exact', \
     filename_final using 1:3 with lines lw 2 title sprintf('t=%.2fs', t_final)

# Plot Pressure (P)
set xlabel 'x (m)'
set ylabel 'Pressure '
set title 'Pressure Comparison'
plot 'data_final/Exact_Sod.dat' using 1:4 with lines lw 2 title 'Exact', \
     filename_final using 1:4 with lines lw 2 title sprintf('t=%.2fs', t_final)
	
# Plot Temperature adimensionnee (-)
set xlabel 'x (m)'
set ylabel 'Temperature (-)'
set title 'T adimensionne Comparison'
plot 'data_final/Exact_Sod.dat' using 1:5 with lines lw 2 title 'Exact', \
     filename_final using 1:5 with lines lw 2 title sprintf('t=%.2fs', t_final)

# Plot Energy (E)
set xlabel 'x (m)'
set ylabel 'Energy '
set title 'Energy Comparison'
plot 'data_final/Exact_Sod.dat' using 1:6 with lines lw 2 title 'Exact', \
     filename_final using 1:6 with lines lw 2 title sprintf('t=%.2fs', t_final)
	 
# Plot Entropy (s)
set xlabel 'x (m)'
set ylabel 's '
set title 'Entropy Comparison'
plot 'data_final/Exact_Sod.dat' using 1:7 with lines lw 2 title 'Exact', \
     filename_final using 1:7 with lines lw 2 title sprintf('t=%.2fs', t_final)

unset multiplot
