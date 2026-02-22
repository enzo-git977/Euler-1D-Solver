set terminal pngcairo size 1200,800 enhanced font "Arial,12"
set output 'data_final/results_t=0.png'

# Enable multiplot mode
set multiplot layout 2,2 title "Initial Solution (t=0)" font ",14"

# Plot Density (rho)
set xlabel "Position x"
set ylabel "Density (rho)"
plot "data_final/results_at_0.00s.dat" using 1:2 with lines lw 2 title "Density (rho)"

# Plot Velocity (u)
set xlabel "Position x"
set ylabel "Velocity (u)"
plot "data_final/results_at_0.00s.dat" using 1:3 with lines lw 2 title "Velocity (u)"

# Plot Pressure (P)
set xlabel "Position x"
set ylabel "Pressure (P)"
plot "data_final/results_at_0.00s.dat" using 1:4 with lines lw 2 title "Pressure (P)"

# Plot Energy (E)
set xlabel "Position x"
set ylabel "Energy (E)"
plot "data_final/results_at_0.00s.dat" using 1:5 with lines lw 2 title "Energy (E)"

# Disable multiplot mode
unset multiplot