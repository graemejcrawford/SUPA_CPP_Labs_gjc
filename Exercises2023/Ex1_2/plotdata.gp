#!/usr/bin/gnuplot

# Get the input data file
data_file = "input2D_float.txt"
data_file2 = "/workspaces/SUPA_CPP_Labs_gjc/Exercises2023/Ex1_2/Outputs/bestFit_output.txt"

set datafile separator ","

# Set plot properties
set title "Best Fit Plot"
set xlabel "X"
set ylabel "Y"
set style data points
set xrange [0:6]
set yrange [0:6]

set terminal pngcairo
set output "Outputs/png/plot_bestFit.png"

# Plot the data from the specified file, skipping the first line

plot data_file every ::1 using 1:2 with points ps 1 pt 7 lc black title "Dummy data", data_file2  w l title "Best Fit"
