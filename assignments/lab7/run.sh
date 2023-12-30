./a.out
gnuplot -p -e "plot 'leibmann_data.gr' using 1:2 with lines" 2> /dev/null &
gnuplot -p -e "plot 'relexations_data.gr' using 1:2 with lines" 2> /dev/null &
gnuplot -p -e "plot 'seidel_data.gr' using 1:2 with lines" 2> /dev/null &
