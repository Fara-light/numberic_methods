#!/bin/bash
./a.out
gnuplot -p -e "plot '.implicitDifferenceScheme' using 1:2 with lines" 2> /dev/null &
gnuplot -p -e "plot '.explicitDifferenceScheme' using 1:2 with lines" 2> /dev/null &
gnuplot -p -e "plot '.crankNicolsonDifferenceScheme' using 1:2 with lines" 2> /dev/null &

