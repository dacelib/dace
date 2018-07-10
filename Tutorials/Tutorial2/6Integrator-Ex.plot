set term pdf
set output 'ex6_1_2.pdf'

set xlabel 'time [s]'
set ylabel 'angle [rad]'
plot 'ex6_1_2.dat' using 1:2 i 0 with lines notitle
set ylabel 'angular velocity [rad/s]'
plot 'ex6_1_2.dat' using 1:3 i 0 with lines notitle
set ylabel 'angle [rad]'
plot 'ex6_1_2.dat' using 1:2 i 1 with lines title 'left boundary', 'ex6_1_2.dat' using 1:3 i 1 with lines title 'right boundary'
set ylabel 'angular velocity [rad/s]'
plot 'ex6_1_2.dat' using 1:4 i 1 with lines title 'left boundary', 'ex6_1_2.dat' using 1:5 i 1 with lines title 'right boundary'


set term pdf
set output 'ex6_1_3.pdf'

unset xlabel
unset ylabel
set key bmargin
set view 0, 0
set view equal xyz
splot for [i=0:6] 'ex6_1_3.dat' using 2:3:(0) index i with lines title 'snapshot '.i


set term pdf
set output 'ex6_1_5.pdf'

unset xlabel
unset ylabel
set size square
set key rmargin
plot [-2:2] [-2:2] \
     'ex6_1_5.dat' index 0 title 'Starting point', \
     'ex6_1_5.dat' index 1 with line title 'Final states'


set term pdf
set output 'ex6_2_2.pdf'

unset xlabel
unset ylabel
set key bmargin
set view 0, 0
set view equal xyz
splot [-2:2] [-2:2] for [i=0:6] 'ex6_2_2.dat' using 2:3:(0) index i with lines title 'snapshot '.i

do for [i=0:6] {
    splot [-2:2] [-2:2] 'ex6_2_2.dat' using 2:3:(0) index i with lines title 'snapshot '.i
}


set term pdf
set output 'ex6_2_4.pdf'

unset xlabel
unset ylabel
set key bmargin
set view 0, 0
set view equal xyz
splot for [i=0:6] 'ex6_2_4.dat' using 2:3:(0) index i with lines title 'snapshot '.i
