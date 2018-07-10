deg(x)=x*180/pi


set term pdf
set output 'ex7_1_2.pdf'

set xlabel 'time [s]'
set ylabel 'angle [rad]'
plot 'ex7_1_2.dat' u 1:(deg($2)) w l t 'angle'
set ylabel 'force [N]'
plot 'ex7_1_2.dat' u 1:3 w l t 'control'


set term pdf
set output 'ex7_1_3.pdf'

set ylabel 'angle [rad]'
plot \
     'ex7_1_3.dat' u 1:(deg($3)):(deg($4)) w filledcurves lc 'gray' not,\
     'ex7_1_3.dat' u 1:(deg($5)) w l t 'angle (left boundary)',\
     'ex7_1_3.dat' u 1:(deg($6)) w l t 'angle (right boundary)',\
     'ex7_1_3.dat' u 1:(deg($2)) w l t 'angle (reference)'
set ylabel 'force [N]'
plot \
     'ex7_1_3.dat' u 1:8:9 w filledcurves lc 'gray' not,\
     'ex7_1_3.dat' u 1:10 w l t 'control (left boundary)',\
     'ex7_1_3.dat' u 1:11 w l t 'control (right boundary)',\
     'ex7_1_3.dat' u 1:7 w l t 'control (reference)'


set term pdf
set output 'ex7_2_2.pdf'

set ylabel 'angle [rad]'
plot \
     'ex7_2_2.dat' u 1:(deg($3)):(deg($4)) w filledcurves lc 'gray' not,\
     'ex7_2_2.dat' u 1:(deg($5)) w l t 'angle (left boundary)',\
     'ex7_2_2.dat' u 1:(deg($6)) w l t 'angle (right boundary)',\
     'ex7_2_2.dat' u 1:(deg($2)) w l t 'angle (reference)'
set ylabel 'force [N]'
plot \
     'ex7_2_2.dat' u 1:8:9 w filledcurves lc 'gray' not,\
     'ex7_2_2.dat' u 1:10 w l t 'control (left boundary)',\
     'ex7_2_2.dat' u 1:11 w l t 'control (right boundary)',\
     'ex7_2_2.dat' u 1:7 w l t 'control (reference)'

