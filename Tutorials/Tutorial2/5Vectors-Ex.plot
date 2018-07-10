set term pdf
set output 'ex5_2_1.pdf'
#set term x11

set pm3d
set hidden3d
set view equal xyz
set title 'Solar flux'
splot 'ex5_2_1.dat' u 1:2:3:4 w l not

#pause -1
