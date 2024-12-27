#!/bin/gnuplot

set terminal qt size 900,800

set xrange [0:2000]
set yrange [0:2000]
set size square

unset key

set object 2 ellipse at 1000,1000 size 200,200 lw 1 fs empty border rgb "blue" front

set logscale cb
#set logscale z

plot "out.txt" matrix with image, "out_2.txt" using ($1 + 1000):($2 + 1000):3 with points lc rgb "black"
#splot "out.txt" matrix u 1:2:3, "out_2.txt" using ($1 + 1000):($2 + 1000):3 with points lc rgb "black"

#splot "out.txt" matrix

pause mouse close
