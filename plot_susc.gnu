set term png
set output "XI_008.png"
set term png size 600, 600 

set title font "Times New Roman,14"
set xlabel font "Times New Roman,14"
set ylabel font "Times New Roman,14"

set title "Gràfica de la susceptibilitat magnètica en funció de la temperatura per partícula"
set xlabel "T*" 
set ylabel "Capacitat calorífica/N"
!set key outside 
set grid xtics
set grid ytics

N= 1024.
set key top left
#set yrange [0:20]
set xrange[1:3.5]

plot "MC-L-008-TEMP-0200.dat" u 2:14 w p t"X"
