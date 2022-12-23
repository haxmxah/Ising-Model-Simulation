set term png
set output "ENERG_008.png"
set term png size 600, 600 

set title font "Times New Roman,14"
set xlabel font "Times New Roman,14"
set ylabel font "Times New Roman,14"

set title "Gràfica de l'energia en funció de la temperatura per particula"
set xlabel "T" 
set ylabel "Energia/N"
!set key outside 
set grid xtics
set grid ytics

set key top left
!set yrange [0:-2]
set xrange[1:3.5]

plot "MC-L-008-TEMP-0200.dat" u 2:4 w p pointtype 7 pointsize 0.5 lt rgb "violet"
#:7 w yerrorbars lt rgb "pink" t "Energia"

