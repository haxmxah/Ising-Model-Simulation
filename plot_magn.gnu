set term png
set output "MAGN_008.png"
set term png size 600, 600 

set title font "Times New Roman,14"
set xlabel font "Times New Roman,14"
set ylabel font "Times New Roman,14"

set title "Gràfica de la magnetitzacio en funció de la temperatura per partícula"
set xlabel "T" 
set ylabel "Magnetitzacio/N"
!set key outside 
set grid xtics
set grid ytics

N= 1024
set key left bottom
#set key outside
!set yrange [0:-2]
set xrange[1:3.5]

plot "MC-L-008-TEMP-0200.dat" u 2:(sqrt($10)) w p lt rgb "forest-green" t "<m^{2}>^{1/2}", "MC-L-008-TEMP-0200.dat" u 2:9:12 w yerrorbars lt rgb "sienna1" t "<|m|>"

