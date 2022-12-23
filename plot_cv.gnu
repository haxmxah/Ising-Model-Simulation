set term png
set output "CV_008.png"
set term png size 600, 600 

set title font "Times New Roman,14"
set xlabel font "Times New Roman,14"
set ylabel font "Times New Roman,14"

set title "Gràfica de la magnetitzacio en funció de la temperatura per particula"
set xlabel "T*" 
set ylabel "Capacitat calorífica/N"
!set key outside 
set grid xtics
set grid ytics

N= 1024.
set key top left
#set yrange [0:2]
set xrange[1:3.5]

plot "MC-L-008-TEMP-0200.dat" u 2:13 w lp pointsize 1 t"Cv", "derivades_008.dat" u 1:2 w l
