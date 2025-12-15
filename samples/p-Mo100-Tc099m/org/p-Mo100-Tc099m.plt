set title "p + Mo100 --> Tc99m" font ",25"
set xlabel "Irradiation time [h]" font ",30"
set ylabel "Yield [Gbq/mA.h]" font ",30" offset -2
set term postscript enhanced font "Times-Roman,16" color
set output "Tc099m.eps"
plot [0:48][0:40] \
 "Tc099m.act" u 1:2 with lines t "Tc-99m"
