reset
set grid xtics
set grid ytics
set grid ztics

set encoding utf8
set terminal pngcairo size 1080,1080 font 'Verdana,13  '
set view map
set ticslevel 0
set title "Spin structure factor Fourier Transfrom Module for h=<Ediff>_N"
set zlabel "Re(S(q,w)_FT)" font "Verdana,14 "
set xlabel "q[in pi units]"
set ylabel "omega[in pi units]" offset -3, 0
set palette model RGB


set key box opaque 
set pm3d map
set pm3d interpolate 0,0
set autoscale fix
set isosamples 2
set colorbox bdefault 
set format z "%.1e"
set format y "%.2f"
set format x "%.2f"


do for[i=12:12:2]{
set dgrid3d ('2*(%d+1)',i), 10*3.14,1.5
set output sprintf('SqwFourierL%dLongTime.png',i)                           
titla = sprintf('%d_d1_spinMapFourier.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 1.0, h = 0.01',i)
}
