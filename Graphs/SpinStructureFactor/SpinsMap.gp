reset
set grid xtics
set grid ytics
set grid ztics
set encoding utf8
system('mkdir -p png')
set terminal pngcairo size 1920,1080 font 'Verdana,13  '
set view map
set ticslevel 0
set title "Spin structure factor"
set zlabel "S(q,w)" font "Verdana,14 "
set xlabel "q[in pi units]"
set ylabel "omega[in pi units]"


set key box opaque 
set pm3d map
set pm3d interpolate 0,0

unset margin
set autoscale fix

do for[i=4:8:2]{
set output sprintf('TinfL%d.png',i)                           
titla = sprintf('T00_L%d_J1_d1.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 1.0, T = inf',i)
}
do for[i=4:8:2]{
set output sprintf('T0L%d.png',i)                           
titla = sprintf('Tt0_L%d_J1_d1.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 1.0,T = 0' ,i)
}
do for[i=4:8:2]{
set output sprintf('T1L%d.png',i)                           
titla = sprintf('Tt1_L%d_J1_d1.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 1.0, T = 1',i)
}

