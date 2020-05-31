cd 'D:\Uni\SEMESTERS\MS\I\Many Body\LAB\Hamiltonian\Graphs\Evolution'
set grid xtics
set grid ytics
set grid ztics
system('mkdir -p png')
set terminal pngcairo size 1920,1080 font 'Verdana,13  '
set title "<psi(t)|phi_T> for Tinit -> inf"
set xlabel "time"
set ylabel "overlapping" offset -1,0 font "Verdana,14 "

do for[i=2:16:2]{
set output sprintf('d0OverlapL%d.png',i)                           
titla = sprintf('%d_d0_timeEvo.dat',i)              
plot titla using 1:2 title sprintf('Re(), L = %d, delta = 0.0',i),titla using 1:3 title sprintf('Im(), L = %d, delta = 0.0',i)
}
