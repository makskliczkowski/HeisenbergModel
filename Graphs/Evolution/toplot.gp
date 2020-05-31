
cd 'D:\Uni\SEMESTERS\MS\I\Many Body\LAB\Hamiltonian\Graphs\Evolution'

set title "Mean energy"
set ylabel "<E>_T" font "Helvetica,14"
set xlabel "Temperature" font "Helvetica,14"
set term gif animate
set output "meanEnergy.gif"
do for[i=2:16:2]{                                   
titla = sprintf('%d_d0_meanEnergy.dat',i)              
plot titla title sprintf('L=%d, delta= 0.0',i)
}
do for[i=2:16:2]{
titla = sprintf('%d_d1_meanEnergy.dat',i)              
plot titla title sprintf('L=%d, delta= 1.0',i)
}
do for[i=2:16:2]{
titla = sprintf('%d_d2_meanEnergy.dat',i)              
plot titla title sprintf('L=%d, delta= 2.0',i)
}

