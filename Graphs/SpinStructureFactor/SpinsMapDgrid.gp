 cd 'D:\Uni\SEMESTERS\MS\I\Many Body\LAB\Hamiltonian\Graphs\SpinStructureFactor' 
f(x)= pi/2*abs(sin(x))
g(x) = pi*abs(sin(x/2))
reset
set grid xtics
set grid ytics
set grid ztics
L = 14

set encoding utf8
set terminal pngcairo size 1080,1080 font 'Verdana,13  '
set view map
set ticslevel 0
set title "Spin structure factor"
#set zlabel "S(q,w)" font "Verdana,14 "
#set xlabel "q[in pi units]"

#set ylabel "omega[in pi units]" offset -3, 0
set palette model RGB
set lmargin at screen 0.2
set rmargin at screen 0.8
set bmargin at screen 0.2
set tmargin at screen 0.8
unset border

set key box opaque 
set pm3d map
set pm3d interpolate 0,0
set autoscale fix
set isosamples 2
set colorbox bdefault 
set xrange [0:2]
set border
set tics out nomirror scale 2
set mxtics 5

set ylabel '{/*1.3 {/Symbol w} in {/Symbol p} units' offset -2.5,0 rotate by 90
set xlabel '{/*1.3 q/{/Symbol p} (a.u.)}' offset 0, 0.3
set cblabel '{/*1.3 S(q,{/Symbol w})}' offset -50,80.5 rotate by 360


do for[i=10:L:2]{

set dgrid3d ('2*(%d+1)',i), 10*3.14,1
set output sprintf('TinfL%d.png',i)                       
titla = sprintf('T00_L%d_J1_d1.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 1.0, T = inf',i)#, f(x) w l lw 2 lc rgb "black" , g(x) w l lw 2 lc rgb "blue"
}
do for[i=10:L:2]{

set dgrid3d ('2*(%d+1)',i), 10*3.14,1
set output sprintf('T0L%d.png',i)      
                 
titla = sprintf('Tt0_L%d_J1_d1.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 1.0,T = 0' ,i)#,  f(x) w l lw 2 lc rgb "black" , g(x) w l lw 2 lc rgb "blue"

}

do for[i=10:L:2]{
set dgrid3d ('2*(%d+1)',i), 10*3.14,1
set output sprintf('T1L%d.png',i)              
    
titla = sprintf('Tt1_L%d_J1_d1.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 1.0, T = 1',i)#, f(x) w l lw 2 lc rgb "black" , g(x) w l lw 2 lc rgb "blue"

}

