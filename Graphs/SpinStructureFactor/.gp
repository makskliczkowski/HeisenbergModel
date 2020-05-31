L = 8

set yrange [0:pi]

unset ylabel
unset xlabel 
set xtics () # clear ytics

set multiplot

unset border
unset xtics; unset ytics;
set lmargin at screen 0.15
set rmargin at screen 0.68
set bmargin at screen 0.15
set tmargin at screen 0.9
set view map
set pm3d map
set pm3d interpolate 0,0
set dgrid3d 2*(L+1),pi*10,1
unset key
plot "fff.png" binary filetype=png w rgbimage


set xrange [0:2*pi]
set border
set tics out nomirror scale 2
set mxtics 5

set ylabel '{/*1.3 {/Symbol w}/J}' offset -2.5,0 rotate by 360
set xlabel '{/*1.3 q/{/Symbol p} (a.u.)}' offset 0, 0.3
set cblabel '{/*1.3 S(q,{/Symbol w})}' offset -50,80.5 rotate by 360



set xtics () # clear ytics
set xtics add ("0" 0)

do for [i=1:L]{
set xtics add (sprintf("%0.2f", 2*(i+0.0)/(L+0.0)) 2*pi/L*i)
}

do for[i=4:L:2]{
set multiplot 
set title sprintf('{/*L = %d .1.2 Spin Structure Factor for T ->inf and {/Symbol D} = %0.0f}\ ',i, 1.0)
set dgrid3d ('2*(%d+1)',i), 10*3.14,1.5
#set output sprintf('TinfL%d.png',i)                           
titla = sprintf('T00_L%d_J1_d1.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 1.0, T = inf',i) with pm3d
unset multiplot
}

do for[i=4:L:2]{
set title sprintf('{/*L = %d .1.2 Spin Structure Factor for T ->0 and {/Symbol D} = %0.0f}\ ',i, 1.0)
set dgrid3d ('2*(%d+1)',i), 10*3.14,1.5
#set output sprintf('T0L%d.png',i)                           
titla = sprintf('Tt0_L%d_J1_d1.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 1.0,T = 0' ,i) with pm3d
}

do for[i=4:L:2]{
set title sprintf('{/*L = %d .1.2 Spin Structure Factor for T =1.0 and {/Symbol D} = %0.0f}\ ',i, 1.0)
set dgrid3d ('2*(%d+1)',i), 10*3.14,1.5
#set output sprintf('T1L%d.png',i)                           
titla = sprintf('Tt1_L%d_J1_d1.dat',i)              
splot titla using 1:2:3 title sprintf('L = %d, delta = 1.0, T = 1',i)  with pm3d
}

set key outside right bottom

plot pi/2*abs(sin(x)) w l lw 2 lc rgb "black" title "{/Symbol p}/2*|sin(q)|", pi*abs(sin(x/2)) w l lw 2 lc rgb "blue" title "{/Symbol p}*|sin(q/2)|"

unset multiplot

