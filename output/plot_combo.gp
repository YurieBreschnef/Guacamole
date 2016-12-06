i = SINDEX
load './output/gnuplot-palettes-master/jet.pal'    


do for [i=SINDEX:EINDEX] {
   unset cbrange
   set xrange [0:LX]
   set yrange [0:LY]

   set size ratio ASPRATIO
   set terminal pngcairo size 1600,800 enhanced font 'Verdana,10'
   set output './output/visual/combo/'.i.'.png'
   set multiplot layout 2,4
   # Axes
   set style line 11 lc rgb '#808080' lt 1
   set border 3 back ls 11
   set tics nomirror out scale 0.75

   set title 'temperature field'
   plot './output/data/temp/'.i.'.temp.dat' using 1:2:3  with image notitle

   set title 'chemical field'
   plot './output/data/chem/'.i.'.chem.dat' using 1:2:3  with image notitle

   set title 'absolute of u '
   plot './output/data/abs_u/'.i.'.abs_u.dat' using 1:2:3  with image notitle

   set title 'z component of u '
   plot './output/data/u/'.i.'.u.dat' using 1:2:4  with image notitle

   set title 'buoyancy field (B_therm*T_z - B_comp * S_z)'
   plot './output/data/buo/'.i.'.buo.dat' using 1:2:3  with image notitle

   set title 'velocity field'
   plot './output/data/u/'.i.'.u.dat' using 1:2:3:4  with vectors

   set title 'vorticity field'
   plot './output/data/vort/'.i.'.vort.dat' using 1:2:3  with image notitle

   unset xrange 
   unset yrange 
   set cbrange[MINSPEC:MAXSPEC]
   set title 'fourier spectrum of chem field'
   plot './output/data/chem_f/'.i.'.chem_f.dat' using 1:2:3 with image notitle
   unset multiplot
}
