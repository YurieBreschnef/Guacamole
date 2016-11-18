i = 0
load './gnuplot-palettes-master/jet.pal'    


do for [i=SINDEX:EINDEX] {
   unset cbrange

   set size ratio ASPRATIO
   set terminal pngcairo size 800,800 enhanced font 'Verdana,10'
   set output './visual/spec_inspect_combo/'.i.'.png'
   set multiplot layout 3,3
   # Axes
   set style line 11 lc rgb '#808080' lt 1
   set border 3 back ls 11
   set tics nomirror out scale 0.75

    unset cbrange
    set xrange [0.0:LX]
    set yrange [0.0:LY]
  	set title ' u_x '
   	plot './data/u/'.i.'.u.dat' using 1:2:3 with image notitle

    set xrange [0.0:XDIM]
    set yrange [0.0:YDIM]
     set cbrange[MINSPEC:MAXSPEC]
  	set title 'ln[Re(F(u_x))] '
   	plot './data/u_f/'.i.'.u_f.dat' using 1:2:3 with image notitle
  	set title 'ln[imag(F(u_x))] '
   	plot './data/u_f/'.i.'.u_f.dat' using 1:2:4 with image notitle

    unset cbrange
    set xrange [0.0:LX]
    set yrange [0.0:LY]
  	set title ' u_y '
   	plot './data/u/'.i.'.u.dat' using 1:2:4 with image notitle

    set xrange [0.0:XDIM]
    set yrange [0.0:YDIM]
    set cbrange[MINSPEC:MAXSPEC]
  	set title 'ln[Re(F(u_y))] '
   	plot './data/u_f/'.i.'.u_f.dat' using 1:2:5 with image notitle
  	set title 'ln[imag(F(u_y))] '
   	plot './data/u_f/'.i.'.u_f.dat' using 1:2:6 with image notitle

    unset cbrange
    set xrange [0.0:LX]
    set yrange [0.0:LY]
  	set title 'absolute magnitude  of velocity field'
   	plot './data/abs_u/'.i.'.abs_u.dat' using 1:2:3 with image notitle

    set xtics
    set ytics
    unset yrange
    set xrange [0.0:64.0]
    set title 'log(energy) vs wavenumber k of first 64 modes[arb]'
    plot './data/k_spec/'.i.'.spec.dat' using 2:3 with lines notitle
    unset multiplot
}
