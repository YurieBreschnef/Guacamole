#load './gnuplot-palettes-master/spectral.pal'    
#load './gnuplot-palettes-master/blues.pal'    
load './gnuplot-palettes-master/jet.pal'    

aspect_ratio = 1
Lx = 8.0*3.14159 
Ly = 8.0*3.14159  
xdim = 256 
ydim = 256 


max_spec_order = 0
min_spec_order = -25

no_of_img = 299 

    ####################STATISTICS###############################
    set terminal pngcairo size 2200,1200 enhanced font 'Verdana,10'

    #set output './visual/stat/Energies.png'
    #set title 'simulation time [arb] vs. E_kin [arb]'
    #plot './data/E_stat/E_stat.dat' using 2:3 title "E_{kin}"   ,\
    #     './data/E_stat/E_stat.dat' using 2:4 title "E_{pot}"   ,\
    #     './data/E_stat/E_stat.dat' using 2:5 title "E_{tot}"  
    #
    #set output './visual/stat/Velocities.png'
    #set title 'simulation [arb] time vs. velocities [arb]'
    #plot './data/u_stat/u_stat.dat' using 2:3 title "v_{max}"   ,\
    #     './data/u_stat/u_stat.dat' using 2:4 title "v_{rms}"  
    #
    #
    #set output './visual/stat/Temp.png'
    #set title 'simulation time vs.Temperature measures [arb]'
    #plot './data/T_stat/T_stat.dat' using 2:3 title "max temp"  , \
    #     './data/T_stat/T_stat.dat' using 2:4 title "mean temp"   
    #
    #set output './visual/stat/Chem.png'
    #set title 'simulation time vs. Chemical field measures [arb]'
    #plot './data/C_stat/C_stat.dat' using 2:3 title "max Chem"  , \
    #     './data/C_stat/C_stat.dat' using 2:4 title "mean Chem"   
    #set output './visual/stat/Chem.png'
    #set title 'simulation time vs. Chemical field measures [arb]'
    #plot './data/C_stat/C_stat.dat' using 2:3 title "max Chem"  , \
    #     './data/C_stat/C_stat.dat' using 2:4 title "mean Chem"   
    #########Multiplot###############
    set output './visual/stat/stat_combo.png'
    set multiplot layout 4,4
    #set xrange [0.0:20.00]

        set title 'simulation time vs.Temperature measures [arb]'
        plot './data/T_stat/T_stat.dat' using 2:3 with lines lw 3 title "max temp"  , \
             './data/T_stat/T_stat.dat' using 2:4 with lines lw 3 title "min temp"  
        set title 'simulation time vs. statistical Temperature measures [arb]'
        plot './data/T_stat/T_stat.dat' using 2:5 with lines lw 3 title "mean temp"    ,\
             './data/T_stat/T_stat.dat' using 2:6 with lines lw 3 title "std-dev temp" , \
             './data/T_stat/T_stat.dat' using 2:7 with lines lw 3 title "variance temp"
        set title 'simulation time vs. Chemical field measures [arb]'
        plot './data/C_stat/C_stat.dat' using 2:3 with lines lw 3 title "max Chem"  , \
             './data/C_stat/C_stat.dat' using 2:4 with lines lw 3 title "min Chem"  
        set title 'simulation time vs. statistical Chemical field measures [arb]'
        plot './data/C_stat/C_stat.dat' using 2:5 with lines lw 3 title "mean Chem" , \
             './data/C_stat/C_stat.dat' using 2:6 with lines lw 3 title "std-dev chem" , \
             './data/C_stat/C_stat.dat' using 2:7 with lines lw 3 title "variance chem"

        set title 'simulation [arb] time vs. vx measures [arb]'
        plot './data/u_stat/u_stat.dat' using 2:3  with lines lw 3 title "vx_{max}"   ,\
             './data/u_stat/u_stat.dat' using 2:4  with lines lw 3 title "vx_{min}"   
        set title 'simulation [arb] time vs. statistical vx measures [arb]'
        plot './data/u_stat/u_stat.dat' using 2:5  with lines lw 3 title "vx_{avg}"   ,\
             './data/u_stat/u_stat.dat' using 2:6  with lines lw 3 title "std-dev vx"   ,\
             './data/u_stat/u_stat.dat' using 2:7  with lines lw 3 title "variance vx" 
        set title 'simulation [arb] time vs. vy measures [arb]'
        plot './data/u_stat/u_stat.dat' using 2:8  with lines lw 3 title "vy_{max}"   ,\
             './data/u_stat/u_stat.dat' using 2:9  with lines lw 3 title "vy_{min}"   
        set title 'simulation [arb] time vs. statistical vy measures [arb]'
        plot './data/u_stat/u_stat.dat' using 2:10  with lines lw 3 title "vy_{avg}"   ,\
             './data/u_stat/u_stat.dat' using 2:11  with lines lw 3 title "std-dev vy"   ,\
             './data/u_stat/u_stat.dat' using 2:12  with lines lw 3 title "variance vy" 

        set title 'simulation [arb] time vs. abs velocities (norm) [arb]'
        plot './data/u_stat/u_stat.dat' using 2:13 with lines lw 3 title "v_{max}"   ,\
             './data/u_stat/u_stat.dat' using 2:14 with lines lw 3 title "v_{min}"   
        set title 'simulation [arb] time vs. statistical abs velocities (norm) [arb]'
        plot './data/u_stat/u_stat.dat' using 2:15 with lines lw 3 title "v_{rms}"   ,\
             './data/u_stat/u_stat.dat' using 2:16 with lines lw 3 title "std-dev v"   ,\
             './data/u_stat/u_stat.dat' using 2:17 with lines lw 3 title "variance v"  
        set title 'simulation time [arb] vs. log(E_kin [arb])'
        plot './data/E_stat/E_stat.dat' using 2:(log($3)) with lines lw 3 title "log(E_{kin})"   ,\
             './data/E_stat/E_stat.dat' using 2:(log($4)) with lines lw 3 title "log(E_{pot})"   ,\
             './data/E_stat/E_stat.dat' using 2:(log($5)) with lines lw 3 title "log(E_{tot})"  
        set title 'simulation time vs. shearstrength '
        plot './data/sys_stat/sys_stat.dat' using 2:4 with lines title "shear strength [arb]"

        set title 'aperiodicity measure in x and y-dir vs. t [arb] '
        plot './data/sys_stat/sys_stat.dat' using 2:8 with lines lw 3 title "y_aperiodicity [arb]", \
             './data/sys_stat/sys_stat.dat' using 2:9 with lines lw 3 title "x_aperiodicity [arb]"
        set title 'simulation time vs. stepwidth dt [arb] '
        plot './data/sys_stat/sys_stat.dat' using 2:5 with lines lw 3 title "dt [arb]"
        set title 'simulation time vs. average vorticity [arb] '
        plot './data/sys_stat/sys_stat.dat' using 2:7 title "average vort [arb]"
        set title 'simulation time vs. average divergence [arb] '
        plot './data/sys_stat/sys_stat.dat' using 2:6 title "max brucker div [arb]" # ,\
 #            './data/sys_stat/sys_stat.dat' using 2:3 title "max div [arb]"


    unset multiplot
    ##########################################################################
    #FU PLOTS
    set terminal pngcairo size 1200,600 enhanced font 'Verdana,10'
    set output './visual/stat/fu_stat_combo.png'
    set multiplot layout 3,1
    #set xrange [0.0:15.00]
        set title 'simulation time [arb] vs. function-terms[arb] (du/dt = adv+diff+buo+shear) '
        plot './data/fu_stat/fu_stat.dat' using 2:5  with lines lw 3 title "avg fu"       ,\
             './data/fu_stat/fu_stat.dat' using 2:10 with lines lw 3 title "avg fu_Nuk (advection)"   ,\
             './data/fu_stat/fu_stat.dat' using 2:15 with lines lw 3 title "avg fu_diff"  ,\
             './data/fu_stat/fu_stat.dat' using 2:20 with lines lw 3 title "avg fu_buo"   ,\
             './data/fu_stat/fu_stat.dat' using 2:25 with lines lw 3 title "avg fu_shear"  

        set title 'simulation time [arb] vs. max function-terms[arb] (du/dt = adv+diff+buo+shear) '
        plot './data/fu_stat/fu_stat.dat' using 2:3  with lines lw 3 title "max fu"       ,\
             './data/fu_stat/fu_stat.dat' using 2:8  with lines lw 3 title "max fu_Nuk (advection)"   ,\
             './data/fu_stat/fu_stat.dat' using 2:13 with lines lw 3 title "max fu_diff"  ,\
             './data/fu_stat/fu_stat.dat' using 2:18 with lines lw 3 title "max fu_buo"   ,\
             './data/fu_stat/fu_stat.dat' using 2:23 with lines lw 3 title "max fu_shear"  

        set title 'simulation time [arb] vs. min function-terms[arb] (du/dt = adv+diff+buo+shear) '
        plot './data/fu_stat/fu_stat.dat' using 2:3  with lines lw 3 title "min fu"       ,\
             './data/fu_stat/fu_stat.dat' using 2:8  with lines lw 3 title "min fu_Nuk (advection)"   ,\
             './data/fu_stat/fu_stat.dat' using 2:13 with lines lw 3 title "min fu_diff"  ,\
             './data/fu_stat/fu_stat.dat' using 2:18 with lines lw 3 title "min fu_buo"   ,\
             './data/fu_stat/fu_stat.dat' using 2:23 with lines lw 3 title "min fu_shear"  

    unset multiplot
    ##########################################################################


do for [i=0:no_of_img] {
    set xrange [0:Lx]
    set yrange [0:Ly]
    set size ratio aspect_ratio 

   set terminal pngcairo size 1200,800 enhanced font 'Verdana,10'
#
#  	set output './visual/temp/'.i.'.png'
#  	set title 'temperature field'
#      	plot './data/temp/'.i.'.temp.dat' using 1:2:3  with image notitle
#
#  	set output './visual/u/'.i.'.png'
#  	set title 'velocity field'
#      	plot './data/u/'.i.'.u.dat' using 1:2:($3*50):($4*50)  with vectors
#
  	set output './visual/abs_u/'.i.'.png'
  	set title 'absolute magnitude  ofvelocity field'
      	plot './data/abs_u/'.i.'.abs_u.dat' using 1:2:3 with image notitle
#
  unset xrange
  unset yrange
 	set output './visual/u_f/'.i.'.png'
  set title 'absolute of fourier spectrum'
 	plot './data/u_f/'.i.'.u_f.dat' using 1:2:3  with image notitle
  set xrange [0:Lx]
  set yrange [0:Ly]
#
#    set output './visual/chem/'.i.'.png'
#  	set title 'chemical field'
#      	plot './data/chem/'.i.'.chem.dat' using 1:2:3  with image notitle

	set output './visual/buo/'.i.'.png'
	set title 'buoyancy field (B_therm*T_z - B_comp * S_z)'
 	plot './data/buo/'.i.'.buo.dat' using 1:2:3  with image notitle

    #MULTIPLOT FOR DIVERGENCE

    set terminal pngcairo size 800,400 enhanced font 'Verdana,10'
  	set output './visual/div/'.i.'.png'
		set multiplot layout 1,2
  	set title 'real part of div (u) (real space)'
      	plot './data/div/'.i.'.div.dat' using 1:2:3 with image notitle
  	set title 'real part brucker divergence kbar times u'
      	plot './data/div/'.i.'.div.dat' using 1:2:4 with image notitle
    unset multiplot


    # Multiplot for better visibility:-----------------------------------
    
    set size ratio aspect_ratio 
    set terminal pngcairo size 1600,800 enhanced font 'Verdana,30'
    set output './visual/combo/'.i.'.png'

    set lmargin 4
    set rmargin 1
    set tmargin 2 
    set bmargin 2
    set cbrange [-1.0:1.0]

    set multiplot layout 1,2
    # Axes
    set style line 11 lc rgb '#808080' lt 1
    set border 3 back ls 11
    set tics nomirror out scale 0.75

    set xrange [0:Lx]
    set yrange [0:Ly]

		set title 'chemical [arb]'
      	plot './data/chem/'.i.'.chem.dat' using 1:2:3  with image notitle
    unset xrange 
    unset yrange 
    unset cbrange
    set lmargin 5
    set xrange [0:xdim]
    set yrange [0:ydim]
    set cbrange[min_spec_order:max_spec_order]
  	set title 'logarithm of Brucker spectrum'
   	plot './data/chem_f/'.i.'.chem_f.dat' using 1:2:3 with image notitle

    #unset parametric
    #set mapping cartesian
    #set view 60,30,1,1
    #set auto
    #set isosamples 60
    #set hidden3d
   	#splot './data/chem_f/'.i.'.chem_f.dat' using 1:2:3 with image notitle

     unset margin
     unset multiplot

    #MULTIPLOT temp/chem fourier ceck--------------------------------------------------
    set terminal pngcairo size 800,800 enhanced font 'Verdana,10'
    set output './visual/spec_inspect_combo/'.i.'.png'
  	set multiplot layout 3,3


    #unset tics
    #unset colorbox
    #unset border 
  	#set title 'temperature field'
    #plot './data/temp/'.i.'.temp.dat' using 1:2:3 with image notitle
  	##set title 'Real part of fourier spectrum of temp field'
   	##plot './data/temp_f/'.i.'.temp_f.dat' using 1:2:3  with image notitle
  	#set title 'Imag part of fourier spectrum of temp field'
   	#plot './data/temp_f/'.i.'.temp_f.dat' using 1:2:4  with image notitle
  	#set title 'real cutof part of fourier spectrum of temp field'
   	#plot './data/temp_f_remap/'.i.'.temp_f_remap.dat' using 1:2:4  with image notitle

  	#set title 'chemical field'
   	#plot './data/chem/'.i.'.chem.dat' using 1:2:3 with image notitle
  	#set title 'Real part of fourier spectrum of chem field'
   	#plot './data/chem_f/'.i.'.chem_f.dat' using 1:2:3 with image notitle
  	##set title 'Imag part of fourier spectrum of chem field'
   	##plot './data/chem_f/'.i.'.chem_f.dat' using 1:2:4 with image notitle
  	#set title 'real cutof part of fourier spectrum of chem  field'
   	#plot './data/chem_f_remap/'.i.'.chem_f_remap.dat' using 1:2:4  with image notitle

    unset cbrange
    set xrange [0.0:Lx]
    set yrange [0.0:Ly]
  	set title ' u_x '
   	plot './data/u/'.i.'.u.dat' using 1:2:3 with image notitle

    set xrange [0.0:xdim]
    set yrange [0.0:ydim]
    set cbrange[min_spec_order:max_spec_order]
  	set title 'ln[Re(F(u_x))] '
   	plot './data/u_f/'.i.'.u_f.dat' using 1:2:3 with image notitle
  	set title 'ln[imag(F(u_x))] '
   	plot './data/u_f/'.i.'.u_f.dat' using 1:2:4 with image notitle

    unset cbrange
    set xrange [0.0:Lx]
    set yrange [0.0:Ly]
  	set title ' u_y '
   	plot './data/u/'.i.'.u.dat' using 1:2:4 with image notitle

    set xrange [0.0:xdim]
    set yrange [0.0:ydim]
    set cbrange[min_spec_order:max_spec_order]
  	set title 'ln[Re(F(u_y))] '
   	plot './data/u_f/'.i.'.u_f.dat' using 1:2:5 with image notitle
  	set title 'ln[imag(F(u_y))] '
   	plot './data/u_f/'.i.'.u_f.dat' using 1:2:6 with image notitle

    unset cbrange
    set xrange [0.0:Lx]
    set yrange [0.0:Ly]
  	set title 'absolute magnitude  of velocity field'
   	plot './data/abs_u/'.i.'.abs_u.dat' using 1:2:3 with image notitle
    #set yrange[0.0:0.0003]

    unset xrange  
    unset yrange
    set xtics
    set ytics

    set xrange [0.0:64.0]
  	set title 'log(energy) vs wavenumber k [arb]'
   	plot './data/k_spec/'.i.'.spec.dat' using 2:3 with lines notitle

  	#set title 'log(energy) vs wavenumber k [arb]'
   	#plot './data/k_spec/'.i.'.spec.dat' using 2:4 with lines notitle

  	#set title 'ln[Re(F(abs(u))] '
   	#plot './data/u_f/'.i.'.u_f.dat' using 1:2:7 with image notitle
  	#set title 'ln[imag(F(abs(u)))] '
   	#plot './data/u_f/'.i.'.u_f.dat' using 1:2:8 with image notitle

  	unset multiplot
    unset size
}
