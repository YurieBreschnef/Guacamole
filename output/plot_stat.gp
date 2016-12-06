    set terminal pngcairo size 1200,800 enhanced font 'Verdana,10'
    set output './visual/stat/stat_combo.png'
    set multiplot layout 4,4

        set title 'simulation time vs.Temperature measures [arb]'
        plot './output/data/T_stat/T_stat.dat' using 2:3 with lines lw 3 title "max temp"  , \
             './output/data/T_stat/T_stat.dat' using 2:4 with lines lw 3 title "min temp"  
        set title 'simulation time vs. statistical Temperature measures [arb]'
        plot './output/data/T_stat/T_stat.dat' using 2:5 with lines lw 3 title "mean temp"    ,\
             './output/data/T_stat/T_stat.dat' using 2:6 with lines lw 3 title "std-dev temp" , \
             './output/data/T_stat/T_stat.dat' using 2:7 with lines lw 3 title "variance temp"
        set title 'simulation time vs. Chemical field measures [arb]'
        plot './output/data/C_stat/C_stat.dat' using 2:3 with lines lw 3 title "max Chem"  , \
             './output/data/C_stat/C_stat.dat' using 2:4 with lines lw 3 title "min Chem"  
        set title 'simulation time vs. statistical Chemical field measures [arb]'
        plot './output/data/C_stat/C_stat.dat' using 2:5 with lines lw 3 title "mean Chem" , \
             './output/data/C_stat/C_stat.dat' using 2:6 with lines lw 3 title "std-dev chem" , \
             './output/data/C_stat/C_stat.dat' using 2:7 with lines lw 3 title "variance chem"

        set title 'simulation [arb] time vs. vx measures [arb]'
        plot './output/data/u_stat/u_stat.dat' using 2:3  with lines lw 3 title "vx_{max}"   ,\
             './output/data/u_stat/u_stat.dat' using 2:4  with lines lw 3 title "vx_{min}"   
        set title 'simulation [arb] time vs. statistical vx measures [arb]'
        plot './output/data/u_stat/u_stat.dat' using 2:5  with lines lw 3 title "vx_{avg}"   ,\
             './output/data/u_stat/u_stat.dat' using 2:6  with lines lw 3 title "std-dev vx"   ,\
             './output/data/u_stat/u_stat.dat' using 2:7  with lines lw 3 title "variance vx" 
        set title 'simulation [arb] time vs. vy measures [arb]'
        plot './output/data/u_stat/u_stat.dat' using 2:8  with lines lw 3 title "vy_{max}"   ,\
             './output/data/u_stat/u_stat.dat' using 2:9  with lines lw 3 title "vy_{min}"   
        set title 'simulation [arb] time vs. statistical vy measures [arb]'
        plot './output/data/u_stat/u_stat.dat' using 2:10  with lines lw 3 title "vy_{avg}"   ,\
             './output/data/u_stat/u_stat.dat' using 2:11  with lines lw 3 title "std-dev vy"   ,\
             './output/data/u_stat/u_stat.dat' using 2:12  with lines lw 3 title "variance vy" 

        set title 'simulation [arb] time vs. abs velocities (norm) [arb]'
        plot './output/data/u_stat/u_stat.dat' using 2:13 with lines lw 3 title "v_{max}"   ,\
             './output/data/u_stat/u_stat.dat' using 2:14 with lines lw 3 title "v_{min}"   
        set title 'simulation [arb] time vs. statistical abs velocities (norm) [arb]'
        plot './output/data/u_stat/u_stat.dat' using 2:15 with lines lw 3 title "v_{rms}"   ,\
             './output/data/u_stat/u_stat.dat' using 2:16 with lines lw 3 title "std-dev v"   ,\
             './output/data/u_stat/u_stat.dat' using 2:17 with lines lw 3 title "variance v"  
        set title 'simulation time [arb] vs. log(E_kin [arb])'
        plot './output/data/E_stat/E_stat.dat' using 2:(log($3)) with lines lw 3 title "log(E_{kin})"   ,\
             './output/data/E_stat/E_stat.dat' using 2:(log($4)) with lines lw 3 title "log(E_{pot})"   ,\
             './output/data/E_stat/E_stat.dat' using 2:(log($5)) with lines lw 3 title "log(E_{tot})"  
        set title 'simulation time vs. shearstrength '
        plot './output/data/sys_stat/sys_stat.dat' using 2:4 with lines title "shear strength [arb]"

        set title 'aperiodicity measure in x and y-dir vs. t [arb] '
        plot './output/data/sys_stat/sys_stat.dat' using 2:8 with lines lw 3 title "y_aperiodicity [arb]", \
             './output/data/sys_stat/sys_stat.dat' using 2:9 with lines lw 3 title "x_aperiodicity [arb]"
        set title 'simulation time vs. stepwidth dt [arb] '
        plot './output/data/sys_stat/sys_stat.dat' using 2:5 with lines lw 3 title "dt [arb]"
        set title 'simulation time vs. average vorticity [arb] '
        plot './output/data/sys_stat/sys_stat.dat' using 2:7 title "average vort [arb]"
	unset output
    unset multiplot

