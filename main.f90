program guacamole 
 !----------------------------------------------------------------------------------------
 !  Copyright (C) 2016  Michael Wenske

 !  This program is free software: you can redistribute it and/or modify
 !  it under the terms of the GNU General Public License as published by
 !  the Free Software Foundation, either version 3 of the License, or
 !  (at your option) any later version.

 !  This program is distributed in the hope that it will be useful,
 !  but WITHOUT ANY WARRANTY; without even the implied warranty of
 !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 !  GNU General Public License for more details.

 !  You should have received a copy of the GNU General Public License
 !  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 !----------------------------------------------------------------------------------------
 ! guacamole : simulation double diffusive convection with homogeneous background shear
 !
 ! main program: the main timestepping loop is implemented here, calling the actual
 ! stepping routines in the timestepping_mod, and writing to files by calling write_all()
 ! from IO_module.
 ! periodic measures are taken according to parameter: measure_every in const_mod

  !use omp_lib
  use iso_c_binding
  use plans
  use sys_state
  use init
  use exit_mod
  use timestepping
  use IO_mod
  use trafo
  use test
  use benchmark
  implicit none
  if(debuglevel .GE. 1) write(*,*) '__________________START____________________________________'
  call init_all()
  call test_all()
  call write_all()      ! write initial conditions to file
  last_written = last_written+write_intervall
  if(debuglevel .GE.1) write(*,*) '__________________TIMESTEPPING_____________________________'
  !do main_stp= 0,steps
  do while (state%t <tmax)
    ! STATISTICS WRITING----------------------------------------------------------------------
    if(benchmarking ==1) bm_step_starttime=  omp_get_wtime()
    if(benchmarking ==1) bm_statwrite_starttime=  omp_get_wtime()
    if(mod(state%step,(measure_every)).EQ.0) then
      call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,sheartime)
      call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,sheartime)
      call transform(state%temp_f%val,state%temp%val,-1,shearing,sheartime)
      call transform(state%chem_f%val,state%chem%val,-1,shearing,sheartime)
      call write_u_stat()
      call write_fu_stat()
      call write_E_stat()
      call write_T_stat()
      call write_C_stat()
      call write_sys_stat()
    end if 
    if(benchmarking ==1) bm_statwrite_endtime=  omp_get_wtime()
    ! FILE WRITING----------------------------------------------------------------------------
    if(benchmarking ==1) bm_filewrite_starttime=  omp_get_wtime()
    if(benchmarking ==1) bm_filewrite_endtime=  omp_get_wtime()
    if(state%t > last_written) then
      if(benchmarking ==1) bm_filewrite_starttime=  omp_get_wtime()
      call write_all()
      call plot_all(int(state%t/write_intervall)-1,int(state%t/write_intervall)-1,1.0,-25.0)
      last_written = last_written+write_intervall
      if(benchmarking ==1) bm_filewrite_endtime=  omp_get_wtime()
      ! check for NAN's occasionally (recourse intensive):
      IF(ANY(IsNaN(real(state%u_f%val))))  then
        write(*,*) 'main: after write_all(): NAN detected in u_f array'
        stop
      end if
    end if 
    ! CONSOLE OUTPUT----------------------------------------------------------------------------
    if((steps>=1000).AND.(mod(state%step,(steps/1000)).EQ.0)) then
        write(*,*) ' permille: ',int(state%step/(steps/1000),2),&
                   '|step:',int(main_stp,2), &
                   '|t:',real(state%t,4),&
                   '| dt:',real(dt,4),&
                   '|shearing:',int(shearing,1),&
                   '|sheartime:(',real(sheartime,4),'/',real(T_rm,4),&
                 ') |steptime:',real(bm_timestepping_endtime-bm_timestepping_starttime,4)
    end if
    ! TIMESTEPPING ---------------------------------------------------------------------------
    if(.TRUE.) bm_timestepping_starttime=  omp_get_wtime()
    !call RK4_adjust_dt()
    !call RK4_step()
    !call euler_step()
    call IF2_step()
    !call div_tester()
    main_stp = main_stp +1
    if(.TRUE.) bm_timestepping_endtime=  omp_get_wtime()
    ! ----------------------------------------------------------------------------------------
    if(benchmarking ==1) bm_step_endtime=  omp_get_wtime()
    !BENCHMARKING------------------------------------
    if(benchmarking ==1) call bm_evaluate(.true.)
  end do

  if(debuglevel <= 1) write(*,*) '__________________END OF TIMESTEPPING______________________'
  if(debuglevel <= 1) write(*,*) '__________________END______________________________________'
  call exit_all()
end program guacamole 



    !if(state%t >1.0_rp*tmax/5.0_rp) then
    !    shearing = 1
    !    shear =  0.1
    !end if

    !if(state%t >2.0_rp*tmax/5.0_rp) then
    !    shearing = 1
    !    shear = 0.05_rp
    !end if

    !if(state%t >3.0_rp*tmax/5.0_rp) then
    !    shearing = 1
    !    shear =0.05+ (state%t-3.0*tmax/5.0_rp)/(tmax/(5.0)) *0.05_rp
    !end if

    !if(state%t >4.0_rp*tmax/5.0_rp) then
    !    shearing = 1
    !    shear =0.10
    !end if
