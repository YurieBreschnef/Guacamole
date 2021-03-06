module benchmark
! measures the execution time of different parts of the program. This output can be switched on by the const
! "benchmarking" in the const_mod 
!
! usage: write variables by adding lines:
!     if(benchmarking ==1) bm_SUBROUTINE_starttime=  omp_get_wtime()
!     if(benchmarking ==1) bm_SUBROUTINE_endtime  =  omp_get_wtime()
! add substraction and output at base of this module
use const
use sys_state

real(kind = rp)               :: bm_step_starttime
real(kind = rp)               :: bm_step_endtime
real(kind = rp)               :: bm_step_time

real(kind = rp)               :: bm_statwrite_starttime
real(kind = rp)               :: bm_statwrite_endtime
real(kind = rp)               :: bm_statwrite_time

real(kind = 8)               :: bm_trafo_starttime
real(kind = 8)               :: bm_trafo_endtime
real(kind = 8)               :: bm_trafo_time

real(kind = 8)               :: bm_remap_starttime
real(kind = 8)               :: bm_remap_endtime
real(kind = 8)               :: bm_remap_time

real(kind = rp)               :: bm_filewrite_starttime
real(kind = rp)               :: bm_filewrite_endtime
real(kind = rp)               :: bm_filewrite_time

real(kind = rp)               :: bm_timestepping_starttime
real(kind = rp)               :: bm_timestepping_endtime
real(kind = rp)               :: bm_timestepping_time

real(kind = rp)               :: bm_set_ik_bar_starttime
real(kind = rp)               :: bm_set_ik_bar_endtime
real(kind = rp)               :: bm_set_ik_bar_time

real(kind = rp)               :: bm_dealiase_starttime
real(kind = rp)               :: bm_dealiase_endtime
real(kind = rp)               :: bm_dealiase_time

real(kind = rp)               :: bm_fu_starttime
real(kind = rp)               :: bm_fu_endtime

real(kind = rp)               :: bm_fu_L_starttime
real(kind = rp)               :: bm_fu_L_endtime
real(kind = rp)               :: bm_fu_L_time

real(kind = rp)               :: bm_fu_N_starttime
real(kind = rp)               :: bm_fu_N_endtime
real(kind = rp)               :: bm_fu_N_time
real(kind = rp)               :: bm_fu_time
real(kind = rp)               :: bm_fu_Nuk_starttime
real(kind = rp)               :: bm_fu_Nuk_endtime
real(kind = rp)               :: bm_fu_Nuk_time
real(kind = rp)               :: bm_fu_buo_starttime
real(kind = rp)               :: bm_fu_buo_endtime
real(kind = rp)               :: bm_fu_buo_time
real(kind = rp)               :: bm_fu_diff_starttime
real(kind = rp)               :: bm_fu_diff_endtime
real(kind = rp)               :: bm_fu_diff_time
real(kind = rp)               :: bm_fu_shear_starttime
real(kind = rp)               :: bm_fu_shear_endtime
real(kind = rp)               :: bm_fu_shear_time

real(kind = rp)               :: bm_ft_starttime
real(kind = rp)               :: bm_ft_endtime
real(kind = rp)               :: bm_ft_N_starttime
real(kind = rp)               :: bm_ft_N_endtime
real(kind = rp)               :: bm_ft_N_time
real(kind = rp)               :: bm_ft_time
real(kind = rp)               :: bm_ft_adv_starttime
real(kind = rp)               :: bm_ft_adv_endtime
real(kind = rp)               :: bm_ft_adv_time
real(kind = rp)               :: bm_ft_diff_starttime
real(kind = rp)               :: bm_ft_diff_endtime
real(kind = rp)               :: bm_ft_diff_time
real(kind = rp)               :: bm_ft_strat_starttime
real(kind = rp)               :: bm_ft_strat_endtime
real(kind = rp)               :: bm_ft_strat_time

real(kind = rp)               :: bm_fc_starttime
real(kind = rp)               :: bm_fc_endtime
real(kind = rp)               :: bm_fc_N_starttime
real(kind = rp)               :: bm_fc_N_endtime
real(kind = rp)               :: bm_fc_N_time
real(kind = rp)               :: bm_fc_time
real(kind = rp)               :: bm_fc_adv_starttime
real(kind = rp)               :: bm_fc_adv_endtime
real(kind = rp)               :: bm_fc_adv_time
real(kind = rp)               :: bm_fc_diff_starttime
real(kind = rp)               :: bm_fc_diff_endtime
real(kind = rp)               :: bm_fc_diff_time
real(kind = rp)               :: bm_fc_strat_starttime
real(kind = rp)               :: bm_fc_strat_endtime
real(kind = rp)               :: bm_fc_strat_time

real(kind = rp)               :: bm_buo_L_starttime
real(kind = rp)               :: bm_buo_L_endtime
real(kind = rp)               :: bm_buo_L_time

real(kind = rp)               :: bm_buo_N_starttime
real(kind = rp)               :: bm_buo_N_endtime
real(kind = rp)               :: bm_buo_N_time

real(kind = rp)               :: bm_buo_adv_starttime
real(kind = rp)               :: bm_buo_adv_endtime
real(kind = rp)               :: bm_buo_adv_time

real(kind = rp)               :: bm_buo_diff_starttime
real(kind = rp)               :: bm_buo_diff_endtime
real(kind = rp)               :: bm_buo_diff_time

real(kind = rp)               :: bm_buo_strat_starttime
real(kind = rp)               :: bm_buo_strat_endtime
real(kind = rp)               :: bm_buo_strat_time



real(kind = rp)               :: bm_percent_unaccounted

contains

subroutine bm_evaluate(write_to_console)
  logical,intent(in)               :: write_to_console



  bm_buo_L_time     = bm_buo_L_endtime  - bm_buo_L_starttime
  bm_buo_N_time     = bm_buo_N_endtime - bm_buo_N_starttime
  bm_buo_adv_time   = bm_buo_adv_endtime -bm_buo_adv_starttime
  bm_buo_diff_time  = bm_buo_diff_endtime  - bm_buo_diff_starttime
  bm_buo_strat_time = bm_buo_strat_endtime - bm_buo_strat_starttime

  bm_step_time  =         bm_step_endtime-bm_step_starttime
  bm_trafo_time=         bm_trafo_endtime-bm_trafo_starttime
  bm_remap_time=         bm_remap_endtime-bm_remap_starttime
  bm_fu_time    =         bm_fu_endtime-bm_fu_starttime
  bm_fu_N_time    =         bm_fu_N_endtime-bm_fu_N_starttime
  bm_fu_L_time    =         bm_fu_L_endtime-bm_fu_L_starttime
  bm_fu_Nuk_time    =     bm_fu_Nuk_endtime-bm_fu_Nuk_starttime
  bm_fu_buo_time    =     bm_fu_buo_endtime-bm_fu_buo_starttime
  bm_fu_diff_time    =     bm_fu_diff_endtime-bm_fu_diff_starttime
  bm_fu_shear_time    =     bm_fu_shear_endtime-bm_fu_shear_starttime
  bm_set_ik_bar_time =     bm_set_ik_bar_endtime-bm_set_ik_bar_starttime
  bm_dealiase_time =     bm_dealiase_endtime-bm_dealiase_starttime

  bm_ft_time    =         bm_ft_endtime-bm_ft_starttime
  bm_ft_N_time    =         bm_ft_N_endtime-bm_ft_N_starttime
  bm_ft_adv_time    =     bm_ft_adv_endtime-bm_ft_adv_starttime
  bm_ft_diff_time   =     bm_ft_diff_endtime-bm_ft_diff_starttime
  bm_ft_strat_time   =     bm_ft_strat_endtime-bm_ft_strat_starttime

  bm_fc_time    =         bm_fc_endtime-bm_fc_starttime
  bm_fc_N_time    =       bm_fc_N_endtime-bm_fc_N_starttime
  bm_fc_adv_time    =     bm_fc_adv_endtime-bm_fc_adv_starttime
  bm_fc_diff_time   =     bm_fc_diff_endtime-bm_fc_diff_starttime
  bm_fc_strat_time   =    bm_fc_strat_endtime-bm_fc_strat_starttime

  bm_timestepping_time =  bm_timestepping_endtime-bm_timestepping_starttime
  bm_statwrite_time =     bm_statwrite_endtime-bm_statwrite_starttime
  bm_filewrite_time =     bm_filewrite_endtime-bm_filewrite_starttime
  bm_timestepping_time=   bm_timestepping_endtime-bm_timestepping_starttime

  if(write_to_console) then
  write(*,*) '______________________________BENCHMARK:_______________________________________'
  write(*,*) 'total step:                :',bm_step_time               ,'sec,'
  write(*,*) '  -statwrite               :  ',bm_statwrite_time          ,'sec,',int(100.0_rp*bm_statwrite_time/bm_step_time),'%'
  write(*,*) '  -filewrite               :  ',bm_filewrite_time,'sec,',int(100.0_rp*bm_filewrite_time/bm_step_time),'%'
  write(*,*) '  -timestepping:           :  ',bm_timestepping_time,'sec,',int(100.0_rp*bm_timestepping_time/bm_step_time),'%'
  write(*,*) '_______________________________IF2____________________________________________'
  write(*,*) '    -function fu_N:        :    ',  bm_fu_N_time,'  sec,',int(100.0_rp*bm_fu_N_time/bm_step_time),'%'
  write(*,*) '       -function fu_Nuk    :      ',bm_fu_Nuk_time,'sec,',int(100.0_rp*bm_fu_Nuk_time/bm_fu_time),'%'
  write(*,*) '       -function fu_buo    :      ',bm_fu_buo_time,'sec,',int(100.0_rp*bm_fu_buo_time/bm_fu_time),'%'
  write(*,*) '       -function fu_shear  :      ',bm_fu_shear_time,'sec,',int(100.0_rp*bm_fu_shear_time/bm_fu_time),'%'
  write(*,*) '    -function fu_L:        :    ',  bm_fu_L_time,'  sec,',int(100.0_rp*bm_fu_L_time/bm_step_time)&
,'%(measured once?!)'
  write(*,*) '       -function fu_diff   :      ',bm_fu_diff_time,'sec,',int(100.0_rp*bm_fu_diff_time/bm_fu_time),'%'
  write(*,*) '    -function buo_N:       :    ',  bm_buo_N_time,'  sec,',int(100.0_rp*bm_buo_N_time/bm_step_time),'%'
  write(*,*) '       -function buo_adv   :      ',bm_buo_adv_time,'sec,',int(100.0_rp*bm_ft_adv_time/bm_buo_N_time),'%'
  write(*,*) '       -function buo_strat :      ',bm_buo_strat_time,'sec,',int(100.0_rp*bm_buo_strat_time/bm_buo_N_time),'%'
  write(*,*) '    -function buo_L:       :    ',  bm_buo_L_time,'  sec,',int(100.0_rp*bm_buo_L_time/bm_step_time)&
,'%(measured once?)'
  write(*,*) '       -function buo_diff  :      ',bm_buo_diff_time,'sec,',int(100.0_rp*bm_buo_diff_time/bm_buo_L_time)&
,'% (measured once?!)'
  write(*,*) '_______________________________Trafo__________________________________________'
  write(*,*) '    -single trafo          :',bm_trafo_time,'sec,',int(100.0_rp*bm_trafo_time/bm_step_time),'%'
  write(*,*) '    -remapping             :',bm_remap_time,'sec,',int(100.0_rp*bm_remap_time/bm_step_time),'%'
  write(*,*) 'percent unnacounted        :',1.0_rp - (bm_statwrite_time+bm_filewrite_time+bm_timestepping_time)/bm_step_time,'%'
  write(*,*) '_______________________________Trafo__________________________________________'
  write(*,*) '    -dealiasing            :',bm_dealiase_time,'sec'
  write(*,*) '    -set_ik_bar		 :',bm_set_ik_bar_time,'sec'
  ! note that ft and fc will take the same ammount of computing time
  end if
 ! call super functions to measure their time to (in IF2 the FU_L and buo_L are not used)
end subroutine

end module
