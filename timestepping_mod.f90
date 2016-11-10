module timestepping
  ! implements time-stepping scemes that are called by e.g. 'call RK4_step()','call euler_step()'
  use sys_state
  use pdgl
  use plans
  use remap
  implicit none

  contains
!------------------------------------------------------------------------------------------
subroutine RK4_step()
  ! performs a timestep with RK4 and stores the new result in u_f,temp_f,chem_f
  !___________ REMAPPING_____________________
  if(remapping==1 .AND.shearing==1.) then
    call remap_stepwise()
  end if
  !_____________________k1_________________________________
  call set_ik_bar(sheartime) 
	state%u_k1%val = fu(state%u_f%val ,state%temp_f%val ,state%chem_f%val,state%t)     
	state%t_k1%val = ft(state%u_f%val ,state%temp_f%val ,state%t)     
	state%c_k1%val = fc(state%u_f%val ,state%chem_f%val ,state%t)     
  !_____________________k2_________________________________
  call set_ik_bar(sheartime) 
	state%u_k2%val = fu(state%u_f%val   +dt_2*state%u_k1%val,&  !f(u_f,temp_f,chem_f,t)
                      state%temp_f%val+dt_2*state%t_k1%val,&
                      state%chem_f%val+dt_2*state%c_k1%val,&
                      state%t+dt_2)     

	state%t_k2%val = ft(state%u_f%val   +dt_2*state%u_k1%val,&  !ft(u_f,temp_f,t)
                      state%temp_f%val+dt_2*state%t_k1%val,&
                      state%t+dt_2)     
	state%c_k2%val = fc(state%u_f%val   +dt_2*state%u_k1%val,&  !ft(u_f,chem_f,t)
                      state%chem_f%val+dt_2*state%c_k1%val,&
                      state%t+dt_2)     
  !_____________________k3_________________________________
  call set_ik_bar(sheartime) 
	state%u_k3%val = fu(state%u_f%val   +dt_2*state%u_k2%val,&
                      state%temp_f%val+dt_2*state%t_k2%val,&
                      state%chem_f%val+dt_2*state%c_k2%val,&
                      state%t+dt_2)     
	state%t_k3%val = ft(state%u_f%val   +dt_2*state%u_k2%val,&
                      state%temp_f%val+dt_2*state%t_k2%val,&
                      state%t+dt_2)     
	state%c_k3%val = fc(state%u_f%val   +dt_2*state%u_k2%val,&  
                      state%chem_f%val+dt_2*state%c_k2%val,&
                      state%t+dt_2)     
  !_____________________k4_________________________________
  call set_ik_bar(sheartime) 
	state%u_k4%val = fu(state%u_f%val   +dt*state%u_k3%val,&
                      state%temp_f%val+dt*state%t_k3%val,&
                      state%chem_f%val+dt*state%c_k3%val,&
                      state%t+dt)     
	state%t_k4%val = ft(state%u_f%val   +dt*state%u_k3%val,&
                      state%temp_f%val+dt*state%t_k3%val,&
                      state%t+dt)     

	state%c_k4%val = fc(state%u_f%val   +dt*state%u_k3%val,&  
                      state%chem_f%val+dt*state%c_k3%val,&
                      state%t+dt)     
  !____________________step______________________________'_
	state%u_f%val     =state%u_f%val     +(dt/6.0_rp)*(      state%u_k1%val&
                                                   +2.0_rp*state%u_k2%val&
                                                   +2.0_rp*state%u_k3%val&
                                                          +state%u_k4%val   )

	state%temp_f%val  =state%temp_f%val  +(dt/6.0_rp)*(       state%t_k1%val&
                                                    +2.0_rp*state%t_k2%val&
                                                    +2.0_rp*state%t_k3%val&
                                                           +state%t_k4%val)

	state%chem_f%val  =state%chem_f%val  +(dt/6.0_rp)*(       state%c_k1%val&
                                                    +2.0_rp*state%c_k2%val&
                                                    +2.0_rp*state%c_k3%val&
                                                           +state%c_k4%val)
  call dealiase_all()
  sheartime = sheartime+dt
  state%t=state%t+dt
  state%step=state%step+1
end subroutine
!------------------------------------------------------------------------------------------
subroutine euler_step()
  !performs a timestep with simple euler and stores the new result in u_f,temp_f,chem_f
  if(debuglevel .GE.3) write(*,*)'euler sub called'
  ! REMAPPING
  if(remapping==1 .AND.shearing==1.) then
    call remap_stepwise()
  end if
  call set_ik_bar(sheartime) 
  state_np1%u_f%val    = state%u_f%val    + dt*fu(state%u_f%val ,state%temp_f%val,state%chem_f%val,sheartime)     
  state_np1%temp_f%val = state%temp_f%val + dt*ft(state%u_f%val ,state%temp_f%val ,sheartime)     
  state_np1%chem_f%val = state%chem_f%val + dt*fc(state%u_f%val ,state%chem_f%val ,sheartime)     
  state%u_f%val     = state_np1%u_f%val
  state%temp_f%val = state_np1%temp_f%val
  state%chem_f%val = state_np1%chem_f%val
  !call dealiase_all()
  sheartime = sheartime+dt
  state%t   = state%t+dt
  state%step= state%step+1
end subroutine
!-------------------------------------------------------------------------------------------
subroutine IF2_step()
  ! perfomrs Integrating factor technique as in Brucker2007 (2nd order estimation of integral)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)  :: u_q
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)  :: u_RHS_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)  :: u_RHS_np1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)  :: u_exp_qh  ! read as exp(q*h)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)  :: u_exp_mqh ! read as exp(q*h)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)      :: t_q
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)      :: t_RHS_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)      :: t_RHS_np1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)      :: t_exp_qh  ! read as exp(q*h)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)      :: t_exp_mqh ! read as exp(q*h)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)      :: c_q
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)      :: c_RHS_n
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)      :: c_RHS_np1
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)      :: c_exp_qh  ! read as exp(q*h)
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1)      :: c_exp_mqh ! read as exp(q*h)
  ! REMAPPING
  if(remapping==1 .AND.shearing==1.) then
    call remap_stepwise()
  end if
  ! set q-values for exponent, note the minus sign in iki_sqr
  u_q(:,:,1) = real(-D_visc  *state%iki_bar_sqr%val(:,:),rp)
  u_q(:,:,2) = real(-D_visc  *state%iki_bar_sqr%val(:,:),rp)
  t_q        = real(-D_therm *state%iki_bar_sqr%val,rp)
  c_q        = real(-D_comp  *state%iki_bar_sqr%val,rp)
  ! calc exponentials for multiplication
  u_exp_qh = real(exp(u_q*dt),rp)
  t_exp_qh = real(exp(t_q*dt),rp)
  c_exp_qh = real(exp(c_q*dt),rp)
  ! with minus in exponent for timestep
  u_exp_mqh = real(exp(-u_q*dt),rp)
  t_exp_mqh = real(exp(-t_q*dt),rp)
  c_exp_mqh = real(exp(-c_q*dt),rp)
  ! calc  RHS_n and RHS_np1 (read as n plus one)
  u_RHS_n   = fu_N(state%u_f%val*u_exp_qh,state%temp_f%val*t_exp_qh,state%chem_f%val*c_exp_qh,sheartime)
  t_RHS_n   = ft_N(state%u_f%val*u_exp_qh,state%temp_f%val*t_exp_qh                          ,sheartime)
  c_RHS_n   = fc_N(state%u_f%val*u_exp_qh                          ,state%chem_f%val*c_exp_qh,sheartime)
  u_RHS_np1 = fu_N(state%u_f%val*u_exp_qh+dt*u_RHS_n      &
                  ,state%temp_f%val*t_exp_qh +dt*t_RHS_n  &
                  ,state%chem_f%val*c_exp_qh +dt*c_RHS_n  &
                  ,sheartime+dt)
  t_RHS_np1 = ft_N(state%u_f%val*u_exp_qh+dt*u_RHS_n      &
                  ,state%temp_f%val*t_exp_qh +dt*t_RHS_n  &
                  ,sheartime+dt)                          
  c_RHS_np1 = fc_N(state%u_f%val*u_exp_qh+dt*u_RHS_n      &
                  ,state%chem_f%val*c_exp_qh +dt*c_RHS_n  &
                  ,sheartime+dt)
  ! make step
  state%u_f%val   =state%u_f%val    *u_exp_mqh + dt_2*(u_RHS_n*u_exp_mqh + u_RHS_np1)
  state%temp_f%val=state%temp_f%val *t_exp_mqh + dt_2*(t_RHS_n*t_exp_mqh + t_RHS_np1)
  state%chem_f%val=state%chem_f%val *c_exp_mqh + dt_2*(c_RHS_n*c_exp_mqh + c_RHS_np1)
  call dealiase_all()
  sheartime = sheartime+dt
  call set_ik_bar(sheartime)
  state%t           = state%t     +dt
  state%step        = state%step  +1
end subroutine

subroutine ABBDF3_step()
  !calculates a timestep according to AB/BDF theme as in [Peyret 2002]
  !starting sceme
  if(state%step==0) then
    write(*,*) 'sub ABBDF3: staring scheme called..'
    state_nm2 = state
    call euler_step()
    state_nm1 = state
    call euler_step()
    write(*,*) 'sub ABBDF3: staring scheme done..'
    return
  else
    !write(*,*) 'sub ABBDF3: calculating first part, current step:',state%step
    state_np1     = state
    state_np1%u_f%val   =(a_1*state%u_f%val    + a_2*state_nm1%u_f%val     + a_3*state_nm2%u_f%val   )
    state_np1%temp_f%val=(a_1*state%temp_f%val + a_2*state_nm1%temp_f%val  + a_3*state_nm2%temp_f%val)
    state_np1%chem_f%val=(a_1*state%chem_f%val + a_2*state_nm1%chem_f%val  + a_3*state_nm2%chem_f%val)

    !write(*,*) 'sub ABBDF3: calculating second part'
    state_np1%u_f%val   =(-state_np1%u_f%val&
          	+dt*( b_0*fu(state%u_f%val    ,state%temp_f%val    ,state%chem_f%val    ,sheartime)&
                       +b_1*fu(state_nm1%u_f%val,state_nm1%temp_f%val,state_nm1%chem_f%val,sheartime-dt)&
                       +b_2*fu(state_nm2%u_f%val,state_nm2%temp_f%val,state_nm2%chem_f%val,sheartime-(2.0_rp*dt))))/a_0
    state_np1%temp_f%val   =(-state_np1%temp_f%val&
          	+dt*( b_0*ft(state%u_f%val    ,state%temp_f%val    ,sheartime)&
                       +b_1*ft(state_nm1%u_f%val,state_nm1%temp_f%val,sheartime-dt)&
                       +b_2*ft(state_nm2%u_f%val,state_nm2%temp_f%val,sheartime-(2.0_rp*dt))))/a_0
    state_np1%chem_f%val   =(-state_np1%chem_f%val&
          	+dt*( b_0*ft(state%u_f%val    ,state%chem_f%val    ,sheartime)&
                       +b_1*ft(state_nm1%u_f%val,state_nm1%chem_f%val,sheartime-dt)&
                       +b_2*ft(state_nm2%u_f%val,state_nm2%chem_f%val,sheartime-(2.0_rp*dt))))/a_0

    !write(*,*) 'sub ABBDF3: resetting states'
    !set new states
    state_nm2 = state_nm1
    state_nm1 = state
    state     = state_np1

    sheartime = sheartime+dt

    if(benchmarking ==1) bm_set_ik_bar_starttime=  omp_get_wtime()
    call set_ik_bar(sheartime)
    if(benchmarking ==1) bm_set_ik_bar_endtime=  omp_get_wtime()
    call set_ik_bar(sheartime)
    state%t           = state%t     +dt
    state%step        = state%step  +1
    !write(*,*) 'sub ABBDF3: make timestep, current step:', state%step
  end if
end subroutine



















































end module
