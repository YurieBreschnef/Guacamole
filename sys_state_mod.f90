module sys_state
  ! state-module: contains all information characterizing a system state at a given time
  ! holds fields, current time and step as well as several helping variables (not super clean!)
  use const
  implicit none

  type sfield
    !scalar field
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1):: val = 0.0_rp
    ! Chem/Temp-field value (i,j), amplitude stored in %val
  end type

  type sfield_stats
    ! statistical measures of a scalar field
    ! index:            meaning
    !
    !   0               max
    !   1               min
    !   2               mean
    !   3               standard deviation 
    !   4               variance
  	real(kind=rp),dimension(0:4):: val = 0.0_rp
  end type

  type vfield
    !vector field
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1,2):: val = 0.0_rp  
    !u(i,j,spatial_dir) acces by state%u%val(i,j, spatial direction whished)
  end type

  type kfield
    !field of k-values for derivaties in spectral domain
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1):: val = 0.0_rp  
  end type

  type system_state
    ! defines the physical system state in its entirety. Many functions and subs
    ! will use this type as return type
    type(vfield)                      :: u
    type(vfield)                      :: dummy
    type(vfield)                      :: dummy_f
    type(vfield)                      :: z_dummy_f
    type(vfield)                      :: u_f
    type(vfield)                      :: u_k1,u_k2,u_k3,u_k4	
    type(vfield)                      :: k_vec

    type(sfield)                      :: temp
    type(sfield)                      :: temp_f
    type(sfield)                      :: t_k1,t_k2,t_k3,t_k4	

    type(sfield)                      :: chem  
    type(sfield)                      :: chem_f
    type(sfield)                      :: c_k1,c_k2,c_k3,c_k4

    type(sfield)                      :: s_dummy
    type(sfield)                      :: s_dummy_f

    type(kfield)                      :: ikx,iky! k's for deriv
    type(kfield)                      :: ikx_sqr,iky_sqr! k's for deriv
    type(kfield)                      :: iki_sqr

    type(kfield)                      :: ikx_bar,iky_bar! k's for deriv
    type(kfield)                      :: ikx_bar_sqr,iky_bar_sqr! k's for deriv
    type(kfield)                      :: iki_bar_sqr

    real(kind = rp)                   :: t    = 0.0_rp	! time variable
    integer	                      :: step = 0       ! acute step of sim
  end type

  !MAIN STATE VARIABLE:
  type(system_state)                                          ::state_np1 !state n+1	 (needed for IF2 timestepping)
  type(system_state)                                          ::state
  type(system_state)                                          ::state_nm1 !state n-1	(needed for ABBDF3)	
  type(system_state)                                          ::state_nm2 !state n-2
contains


  function measure_sfield_stats(in_field)
    ! measures statistics of a given scalar field. see typedef for indice-explanation
    type(sfield),intent(in)           ::in_field
    type(sfield_stats)                ::measure_sfield_stats

    !minimum/maximum / average
    measure_sfield_stats%val = 0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
        !max
        if(measure_sfield_stats%val(0) < real(in_field%val(i,j),rp)) measure_sfield_stats%val(0)= real(in_field%val(i,j),rp)
        !min
        if(measure_sfield_stats%val(1) > real(in_field%val(i,j),rp)) measure_sfield_stats%val(1)= real(in_field%val(i,j),rp)
        !sum for average
        measure_sfield_stats%val(2) = measure_sfield_stats%val(2) + real(in_field%val(i,j),rp)
    end do
    end do
    !norm average
    measure_sfield_stats%val(2) = measure_sfield_stats%val(2)/real(xdim*ydim,rp) 
    !std-dev and variance
    do i=0,xdim-1
    do j=0,ydim-1
        !sum up deviations from mean for std-dev
        measure_sfield_stats%val(3) = measure_sfield_stats%val(3) + abs(real(in_field%val(i,j),rp)-measure_sfield_stats%val(2))
        !sum up squared deviations from mean for variance
        measure_sfield_stats%val(4) = measure_sfield_stats%val(4) + abs(real(in_field%val(i,j),rp)-measure_sfield_stats%val(2))**2
    end do
    end do
    measure_sfield_stats%val(3) = measure_sfield_stats%val(3)/real(xdim*ydim,rp)
    measure_sfield_stats%val(4) = measure_sfield_stats%val(4)/real(xdim*ydim,rp) 
    
  end function 

  function measure_Ekin()
    !measures a absolute value for kinetic energy within system.
    !TODO - devide by relative density
    real(kind=real_outp_precision)            ::measure_Ekin
    real(kind=rp)                             ::Ekin
    Ekin = 0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
      Ekin = Ekin+0.5_rp*(real(state%u%val(i,j,1))**2 &
                  +real(state%u%val(i,j,2))**2)
    end do
    end do
    measure_Ekin = real(Ekin,real_outp_precision)/real(xdim*ydim,real_outp_precision)
  end function

  function measure_vmax()
    ! returns the highest absolute velocitiy within the system
    real(kind=real_outp_precision)            ::measure_vmax
    real(kind=rp)                             ::vmax
    vmax =0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
    do l=1,2
        if(vmax < real(state%u%val(i,j,l),rp)) vmax = real(state%u%val(i,j,l),rp)
    end do
    end do
    end do
    measure_vmax = vmax
  end function

  function measure_Epot()
    !measures a absolute value for kinetic energy within system.
    !TODO - devide by relative density
    real(kind=real_outp_precision)            ::measure_Epot
    real(kind=rp)                             ::Epot
    Epot= 0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
      Epot= Epot+abs(real((B_therm*state%temp%val(i,j) - B_comp*state%chem%val(i,j)),rp))*(real(j,rp)/real(ydim,rp)*Ly)
    end do
    end do
    measure_Epot = real(Epot,real_outp_precision)/real(xdim*ydim,real_outp_precision)
  end function

  function measure_av(input)
    !measures a relative average value for input sfield .
  	complex(kind=rp),dimension(0:xdim-1,0:ydim-1) :: input
    real(kind=real_outp_precision)                ::measure_av
    real(kind=rp)                                 ::tmp
    tmp= 0.0_rp
    do i=0,xdim-1
    do j=0,ydim-1
      tmp= tmp+(input(i,j)) 
    end do
    end do
    measure_av= real(tmp,real_outp_precision)/real(xdim*ydim,real_outp_precision)
  end function

  function measure_u_rms()
    !measures a relative value for temp  within system.
    real(kind=real_outp_precision)            ::measure_u_rms
    real(kind=rp)                             ::u_rms
    u_rms= 0.0_rp
    do j=0,ydim-1
    do i=0,xdim-1
    do k=1,2
      u_rms= u_rms+(state%u%val(i,j,1)**2 + state%u%val(i,j,2)**2)
    end do
    end do
    end do
    measure_u_rms= sqrt(real((u_rms),real_outp_precision)/real(xdim*ydim,real_outp_precision))
  end function


  subroutine set_ik_bar(ktime)
    ! resets the ik_bar wave vectors for a given time
    ! note how the x component is not really used, but still there for generality and possible
    ! future changes
    real(kind = rp),intent(in)                   :: ktime
    !TODO very ineffective to reset k's if shearing is off.

     !IF(ALL((state%ikx%val ==0.0_rp).OR.ALL(state%iky%val ==0.0_rp)))  then
     !  write(*,*) 'sub set_ik_bar(): ALL ikx or iky are ZERO. BAD. VERY BAD.'
     !  stop
     !end if

    if(shearing ==1) then
      !$omp parallel &
      !$omp private (i,j)
      !$omp do
        do j=0,ydim-1
          do i=0,xdim-1
          state%ikx_bar%val(i,j) = state%ikx%val(i,j) 
          state%iky_bar%val(i,j) = state%iky%val(i,j)-shear*ktime*state%ikx%val(i,j)

          state%ikx_bar_sqr%val(i,j) = state%ikx_bar%val(i,j)**2
          state%iky_bar_sqr%val(i,j) = state%iky_bar%val(i,j)**2
          state%iki_bar_sqr%val(i,j) = state%ikx_bar%val(i,j)**2 + state%iky_bar%val(i,j)**2
          state%k_vec%val(i,j,1) = real(imag*state%ikx_bar%val(i,j),rp)
          state%k_vec%val(i,j,2) = real(imag*state%iky_bar%val(i,j),rp)
          end do
        end do
      !$omp end do
      !$omp end parallel
    end if

  end subroutine

end module
