module IO_mod
  ! main input/output module. add output routine either to write_all() which is invoked every writestep, or invoke
  ! manually. some outputs are very recource intensive.
  use sys_state
  use plans
  use nabla
  use trafo
  use pdgl
  implicit none
  contains

  subroutine write_all()
    !if(debuglevel <= 2) write(*,*) '-calling write_all()'
    
    call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,shearing,sheartime)
    call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,shearing,sheartime)
    call transform(state%temp_f%val,state%temp%val,-1,shearing,sheartime)
    call transform(state%chem_f%val,state%chem%val,-1,shearing,sheartime)

    ! ESSENTIAL
    call write_u()
    call write_temp()
    call write_chem()

    ! OPTIONAL 
    call write_abs_u()
    call write_buo()
    call write_div()
    call write_vort()
    call write_u_f()
    call write_k_spec()
    call write_chem_f()
    call write_chem_f_remap()
    call write_temp_f()
    call write_temp_f_remap()

    !STATISTICS
    call write_u_stat()     ! u-relatetd measures
    call write_fu_stat()    ! pdgl-relatetd measures (influences of terms)
    call write_E_stat()     ! Energy related measures
    call write_T_stat()     ! Energy related measures
    call write_C_stat()     ! Energy related measures
    call write_sys_stat()   ! System wide measures

    !if(debuglevel <= 2) write(*,*) '-done with write_all.'
  end subroutine
  subroutine plot_all(sindex,eindex,maxspec,minspec)
    integer		::sindex,eindex
    real        	::maxspec,minspec
    !------------------------------------------
    character(len=1000)  :: command
    character(len=100)  :: x_dim,L_x
    character(len=100)  :: y_dim,L_y
    character(len=100)  :: aspect
    character(len=100)  :: start_index,end_index
    character(len=100)  :: max_spec,min_spec

    write(x_dim,*)   xdim
    write(y_dim,*)   ydim
    write(L_x,*)   Lx
    write(L_y,*)   Ly
    write(aspect,*)   int(Ly/Lx)
    write(start_index,*)  sindex 
    write(end_index,*)  eindex 
    write(max_spec,*) maxspec 
    write(min_spec,*) minspec 

    command =  'bash ./output/plot_all.sh'//' '&
                      //trim(adjustl(x_dim))//' '&
                      //trim(adjustl(y_dim))//' '&
                      //trim(adjustl(L_x))//' '&
                      //trim(adjustl(L_y))//' '&
                      //trim(adjustl(aspect))//' '&
                      //trim(adjustl(start_index))//' '&
                      //trim(adjustl(end_index))//' '&
                      //trim(adjustl(max_spec))//' '&
                      //trim(adjustl(min_spec))	
write(*,*) 'command: ' , command
 call system(command)
  end subroutine
  subroutine write_vort()
    !write chemical field to file. Note the path variable specifying the destination
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=19),parameter					:: path ='./output/data/vort/'
		!write(suffix,"(I5,A9)") state%step/(steps/maxfiles), ".vort.dat"
		write(suffix,"(I5,A9)") int(state%t/write_intervall), ".vort.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)

    !calc vorticity in fourier space
    state%s_dummy_f%val(:,:) =( state%ikx_bar%val(:,:)*state%u_f%val(:,:,2) &
                             -state%iky_bar%val(:,:)*state%u_f%val(:,:,1))

    call transform(state%s_dummy_f%val,state%s_dummy%val,-1,shearing,sheartime)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_vort!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
                      ,real(state%s_dummy%val(i,j),real_outp_precision)
	  			           
			end do
		end do
    close(20)
  
  end subroutine

  
  subroutine write_u()
    !write velocity field to file in specified realtive path
    integer                             :: io_error = 0
    real(kind = real_outp_precision)    :: out_x,out_y,out_z
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=16),parameter					:: path ='./output/data/u/'
		write(suffix,"(I5,A6)") int(state%t/write_intervall) , ".u.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_u!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
          out_x = real(state%u%val(i,j,1),real_outp_precision)
          out_y = real(state%u%val(i,j,2),real_outp_precision)
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
                      ,out_x,out_y

			end do
		end do
    close(20)
  end subroutine

  subroutine write_chem()
    !write chemical field to file. Note the path variable specifying the destination
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=19),parameter					:: path ='./output/data/chem/'
		write(suffix,"(I5,A9)")  int(state%t/write_intervall), ".chem.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_chem!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
	  			            ,real(state%chem%val(i,j),real_outp_precision)
			end do
		end do
    close(20)
  end subroutine

  subroutine write_abs_u()
    !write chemical field to file. Note the path variable specifying the destination
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=20),parameter					:: path ='./output/data/abs_u/'
		write(suffix,"(I5,A10)")  int(state%t/write_intervall), ".abs_u.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_chem!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
	  			            ,real(sqrt(real(state%u%val(i,j,1))**2+real(state%u%val(i,j,2))**2),real_outp_precision)
			end do
		end do
    close(20)
  end subroutine
  subroutine write_temp()
    !write temp-field to file. analogous to write_chem
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=19),parameter					:: path ='./output/data/temp/'
		write(suffix,"(I5,A9)")  int(state%t/write_intervall), ".temp.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_temp!'
		  do i=0,xdim-1
        do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
                  ,real(state%temp%val(i,j),real_outp_precision)
			  !end do
			end do
		end do
    close(20)
  end subroutine

  subroutine write_div()
    ! writes divergence field to file
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=18),parameter					:: path ='./output/data/div/'
		write(suffix,"(I5,A8)")  int(state%t/write_intervall), ".div.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_div!'
    !write div u  (k_i * u_i ==0)condition to dummy
    call transform(state%u_f%val(:,:,1),state%u%val(:,:,1),-1,1,sheartime)
    call transform(state%u_f%val(:,:,2),state%u%val(:,:,2),-1,1,sheartime)
    call transform(state%u%val(:,:,1),state%dummy_f%val(:,:,1),-1,1,sheartime)
    call transform(state%u%val(:,:,2),state%dummy_f%val(:,:,2),-1,1,sheartime)

    ! normal div
    state%dummy_f%val(:,:,1) = state%ikx%val(:,:)*state%dummy_f%val(:,:,1) &
                              +state%iky%val(:,:)*state%dummy_f%val(:,:,2) 

    !write brucker  (k_bar_i * u_i ==0)condition to other dummy
    state%dummy_f%val(:,:,2) = state%ikx_bar%val(:,:)*state%u_f%val(:,:,1) &
                              +state%iky_bar%val(:,:)*state%u_f%val(:,:,2) 
    call transform(state%dummy_f%val(:,:,1),state%dummy%val(:,:,1),1,1,sheartime)
    call transform(state%dummy_f%val(:,:,2),state%dummy%val(:,:,2),1,1,sheartime)
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
	  			           ,real(state%dummy%val(i,j,1),real_outp_precision)&
                     ,real(state%dummy%val(i,j,2),real_outp_precision)
			end do
		end do
    close(20)
  end subroutine

  subroutine write_buo()
    !write buoyancy field to file
    integer                             :: io_error = 0
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=18),parameter					:: path ='./output/data/buo/'
		write(suffix,"(I5,A8)")  int(state%t/write_intervall), ".buo.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_buo!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) real(i)*(Lx/real(xdim)),real(j)*(Ly/real(ydim))&
	  			,real(state%temp%val(i,j)*B_therm &
              - B_comp*state%chem%val(i,j),real_outp_precision)
			end do
		end do
    close(20)
  end subroutine
  
  subroutine write_k_spec()
    use const
    ! writes the energy spectrum to file 
    integer                             :: io_error = 0
    integer,parameter                   ::  bins = 100
    real(kind = real_outp_precision)    :: out_x,out_y,out_z
    real(kind = rp)                     :: k_abs,r_dummy,dk         ! maximum k-norm
    real(kind = rp),dimension(bins)     :: bin ! maximum k-norm
    real(kind = rp),dimension(bins)     :: bin_area ! maximum k-norm
    type(sfield)                        :: kx_dummy,ky_dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/k_spec/'
		write(suffix,"(I5,A9)")  int(state%t/write_intervall), ".spec.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_k_spec!'
    bin = 0.0_rp
    call set_ik_bar(sheartime)
    ! calculate maximal distance from center (note that outer edges are not properly resolved)
    k_abs = sqrt(aimag(state%ikx%val(xdim/2,1))**2 + aimag(state%iky%val(1,ydim/2))**2)
    ! bin-width dk
    dk = k_abs / real(bins,rp)
    ! calculate bin weight by area
    do k = 0,bins-1
      bin_area(k) = (pi*(real(k+1,rp)**2))-(pi*(real(k,rp)**2))
    end do
    r_dummy = 0.0_rp
		do i=0,xdim-1
	  	do j=0,ydim-1
         do k = 0,bins-1
            ! how far is k-vector from zero mode?
            r_dummy = real(sqrt(aimag(state%ikx_bar%val(i,j))**2 + aimag(state%iky_bar%val(i,j))**2),rp)
            if(r_dummy >= (real(k,rp)*dk)) then
              if(r_dummy <= (real(k+1,rp)*dk)) then
                  ! sort energy in this mode into proper bin
                  bin(k)  = bin(k) +(sqrt(real(state%u_f%val(i,j,1))**2 + aimag(state%u_f%val(i,j,1))**2 ))!&
                                   !+(sqrt(real(state%u_f%val(i,j,2))**2 + aimag(state%u_f%val(i,j,2))**2 ))
              end if
            end if
          end do
			end do
		end do
		do k=0,bins-1
	  			write(20,*) k,&
                      real(k,rp)*dk,&
                      real((bin(k)/bin_area(k)),real_outp_precision),&
                      log(real((bin(k)/bin_area(k)),real_outp_precision))

		end do
    close(20)
  end subroutine



  subroutine write_u_f()
    ! writes the absolute amplitude of fourier coefficients at specified loc  to file
    integer                             :: io_error = 0
    real(kind = real_outp_precision)    :: out_x,out_y,out_z
		character(len=1024) 		  					:: filename
    type(sfield)                        :: x_r_dummy
    type(sfield)                        :: x_i_dummy
    type(sfield)                        :: y_r_dummy
    type(sfield)                        :: y_i_dummy
    type(sfield)                        :: abs_r_dummy
    type(sfield)                        :: abs_i_dummy
    type(vfield)                        :: dummy
    type(sfield)                        :: abs_dummy
		character(len=50) 		  						:: suffix
		character(len=18),parameter					:: path ='./output/data/u_f/'
		write(suffix,"(I5,A8)")  int(state%t/write_intervall), ".u_f.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    !dummy = state%u_f
		!do i=0,xdim-1
	  !	do j=0,ydim-1
	  ! 		dummy%val(i,j,1) = dummy%val(i,j,1)*exp(shear*state%t*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
	  ! 		dummy%val(i,j,2) = dummy%val(i,j,2)*exp(shear*state%t*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
    !  end do
    !end do
    !dummy = rearrange_2Dspectrum(deal_mask(dummy))
    x_r_dummy%val(:,:) = log(abs(real(state%u_f%val(:,:,1),real_outp_precision)))
    x_i_dummy%val(:,:) = log(abs(real(aimag(state%u_f%val(:,:,1)),real_outp_precision)))
    y_r_dummy%val(:,:) = log(abs(real(state%u_f%val(:,:,2),real_outp_precision)))
    y_i_dummy%val(:,:) = log(abs(real(aimag(state%u_f%val(:,:,2)),real_outp_precision)))

    abs_dummy%val(:,:) = sqrt(real(state%u%val(:,:,1),rp)**2 + real(state%u%val(:,:,2),rp)**2)
    call transform(abs_dummy%val,abs_dummy%val,1,shearing,sheartime)
    abs_r_dummy%val(:,:) = log(abs(real(abs_dummy%val(:,:),real_outp_precision)))
    abs_i_dummy%val(:,:) = log(abs(real(aimag(abs_dummy%val(:,:)),real_outp_precision)))

    x_r_dummy = rearrange_2Dspectrum(deal_mask(x_r_dummy))
    x_i_dummy = rearrange_2Dspectrum(deal_mask(x_i_dummy))
    y_r_dummy = rearrange_2Dspectrum(deal_mask(y_r_dummy))
    y_i_dummy = rearrange_2Dspectrum(deal_mask(y_i_dummy))
    abs_r_dummy = rearrange_2Dspectrum(deal_mask(abs_r_dummy))
    abs_i_dummy = rearrange_2Dspectrum(deal_mask(abs_i_dummy))

		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_u_f!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) i,j,real(x_r_dummy%val(i,j),real_outp_precision),&      !3: real(u_x)
                          real(x_i_dummy%val(i,j),real_outp_precision),&      !4: imag(u_x)
                          real(y_r_dummy%val(i,j),real_outp_precision),&      !5: real(u_y)
                          real(y_i_dummy%val(i,j),real_outp_precision),&      !6: imag(u_y)
                          real(abs_r_dummy%val(i,j),real_outp_precision),&    !7: real(abs_u)
                          real(abs_i_dummy%val(i,j),real_outp_precision)      !8: imag(abs_u)
			end do
		end do
    close(20)
  end subroutine

  subroutine write_chem_f()
    !write fourier chem-field to file. Note path variable
    integer                             :: io_error = 0,remap
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/chem_f/'
		write(suffix,"(I5,A11)")  int(state%t/write_intervall), ".chem_f.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)

    dummy = state%chem_f
		!do i=0,xdim-1
	  !	do j=0,ydim-1
	  ! 		dummy%val(i,j) = dummy%val(i,j)*exp(shear*state%t*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
    !  end do
    !end do
    dummy = rearrange_2Dspectrum(deal_mask(dummy))
    !dummy = deal_mask(state%chem_f)

		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_chem_f!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) i,j,log(abs(real(dummy%val(i,j),real_outp_precision))),&
                          log(abs(real(aimag(dummy%val(i,j)),real_outp_precision)))
			end do
		end do
    close(20)
  end subroutine

  subroutine write_chem_f_remap()
    !write fourier chem-field to file. Note path variable
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=27),parameter					:: path ='./output/data/chem_f_remap/'
		write(suffix,"(I5,A17)")  int(state%t/write_intervall), ".chem_f_remap.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)

    dummy%val = cmplx(1.0_rp,1.0_rp,rp)
    dummy%val = dealiase_field(dummy%val)

    do i =0,xdim-1
      do j =0,ydim-1
        if(dummy%val(i,j) .NE.cmplx(0.0_rp,0.0_rp,rp)) then
            dummy%val(i,j) = cmplx(0.0,0.0)
          else
          dummy%val(i,j) = state%chem_f%val(i,j)
        end if
      end do
    end do

    dummy = rearrange_2Dspectrum((dummy))
    !dummy = deal_mask(state%chem_f)
    !dummy = rearrange_2Dspectrum(deal_mask(dummy))

    !do i =0,xdim-1
    !  do j =0,ydim-1
    !    if((i>xdim/6 .AND.i<(5*xdim/6).AND.(j>ydim/6 .AND.j<(5*ydim/6)))) then
    !    dummy%val(i,j) =cmplx(0.0,0.0)
    !    end if
    !  end do
    !end do


		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_chem_f_remap!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) i,j,abs(real(dummy%val(i,j),real_outp_precision)),&
                          abs(real(aimag(dummy%val(i,j)),real_outp_precision))
			end do
		end do
    close(20)
  end subroutine

  subroutine write_temp_f_remap()
    !write fourier chem-field to file. Note path variable
    integer                             :: io_error = 0,remap
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=27),parameter					:: path ='./output/data/temp_f_remap/'
		write(suffix,"(I5,A17)")  int(state%t/write_intervall), ".temp_f_remap.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)

    dummy%val(:,:) = cmplx(1.0_rp,1.0_rp,rp)
    dummy%val = dealiase_field(dummy%val)

    do i =0,xdim-1
      do j =0,ydim-1
        if(dummy%val(i,j).NE.cmplx(0.0_rp,0.0_rp,rp)) then
            dummy%val(i,j) = cmplx(0.0,0.0)
          else
            dummy%val(i,j) = state%temp_f%val(i,j)
        end if
      end do
    end do

    dummy = rearrange_2Dspectrum((dummy))
    !dummy = deal_mask(dummy)
    !dummy = rearrange_2Dspectrum(deal_mask(dummy))

    !do i =0,xdim-1
    !  do j =0,ydim-1
    !    if((i>xdim/6 .AND.i<(5*xdim/6).AND.(j>ydim/6 .AND.j<(5*ydim/6)))) then
    !    dummy%val(i,j) =cmplx(0.0,0.0)
    !    end if
    !  end do
    !end do


		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_temp_f_remap!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) i,j,(abs(real(dummy%val(i,j),real_outp_precision))),&
                          (abs(real(aimag(dummy%val(i,j)),real_outp_precision)))
			end do
		end do
    close(20)
  end subroutine

  subroutine write_temp_f()
    !write fourier temp-field to file. analogous to write_chem
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/temp_f/'
		write(suffix,"(I5,A11)")  int(state%t/write_intervall), ".temp_f.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
		!do i=0,xdim-1
	  !	do j=0,ydim-1
	  ! 		dummy%val(i,j) = dummy%val(i,j)*exp(shear*state%t*state%ikx%val(i,j)*(real(j,rp)/real(ydim,rp))*Ly)
    !  end do
    !end do

    dummy = state%temp_f
    dummy = rearrange_2Dspectrum(deal_mask(dummy))

		open(unit=20,file=filename,status='replace',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_temp_f!'
		  do i=0,xdim-1
	    	do j=0,ydim-1
	  			write(20,*) i,j,log(abs(real(dummy%val(i,j),real_outp_precision))),&
                          log(abs(real(aimag(dummy%val(i,j)),real_outp_precision)))
			end do
		end do
    close(20)
  end subroutine
  subroutine write_fu_stat()
    !write measured u- diagnostics to file.
    integer                             :: io_error = 0
    type(vfield)                        :: vdummy

    type(sfield)                        :: fu_dummy
    type(sfield)                        :: fu_Nuk_dummy
    type(sfield)                        :: fu_diff_dummy
    type(sfield)                        :: fu_buo_dummy
    type(sfield)                        :: fu_shear_dummy

    type(sfield_stats)                  :: fu_stats
    type(sfield_stats)                  :: fu_Nuk_stats
    type(sfield_stats)                  :: fu_diff_stats
    type(sfield_stats)                  :: fu_buo_stats
    type(sfield_stats)                  :: fu_shear_stats

		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=22),parameter					:: path ='./output/data/fu_stat/'


    !fu
    vdummy%val = fu(state%u_f%val,state%u%val,state%temp_f%val,state%chem_f%val,sheartime)
    fu_dummy%val = cmplx((real(vdummy%val(:,:,1),rp)**2+(real(vdummy%val(:,:,2),rp)**2)),0.0_rp)
    fu_stats = measure_sfield_stats(fu_dummy)
    !fu_Nuk
    vdummy%val = fu_Nuk(state%u_f%val,state%u%val,sheartime)
    fu_Nuk_dummy%val = cmplx(sqrt(real(vdummy%val(:,:,1),rp)**2+(real(vdummy%val(:,:,2),rp)**2)),0.0_rp)
    fu_Nuk_stats = measure_sfield_stats(fu_Nuk_dummy)
    !fu_diff
    vdummy%val = fu_diff(state%u_f%val,sheartime)
    fu_diff_dummy%val = cmplx(sqrt(real(vdummy%val(:,:,1),rp)**2+(real(vdummy%val(:,:,2),rp)**2)),0.0_rp)
    fu_diff_stats= measure_sfield_stats(fu_diff_dummy)
    !fu_buo
    vdummy%val = fu_buo(state%u_f%val,state%temp_f%val,state%chem_f%val,sheartime)
    fu_buo_dummy%val = cmplx(sqrt(real(vdummy%val(:,:,1),rp)**2+(real(vdummy%val(:,:,2),rp)**2)),0.0_rp)
    fu_buo_stats= measure_sfield_stats(fu_buo_dummy)
    !fu_shear
    vdummy%val = fu_shear(state%u_f%val,sheartime)
    fu_shear_dummy%val = cmplx(sqrt(real(vdummy%val(:,:,1),rp)**2+(real(vdummy%val(:,:,2),rp)**2)),0.0_rp)
    fu_shear_stats= measure_sfield_stats(fu_shear_dummy)

		write(suffix,"(A11)") "fu_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    if(state%step>=1) then
		open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_fu_stat!'
	 	write(20,*)               state%step,                           & !1
                              state%t,                              & !2

                              fu_stats%val(0),                      & !3  max
                              fu_stats%val(1),                      & !4  min
                              fu_stats%val(2),                      & !5  mean
                              fu_stats%val(3),                      & !6  std-dev
                              fu_stats%val(4),                      & !7  variance

                              fu_Nuk_stats%val(0),                  & !8  max
                              fu_Nuk_stats%val(1),                  & !9  min
                              fu_Nuk_stats%val(2),                  & !10 mean
                              fu_Nuk_stats%val(3),                  & !11 std-dev
                              fu_Nuk_stats%val(4),                  & !12 variance

                              fu_diff_stats%val(0),                 & !13 max
                              fu_diff_stats%val(1),                 & !14 min
                              fu_diff_stats%val(2),                 & !15 mean
                              fu_diff_stats%val(3),                 & !16 std-dev
                              fu_diff_stats%val(4),                 & !17 variance

                              fu_buo_stats%val(0),                  & !18 max
                              fu_buo_stats%val(1),                  & !19 min
                              fu_buo_stats%val(2),                  & !20 mean
                              fu_buo_stats%val(3),                  & !21 std-dev
                              fu_buo_stats%val(4),                  & !22 variance

                              fu_shear_stats%val(0),                & !23 max
                              fu_shear_stats%val(1),                & !24 min
                              fu_shear_stats%val(2),                & !25 mean
                              fu_shear_stats%val(3),                & !26 std-dev
                              fu_shear_stats%val(4)                   !27 variance

    close(20)
    end if
  end subroutine
  subroutine write_u_stat()
    !write measured u- diagnostics to file.
    integer                             :: io_error = 0
    type(sfield)                        :: vx_dummy
    type(sfield)                        :: vy_dummy
    type(sfield)                        :: v_abs_dummy
    type(sfield_stats)                  :: vx_stats
    type(sfield_stats)                  :: vy_stats
    type(sfield_stats)                  :: v_abs_stats
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/u_stat/'
		write(suffix,"(A11)") "u_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    if(state%step>=1) then

    vx_dummy%val = state%u%val(:,:,1)
    vy_dummy%val = state%u%val(:,:,2)
    v_abs_dummy%val = sqrt(real(state%u%val(:,:,2),rp)**2 +real(state%u%val(:,:,1),rp)**2)
    vx_stats = measure_sfield_stats(vx_dummy)
    vy_stats = measure_sfield_stats(vy_dummy)
    v_abs_stats = measure_sfield_stats(v_abs_dummy)


		open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error) 
    if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_u_stat!'
	 	write(20,*)               state%step,                           & !1
                              state%t,                              & !2
                              vx_stats%val(0),                      & !3  max
                              vx_stats%val(1),                      & !4  min
                              vx_stats%val(2),                      & !5  mean
                              vx_stats%val(3),                      & !6  std-dev
                              vx_stats%val(4),                      & !7  variance

                              vy_stats%val(0),                      & !8  max
                              vy_stats%val(1),                      & !9  min
                              vy_stats%val(2),                      & !10 mean
                              vy_stats%val(3),                      & !11 std-dev
                              vy_stats%val(4),                      & !12 variance

                              v_abs_stats%val(0),                   & !13 max
                              v_abs_stats%val(1),                   & !14 min
                              v_abs_stats%val(2),                   & !15 mean
                              v_abs_stats%val(3),                   & !16 std-dev
                              v_abs_stats%val(4)                      !17 variance

    close(20)
    end if
  end subroutine
  subroutine write_E_stat()
    !write measured E-diagnostics to file.
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/E_stat/'
		write(suffix,"(A11)") "E_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    if(state%step>=1) then
		  open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error)
      if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_E_stat!'
	 	  write(20,*) state%step,state%t,measure_Ekin(),measure_Epot(),measure_Ekin()+measure_Epot()
      close(20)
    end if
  end subroutine

  subroutine write_sys_stat()
    !write measured system wide diagnostics to file.
    integer                             :: io_error = 0
    type(sfield)                        :: int_dummy_f
    type(sfield)                        :: int_dummy
    type(sfield)                        :: int1_dummy_f
    type(sfield)                        :: int1_dummy
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=23),parameter					:: path ='./output/data/sys_stat/'
    real(kind = real_outp_precision)    :: y_aperiodicity_measure      ! how periodic is the field in y-dir
    real(kind = real_outp_precision)    :: x_aperiodicity_measure      ! how periodic is the field in y-dir
		write(suffix,"(A12)") "sys_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    if(state%step>=1) then
		  open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error)
      if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_sys_stat!'

      !call dfftw_execute_dft(ifull2D,int_dummy_f%val,int_dummy%val)
      !divergence physical
      call transform( state%u_f%val, int_dummy%val,-1,1,sheartime) 
      call transform(int_dummy%val,int_dummy_f%val,1,0,sheartime) 
      int_dummy_f%val(:,:) = state%ikx%val(:,:)*state%u_f%val(:,:,1) &
                                +state%iky%val(:,:)*state%u_f%val(:,:,2) 
      !divergence brucker
      int1_dummy_f%val(:,:) = state%ikx_bar%val(:,:)*state%u_f%val(:,:,1) &
                               +state%iky_bar%val(:,:)*state%u_f%val(:,:,2) 

      ! measure how similar the upper and lower end of sim-box are
      y_aperiodicity_measure = 0.0_rp
      x_aperiodicity_measure = 0.0_rp
      do i = 0,xdim-1
        y_aperiodicity_measure = y_aperiodicity_measure&
                    + real(abs(state%u%val(i,ydim-1,1)-state%u%val(i,0,1)),real_outp_precision)
        y_aperiodicity_measure = y_aperiodicity_measure &
                    + real(abs(state%u%val(i,ydim-1,2)-state%u%val(i,0,2)),real_outp_precision)
      end do
      do j = 0,ydim-1
        x_aperiodicity_measure = x_aperiodicity_measure&
                    + real(abs(state%u%val(0,j,1)-state%u%val(xdim-1,j,1)),real_outp_precision)
        x_aperiodicity_measure = x_aperiodicity_measure &
                    + real(abs(state%u%val(0,j,2)-state%u%val(xdim-1,j,2)),real_outp_precision)
      end do
     y_aperiodicity_measure = 1.0_rp/(y_aperiodicity_measure/(measure_u_rms()*real(xdim)))
     x_aperiodicity_measure = 1.0_rp/(x_aperiodicity_measure/(measure_u_rms()*real(ydim)))

	 	  write(20,*) state%step,                                               & !1
                  state%t,                                                  & !2
                  maxval(real(int_dummy_f%val,real_outp_precision)),        & !3 maximum of divergence
                  shear,                                                    & !4 strength of shear
                  dt,                                                       & !5
                  maxval(real(int1_dummy_f%val,real_outp_precision)),         & !6 max brucker div 
                  measure_av(state%s_dummy_f%val),                          & !7 average vorticity 
                  y_aperiodicity_measure,         &                           !8 y_aperiodicity
                  x_aperiodicity_measure                                      !8 x_aperiodicity

      close(20)
    end if
  end subroutine

  subroutine write_T_stat()
    !write measured T-diagnostics to file.
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
    type(sfield_stats)                  :: t_stats
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/T_stat/'
		write(suffix,"(A10)") "T_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    !write(*,*) 'sub write_T_stat() has been called'
    if(state%step>=1) then

    t_stats = measure_sfield_stats(state%temp)
		  open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error)
      if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_T_stat!'
	 	  write(20,*) state%step,                                       &! 1
                  state%t,                                          &! 2  
                  t_stats%val(0),                                   &! 3   
                  t_stats%val(1),                                   &! 4  
                  t_stats%val(2),                                   &! 5  
                  t_stats%val(3),                                   &! 6  
                  t_stats%val(4)                                     ! 7
      close(20)
    end if
  end subroutine

  subroutine write_C_stat()
    !write measured T-diagnostics to file.
    integer                             :: io_error = 0
    type(sfield)                        :: dummy
    type(sfield_stats)                  :: c_stats
		character(len=1024) 		  					:: filename
		character(len=50) 		  						:: suffix
		character(len=21),parameter					:: path ='./output/data/C_stat/'
		write(suffix,"(A11)") "C_stat.dat"
    suffix = trim(adjustl(suffix))
    filename = path //suffix
		filename = adjustl(filename)
		filename = trim(filename)
    if(state%step>=1) then
    c_stats = measure_sfield_stats(state%chem)
		  open(unit=20,file=filename,status='unknown',position='append',action='write',iostat=io_error)
      if(io_error .NE. 0) write(*,*) 'ERROR: could not open file in sub write_C_stat!'
	 	  write(20,*) state%step,                                               &! 1
                  state%t,                                                  &! 2
                  c_stats%val(0),                                           &! 3
                  c_stats%val(1),                                           &! 4
                  c_stats%val(2),                                           &! 5
                  c_stats%val(3),                                           &! 6
                  c_stats%val(4)                                             ! 7
      close(20)
    end if
  end subroutine

  function deal_mask(u_f)
  	! show the cutoff frequencies highlighted in output array. only for visualisation
  	type(sfield)				:: deal_mask
  	type(sfield)				:: u_f
  	complex(kind = rp)											:: max_abs
  	max_abs = cmplx(maxval(abs(real(u_f%val))),maxval(abs(real((aimag(u_f%val))))),rp)
  	deal_mask%val = u_f%val

  	deal_mask%val(xdim/3 ,:           ) 	= max_abs			
  	deal_mask%val(2*xdim/3+1,:        ) 	= max_abs			
  	deal_mask%val(:       ,2*ydim/3+1 ) 	= max_abs			
  	deal_mask%val(:       ,ydim/3     ) 	= max_abs
  end function 

  function rearrange_2Dspectrum(arr_f)
	  type(sfield)				:: arr_f
	  type(sfield)				:: rearrange_spec
	  type(sfield)				:: rearrange_2Dspectrum
  
	  do i =0,xdim/2
	  	do j =0,ydim/2
	  		rearrange_spec%val(i,j) = arr_f%val(i+xdim/2-1,j+ydim/2-1)  !lower left
	  	end do
	  end do
	  do i =xdim/2+1,xdim-1
	  	do j =ydim/2+1,ydim-1
	  		rearrange_spec%val(i,j) = arr_f%val(i-xdim/2-1,j-ydim/2-1) 	!upper right
	  	end do
	  end do
	  do i =0,xdim/2
	  	do j =ydim/2+1,ydim-1
	  		rearrange_spec%val(i,j) = arr_f%val(i+xdim/2-1,j-ydim/2-1) ! upper left
	  	end do
	  end do
	  do i =xdim/2+1,xdim-1
	  	do j =0,ydim/2
	  		rearrange_spec%val(i,j) = arr_f%val(i-xdim/2-1,j+ydim/2-1) ! lower right
	  	end do
	  end do
    rearrange_2Dspectrum = rearrange_spec
  end function

end module
