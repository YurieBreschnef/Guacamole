module const
! constants-module: contains simulation wide constants of parameter type
! if variables are not flagged parameter, they may be temporarily changed (e.g. stepwidth)
! numerical parameters are also given here.
  use ISO_C_BINDING 
  use omp_lib
  implicit none
include "fftw3.f03"
integer,parameter				::rp			= 8			!real-precision
integer,parameter				::real_outp_precision	= 4			!output number precision
integer,parameter				::ip			= 4			!integer-precision

	integer,parameter			::fftw_plan_thoroughness= FFTW_MEASURE
	! possible also FFTW_MEASURE
	integer(kind=ip),parameter		::xdim			= 256 
	integer(kind=ip),parameter		::ydim			= 256   

	integer(kind = ip),parameter		::seed			= 111	! seed for random init
	integer(kind = ip),parameter		::maxfiles		= 300 ! maximum no of output files per type
	integer(kind = ip),parameter		::measure_every		= 50 ! measure diagnostics every X steps
	integer(kind = ip),parameter		::debuglevel		= 1			
  ! level 0: no output, level 1: short, level 2: extensive
	real(kind = rp),parameter		::pi		= 3.1415926535897932384626433833_rp
	
	real(kind = rp),parameter		::Lx		  = 2.0_rp *pi !50.0_rp
	real(kind = rp),parameter		::Ly		  = 2.0_rp *pi !50.0_rp

	complex(kind = rp),parameter	:: imag			= (0.0_rp,1.0_rp)


	integer(kind = ip)				:: steps		
	integer(kind = ip)				:: i,j,k,l,main_stp      !used for all kinds of loops

	real(kind = rp),parameter			      :: tmax                      = 100.00_rp
	real(kind = rp)					      :: dt			   = 5.0e-3_rp

	real(kind = rp)					      :: dt_max                    = 1.0e-3_rp
	real(kind = rp)					      :: dt_min                    = 1.0e-6_rp
	real(kind = rp)					      :: max_step_error            = 1.0e-6_rp
	real(kind = rp)					      :: stepwidth_adjustment_rate = 0.05_rp
	real(kind = rp)					      :: write_intervall           = tmax/real(maxfiles,rp)
	real(kind = rp)					      :: last_written              = 0.0_rp

	real(kind = rp)					      :: dt_2           !(1/2) *dt
	real(kind = rp)					      :: dt_3           !(1/3) *dt
	real(kind = rp)					      :: dt_4           !(1/3) *dt
	real(kind = rp)					      :: dt_8           !(1/3) *dt
	real(kind = rp)					      :: dt_34          !(3/4) * dt
	real(kind = rp)					      :: dt_29          !(2/9) * dt
	real(kind = rp)					      :: dt_49          !(4/9) * dt
	real(kind = rp)					      :: dt_724         !(7/24) * dt

	real(kind = rp),parameter			      :: a_0=  11.0_rp / 6.0_rp ! used for ABBDF3
	real(kind = rp),parameter			      :: a_1=  -3.0_rp
	real(kind = rp),parameter			      :: a_2=   3.0_rp / 2.0_rp
	real(kind = rp),parameter			      :: a_3=  -1.0_rp / 3.0_rp

	real(kind = rp),parameter			      :: b_0=   3.0_rp 
	real(kind = rp),parameter			      :: b_1=  -3.0_rp
	real(kind = rp),parameter			      :: b_2=   1.0_rp

        !shear related parameters
	real(kind = rp)               :: shear    = 0.050_rp
        integer(kind = ip)            :: shearing = 1
	real(kind = rp)               :: sheartime= 0.0_rp
	real(kind = rp)               :: T_rm 
        !TODO. store ky_bar_max in array so it is not recalculated every time and reset on set_ik_bar
	real(kind = rp)					      :: ky_max 
	real(kind = rp)					      :: ky_min
	real(kind = rp)					      :: kx_max 
	real(kind = rp)					      :: kx_min
  
  integer(kind = ip)            :: benchmarking = 0

  integer(kind = ip)             :: remapping = 1
  integer(kind = ip)             :: remapping_rate = 10 


  integer(kind = ip)             :: threads   = 0 
  integer(kind = ip)             :: my_thread_id= 0   !thread specific values

  integer(kind = ip)             :: my_x_start
  integer(kind = ip)             :: my_x_end
  integer(kind = ip)             :: my_y_start
  integer(kind = ip)             :: my_y_end

  real(kind = rp),parameter                       :: D_visc   = 1.0_rp*0.100_rp 
  real(kind = rp),parameter			:: D_therm  = 1.0_rp*0.01000_rp
  real(kind = rp),parameter			:: D_comp   = 1.0_rp*0.003_rp
  
  real(kind = rp),parameter			:: B_therm  = 1.1_rp
  real(kind = rp),parameter			:: B_comp   = 1.0_rp
  
  real(kind = rp),parameter			:: S_therm  = 1.0_rp  
  real(kind = rp),parameter			:: S_comp   = 1.0_rp 

	!real(kind = rp),parameter                       :: D_visc   = 1.0_rp*0.100_rp 
	!real(kind = rp),parameter			:: D_therm  = 1.0_rp*0.01000_rp
	!real(kind = rp),parameter			:: D_comp   = 1.0_rp*0.003_rp

	!real(kind = rp),parameter			:: B_therm  = 1.1_rp
	!real(kind = rp),parameter			:: B_comp   = 2.0_rp

	!real(kind = rp),parameter			:: S_therm  = 2.0_rp  
	!real(kind = rp),parameter			:: S_comp   = 1.0_rp 

end module
