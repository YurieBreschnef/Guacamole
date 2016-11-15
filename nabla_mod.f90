module nabla
  ! place for vector operations
  use sys_state
  use const
  use trafo
  implicit none

  contains

  function crossp(arr1,arr2)
  ! conducts standard 3D crossproduct of given input arguments
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            ::arr1 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            ::arr2 
  complex(kind=rp),dimension(0:xdim-1,0:ydim-1,1:2)            ::crossp
  crossp(:,:,1) = (arr1(:,:,1))*(arr2(:,:,2)) &
                       -(arr1(:,:,2))*(arr2(:,:,1))
  crossp(:,:,2) = (arr1(:,:,2))*(arr2(:,:,1)) &
                       -(arr1(:,:,1))*(arr2(:,:,2))

  !IF(ANY(IsNaN(real(crossp))))  then
  !  write(*,*) 'func crossp(): NAN detected in output array'
  !  stop
  !end if
  end function 

end module
