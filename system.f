ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module grids_module
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 :: p1,p2,p3
       integer :: np1,np2,np
       real*8,allocatable,dimension(:):: pp,pw
       real*8,allocatable,target,dimension(:,:,:) :: fl
       real*8,allocatable,dimension(:) :: rr,rrw
       real*8 :: rmax ! used for the Fourier transform of potential
       integer :: nnu
      end module
c----------------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module system
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 :: mu
       integer ::l
      end module
c----------------------------------------------------------------------