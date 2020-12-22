      module bound2b
       use grids_module
       use constants
       use pot, only : vpp_store
       implicit none
       real*8,allocatable,dimension(:,:) :: A
       private :: A

      contains

      subroutine bound()
      use grids_module
      use system
      implicit none
      real*8,dimension(1:np) :: wf,wf1,wf2

      real*8 :: lambda1, lambda2, lambda
      real*8 :: E1,E2,E
      real*8 :: sumdg,sum1,fz
      integer :: i

      if (.not. allocated(A)) allocate(A(1:np,1:np))


       fz = 1.d0 ! use for sign of wave function

       E1=-1.0/hbarc ! 给定secant方法的初始值
!!!!secant method to find E1, iteration method to solve lambda

cccccccccccccccccc完善secant method寻找正确的束缚态能量cccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            E2=-2.0/hbarc
            do while(abs(E2-E1)>1.E-6)
            call eigenvalue(E1,lambda1,wf1)
            call eigenvalue(E2,lambda2,wf2)
            E=E2-(E1-E2)*(lambda2-1)/(lambda1-lambda2)
            E1=E2
            E2=E
            end do
        write(*,*) "bound state energy is", E*hbarc

       sumdg=0.d0
       do i=1,np
       sumdg=sumdg+wf2(i)*wf2(i)
       end do

       write (*,*)
       write (*,*) ' norm of eigenvector from dgeev :',sumdg

!       normalize wavefunction
      sum1 = 0.d0
       do i=1,np
       sum1=sum1+wf(i)*wf(i)*pw(i)*pp(i)*pp(i)
       enddo

      sum1=1.d0/sqrt(sum1)

      do i=1,np
      wf(i)=wf(i)*sum1
      end do

c        sign of wave function

      if (wf(3).lt.0.d0) fz=-1.d0

      do i=1,np
      wf(i)=wf(i)*fz
      write(13,*) pp(i), wf(i)
      end do


     

      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine eigenvalue(E,lambda,wf)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use grids_module
      use system
      implicit none
      real*8 :: E,lambda
      integer :: n,i
      integer :: info,lwork,loc(1)
      real*8, dimension(1:np) :: wr,wi,wf
      real*8, dimension(1:4*np) :: work
      real*8,dimension(1:np,1:np) :: vl,vr
   

      call Amat_bound(E)
! prepare calling the LAPACK routine DGEEV
      lwork=4*np !  the value is taken from the documentation of DGEEV
! the routine will destroy the original array A
! and return the real parts of the eigenvalues in WR and the imaginary parts in WI
! VR(:,i) contain the eigenstate (=wave function) of the right eigenvector of eigenvalue i
! VL(:,i) contain the eigenstate of the left eigenvector of eigenvalue i


      call DGEEV('N','V',NP,A,NP,WR,WI,VL,NP,VR,NP,WORK,LWORK,INFO)
                                                  

      if (info/=0) then
         write(*,*) "problem happened when calling DGEEV"
         stop
      end if

! locate maximal eigenvector
! here assume that this eigenvector is real
! then return this eigenvalue in ETA
! and the eigenvector in the global variable WF (which should be of dimension NP)

C     do n=1,nex-1
C      loc=MAXLOC(wr) ! return value is type of array, must defined with array type
C      wr(loc(1))=0
C     end do
       loc=MAXLOC(wr)

       lambda=wr(loc(1))
       wf=vr(:,loc(1))

      end subroutine

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Amat_bound(E)
       use system
       use grids_module
       use precision
       implicit none
       real*8 :: g0,E
       integer :: ip,jp
       A=0.0d0! 用来储存计算后的A矩阵，A(1:np,1:np)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC请计算A矩阵CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC并画出对角线部分CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do jp=1,np
        do ip=1,np
         A(ip,jp)=pw(jp)*vpp_store(ip,jp)*pp(jp)**2/(E-pp(ip)**2/(2*mu))
        end do
      end do
      !  output A矩阵对角线元素
      !   do ip=1,np
      !      write(999,*)A(ip,ip)
      !   end do
      end subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      end module
