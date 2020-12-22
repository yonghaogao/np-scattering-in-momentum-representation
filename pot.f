!     module for potential 
      module pot 
      implicit none 
      real*8,allocatable, dimension(:,:) :: vpp_store ! potential in momentum space 
      real*8,allocatable,dimension(:) :: vr ! potential in coordinate space 
      contains
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
      subroutine potsys()
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use precision
       use system
       use grids_module
       implicit none
       
       integer :: ip,jp ! index
       
       ! calculate bessel function which used for Fourier transformation
       if (.not. allocated(fl)) then
        allocate(fl(1:nnu,1:np,0:l))  ! bessel function
        call bessel() ! calculate bessel function 
       end if
       if (.not. allocated(vpp_store)) allocate(vpp_store(1:np,1:np))
       
       
        if (.not. allocated(vr)) allocate(vr(1:nnu))
         call vpotr() ! calculate potential in R-sapce
         
         ! Fourier transform to P-space 
         do ip=1,np
           do jp=1,np
             call vpotpp(ip,jp,vpp_store(ip,jp))
           end do       
         end do 
       
       do ip=1,np 
          write(12,*) pp(ip) , vpp_store(ip,ip)
       end do 
       write(12,*) "&"
 
      end subroutine

      
    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  vpotpp(ip,jp,vpp)
       use constants
       use grids_module
       use system
       use precision
       implicit none
       integer :: ip,jp
       real*8 :: vpp
       integer :: ir
       real*8, pointer,dimension(:) :: flkk,flk
       real*8:: temp
       vpp=0.0d0 ! 来储存傅里叶变换后的势，vpp(1:np,1:np)
       flkk=>fl(:,ip,l) ! Bessel function 
       flk =>fl(:,jp,l) ! Bessel function 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCC请写出傅里叶变换的程序并得到CCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCC动量表象下的势并画出对角线部分CCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC即pp(ip),vpp(ip,ip)CCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       temp=0.
       do ir=1,nnu
       temp=temp+flkk(ir)*vr(ir)*flk(ir)*rrw(ir)
       enddo
       vpp=2./pi/pp(ip)/pp(jp)*temp/hbarc
       
      end subroutine
      
c----------------------------------------------------------------------
    
       subroutine vpotr()
       use grids_module
       implicit none
       integer :: ir 
       real*8 :: r,v0

       vr=0.0d0
       V0=-72.167 ! 给定任意初始势能深度
       do ir=1, nnu
         r=rr(ir) 
         vr(ir) = gausspot(r,v0,0.0d0,1.484d0)  ! 注意此处没有加库伦力
         write(11,*)rr(ir), vr(ir)
      end do
      write(11,*) "&"

      end subroutine
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  bessel()
       use coulfunc
       use grids_module
       use system
       implicit none
       real*8,dimension(0:l) :: GC,FCP,GCP
       real*8,pointer,dimension(:) :: FC
       real*8 :: rho
       integer :: ir, ip,ifail
       ifail=0
       do ip=1,np
         do ir=1,nnu
            rho=pp(ip)*rr(ir)
            fc=>fl(ir,ip,0:l)
            call coul90(rho,0.0d0,0.0d0,l,fc,gc,fcp,gcp,0,ifail)
            if (ifail/=0) then
               write(*,*) 'coul90: ifail=',ifail; stop
            endif
         end do
! problem happened when l is larger than 70, when rho is less than 0.01
! the value of fc will be NAN
!         fc=>fl(1,1,0:l)
!         rho=kk(1)*rr(1)
!         call coul90(0.001d0,0.0d0,0.0d0,l,fc,gc,fcp,gcp,0,ifail)
!         write(*,*) "fl=", fl(1,1,0)
       end do

      end subroutine  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             functions of different types of potential
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *** Woods-Saxon (volume)
       function ws(r,v0,r0,a)
       implicit none
       real*8 r,v0,r0,a,ws,aux
       ws=0d0
        if (abs(v0).lt.1e-6) return
        if (a>1e-6) then
           aux=exp(-(r-r0)/a)
         ws=v0/(1d0+1d0/aux)
        else
         write(0,*)'WS: a too small!'
        end if
        return
        end function
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c *** Spin-orbit with WS derivative
c     This form factor will be then multiplied by l.s
c     We use fresco definiton for spin-orbit
      function wsso(r,v0,r0,a)
        implicit none
        real*8 r,v0,r0,a,wsso,aux,conls
        parameter(conls=2.0)
        wsso=0.0d0
         if (r<1e-6) r=1e-6
         if (a>1e-6) then
          aux=exp(-(r-r0)/a)
          wsso=-conls*v0*aux/(1d0+aux)**2/a/r
         else
          write(0,*)'WS spin-orbit : a too small!'
         endif
         return
      end function


c *** WS derivative
      function dws(r,v0,r0,a)
      implicit none
      real*8 r,v0,r0,a,dws,aux
        if (r<1e-6) r=1e-6
           if (a>1e-6) then
             aux=exp(-(r-r0)/a)
             dws=-v0*aux/(1d0+aux)**2/a
           else
             write(0,*)'derivative WS: a too small!'
             write(*,*)"a=",a
             stop
         endif
         return
      end function
c *** Gaussian
      function gausspot(r,v0,r0,a)
      implicit none
       real*8 r,v0,r0,gausspot,a
         if (a > 1e-6) then
           gausspot=V0*exp(-(r-r0)**2/a**2)
             else
               write(*,*) "a=",a
               write(*,*)'a too small in gausspot!'
               stop
         endif
         return
      end function

c-----------------------------------------------------------------------
c *** Spin-orbit with Gauss derivative
c     This form factor will be then multiplied by l.s
c     We use fresco definiton for spin-orbit
        function gausder(r,v0,r0,a)
      implicit none
      real*8 r,v0,r0,a,gausder,conls,rh
      parameter(conls=2.0)
      gausder=0.0d0
       if (r<1e-6) r=1e-6
         if (a>1e-6) then
            rh=(r-r0)/a
            gausder=-exp(-rh**2)**rh*conls*v0/(r*a)
         else
           write(0,*)'Gauss spin-orbit : a too small!'
       endif
       return
      end function
c-----------------------------------------------------------------------
c Coulomb potential
      FUNCTION VCOUL(R,z12,Rc)
          use constants
          implicit none
          real*8 r,rc,rc2,aux,vcoul,z12

          RC2=RC*2d0
          aux=e2*Z12
          vcoul=0
          if (z12.lt.1e-4) return
          if (rc.lt.1e-6) rc=1e-6

          IF(R.GT.RC)GO TO 1
          VCOUL=AUX*(3.-(R/RC)**2)/RC2
          RETURN
1         VCOUL=AUX/R
          RETURN
        END

c------------------------------------------      

      
      
      
      
      end module 