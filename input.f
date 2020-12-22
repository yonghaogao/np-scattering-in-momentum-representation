      module inputfile 
      
      contains
!***********************************************************************
        subroutine input()
!   defines input for potential, grid, and t-matrix
!***********************************************************************
        use grids_module
        use precision
        use system
        use gauss
        use pot
        use constants
        implicit none
        real*8 :: massp, massn 
        integer :: ir 

        namelist /global/ l,massp,massn
        namelist /grids/ p1,p2,p3,np1,np2,rmax,nnu

        l=0;massp=1.0078;massn=1.0087
      
        read (5,nml=global)
        write (6,nml=global)
        mu=massp*massn/(massp+massn)*amu/hbarc

        p1=2.0_dpreal;p2=5.0_dpreal;p3=10.0_dpreal
        np1=30;np2=30;rmax=15;nnu=500
        read (5,nml=grids)
        write (6,nml=grids)

        write (6,*)
        write (6,*) '============  Output  ==================='
        write (6,*)


        ! define momentum grids
        np=np1+np2
        allocate(pp(1:np),pw(1:np))
        ! call trns map
        !  ! get a grid of mometum points
        ! NP1/2 points will be in [0,P1]
        ! NP1/2 points will be in [P1,P2]
        ! NP2   points will be in [P2,P3]
        call trns(np1,np2,np,p1,p2,p3,pp(1:np),pw(1:np))

        ! define r-space grids
        allocate(rr(1:nnu),rrw(1:nnu))
        call gauleg(nnu,0.0_dpreal,rmax,rr,rrw)
C       call simpson(nnu,rmax,rr,rrw)  
          
C         do ir=1, np
C          write(99,*) pp(ir), 0
C        end do 

        end subroutine
        
        
      end module 