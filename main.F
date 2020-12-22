!**********************************************************************
      program main
!**********************************************************************
      use grids_module  !  defines momentum grid
      use bound2b  !  solve bound state (uses LAPACK)
      use constants  !  defines constants used in the calculation
      use pot        !define the potential in momentum space
      use inputfile ! for input 
      implicit none
    
     

      call input() ! define the input parameters
      call potsys() ! get the potential in momentum space 
      call bound() ! calculate the bound state 
      call fkind()
      
      end program
!-----------------------


!***********************************************************************
       subroutine fkind()
cWrite output files units and names
!***********************************************************************
       use system
       character*40 flkind(9999)
       integer writf(9999),nwrit,i
       logical,dimension(1:9999) ::  written
       flkind(:) = ' '

       flkind(11)='pot in r-space'
       written(11)=.TRUE.
       flkind(12)='diagonal part of pot in p-space'
       written(12)=.TRUE.
       flkind(13)='wave-function in p-space'
       written(13)=.TRUE.



       nwrit = 0
       do i=1,9999
        if(written(i)) then
        flkind(i) = trim(flkind(i))//'.'
        nwrit = nwrit+1
        writf(nwrit) = i
        endif
       enddo
       write(*,990) (writf(i),flkind(writf(i)),i=1,nwrit)
990    format(/'  The following files have been created:',
     X  /(2x,2(i3,':',1x,a40)))
       return
       end subroutine


