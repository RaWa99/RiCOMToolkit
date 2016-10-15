 
      Program zoom_sigma

!     This program modifies a uniform sigma coordinate profile with zooming
!     in the upper and/or lower part.

      implicit none

      integer :: k,npv
      real :: du,dl,sig(1000)
      
! *** enter data
      write(*,*) ' Enter number of vertical points'
      read(*,*) npv
      write(*,*) ' Enter UPPER zoom factor (>=0)'
      read(*,*) du
      write(*,*) ' Enter LOWER zoom factor (>=0)'
      read(*,*) dl
      
! *** calculate sigma profile ordered from top to bottom     
      do k=1,npv
        sig(npv-k+1) = (tanh((du+dl)*(float(k-1)/float(npv-1))-dl)+tanh(dl))/(tanh(du)+tanh(dl)) - 1.
      enddo
      
! *** write profile for input to RiCOM, ordered from top to bottom     
      open(unit=33,file='sprofile.dat',status='unknown')
      do k=1,npv
        write(33,'(1x,f12.6)') sig(k)
      enddo
      
      
      stop
      end

