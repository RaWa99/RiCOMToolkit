

      Program linecheck
     
      implicit none
      
      integer j,nskip, strlen
      logical ResOK, openfile
      character*128 charline
      CHARACTER FNAME*256

      nskip = 20730477 


      ResOK = OpenFile(1,'Open input File',fname, &
     &    "Input file(*.dat),*.dat;All files(*.*),*.*;")

      if(.not.resOK) then
        write(*,*), 'Unable to open input file'
        call exit(71)
      endif

      write(*,*) ' Enter nskip'
      read(*,*) nskip

      do j=1,nskip
        read(1,*)
      enddo

      do
        read(1,'(a)') charline
        strlen = len_trim(charline)
        write(*,'(a)') charline(1:strlen)
        read(*,*)
      enddo

      stop
      end
