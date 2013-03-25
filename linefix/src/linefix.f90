

      Program linecheck
     
      implicit none
      
      integer j,nskip, strlen, nfix,istat
      logical ResOK, openfile
      character*128 charline,newline
      CHARACTER FNAME*256

      nskip = 20730477 


      ResOK = OpenFile(1,'Open input File',fname, &
     &    "Input file(*.dat),*.dat;All files(*.*),*.*;")

      if(.not.resOK) then
        write(*,*), 'Unable to open input file'
        call exit(71)
      endif

      open(21,file='linefix.dat',status='unknown')

      write(*,*) ' Enter nskip'
      read(*,*) nskip

      do j=1,nskip
        read(1,'(a)') charline
        strlen = len_trim(charline)
        write(21,'(a)') charline(1:strlen)
      enddo

      do
        read(1,'(a)') charline
        strlen = len_trim(charline)
        write(*,'(a)') charline(1:strlen)
        write(*,'(a)') ' enter 1 to fix, 0 for next, -1 to continue'
        read(*,*) nfix
        if(nfix.eq.1) then
          write(*,'(a)') ' enter corrected line'
          read(*,'(a)') newline
          strlen = len_trim(newline)
          write(21,'(a)') newline(1:strlen)
        elseif(nfix.eq.0) then
          write(21,'(a)') charline(1:strlen)
        else
          write(21,'(a)') charline(1:strlen)
          exit
        endif
      enddo

      do
        read(1,'(a)',iostat=istat) charline
        if(istat.eq.0) then
          strlen = len_trim(charline)
          write(21,'(a)') charline(1:strlen)
        elseif(istat.lt.0) then !end of file
          close(1)
          close(21)
          exit
        elseif(istat.gt.0) then !read error
          write(*,'(a)') ' read error on line'
          write(*,'(a)') charline
        endif
      enddo

      stop
      end
