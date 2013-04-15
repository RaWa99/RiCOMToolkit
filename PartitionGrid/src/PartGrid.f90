

      program partgrid
      
      use MainArrays
      
      implicit none
      
      logical quit
      integer i,j,k
      integer ncommon, nparts, etype, numflag, edgecut, istat
      integer options(5), objval
      integer, allocatable ::  eptr(:), eind(:), vwgt(:), vsize(:) 
      integer, allocatable :: epart(:), npart(:)
      real, allocatable :: tpwgts(:)
      
      call ReadGridFile(quit)
      
!  form element vectors for decomp routines
      
      write(*,*) ' Enter the number of partitions'
      read(*,*) nparts

!  call metis routines - v4
      
      numflag = 1
      etype = 1

      ALLOCATE (eind(4*ne),epart(ne),npart(np), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate metis arrays'
        call exit(71)
      endif

      do j=1,ne
        do k=1,3
          eind(3*(j-1)+k) = nen(k,j)
        enddo
      enddo

      call metis_partmeshdual(ne,np,eind,etype,numflag,nparts,&
                              edgecut,epart,npart)

                              
!  call metis routines - v5
      
!      ncommon = 2

!      ALLOCATE (eptr(ne+1),eind(4*ne),epart(ne),npart(np),&
!                vwgt(ne), vsize(ne),tpwgts(ne), STAT = istat )
!      if(istat.ne.0) then
!        write(*,*) 'FATAL ERROR: Cannot allocate metis arrays'
!        call exit(71)
!      endif

!      call metis_partmeshdual(ne,np,eptr,eind,vwgt,vsize,ncommon,nparts,&
!                              tpwgts,options,objval,epart,npart)

!  write output file of partitions(element number)
      open(21,file='epart.dat', status='unknown')
      do j=1,ne
        write(21,*) epart(j)
      enddo
      
      end
