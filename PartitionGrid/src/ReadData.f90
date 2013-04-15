  !***********************************************************************

! --------------------------------------------------------------------------*
!     Routines to read input files. 
! --------------------------------------------------------------------------*

      SUBROUTINE ReadGridFile( Quit )

! Purpose: Prompts user for the input grid file.
! Givens : None
! Returns: None
! Effects: User is repeatedly prompted for a grid file until a valid
!          filename is entered.  Grid file is opened, points are read,
!          and grid file is closed.

      use MainArrays

      implicit none

      LOGICAL Quit

      integer i,ii,j
      integer Fnlen, nunit, istat
      CHARACTER*256 fname
      character(80) Firstline

!------------------BEGIN------------------

      Quit = .FALSE.
      nunit = 8

      write(*,*) 'Enter Grid File name (*.ngh, *.xye)'
      read(*,'(a)') FNAME

      fnlen = len_trim( fname )
      OPEN(UNIT=nunit,file=FNAME(1:fnlen),status='old')

      READ(nunit,'(a)', err=9999, end=99999) Firstline

      if(firstline(1:4).eq."#NOD") then  !node file, new format
        write(*,*) 'Node file.. Wrong format for grid.','ReadGrid'
        Quit = .true.
        Return
      elseif(firstline(1:4).eq."#XYE") then  !xyz and element grid file, new format
        call ReadXYEData (nunit, Quit)
        return
      elseif(firstline(1:4).eq."#NGH") then  !neigh grid file, new format
        do
          READ(nunit,'(a)', err=9999, end=99999) Firstline
          if(firstline(1:1).ne."#") then    !comment lines
!           following line is internal read of firstline          
! - read offsets, scale factors, coordinate type
            READ( firstline, *, err = 9999, end=99999 ) x0off, y0off, scaleX, scaleY, igridtype
            exit
          endif
        enddo

! - read max number of nodes in this grid
        read(nunit,*,err=9999, end=99999) np

! - read max number of neighbours required for this grid
        READ( nunit, *, err = 9999, end=99999 ) NNBR

      else  ! old format
! - read max number of nodes in this grid
        read(firstline,*,err=9999, end=99999) np
! - read max number of neighbours required for this grid
        READ( nunit, *, err = 9999, end=99999 ) NNBR
! - read offsets, scale factors, coordinate type
        READ( nunit, *, err = 9999, end=99999 ) x0off, y0off, scaleX, scaleY
        igridtype = 0
      endif  ! if new or old NGH format

!  allocate arrays
!  Main arrays for grid
      ne = 2*np
      ALLOCATE (x(np),y(np),z(np),ncode(np),NL(NNBR,np),nen(4,ne),Ecode(ne), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate main storage arrays'
        call exit(71)
      endif
   
      do i=1,np  
        READ(nunit,*,end=999,err=9999) ii,X(i),Y(i),NCODE(i),z(i),( NL(j,i),j = 1, NNBR )
      enddo
      
      close(nunit)
!      if(.not.quit) call DoCheckEdges()

!  OK, now read element list
      nunit = 18
      write(*,*) 'Enter Element File name (*.el*)'
      read(*,'(a)') FNAME

      fnlen = len_trim( fname )
      OPEN(UNIT=nunit,file=FNAME(1:fnlen),status='old')
      
      do j=1,ne
        read(nunit,*) (nen(i,j),i=1,4),Ecode(j)
      enddo
      close(nunit)

      RETURN

      999   continue
      write(*,*)'WARNING: Premature end of file.'
      write(*,'(a,i6,a,i6,a)') 'Only ',ii,' of expected ',np,' nodes read in'
      Quit = .TRUE.
      return
   
9999  continue  
      write(*,*) 'ERROR reading grid file: Most likely a format error'
      Quit = .TRUE.
      return
      
99999 continue  
      write(*,*) 'ERROR reading grid file: Premature end of file in header.'
      Quit = .TRUE.
      return
      END

! --------------------------------------------------------------------------*

      SUBROUTINE ReadXYEData (nunit, quit)
   
! Purpose : To read grid data into memory

      use MainArrays

!     - PASSED VARIABLES
      LOGICAL Quit

!     - LOCAL VARIABLES
      INTEGER i , j, nunit
      integer TotBndys,TotIntBndys,PtsThisBnd(100)
!        - counters
      INTEGER   istat

!------------------BEGIN-------------------------------

      QUIT = .FALSE.
      
! - read max number of nodes and elements in this grid
      read(nunit,*,IOSTAT=istat) np, ne
      if(istat.lt.0) then
        write(*,*) 'ERROR end of file at line 2','ReadGrid'
        Quit = .TRUE.
        return
      elseif(istat.gt.0) then
        write(*,*) 'ERROR reading XYE file at line 2','ReadGrid'
        Quit = .TRUE.
        return
      endif

!  allocate arrays
!  Main arrays for grid
      ALLOCATE (x(np),y(np),z(np),ncode(np),nen(4,ne),Ecode(ne), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate main storage arrays'
        call exit(71)
      endif

      do i=1,np
        read(nunit,*,IOSTAT=istat) x(i),y(i),z(i)
        if(istat.lt.0) then
          write(*,*) 'ERROR premature end of file','ReadGrid'
          Quit = .TRUE.
          return
        elseif(istat.gt.0) then
          write(*,*) 'ERROR reading XYE file: Most likely a format error','ReadGrid'
          Quit = .TRUE.
          return
        endif 
      enddo

      ncode = 0
      nen = 0
      ECode = 1

      do i=1,ne
        read(nunit,*,IOSTAT=istat) (nen(j,i),j=1,3)
        if(istat.lt.0) then
          write(*,*) 'ERROR premature end of file','ReadGrid'
          Quit = .TRUE.
          return
        elseif(istat.gt.0) then
          write(*,*) 'ERROR reading XYE file: Most likely a format error','ReadGrid'
          Quit = .TRUE.
          return
        endif
      enddo

! *** generate neighbor list
      nindx = np
!      CALL ALTER (nindx,np,ncode,nnbr,nnbr,NL,ne,nen,ECode,&
!                     TotBndys,TotIntBndys,PtsThisBnd)

      return
      END

! --------------------------------------------------------------------------*

      SUBROUTINE DoCheckEdges()

!  Purpose:
!  1.  To check that all edges defined by and node's neighbour list are
!      also defined by the neighbouring node's neighbour list. This check
!      is necessary because the definition of each edge is defined in two
!      places: once in each nodes neighbour list. This means it is easy
!      for one of the lists to be modified and the other not modified.
!      For any edges defined only at one end, we add the other end's 
!      definition.
!  2.  To check that the neighbours in any node's neighbour list are
!      not listed more than once. Duplicates are deleted.
!  3.  To check that the neighbours exist. If they do not exist, they are
!      deleted from the neighbour list.
!  4.  To ensure that a node is not connected to itself. If it is, this
!      connection is deleted.
!  Method:
!      For each node in turn, check 2 is implemented first, and duplicates
!      removed without confirming with the user, but with information.
!      Then for each neighbour, the neighbours neighbour list is checked
!      to ensure that it contains this node. If not, it is added, subject
!      to user confirmation.

! dummy arguments
!      integer Nrec
! include files

      use MainArrays

! local variables
      integer i, j, k, nei, last
      character*80 cstr
! code
      do i=1,np
!          if(exist(i).or.code(i).ge.0)then
          if(ncode(i).ge.0)then
! test 3: remove non-existent neighbours
              do j=1,nnbr
                  if(NL(j,i).gt.0)then
!                      if(.not.exist(NL(j,i)).or.code(i).lt.0) then
                      if(ncode(i).lt.0) then
!                         pointing to a non-existent node
!                          NL(j,i) = 0
                          write(*,'(a,i7,a,i7)')'Removing connection from node ',i,' to non-existent (deleted) node',NL(j,i)
                          NL(j,i) = 0
                      else if(NL(j,i).gt.np) then
!                         pointing to a node > last node
                          write(*,'(a,i7,a,i7)')'Removing connection from node ',i,' to non-existent ( >Nrec) node',NL(j,i)
                          NL(j,i) = 0
                      else if(NL(j,i).eq.i) then
! test 4: cannot be joined to self
                          write(*,'(a,i7,a)')'Removing connection from node ',i,' to itself!'
                          NL(j,i) = 0
                      endif
                  endif
              enddo
! test 2: compress neighbour list
!             step 1: replace duplicates with zeros
              do j=1,nnbr
                  if(NL(j,i).gt.0) then
                    do k=j+1,nnbr
                      if(NL(j,i).eq.NL(k,i)) then
!                         duplicate found
                          write(*,'(a,i6,a,i6,a,i6)')'Removing duplicated neighbour ',k,'(',NL(k,i),') from node',i
                          NL(k,i) = 0
                      endif
                    enddo
                  endif
              enddo
!             step 2: remove zero neighbours, move others to left in list
              last = nnbr
              do while  ((NL(last,i).eq.0).and.(last.gt.1))
                  last = last - 1
              enddo
              do j=1,last
                  if(NL(j,i).eq.0) then
                      do k=j,last
                          NL(k,i) = NL(k+1,i)
                      enddo
                  endif
              enddo
! test 1: check neighbour list of all neighbours includes this node
              do j=1,nnbr
                  nei = NL(j,i)
                  if(nei.ne.0) then
                      k=1
                      do while(    (NL(k,nei).ne.i).and.(k.lt.nnbr))
                          k = k+1
                      enddo
                      if(NL(k,nei).ne.i) then
!                         this node not found in neighbour's neighbour list
                          write (*,'(a,i7,a,i7)') 'Adding node ',i,' to neighbour list of node ',nei
                          k = 1
                          do while(    (NL(k,nei).ne.0).and.(k.lt.nnbr))
                              k = k + 1
                          enddo
                          if(k.le.nnbr)then
                              NL(k,nei) = i
                          else
                              write(*,*)'Too many neighbours already to add new one.'
                          endif
                      endif
                  endif
              enddo
          endif
      enddo

      end

! *********************************************************************

    SUBROUTINE ALTER(nindx,np,ncode,maxngh,numngh,nbrs,ne,nen,IECode,&
                     numbnd,NumIntBnd,numbndpts)

! ***********************************************************************
! This routine converts from triangle list (TRIANG 
! format) and node file (NODE format) to NEIGH format
! ******************************************************

    IMPLICIT NONE

! *** passed VARIABLES ***
      integer nindx
      integer np,ncode(np),maxngh,numngh,nbrs(maxngh,np)
      integer ne,nen(4,ne),IECode(ne)
      integer numbnd,NumIntBnd,NumBndPts(NumBnd)

! *** LOCAL VARIABLES ***

      integer nonbrs(nindx)
!   CUVX - current vertex being considered
      INTEGER CUVX
!   CUNBR - current neighbour being considered
    INTEGER CUNBR
    INTEGER II,JJ,KK,LL,MM
    LOGICAL pass1, NEWNBR
    character cstr*80

!   Starting neighbor list

!   Set count of nbrs to zero for all nodes

    DO 10 KK = 1, np
!      Exist(KK) = .TRUE.
        NONBRS(KK) = 0
!   and set all neighbour arrays to zero
          DO 15 JJ = 1, maxngh
          nbrs(JJ,KK) = 0
15        CONTINUE
10      CONTINUE


! Nodes in output NEIGH file are in same order as in input NODE file

! *** Check each triangle and check that each vertex is in each of 
! *** the other two vertices' neighbour lists

        pass1 = .true.
        DO 101 JJ = 1, ne
            IECode(JJ) = 1 
! *** Check each vertex in current triangle
          DO 102 II = 1, 3
! *** Choose current vertex
            CUVX = nen(II,JJ)
! *** Take other two vertices in turn
            DO 103 LL = 1,3
              IF(LL.eq.II) GO TO 98
!              Choose current neighbour of chosen vertex 
              CUNBR = nen(LL,JJ)
!              Check if CUNBR is already in neighbour list of CUVX
              NEWNBR = .TRUE.
              IF(NONBRS(CUVX).eq.0) GO TO 99
              DO 104 MM = 1, NONBRS(CUVX)
                IF(CUNBR.eq.nbrs(MM,CUVX)) NEWNBR = .FALSE.
104           CONTINUE
99            CONTINUE
!              If CUNBR is new neighbour of CUVX, add to list
              IF(NEWNBR) THEN     
                if(nonbrs(cuvx).ge.maxngh) then
                  if(pass1) then
                    write(*,*)' Too many neighbor points - truncating:'
!                   stop
                    pass1 = .false.
                  endif
                else
                  NONBRS(CUVX) = NONBRS(CUVX) + 1
                  nbrs(NONBRS(CUVX),CUVX) = CUNBR
                endif
              ENDIF
98            CONTINUE
103         CONTINUE
102       CONTINUE
101     CONTINUE

!       Find max number of neighbours
        numngh = 0
        DO 105 II = 1, np
!          if(NONBRS(II).le.0) ncode(II) = -9
          IF(NONBRS(II).gt.numngh) numngh = NONBRS(II)
105     CONTINUE

!       Set up computational codes : 1 everywhere on outer boundary
!                                    2 on island boundaries
!                                    90 at line constraints
!                                    0 at interior nodes
        CUVX = 0
        DO 108 JJ = 1, NumBnd
          DO 109 KK = 1, NumBndPts(JJ)
             CUVX = CUVX + 1
             IF (JJ.eq.1) THEN
                nCODE(CUVX) = 1
             ELSEIF (jj.le.NumBnd-NumIntBnd) then
                nCODE(CUVX) = 2
             ELSE
                nCODE(CUVX) = 90
             ENDIF
109       CONTINUE
108     CONTINUE

        do jj=cuvx+1,np
          ncode(jj) = 0
        enddo

    RETURN
    END

! *****************************************************************
! -------------------------------------------------------------------
