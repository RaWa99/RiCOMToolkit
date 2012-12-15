      PROGRAM REFINEx2
C------------------------------------------------------------------------
C   Purpose :
c	This program takes the neigh [NEIGH] and triangle [TRIANG] files
c	and refines by 2 and creates a file in Trigrid NEIGH format. 
c
c	Two files are required on input:
c		device 18 contains the neigh file
c		device 11 contains the triangle file
c	The output file is device 2.
C-----------------------------------------------------------------------
C
c      include '..\includes\master2.par'

      parameter (maxnp=500000, maxne=2*maxnp, nbtot=50 )

      common /stuff/ NEN(maxne,6),X(maxnp),y(maxnp),dep(maxnp),
     *        NTMP(maxne,3),ic(maxnp),ngh(maxnp,50),nongh(maxnp)
      character*80 fname
      CHARACTER ANS*1
      LOGICAL EXIS
	logical ResOK, openfile
      data maxngh2/50/
      DATA VOID/-1.E30/

!***********************************************************************
! Program info displayed to user at run time
      PRINT *,'-------------------------------------------------------'
      PRINT *,'Program REFINEx2 reads a NEIGH and a TRIANG format'
      PRINT *,'file, refines the grid by subdividing the triangles'
      PRINT *,'at the midsides, and greates a new NEIGH format file'
      PRINT *,' '
      PRINT *,'Array size settings for input and output grids: '
      PRINT *,'  Max. number of nodes               = ',maxnp  
      PRINT *,'  Max. number of neighbours per node = ',NBTOT
      PRINT *,'  Max. number of elements            = ',maxne  
      PRINT *,'------------------------------------------------------'
      PRINT *,' '
! -----------------------------------------------------------------------

      izero = 0
      ione = 1
      DO J=1,maxne
        do k=1,3
          ntmp(j,k)=0
        enddo
        DO K=1,6
          NEN(J,K)=0
	  enddo
	enddo

      DO J=1,maxnp
        X(J)=VOID
        y(J)=VOID
        dep(J)=VOID
        ic(j) = 0
        nongh(j) = 0
        do k=1,50
          ngh(j,k) = 0
        enddo
	enddo

! *** read input data
!340   continue
!      write(*,*) ' enter filename for neigh input file'
!      read(*,'(a80)') fname
!      INQUIRE (FILE=FName, EXIST=Exis)
!      IF (Exis) THEN
!      else!
!	PRINT *,'FILE DOES NOT EXIST - QUIT(Y/N)?'
!	READ(5,'(A)') ANS
!	IF (ANS.EQ.'Y'.or.ans.eq.'y') then
!          stop
!        else
!	  GOTO 340
!        endif
!      ENDIF
!      OPEN (UNIT=18,file=fname,STATUS='old')

      ResOK = OpenFile(18,'Open NGH Input File',fname,                  &
     &    "Input file(*.ngh),*.ngh;All files(*.*),*.*;")

      if(.not.resOK) then
	  stop
	endif

	read(18,*) np
      if(np.gt.maxnp) then
        write(*,*) ' np too large..increase dimensions'
        stop
      endif
	read(18,*) nlim
	nneigh=nlim+1
	read(18,*) xmax,ymax,xmin,ymin
	do k=1,np
	  no=nongh(k)
	  read(18,*) kk,x(k),y(k),ic(k),dep(k),(ngh(k,l),l=2,nneigh)
	enddo
      close (18)

!360   continue
!      write(*,*) ' enter filename for triang input file'
!      read(*,'(a80)') fname
!      INQUIRE (FILE=FName, EXIST=Exis)
!      IF (Exis) THEN
!      else
!	PRINT *,'FILE DOES NOT EXIST - QUIT(Y/N)?'
!	READ(5,'(A)') ANS
!	IF (ANS.EQ.'Y'.or.ans.eq.'y') then
!          stop
!        else
!          GOTO 360
!        endif
!      ENDIF
!      open(unit=11,file=FNAME,status='old')

      ResOK = OpenFile(11,'Open Element Input File',fname,              &
     &    "Input file(*.ele),*.el*;All files(*.*),*.*;")

      if(.not.resOK) then
	  stop
	endif

      ii = 1
      do j=1,maxne
! *** IOS
!       read(11,*,end=207) jj,(ntmp(j,k),k=1,3)
! *** TRI format
        read(11,*,end=207) (ntmp(j,k),k=1,3)
        ii = ii + 1
      enddo
207   ne = ii-1
      CLOSE(UNIT=11)
      write(*,*) ' ne,np=',ne,np

      do j=1,ne
        do k=1,3
          kk=2*k-1
          nen(j,kk)=ntmp(j,k)
        enddo
      enddo

      NCN=6
      NPN=NP

      DO 350 L=1,NE
        DO 300 M=1,NCN,2
          N1=IABS(NEN(L,M))
          N2=IABS(NEN(L,M+1))
          N3=M+2
          IF(N3.GT.NCN) N3=1
          N3=IABS(NEN(L,N3))
          IF(N2.NE.0) GO TO 300
          NPN=NPN+1
          NEN(L,M+1)=NPN
      if(ic(n1).le.0.and.ic(n3).le.0.) then
	  
      elseif(ic(n1).eq.1.and.ic(n3).eq.1) then
        ic(npn) = 1
      elseif(ic(n1).eq.2.and.ic(n3).eq.2) then
        ic(npn) = 2
      elseif((ic(n1).eq.5.or.ic(n1).eq.6)
     *  .and.(ic(n3).eq.5.or.ic(n3).eq.6)) then
        ic(npn) = 5
      elseif((ic(n1).eq.5.or.ic(n1).eq.6)
     *  .and.(ic(n3).eq.10.or.ic(n3).eq.11)) then
        ic(npn) = 11
      elseif((ic(n1).eq.10.or.ic(n1).eq.11)
     *  .and.(ic(n3).eq.5.or.ic(n3).eq.6)) then
        ic(npn) = 11
      elseif((ic(n1).eq.5.or.ic(n1).eq.6)
     *  .and.(ic(n3).eq.1)) then
        ic(npn) = 1
      elseif((ic(n1).eq.5.or.ic(n1).eq.6)
     *  .and.(ic(n3).eq.2)) then
        ic(npn) = 2
      elseif((ic(n1).eq.1)
     *  .and.(ic(n3).eq.5.or.ic(n3).eq.6)) then
        ic(npn) = 1
      elseif((ic(n1).eq.2)
     *  .and.(ic(n3).eq.5.or.ic(n3).eq.6)) then
        ic(npn) = 2
      elseif((ic(n1).eq.7.or.ic(n1).eq.8)
     *  .and.(ic(n3).eq.7.or.ic(n3).eq.8)) then
        ic(npn) = 7
      elseif((ic(n1).eq.10.or.ic(n1).eq.11)
     *  .and.(ic(n3).eq.1)) then
        ic(npn) = 1
      elseif((ic(n1).eq.10.or.ic(n1).eq.11)
     *  .and.(ic(n3).eq.2)) then
        ic(npn) = 2
      elseif((ic(n1).eq.1)
     *  .and.(ic(n3).eq.10.or.ic(n3).eq.11)) then
        ic(npn) = 1
      elseif((ic(n1).eq.2)
     *  .and.(ic(n3).eq.10.or.ic(n3).eq.11)) then
        ic(npn) = 2
      elseif((ic(n1).eq.10.or.ic(n1).eq.11)
     *  .and.(ic(n3).eq.10.or.ic(n3).eq.11)) then
        ic(npn) = 11
      elseif(ic(n1).eq.90.and.ic(n3).eq.90) then
        ic(npn) = 90
      else
        ic(npn) = 0
      endif
      X(NPN)=(X(N1)+X(N3))*.5
      Y(NPN)=(Y(N1)+Y(N3))*.5
      dep(NPN)=(dep(N1)+dep(N3))*.5
!-.....LOOK FOR THE SAME SIDE
210   ICNT=0
      DO 250 J=1,NE
      IF(J.EQ.L) GO TO 250
        DO 200 K=1,NCN,2
          NN1=IABS(NEN(J,K))
          IF(NN1.NE.N3) GO TO 200
          NN3=K+2
          IF(NN3.GE.NCN) NN3=1
          NN3=IABS(NEN(J,NN3))
          IF(NN3.NE.N1) GO TO 200
            NN2=IABS(NEN(J,K+1))
            IF(NN2.NE.0) WRITE(6,6090) J,K,NN1,NN2,NN3,L,M,N1,NPN,N3
            NEN(J,K+1)=NPN
            ICNT=ICNT+1
            IF(ICNT.GT.1) WRITE(6,6095) ICNT,L,M,J,K
            if(ic(npn).ne.90) ic(npn)=0
          go to 300
200     CONTINUE
250   CONTINUE

300     CONTINUE
350   CONTINUE

      DO 100 J=1,NE
        DO 50 K=1,4
          N1= MOD(2*K+3,6) +1
          N2=2*K-1
          N3=N2+1
          IF(K.EQ.4) N2=2
          IF(K.EQ.4) N3=4
          IE=4*(J-1)+K
          NTMP(IE,1)=IABS(nen(J,N1))
          NTMP(IE,2)=IABS(nen(J,N2))
          NTMP(IE,3)=IABS(nen(J,N3))
50      CONTINUE
100   CONTINUE

760   NE=4*NE
      np = npn

      write(*,*) ' After refine, ne,np=',ne,np

! *** write output triangle data
!840   continue
!      write(*,*) ' enter filename for triangle output file'
!      read(*,'(a80)') fname
!      INQUIRE (FILE=FName, EXIST=Exis)
!      IF (Exis) THEN
!	PRINT *,'FILE ALREADY EXISTS - RENAME(Y/N)?'
!	READ(5,'(A)') ANS
!	IF (ANS.EQ.'Y'.or.ans.eq.'y') then
!	  GOTO 840
 !       endif
!      ENDIF
!      OPEN (UNIT=18,file=fname,STATUS='new')

      ResOK = OpenFile(18,'Open Element OUTput File',fname,             &
     &    "OUTput file(*.ele),*.el*;All files(*.*),*.*;")

      if(.not.resOK) then
	  stop
	endif
      do j=1,ne
        write(18,1855) (ntmp(j,k),k=1,3),izero, ione
      enddo

1855  format(4(1x,i7),i4)
      close (18)

!	for each node check to see what triangles it is in
!	and collect the neighbours
        nlim=0
	do 80 k=1,np
	  ngh(k,1)=k
	  nongh(k)=1
	  do 81 l=1,ne
	    do 82 ll=1,3
	      if(ntmp(l,ll).eq.k) go to 83
82	    continue
	    go to 81
83	    continue
!	found a element that contains node k
!	if the new nodes are not already in the list, add them
	    no=nongh(k)
	    do 84 ll=1,3
	      do 851 kk=1,no
	        if(ntmp(l,ll).eq.ngh(k,kk)) go to 84
851	      continue
!	new neighbour
	      no=no+1
            if(no.le.maxngh2) then
              ngh(k,no)=ntmp(l,ll)
              nongh(k)=no
            endif
	      if(no.gt.nlim+1) nlim=no-1
84	    continue
81	  continue
80	continue
      write(*,*) ' max no of ngh=',nlim
      if(nlim+1.gt.maxngh2) then
        write(*,*)
     +        ' *** ERROR -too many neighbors,increase dimensions'
        write(*,*) ' maxngh=',maxngh2,' n0=',nlim+1
      endif
! *** remove dead nodes
      write(*,*) ' *** removing dead nodes'
      nend=np
      j=0
911   j=j+1
914   indx=0
      if(ngh(j,2).le.0) then
! *** roll down arrays
         write(*,*) ' *** rolling ',j
         indx=1
         nend=nend-1
         do 921 k=j,nend
           ngh(k,1)=k
           x(k)=x(k+1)
           y(k)=y(k+1)
           ic(k)=ic(k+1)
           dep(k)=dep(k+1)
           do 917 L=2,16
              ngh(k,L)=ngh(k+1,L)
917        continue
921      continue
         do 925 k=1,nend
         do 925 L=2,16
           if(ngh(k,L).gt.j) ngh(k,L)=ngh(k,L)-1
925      continue
      endif
      if(j.lt.nend) then
         if(indx.eq.1) then
           go to 914
         else
           go to 911
         endif
      endif
      np = nend

! *** write output data
!      write(*,*) ' enter filename for neigh output file'
!      read(*,'(a80)') fname
!      OPEN (UNIT=18,file=fname,STATUS='new')

      ResOK = OpenFile(18,'Open NGH OUTput File',fname,                 &
     &    "OUTput file(*.ngh),*.ngh;All files(*.*),*.*;")

      if(.not.resOK) then
	  stop
	endif

	write(18,781) np
	nneigh=nlim+1
	write(18,781) nlim
781	format(i10)
	write(18,782) xmax,ymax,xmin,ymin
782	format(4(1x,e15.7))
	write(*,*) ' wait- writing out the new file'
      if(np.lt.10000) then
	  do 92 k=1,np
          if(ic(k).eq.11) ic(k)=10
	    no=nongh(k)
	    write(18,91) k,x(k),y(k),ic(k),dep(k),(ngh(k,l),l=2,nneigh)
92	  continue
91	  format(i8,2(1x,e15.7),1x,i4,1x,e15.7,15i8)
      elseif(np.lt.100000) then
	  do 892 k=1,np
          if(ic(k).eq.11) ic(k)=10
	    no=nongh(k)
	    write(18,91) k,x(k),y(k),ic(k),dep(k),(ngh(k,l),l=2,nneigh)
892	  continue
891	  format(i6,2(1x,f12.5),1x,i2,1x,f9.2,15(1x,i5))
      else
	  do k=1,np
          if(ic(k).eq.11) ic(k)=10
	    no=nongh(k)
	    write(18,91) k,x(k),y(k),ic(k),dep(k),(ngh(k,l),l=2,nneigh)
	  enddo
895	  format(i7,2(1x,f12.5),1x,i2,1x,f9.2,15(1x,i7))
      endif
	close(18)

      STOP

!-.....INPUT DATA CARD FORMATS.....
5090  FORMAT(A)
6090  FORMAT(5X,'SIDE ERROR- ',10I6)
6095  FORMAT(5X,I3,' DUPLICATE SIDES, L,M=',2I5,'  J,K=',2I5)

      END
