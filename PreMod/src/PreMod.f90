!*****************************************************************

      Module PMMod

      SAVE

! Scalar data
      integer ne,np,npr,ntype,ncnmax,ncn,nsides,limngh,lp
      integer izup,itest,ireorder,isideopt,iprt,nopt
      integer ifront,maxrow,maxsto

! 1D data
      integer, allocatable ::  nbc(:), matl(:)
      real, allocatable ::  Area(:)
      real, allocatable ::  sdep(:),slen(:),refdep(:),sdx(:),sdy(:)
      real, allocatable ::  dlinv(:),gamma(:)
      real, allocatable ::  SPEC(:),alfa(:)

! 2D data
      integer, allocatable ::  nbrs(:,:), ies(:,:)
      integer, allocatable ::  nen(:,:), ieadj(:,:)
      integer, allocatable ::  numside(:,:), iside(:,:), iends(:,:), isidebc(:,:)
      real, allocatable ::  sxy(:,:), xyz(:,:)

      character(256) OutputFileName
      
      real*8 xnp,xnn,dnan,dnbn,weight
      integer ngp,ncn3,jordr
!      integer nbrs(maxngh,maxnp)
      common /ebasis/ xnn(6,7),xnp(6,3),jordr,ngp,ncn3

!  Following needed for netCDF
      integer	::	nxyz=3		! 1st dim of xyz 
      integer	::	n_numside=4	! 2nd dim of numside 
      integer	::	nieadj=5	! 1st dim of ieadj 
      integer	::	ntwo=2		! 1st dim of sxy & iside 


      END MODULE

!*****************************************************************

      PROGRAM PreModTQ

      USE PMMod

! *** read network input data, perform tests, and write binary output file
! *** Modified for p methods. Expand bases to order = jordr

!	To make this program work properly with Visual Fortran, go to Settings configuration and do this:
!	1.  In Fortran/Compatibility click "Filenames from Command Line"
!      This will allow file browsing. Otherwise use command line option in subroutine OpenFileCMD.f90

      real*8 etime, ticktock

! *** START

      write(*,*) ' PreModTQ Version 03.4.15'
      write(*,*) ' A preprocessor for RAW models, triangles and quads'
      write(*,*) ' Created by Roy Walters, r.walters@niwa.cri.nz'
      write(*,*) ' Program size: Allocated arrays'
      write(*,*) ' '

      call DataInput
      
      if(ifront.gt.0) then
        write(*,*) ' Calling PreFront'
        write(lp,*) ' Calling PreFront'
        CALL PreFront
      endif

      if(NOPT.eq.1) then
        call OutputBinary
      elseif(NOPT.eq.2) then
        call Output_netCDF
      else
        write(*,*)'No output file'
      endif

      stop
      end

!*****************************************************************

      Subroutine DataInput

      USE PMMod

! *** read network input data, perform tests, and write binary output file
! *** Modified for p methods. Expand bases to order = jordr

!
!	Input deck looks like this:
!	xscale,yscale,zscale are the scale factors to apply to x, y and z to get the correct units (m).
!	zoffset	= the amount by which to displace z so z=z+zoffset
!	zlimit	= the minimum depth if izup is 0 (z is positive downward), all depths will be greater than this
!			= the maximum elevation if izup is 1 (z is positive upward).
!	izup	= 1 => z is positive upwards (RiCom)
!			= 0 => z is positive downwards, ie depth (TIDE2D, WEQ2D, etc)
!	ifulltest	not implemented
!	ireorder	not implemented
!	iord		element basis polynomial order (valid for isideopt=0, not for RiCOM)
!	ifmt	=0	=> triangle list, 3 nodes per line
!			=1	=> element list + material type (4 nodes + element code per line)
!			=2	=> element no, element list (4 integers per line)
!			=3	=> element list (4 integers per line) from Tecplot for mixed tri/quad
!			=4	=> xyz      and element list + material type (4 integers per line from Tecplot for mixed tri/quad + 1 integer)
!			=5	=> xyz+code and element list + material type (4 integers per line from Tecplot for mixed tri/quad + 1 integer)
!	ifront	= 0 => no frontal solver information generated
!			= 1 => frontal solver information generated (for TIDE2D)
!	isideopt = 0	=> do not generate side info 
!			 = 1	=> generate side info (for RiCOM)
!
!	Following lines are:
!	filename of .ngh file
!	filename of .tri file
!	ialfa = no of nodes where boundary angle is specified
!	nbx = no of nodes where boundary node codes are specified
!	iptr=0 => don't write text file of results
!	nopt	= 0 => ascii output
!			= 1 => binary output
!			= 2 => netCDF output
!	filename of output file

      integer istat, igridtype
      real x0off, y0off, scaleX, scaleY
      real*8 etime, ticktock
      character*80 title, ans*1
      character*256 FNAME
      character(80) Firstline
      logical exis
      logical(4) resOK, openfile
      DATA VOID/-1.E20/

! *** START


      NCNmax=4
!  set input and output unit numbers
      in = 9
      lp = 10

      ResOK = OpenFile(in,'Open Control File',fname, &
     &  "Control file(*.pmd),*.pmd;All files(*.*),*.*;")

      if(.not.resOK) then
        stop
      endif

      ResOK = OpenFile(lp,'Open Output File',fname, &
     &  "Output file(*.pre),*.pre;All files(*.*),*.*;")

      if(.not.resOK) then
        stop
      endif

      etime = ticktock()

!-.....READ AND PRINT CONTROL PARAMETERS
      READ(in,'(a)') TITLE
      WRITE(lp,*)  TITLE
      READ(in,*) XSCALE,YSCALE,ZSCALE,zoffset,zlimit,izup
      read(in,*) itest,ireorder,ifront,iord,ifmt,isideopt

!-.....READ NODAL COORDINATES

      write(*,*) ' Reading ngh file....'
      read(in,'(a)') FNAME
      if(ifmt.lt.10) then
        OPEN(UNIT=1,file=FNAME,status='OLD')
        READ(1,'(a)') Firstline
        if(firstline(1:4).eq."#NGH") then  !neigh grid file
          do
            READ(1,'(a)') Firstline
            if(firstline(1:1).ne."#") then    !comment lines
!           following line is internal read of firstline          
! - read offsets, scale factors, coordinate type
              READ( firstline, * ) x0off, y0off, scaleX, scaleY, igridtype
              exit
            endif
          enddo
          read(1,*) np
          read(1,*) limngh
        else
          read(firstline,*) np
          read(1,*) limngh
          read(1,*) x0off, y0off, scaleX, scaleY
        endif
      else
        OPEN(UNIT=1,file=FNAME,form='unformatted',status='OLD')
        read(1) np
        read(1) limngh
        read(1) x0off, y0off, scaleX, scaleY
      endif

!  Main arrays for grid
      write(*,*) '   Allocating main arrays'
      ne = 2*np
      ALLOCATE (nen(ncnmax,ne),matl(ne),xyz(3,np),area(ne),nbc(np),alfa(np), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate main storage arrays'
        stop
      endif

      nsides = 0
      alfa = 0.
      xyz = 0.
      nbc = 0
      nen = 0
      matl = 1
      area = 0.

      ALLOCATE (nbrs(limngh,np), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate neighbor array'
        stop
      endif
      nbrs = 0

      if(ifmt.lt.10) then

        if(ifmt.eq.3.or.ifmt.eq.4) then
          do j=1,np
            READ(1,*,end=89,err=989) (XYZ(k,J),K=1,3)
          enddo
        elseif(ifmt.eq.5) then
          do j=1,np
            READ(1,*,end=89,err=989) (XYZ(k,J),K=1,3), nbc(j)
          enddo
        else
          do j=1,np
            READ(1,*,end=89,err=989) jj,(XYZ(k,J),K=1,2),nbc(j),xyz(3,j),(nbrs(l,j),l=1,limngh)
          enddo
        endif
      else
        read(1)
        if(ifmt.eq.13.or.ifmt.eq.14) then
          do j=1,np
            READ(1,end=89,err=989) (XYZ(k,J),K=1,3)
          enddo
        elseif(ifmt.eq.15) then
          do j=1,np
            READ(1,end=89,err=989) (XYZ(k,J),K=1,3), nbc(j)
          enddo
        else
          do j=1,np
            READ(1,end=89,err=989) jj,(XYZ(k,J),K=1,2),nbc(j),xyz(3,j),(nbrs(l,j),l=1,limngh)
          enddo
        endif
      endif
89    CLOSE(unit=1)

      DO J = 1,NP
        XYZ(1,J) = XYZ(1,J)*XSCALE
        XYZ(2,J) = XYZ(2,J)*YSCALE
        XYZ(3,J) = XYZ(3,J)*ZSCALE +zoffset
          if(izup.gt.0.) then
            if(XYZ(3,J).gt.zlimit) XYZ(3,J) = zlimit
          else
            if(XYZ(3,J).lt.zlimit) XYZ(3,J) = zlimit
          endif
      enddo

!-.....FIND MAX NODE AND ELEMENT NUMBERS.....
      read(in,'(a)') FNAME
      OPEN(UNIT=1,file=FNAME,status='OLD')

      write(*,*) ' Reading element file....'
      ncn = 3
      j=0
67    j=j+1
      ifmt = mod(ifmt,10)
      if(ifmt.eq.0) then
        READ(1,*,end=69,err=969) (NEN(k,J),K=1,3)
      elseif(ifmt.eq.1) then
        READ(1,*,end=69,err=969) (NEN(k,J),K=1,4),matl(j)
        if(NEN(4,J).ne.0) ncn=4
      elseif(ifmt.eq.2) then
        READ(1,*,end=69,err=969) JJ,(NEN(k,J),K=1,3)
      elseif(ifmt.eq.3) then  !tecplot write file
        READ(1,*,end=69,err=969) (NEN(k,J),K=1,4)
        if(NEN(4,J).eq.NEN(3,J)) NEN(4,J)=0
        if(NEN(4,J).ne.0) ncn=4
       elseif(ifmt.eq.4.or.ifmt.eq.5) then  !tecplot write file +  matl
        READ(1,*,end=69,err=969) (NEN(k,J),K=1,4),matl(j)
        if(NEN(4,J).eq.NEN(3,J)) NEN(4,J)=0
        if(NEN(4,J).ne.0) ncn=4
     else
        write(*,*) ' Unknown triangle format'
        write(lp,*) ' Unknown triangle format'
        stop
      endif
      go to 67

69    CLOSE(unit=1)
      ne=j-1
      nump = 0
      ntype = 1
      DO J =1, ne
        if(matl(j).gt.ntype) ntype=matl(j)
        DO K = 1, NCN
          IF( NEN(k,J) .GT. nump ) nump = NEN(k,J)
        enddo
      enddo

      if(np.ne.nump) then
        write(lp,*) np, nump,' np in NEI file .ne. nump in ELE file'
      endif

      write(*, *) ' NE=',ne,', NP=',np,', ntype=',ntype 
      write(lp, *) ' NE=',ne,', NP=',np,', ntype=',ntype

      if(ireorder.gt.0) then
        write(lp,*) ' Reordering not implemented yet'
      endif

      write(*,*) ' Nodal and element data complete'
      write(lp,*) ' Nodal and element data complete'

! ** translate TRIGRID codes to model boundary codes
! ** code=1, land boundary, slip condition
! ** code=2, island boundary, slip condition
! ** code=3, non-zero normal velocity, no translation possible
! ** code=4, geostrophic condition, no translation possible
! ** code=5, open boundary, sea level condition
! ** code=6, join open and land boundary, slip condition
! ** code=7, open boundary, radiation condition
! ** code=8, join 1 and 7
! ** code=9, flux condition
! ** code=10, land boundary, u=v=0
! ** code=90, interior node, fixed location
! ** code=97, moving bottom
      DO J=1,NP
        if(nbc(j).eq.1) then
           nbc(j)=10
        elseif(nbc(j).eq.2) then
           nbc(j)=10
        elseif(nbc(j).eq.5) then
           nbc(j)= 2
        elseif(nbc(j).eq.6) then
           nbc(j)=12
        elseif(nbc(j).eq.7) then
           nbc(j)=7
        elseif(nbc(j).eq.8) then
           nbc(j)=17
        elseif(nbc(j).eq.9) then
           nbc(j)=9
        elseif(nbc(j).eq.10) then
           nbc(j)=19
        elseif(nbc(j).eq.11) then
           nbc(j)=110
        elseif(nbc(j).eq.97) then
          !leave alone
        else
           nbc(j)=0
        endif
      enddo

      if(itest.eq.0) then
        write(*,*) ' NO testing, do NOT generate boundary slopes'
        write(lp,*) ' NO testing, do NOT generate boundary slopes'
      elseif(itest.eq.1) then
        write(*,*) ' Quick test and generate boundary slopes'
        write(lp,*) ' Quick test and generate boundary slopes'
        Call QuickTest(nes)
		Call BoundarySlope(nes)
	  elseif(itest.eq.2) then
        write(*,*) ' Full test and generate boundary slopes'
        write(lp,*) ' Full test and generate boundary slopes'
        Call FullTest(nes)
		Call BoundarySlope(nes)
	  else
        write(*,*) ' Test option not valid'	  
        write(lp,*) ' Test option not valid'	  
	  endif

!-.....READ SPECIFIED NODAL SLOPES
      read(in,*) ialfa
      if(ialfa.gt.0) then
        do j=1,ialfa
          READ(in,*) k,ALFA(k)
        enddo
      endif

!-.....READ OTHER BOUNDARY CONDITIONS
  279 read(in,*) nbx
      if(nbx.gt.0) then
        do j=1,nbx
          READ(in,*) k,nbc(k)
        enddo
      endif

      DO 390 J=1,NP
        NFX=NBC(J)/10
        IF(NFX.NE.1) then
          ALFA(J)=0.
        else
          ALFA(J)=ATAN(ALFA(J))
        endif
  390 CONTINUE

      write(*,*) ' Finished with slopes'
      write(lp,*) ' Finished with slopes'

! *** expand bases for degree jordr
      npr = np

      if(isideopt.gt.0) then
        call Setup
      else
! *** compute area
        write(*,*) 'Finding areas'
        DO nn=1,NE
          I1=Nen(1,nn)
          I2=Nen(2,nn)
          I3=Nen(3,nn)
          Area(nn)=0.5*(xyz(1,I1)*(xyz(2,I2)-xyz(2,I3)) &
     &               +xyz(1,I2)*(xyz(2,I3)-xyz(2,I1)) &
     &               +xyz(1,I3)*(xyz(2,I1)-xyz(2,I2)))
          if(ncn.eq.4) then
            I4=Nen(ncnmax,nn)
            if(I4.gt.0) then
              Area(nn)=Area(nn) &
     &              +0.5*(xyz(1,I3)*(xyz(2,I4)-xyz(2,I1)) &
     &                   +xyz(1,I4)*(xyz(2,I1)-xyz(2,I3)) &
     &                   +xyz(1,I1)*(xyz(2,I3)-xyz(2,I4)))
            endif
          endif
        enddo

        if(iord.gt.1) then
! *** load area coordinates at nodes
          iflg = iord
          jordr = iord
          iord1 = 1
          npts = (iord+1)*(iord+2)/2
          call shape(iflg,iord1,nipt)
          do k=1,3
            do j=1,npts
              xnp(j,k) = xnn(k,j)
            enddo
          enddo
! *** expand elements for higher order bases
          call expand_basis(jordr)
        endif
		
      endif

      DO nn=1,NE
        if(area(nn).le.0.) then
          write(*,*) ' WARNING: Area le 0 at', nn
          write(lp,*) ' WARNING: Area le 0 at', nn
        endif
      enddo

      read(in,*) iprt
      IF(IPRT.LE.0)  GO TO 233
      write(*,*) ' Writing printer output'
      write(lp,*) ' Writing printer output'

!-.....PRINT NODAL CHARACTERISTICS
      JNT= (NP+2)/3
      WRITE(lp,6040)
      DO 231 J=1,JNT
        WRITE(lp,6050)(N,NBC(N),(XYZ(k,N),K=1,3),ALFA(N),N=J,NP,JNT)
  231 CONTINUE

!-.....COMPUTE AREA AND PRINT ELEMENT CHARACTERISTICS
      jnt = (ne+2)/3
      WRITE(lp,6060)
      DO J=1,JNT
        WRITE(lp,6070)(N,(NEN(m,N),M=1,4),AREA(N),N=J,NE,JNT)
      enddo
  233 CONTINUE

      CLOSE(unit=2)

      etime = ticktock() - etime
      write(*,*) ' End input routine, elapsed time=', etime
      write(lp,*) ' End input routine, elapsed time=', etime

      read(in,*) nopt
      IF(NOPT.gt.0) then
        read(in,'(a)') OutputFileName
      endif

      return

! *** Error traps
969   CLOSE(unit=1)
      write(lp,*) ' ERROR reading triangle file'
      stop

989   CLOSE(unit=1)
      write(lp,*) ' ERROR reading node file'
      stop

!-.....INPUT DATA CARD FORMATS.....
 6040 FORMAT( // 2X,'NODAL COORDINATES,DEPTHS BELOW MLLW, AND SLOPES'//&
           &3('  NODE  B     X        Y     Z     BNDRY') /&
           &3(' NUMBER C    (M)      (M)   (M)    SLOPE'))
 6050 FORMAT(3(I6,I4,2E9.3,F6.1,F7.3))
 6060 FORMAT( // 2X, 'ELEMENT CHARACTERISTICS' //&
        &1X,3(' ELEM   NODES(CCW)   TYPE AREA(M2) '))
 6062 FORMAT( // 2X, 'ELEMENT CHARACTERISTICS' //&
        &1X,1(' ELEM   NODES(CCW)   TYPE AREA(M2) '))
 6070 FORMAT(1X,3(5I6,(1PE10.2)))
 6072 FORMAT(1X,8I6,(1PE10.2))
 6095 FORMAT(5X,I5,' DUPLICATE SIDES, L,M=',2I5,' J,K=',2I5)
 6201 FORMAT(//2X,'***SLOPE ERROR, ELE ',I5,' NODES',3I5,3F9.4)
 6210 FORMAT( // 2X,'ELEMENT INPUT COMPLETE'/5X,'MAX ELEMENT NUMBER =',&
        &I8, ' .. MAX NODE NUMBER =', I8 )

      END

!*****************************************************************

      Subroutine QuickTest(nes)

      USE PMMod

      
      nes = 0
      maxsides = ne

      ALLOCATE (ies(3,maxsides), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate boundary side arrays'
        stop
      endif
	  
      ies = 0
	  
      DO L=1,NE
        ncn2 = ncnmax - 1  + min(1,NEN(ncnmax,L))
        M300: DO M=1,NCN2
          N1=IABS(NEN(M,L))
          if(nbc(n1).eq.0) cycle  !go to 300
          N3=M+1
          IF(N3.GT.NCN2) N3=1
          N3=IABS(NEN(N3,L))
          if(nbc(n3).eq.0) cycle  !go to 300
!.....LOOK FOR THE SAME SIDE
          ICNT=0
          DO J=1,NE
            IF(J.EQ.L) cycle !GO TO 250
            ncn3 = ncnmax - 1  + min(1,NEN(ncnmax,J))
            DO K=1,NCN3
              NN1=IABS(NEN(K,J))
              IF(NN1.NE.N3) cycle  !GO TO 200
              NN3=K+1
              IF(NN3.GT.NCN3) NN3=1
              NN3=IABS(NEN(NN3,J))
              IF(NN3.NE.N1) cycle  !GO TO 200
!.....FOUND ONE
              cycle M300 !to 300 outer  
            enddo
          enddo
          NES=NES+1
          if(nbc(n1).eq.0) NBC(N1)=10
          if(nbc(n3).eq.0) NBC(N3)=10
          if(nes.le.maxsides) then
            IES(1,NES)=N1
            IES(2,NES)=N3
          endif
        enddo M300
!300       CONTINUE
      enddo

      if(nes.gt.maxsides)write(lp,*)' **WARNING** NES exceeds Max boundary sides'
      write(lp,*) '  NES =',NES,'  maxBSides =',maxsides


      END

!*****************************************************************

      Subroutine FullTest(nes)

      USE PMMod

	  nes = 0
	  maxsides = ne

      ALLOCATE (ies(3,maxsides), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate boundary side arrays'
        stop
      endif
	  
	  ies = 0

	  DO L=1,NE
        ncn2 = ncnmax - 1  + min(1,NEN(ncnmax,L))
        DO M=1,NCN2
          N1=IABS(NEN(M,L))
          if(nbc(n1).eq.0) cycle  !go to 300
          N3=M+1
          IF(N3.GT.NCN2) N3=1
          N3=IABS(NEN(N3,L))
          if(nbc(n3).eq.0) cycle  !go to 300
!.....LOOK FOR THE SAME SIDE
          ICNT=0
          DO J=1,NE
            IF(J.EQ.L) cycle !GO TO 250
            ncn3 = ncnmax - 1  + min(1,NEN(ncnmax,J))
            DO K=1,NCN3
              NN1=IABS(NEN(K,J))
              IF(NN1.NE.N3) cycle  !GO TO 200
              NN3=K+1
              IF(NN3.GT.NCN3) NN3=1
              NN3=IABS(NEN(NN3,J))
              IF(NN3.NE.N1) cycle  !GO TO 200
!.....FOUND ONE
              ICNT=ICNT+1
              IF(ICNT.GT.1) WRITE(lp,6095) ICNT,L,M,J,K
            enddo
		  enddo
          IF(ICNT.NE.0) cycle  !GO TO 540
          if(nbc(n1).eq.0) NBC(N1)=10
          if(nbc(n3).eq.0) NBC(N3)=10
		  NES=NES+1
          if(nes.le.maxsides) then
            IES(1,NES)=N1
            IES(2,NES)=N3
          endif
        enddo
300       CONTINUE
      enddo

      if(nes.gt.maxsides)write(lp,*)' **WARNING** NES exceeds Max boundary sides'
      write(lp,*) '  NES =',NES,'  maxBSides =',maxsides

 6095 FORMAT(5X,I5,' DUPLICATE SIDES, L,M=',2I5,' J,K=',2I5)

      END

!*****************************************************************

      Subroutine BoundarySlope(nes)

      USE PMMod

! .....SET UP BOUNDRY SLOPES

      ICNT=0
      N1=IES(1,1)
      NN1=IES(2,1)
      ies(3,1) = 0

360   continue
!      write(*,*) ' Icnt=',Icnt
      DO 375 L=1,NES
!      write(*,*) ' Working on side',L
        IF(IES(1,L).NE.NN1) then
          GO TO 375
        elseif(ies(3,L).eq.0) then
          go to 374
        endif
! *** find new starting segment
        do 373 ll=2,nes
          if(ies(3,ll).ne.0) go to 373
            n1 = ies(1,ll)
            nn1 = ies(2,ll)
            write(lp,*)  ' New boundary segment, L,LL,Icnt=',L,LL,icnt,n1,nn1
            go to 360
373     continue
        write(*,*) 'End of boundary sides, L=',L,' icnt=',icnt
        write(lp,*) 'End of boundary sides, L=',L,' icnt=',icnt
        go to 380
374     continue
        N3=IES(2,L)
!      AL=SQRT((2.*XYZ(NN1,3)+XYZ(N3,3))/(2.*XYZ(NN1,3)+XYZ(N1,3)))
!      AI=1./AL
!      DY=AL*(XYZ(N3,2)-XYZ(N1,2))+(AL-AI)*(XYZ(N1,2)-XYZ(NN1,2))
!      DX=AI*(XYZ(N3,1)-XYZ(N1,1))+(AL-AI)*(XYZ(N3,1)-XYZ(NN1,1))
        DY = xyz(2,n3) - xyz(2,n1)
        DX = xyz(1,n3) - xyz(1,n1)
        if(nbc(nn1).eq.12.or.nbc(nn1).eq.17.or.nbc(nn1).eq.19) then
          if(nbc(n1).eq.10) then
            dy = xyz(2,nn1) - xyz(2,n1)
            dx = xyz(1,nn1) - xyz(1,n1)
          elseif((nbc(n3).eq.10).or.(nbc(n3).eq.100)) then
            dy = xyz(2,n3) - xyz(2,nn1)
            dx = xyz(1,n3) - xyz(1,nn1)
          endif
        endif
        IF(ABS(DX).LT.1.e-6) then
          ibc = mod(nbc(nn1),10)
          if(nbc(nn1)/10.eq.1) then
            NBC(NN1)=100 + ibc
            ALFA(NN1)=0.
          endif
        else
!        write(*,*) ' Calculating slope at node',nn1,dy,dx
          ALFA(NN1)= DY/DX
        endif
        N1=NN1
        NN1=N3
        ICNT=ICNT+1
        ies(3,L) = 1
!      write(*,*) ' Found match,L,Icnt=',L,icnt,n1,nn1
        GO TO 360
  375 CONTINUE
      write(lp,*) 'ERROR: Cant find matching node, L=',L,' nn1=',nn1

  380 CONTINUE


      END

!*****************************************************************

      Subroutine PreFront

      USE PMMod

      integer, allocatable :: NK(:),LDEST(:),LHED(:),neq(:,:)

! *** start
      ALLOCATE (nk(np),ldest(np),lhed(np),neq(np,1), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate PreFront arrays'
        stop
      endif

      jdf = 1
      nmax = np
      NEQ = 0

!....   FIND LAST APPEARANCE OF EACH NODE
      DO 12 I=1,NP
        NLAST=0
        DO 8 N=1,NE
          DO 4 L=1,NCN
            IF(NEN(L,N).NE.I)GO TO 4
            IF(NLAST.NE.N)GO TO 3
            WRITE(lp,6416) N
    3       CONTINUE
            NLAST=N
            L1=L
    4     CONTINUE
    8   CONTINUE
        IF(NLAST.EQ.0) GO TO 12
        J = NEN(L1,NLAST)
        NEN(L1,NLAST)= -J
        DO 9 K=1,jDF
    9     NEQ(J,K) = 1
          NLAST=0
   12 CONTINUE

!-.....SET UP NEQ ARRAY
      DO 110 N=1,NP
        NLAST = NBC(N)
        DO 105 M=1,jDF
          NA = 10**(jDF-M)
          NB = NLAST/NA
          NB=MOD(NB,10)
          IF(NB.gt.0) NEQ(N,M) = 0
  105   CONTINUE
  110 CONTINUE

!-.....NUMBER EQUATIONS
      NEQS = 0
      DO 150 N = 1, NP
        DO 150 M=1,jDF
          IF(NEQ(N,M).EQ.0) GO TO 150
          NEQS = NEQS +1
          NEQ(N,M)=NEQS
  150 CONTINUE

!-.....DETERMINE SIZE OF THE PROBLEM
      NDFE = NCN*jdf
      IEQ   = 0
      LEND  = 0
      LCMAX = 0
      NELL  = 0
      LCOL  = 0
   18 NELL=NELL+1
      IA=0
      DO 22 J=1,NCN
        IB=NEN(j,NELL)
        IC=IABS(IB)
        DO 22 K=1,jDF
          IA=IA+1
          NK(IA)=NEQ(IC,K)
          IF(NK(IA).LT.0) NK(IA)=0
          IF(IB.LT.0) NK(IA)=-NK(IA)
   22 CONTINUE
   23 CONTINUE

!....   SET UP HEADING VECTORS
      DO 52 LK=1,NDFE
        NODE=NK(LK)
        IF(NODE.EQ.0) GO TO 52
        IF(LCOL.EQ.0)GO TO 28
        DO 24 L=1,LCOL
          LL=L
          IF(IABS(NODE).EQ.IABS(LHED(L)))GO TO 32
   24   CONTINUE
   28   LCOL=LCOL+1
        LDEST(LK)=LCOL
        LHED(LCOL)=NODE
        GO TO 52
   32   LDEST(LK)=LL
        LHED(LL)=NODE
   52 CONTINUE
      KROW=LCOL
      IF(LCOL.LE.NMAX)GO TO 54
      WRITE(lp,6417)
      STOP
   54 CONTINUE

!....   FIND OUT WHICH MATRIX ELEMENTS ARE FULLY SUMMED
   60 JPIV=0
      DO 64 L=1,LCOL
        IF(LHED(L).GE.0)GO TO 64
        JPIV=L
        GO TO 66
   64 CONTINUE
      IF(NELL.LT.NE) GO TO 18
      WRITE(lp,6418) LCOL
      go to 136
   66 CONTINUE
      IEQ=IEQ+1
      LEND=lend+LCOL
      IF(LCOL.GT.LCMAX) LCMAX=LCOL
!     REARRANGE HEADING VECTORS
      LCOL=LCOL-1
      KROW=LCOL
      IF(JPIV.EQ.LCOL+1)GO TO 136
      DO 132 L=JPIV,LCOL
      LHED(L)=LHED(L+1)
 132  CONTINUE
 136  CONTINUE

!-.....DETERMINE WHETHER TO CONTINUE OR BACKSUBSTITUTE
      IF(NELL.LT.NE.OR.KROW.GT.1)GO TO 60
      IEQ=IEQ+1
      WRITE(lp,6490) IEQ,LCMAX,LEND
      maxrow = lcmax
      maxsto = lend

      RETURN

! *** formats
 6416 FORMAT(2X,'ELEMENT',I5,' HAS A DUPLICATE NODE NUMBER')
 6417 FORMAT(2X,'NMAX IS TO SMALL FOR FULL ASSEMBLY')
 6418 FORMAT(2X,'THERE ARE NO ROWS FULLY SUMMED. LCOL=',I5)
 6490 FORMAT( / 2X,'PREFRONT COMPLETE' / &
     &       5X,'THERE ARE ',I7,' ROWS, ',I5,' COLUMNS, AND STORAGE'&
     &    ,' LENGTH (complex) IS ',I10)

      END

!*****************************************************************

      Subroutine OutputBinary

      USE PMMod

      character ans*1
      logical exis

!-.....WRITE NODAL BASED NETWORK FILE
      INQUIRE(File=OutputFileName, Exist=exis)

      do while (exis)
        write(*,*) ' Output file already exists. Change name(y/n)? '
        read(*,'(a)') ans(1:1)
        if(ans(1:1).eq.'Y'.or.ans(1:1).eq.'y') then
          write(*,*) ' Enter new filename: '
          read(*,'(a)') OutputFileName
          INQUIRE(File=OutputFileName, Exist=exis)
        elseif(ans(1:1).eq.'N'.or.ans(1:1).eq.'n') then
          exis = .false.
        else
          write(*,*) ' Answer y or n'
        endif
      enddo

      write(*,*) ' Writing output file: np=',np,' ne=',ne

      OPEN(UNIT=21,file=OutputFileName,status='unknown',FORM='UNFORMATTED')
      WRITE(21) NE,NTYPE,NP,NPR,NCN,nsides,limngh,izup,ifront,maxrow,maxsto
      WRITE(21) ((NEN(k,J),K=1,NCN),J=1,NE),(MATL(J),J=1,NE), &
                ((XYZ(k,J),J=1,NP),K=1,3),(ALFA(J),J=1,NP), &
                 (NBC(J),J=1,NP),((nbrs(k,j),k=1,limngh),j=1,np)
      if(nsides.gt.0) then
        WRITE(21) (area(j),j=1,ne),((ieadj(k,j),k=1,5),j=1,ne),((numside(j,k),j=1,ne),k=1,4) 
        WRITE(21) ((iside(k,j),k=1,2),j=1,nsides),((sxy(k,j),k=1,2),j=1,nsides)
        WRITE(21) (refdep(j),j=1,nsides),(slen(j),j=1,nsides)
        WRITE(21) (sdx(j),j=1,nsides),(sdy(j),j=1,nsides),(dlinv(j),j=1,nsides)
        WRITE(21) ((iends(k,j),k=1,2),j=1,nsides)
      endif

      CLOSE(unit=21)

      return
      end

!*****************************************************************
!*****************************************************************

      SUBROUTINE Setup

      USE PMMod

      ALLOCATE (ieadj(5,ne), numside(ne,4), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate side adj arrays'
        stop
      endif

      write(*,*) 'Starting Setup for sides'
      nsides = 0
      numside = 0
      ieadj = 0
      DO L=1,NE
        ieadj(1,L) = L
        ncn2 = ncnmax - 1  + min(1,NEN(ncnmax,L))
        mloop: DO M=1,ncn2 !3
          if(numside(L,M).eq.0) then
            nsides=nsides+1
            numside(L,M) = nsides
          endif
          if(ieadj(M+1,L).eq.0) then
            N1=NEN(m,L)
            N2=NEN(mod(m,ncn2)+1,L)
! *** LOOK FOR THE SAME SIDE
            DO J=L+1,NE
                ncn3 = ncnmax - 1  + min(1,NEN(ncnmax,J))
                DO K=1,ncn3 !3
                  NN1=NEN(k,J)
                  NN2=NEN(mod(K,ncn3)+1,J)
                  IF(NN1.eq.N2.and.nn2.eq.n1) then
! *** FOUND ONE
                    ieadj(M+1,L) = J
                    ieadj(K+1,J) = L
                    numside(J,K) = nsides
                    if(mod(nsides,100000).eq.0) write(*,*) nsides,' sides done'
                    cycle mloop
                  endif
                enddo
            enddo
          endif
        enddo mloop
      enddo

!  Side-based arrays  !,sv(nsides)
      ALLOCATE (slen(nsides),refdep(nsides),sdx(nsides),sdy(nsides), &
      sxy(2,nsides),dlinv(nsides),iside(2,nsides),iends(2,nsides),STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate side-based storage arrays'
        stop
      endif
      iside = 0
      iends = 0
      dlinv = 0.

! *** Fill side arrays
      write(*,*) 'Finding areas'
      DO nn=1,NE
        I1=Nen(1,nn)
        I2=Nen(2,nn)
        I3=Nen(3,nn)
        Area(nn)=0.5*(xyz(1,I1)*(xyz(2,I2)-xyz(2,I3)) &
     &               +xyz(1,I2)*(xyz(2,I3)-xyz(2,I1)) &
     &               +xyz(1,I3)*(xyz(2,I1)-xyz(2,I2)))
        if(ncn.eq.4) then
          I4=Nen(ncnmax,nn)
          if(I4.gt.0) then
            Area(nn)=Area(nn) &
     &              +0.5*(xyz(1,I3)*(xyz(2,I4)-xyz(2,I1)) &
     &                   +xyz(1,I4)*(xyz(2,I1)-xyz(2,I3)) &
     &                   +xyz(1,I1)*(xyz(2,I3)-xyz(2,I4)))
          endif
        endif
      enddo

! *** Fill side arrays
      write(*,*) 'Filling side arrays'
      DO L=1,NE
        ncn2 = ncnmax - 1  + min(1,NEN(ncnmax,L))
        DO M=1,ncn2
          idxside = numside(L,M)
          if(iside(1,idxside).eq.0) then
            iside(1,idxside) = L
            N1=NEN(m,L)
            N2=NEN(mod(M,ncn2)+1,L)
            iends(1,idxside) = n1
            iends(2,idxside) = n2
            sxy(1,idxside)  = 0.5*(xyz(1,n2)+xyz(1,n1))
            sxy(2,idxside)  = 0.5*(xyz(2,n2)+xyz(2,n1))
            refdep(idxside) = 0.5*(xyz(3,n1)+xyz(3,n2))
            slen(idxside)=sqrt((xyz(1,n1)-xyz(1,n2))**2 +(xyz(2,n1)-xyz(2,n2))**2)
            sdx(idxside) = (xyz(1,n2)-xyz(1,n1))/slen(idxside)
            sdy(idxside) = (xyz(2,n2)-xyz(2,n1))/slen(idxside)
            if(ncn2.eq.3) then
              dlinv(idxside) = 2.*area(L)/(3.*slen(idxside))
            elseif(ncn2.eq.4) then
              call GetDist2cm(L,M,cmdist)
              dlinv(idxside) = cmdist  !0.5*slen(idxside)
            else
              write(*,*) "Number of element nodes not 3 or 4",L
              call exit(72)
            endif
          elseif(iside(2,idxside).eq.0) then
            iside(2,idxside) = L
            if(ncn2.eq.3) then
              dlinv(idxside) = dlinv(idxside) + 2.*area(L)/(3.*slen(idxside))
            elseif(ncn2.eq.4) then
              call GetDist2cm(L,M,cmdist)
              dlinv(idxside) = dlinv(idxside) + cmdist   !0.5*slen(idxside)
            else
              write(*,*) "Number of element nodes not 3 or 4",L
              stop
            endif
          endif
        enddo
      enddo

      do j=1,nsides
        dlinv(j) = 1./dlinv(j)
      enddo

! *** FIND LAST APPEARANCE OF EACH NODE
      write(*,*) 'Set front solution indexing'
      maxrow = 0
      maxsto = 0
      if(ifront.gt.0) then
        DO I=1,Ne
          NLAST=0
          eloop: DO N=NE,1,-1
            DO L=ncn+1,1,-1
              IF(ieadj(L,N).eq.I) then
                NLAST=N
                L1=L
                exit eloop
              endif
            enddo
          enddo eloop
          IF(NLAST.ne.0) then
            J = ieadj(L1,NLAST)
            ieadj(L1,NLAST)= -J
            NLAST=0
          else
            write(*,*) '  BIG trouble in front indexing'
            stop
          endif
        enddo
      endif

      write(*,*) ' Nsides=',nsides

      return

      open(97,file='chkarray.dat',status='unknown')
      write(97,*) 'Writing ne, area, ieadj(1:5), numside(1:4):'
      do j=1,ne
        write(97,'(i5,e10.3,9i5)') j,area(j),(ieadj(k,J),k=1,5),&
     &                        (numside(j,k),k=1,4)
      enddo
      write(97,*) 'Writing ns,iside(1:2),slen,dist,dij,refdep,sxy(1:2),sdx,sdy:'
      do j=1,nsides
        write(97,'(3i5,10e12.4)') j,(iside(k,j),k=1,2),slen(j),1./dlinv(j),&
     &                       dlinv(j),refdep(j),sxy(1,j),sxy(2,j),sdx(j),sdy(j)
      enddo
      close(97)
      stop

    end

!******************************************************************************

      SUBROUTINE GetDist2cm (nn,M,cmdist)

      USE PMMOD
      
      implicit none
      
! *** passed variables
      integer :: nn,M
      real :: cmdist

! *** local variables
      integer :: i1,i2,i3,i4,nm,n1,n2
      real :: area1,area2,asum
      real :: xc1,yc1,xc2,yc2,xc,yc,asl,bsl,csl,dsl

      I1=Nen(1,nn)
      I2=Nen(2,nn)
      I3=Nen(3,nn)
      I4=Nen(4,nn)

!      Area1 = 0.5*(mf1(I2)*(xyz(I2,1)-xyz(I1,1))*(xyz(I3,2)-xyz(I1,2)) 
!                    +mf1(I3)*(xyz(I3,1)-xyz(I1,1))*(xyz(I1,2)-xyz(I2,2)))
      Area1 = 0.5*((xyz(1,I2)-xyz(1,I1))*(xyz(2,I3)-xyz(2,I1)) &
     &               +(xyz(1,I3)-xyz(1,I1))*(xyz(2,I1)-xyz(2,I2)))

      xc1 = (xyz(1,I1) + xyz(1,I2) + xyz(1,I3))/3.
      yc1 = (xyz(2,I1) + xyz(2,I2) + xyz(2,I3))/3. 

!      Area2 = +0.5*(mf1(I3)*(xyz(I3,1)-xyz(I1,1))*(xyz(I4,2)-xyz(I1,2)) &
!     &                   +mf1(I4)*(xyz(I4,1)-xyz(I1,1))*(xyz(I1,2)-xyz(I3,2)))
      Area2 = +0.5*((xyz(1,I3)-xyz(1,I1))*(xyz(2,I4)-xyz(2,I1)) &
     &                   +(xyz(1,I4)-xyz(1,I1))*(xyz(2,I1)-xyz(2,I3)))

      xc2 = (xyz(1,I1) + xyz(1,I3) + xyz(1,I4))/3.
      yc2 = (xyz(2,I1) + xyz(2,I3) + xyz(2,I4))/3.
  
      Asum = Area1 + Area2
      xc = (Area1*xc1 + Area2*xc2)/Asum 
      yc = (Area1*yc1 + Area2*yc2)/Asum 
      
      nm = numside(nn,M)
      n2 = nen(mod(M,4)+1,nn)
      n1 = nen(M,nn)
      asl = -(xyz(2,n2) - xyz(2,n1))
      bsl =  (xyz(1,n2) - xyz(1,n1))   !*mf1(n2)
      csl = -asl*xyz(1,n1) -bsl*xyz(2,n1)
      dsl =  sqrt(asl*asl+bsl*bsl)

      cmdist = abs((asl*xc + bsl*yc + csl)/dsl)
  
      return
      end

!***************************************************************

        real*8 function ticktock()

!        USE DFPORT

        real*4 a(2)

        ticktock = 0.
!       for sun      
        ticktock = etime(a)
       
        return
        end

!***************************************************************
        subroutine shape(iflg,iord,nipt)

        USE PMMod

        implicit double precision (a-h,o-z)      

! this subroutine evaluates the shape functions at a set of integration
! points

! iflg -- =  0 to eval at integ pts 
!         =  > 0 to eval at the nodes of degree iflg
! iord -- = degree of poly basis
! nipt -- = total number of integ points, if iflg > 0 this value
!           is calculated by the routine
      
!      include 'param.inc'
!      include 'prem2com.h'

        parameter (izr=2,iza=(izr+1)*(izr+2)/2,iz6=10)
!      parameter (maxncn=50)
!      parameter (maxni=100)
 
!       this block includes the shape functions and derivs evaluated
!       at the integ pts        
        real*8 psit(iza,iz6),dpsic(iza,iz6),dpsie(iza,iz6)

!       these are the integ pts and weights         
        real*8 xi(iz6,3),wxi(iz6),iptloc(100)

!     common /ebasis/ xnn(maxncn,maxni),weight(maxni),
!    *         dnan(maxncn,maxni),dnbn(maxncn,maxni),
!    *         xnh(maxncn,maxni),wth(maxni),xnp(maxncn,3),
!    *         dnah(maxncn,maxni),dnbh(maxncn,maxni),
!    *         xm3(maxncn,maxni),dma3(maxncn,maxni),dmb3(maxncn,maxni),
!    *         jordr,ngp,ncn3

!       real*8 xnh,dnah,dnbh,wth
!       real*8 xnp,xnn,dnan,dnbn,weight,xm3,dma3,dmb3
!       integer ngp,ncn3,jordr

        integer iflg,iord,nipt
        real*8 psi(iza,iz6),dpsi(iza,2,iz6)

!       area coords and partials       
        real*8 a1,a2,a3
        real*8 da1dcassi,da2dcassi,da3dcassi,da1deta,da2deta,da3deta

!       evaluation points
        real*8 xp(iz6,3)

!       permutation index to get the ordering we want
        integer iperm(iza), iperm2(iza)


!       compute the number of nodes 
        numnod = (iord+1)*(iord+2)/2
        if(numnod .gt. iza)then
          write(lp,*)' error shape function too high a degree '
          write(lp,*)' nodes: ',numnod,iza
          stop
        endif

        if(nipt .gt. iz6)then
          write(lp,*)' attempted to evaluate shapes at too many pts '
          write(lp,*)' points: ',nipt,iz6
        endif

!       check the integration rule
        if(iflg .eq. 0)then
          if(iptloc(nipt) .eq. 0)then
            write(lp,*)' error integ rule not supported ',nipt
            stop
          endif
        endif


!       load the nodal evaluation points
        if(iflg .gt. 0)then
          nipt = (iflg+1)*(iflg+2)/2
          icnt = 0
          do j = 0,iflg
            do i = 0,iflg
              if(i+j .le. iflg)then
                 
                icnt = icnt + 1 
                cassi = dfloat(i)/dfloat(iflg)
                eta = dfloat(j)/dfloat(iflg)
                side =  dfloat(iflg - i - j)/dfloat(iflg)
                xp(icnt,1) = side
                xp(icnt,2) = cassi
                xp(icnt,3) = eta
                  
              endif
            enddo
          enddo

!       now permute to the ordering we want
!       load permutation index:
        isid1 = 0
        isid2 = 0
        isid3 = 0
        iint = 0
        ish = 0
!       loop on current ordering        
          do j = 0,iflg
          do i = 0,iflg

!           note constraint  0 \leq i+j \leq iflg               
            if(i+j .le. iflg)then

!             increment the shape function index
              ish = ish + 1

!             compute k value                 
              k = iflg - i - j

!             corner indices
              if(i .eq. 0 .and. j .eq. 0)then
                 
                iperm2(ish) = 1
                
              elseif(i .eq. iflg .and. j .eq. 0)then
                 
                iperm2(ish) = iflg + 1
                
              elseif(i .eq. 0 .and. j .eq. iflg)then
                 
                iperm2(ish) = 2*iflg + 1
                
!             side indices                
              elseif(i .gt. 0 .and. i .lt. iflg .and. j .eq. 0)then
                 
                isid1 = isid1 + 1
                iperm2(ish) = 1 + isid1
                
              elseif(j .gt. 0 .and. j .lt. iflg .and. i .eq. 0)then
                 
                isid3 = isid3 + 1
                iperm2(ish) = 3*iflg + 1 - isid3
                
              elseif(i .lt. iflg .and. j .lt. iflg .and. k .eq. 0)then
                 
                isid2 = isid2 + 1
                iperm2(ish) = iflg + 1 + isid2

!             interior indices                
              elseif(i .gt. 0 .and. j .gt. 0 .and. k .gt. 0)then
                 
                iint = iint + 1
                iperm2(ish) = 3*iflg + iint
                
              else
                 
                write(lp,*)' error in iflg permuattion index ',i,j,k
                stop
                
              endif
            endif
          enddo
          enddo

        elseif(iflg .eq. 0)then

          do in=1,nipt
            xp(in,1) = xi(in,1)
            xp(in,2) = xi(in,2)
            xp(in,3) = xi(in,3)            
          enddo
          
        endif

        
!       set order parameter
        r = dfloat(iord)

!       set some partial derivs
        da1dcassi = -1.d0
        da2dcassi = 1.d0
        da3dcassi = 0.d0

        da1deta = -1.d0
        da2deta = 0.d0
        da3deta = 1.d0


!       evaluate shape funcs at integs;
!       loop on integs        
        do in=1,nipt

!         dump out where we evaluate; integs
          if(iflg .eq. 0)then

!           write(6,*)' evaluating shapes at integ pts basis '

!         or nodes 
          elseif(iflg .gt. 0)then

!           write(6,*)' evaluating shapes at the nodes '

          else

            write(lp,*)' cant eval shapes for this iflg ',iflg

          endif

!         well, copy them over - these are area coords
!         a1 - side, a2 - cassi, a3 - eta          
          a1 = xp(in,1)
          a2 = xp(in,2)
          a3 = xp(in,3)

!         zero the shape function index
          ish = 0


!         loop on area coord indices respectively           
          do j = 0,iord
            do i = 0,iord

!             note constraint  0 \leq i+j \leq iord               
              if(i+j .le. iord)then

!               increment the shape function index
                ish = ish + 1

!               compute k value                 
                k = iord - i - j              

!               compute the three product terms and deriv terms:
!               cassi
                cprod = 1.d0
                do l = 0,i-1
                  cprod = cprod * (dfloat(l)-r*a2)/(dfloat(l)-dfloat(i))
                enddo
! ROY
                cprodp = 0.d0
                do ip = 0,i-1
                  aprod = 1.d0
                  do l = 0,i-1
                    if(l .ne. ip)then
                      aprod=aprod*(dfloat(l)-r*a2)/(dfloat(l)-dfloat(i))
                    endif
                  enddo
                  cprodp = cprodp + aprod*(-r)/(dfloat(ip)-dfloat(i))
                enddo


                
!               eta
                eprod = 1.d0
                do m = 0,j-1
                  eprod = eprod * (dfloat(m)-r*a3)/(dfloat(m)-dfloat(j))
                enddo
! ROY
                eprodp = 0.d0
                do ip = 0,j-1
                  aprod = 1.d0
                  do m = 0,j-1
                    if(m .ne. ip)then
                      aprod=aprod*(dfloat(m)-r*a3)/(dfloat(m)-dfloat(j))
                    endif
                  enddo
                  eprodp = eprodp + aprod*(-r)/(dfloat(ip)-dfloat(j))
                enddo



                
!               side
                sprod = 1.d0
                do n = 0,k-1
                  sprod = sprod * (dfloat(n)-r*a1)/(dfloat(n)-dfloat(k))
                enddo
! ROY
                sprodp = 0.d0
                do ip = 0,k-1
                  aprod = 1.d0
                  do n = 0,k-1
                    if(n .ne. ip)then
                      aprod=aprod*(dfloat(n)-r*a1)/(dfloat(n)-dfloat(k))
                    endif
                  enddo
                  sprodp = sprodp + aprod*(-r)/(dfloat(ip)-dfloat(k))
                enddo
                

!               now the value of the shape function 
                psi(ish,in) = cprod*eprod*sprod

!               derivatives w/res to area coords, a1, a2, a3                
                dpsida1 = cprod*eprod*sprodp
                dpsida2 = cprodp*eprod*sprod
                dpsida3 = cprod*eprodp*sprod

!               derivs w/res to cassi, eta                
                dpsi(ish,1,in) = dpsida1*da1dcassi + dpsida2*da2dcassi + dpsida3*da3dcassi                 
                
                dpsi(ish,2,in) = dpsida1*da1deta + dpsida2*da2deta + dpsida3*da3deta                 


              endif
            enddo
          enddo
        enddo


!       now permute to the ordering we want
!       load permutation index:
        isid1 = 0
        isid2 = 0
        isid3 = 0
        iint = 0
        ish = 0
!       loop on current ordering        
        do j = 0,iord
          do i = 0,iord

!           note constraint  0 \leq i+j \leq iord               
            if(i+j .le. iord)then

!             increment the shape function index
              ish = ish + 1

!             compute k value                 
              k = iord - i - j

!             corner indices
              if(i .eq. 0 .and. j .eq. 0)then
                 
                iperm(ish) = 1
                
              elseif(i .eq. iord .and. j .eq. 0)then
                 
                iperm(ish) = iord + 1
                
              elseif(i .eq. 0 .and. j .eq. iord)then
                 
                iperm(ish) = 2*iord + 1
                
!             side indices                
              elseif(i .gt. 0 .and. i .lt. iord .and. j .eq. 0)then
                 
                isid1 = isid1 + 1
                iperm(ish) = 1 + isid1
                
              elseif(j .gt. 0 .and. j .lt. iord .and. i .eq. 0)then
                 
                isid3 = isid3 + 1
                iperm(ish) = 3*iord + 1 - isid3
                
              elseif(i .lt. iord .and. j .lt. iord .and. k .eq. 0)then
                 
                isid2 = isid2 + 1
                iperm(ish) = iord + 1 + isid2

!             interior indices                
              elseif(i .gt. 0 .and. j .gt. 0 .and. k .gt. 0)then
                 
                iint = iint + 1
                iperm(ish) = 3*iord + iint
                
              else
                 
                write(lp,*)' error in permuattion index ',i,j,k
                stop
                
              endif
            endif
          enddo
        enddo

!       now copy them over
        do in=1,nipt
!          weight(in) = wxi(in)
          do ish=1,numnod
            psit(iperm(ish),in) = psi(ish,in)
            dpsic(iperm(ish),in) = dpsi(ish,1,in)
            dpsie(iperm(ish),in) = dpsi(ish,2,in)
            if(iflg.gt.0) then
              xnn(iperm(ish),iperm2(in)) = psi(ish,in)
!              dnan(iperm(ish),iperm2(in)) = dpsi(ish,1,in)
!              dnbn(iperm(ish),iperm2(in)) = dpsi(ish,2,in)
            else
              xnn(iperm(ish),in) = psi(ish,in)
!              dnan(iperm(ish),in) = dpsi(ish,1,in)
!              dnbn(iperm(ish),in) = dpsi(ish,2,in)
            endif
          enddo
        enddo

! this routine stores the values of the shape functions at the requested points
! in array psit,
!          psit(shape_function_index,evaluation_point_index)
! where the evaluation_point_index is typically the integration point index
! the shape functions are ordered as follows
! first the three corners,
!       corner 1: cassi = 0, eta = 0
!       corner 2: cassi = 1, eta = 0
!       corner 3: cassi = 0, eta = 1

! next the three sides;
!       side 1: eta = 0, increasing shape function index with increasing cassi
!       side 2: cassi+eta = 1, increasing shape function index with increasing eta
!                              and decreasing cassi
!       side 3: cassi = 0 increasing shape function index with increasing eta

! next the interior nodes; ordering is probably irrelevant

! cassi is the horizontal axis
! eta is the vertical axis
! derivatives w/res to cassi are in dpsic(,); derivs w/res to eta are in dpsie(,)
!        

        
!       that's all folks
        return
        end

!***************************************************************
      subroutine expand_basis(iord)

      USE PMMod

!	This program takes the neigh [NEIGH] and triangle [TRIANG] files
!	and refines and expands the appropriate arrays.

      integer, allocatable :: nodeb(:)
      integer nbn, stat
      integer, parameter :: maxncn=6
      logical CalcAlfa

      write(lp,*) ' Start expand'
      maxne = 4*ne
      DO J=1,maxne
        DO K=4,maxncn
          NEN(J,K)=0
        enddo
        nen(j,2*iord+1)=nen(j,3)
        nen(j,3)=0
        nen(j,iord+1)=nen(j,2)
        nen(j,2)=0
      enddo

      open(83,file='bpoints.dat',status='old')
      read(83,*) maxnp

      ALLOCATE (nodeb(maxnp), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate nodeb to expand bases'
        stop
      endif
      
      nbn = 0
      nb1 = 1
      do j=1,maxnp
        read(83,*,end=19) n1
        if(n1.gt.0) then
          nbn = nbn + 1
          nodeb(nbn) = n1
        elseif((nbn.gt.nb1).and.(n1.eq.0)) then
          nbn = nbn + 1
          nodeb(nbn) = nodeb(nb1)
          nb1 = nbn + 1
        endif
      enddo
19    close(83)
      nodeb(nbn+1) = nodeb(nb1)
      write(lp,*) ' There are',nbn,' boundary points'
!      write(6,*) (nodeb(j),j=1,nbn)

      NCN=(iord+1)*(iord+2)/2
      nsn = 3*iord
!      ncn3 = 3
      NPN=NP
      DO 350 L=1,NE
        DO 300 M=1,NsN,iord
          N1=IABS(NEN(L,M))
          N3=M+iord
          IF(N3.GT.NsN) N3=1
          N3=IABS(NEN(L,N3))
! Is this a boundary side?
          iboundary=0
          do jj=1,nbn
            if(nodeb(jj).eq.n1) then
              if(nodeb(jj+1).eq.n3) then
                write(lp,*) ' Boundary between',n1,n3
                iboundary=1
                go to 21
              endif
            endif
          enddo
21        continue
! *** loop over mid-side nodes
          do kord = 1,iord-1
            CalcAlfa = .false.
            fct = float(kord)/float(iord)
            N2=IABS(NEN(L,M+kord))
            IF(N2.NE.0) GO TO 300
            NPN=NPN+1
            if(npn.gt.maxnp) then
              write(lp,*) ' npn exceeds maxnp',npn,maxnp
              stop
            endif
            alfa(npn) = 0.
            NEN(L,M+kord)=NPN
            if(iboundary.eq.1) then

            ibc1 = nbc(n1)
            ibc3 = nbc(n3)
            ih1 = mod(ibc1,10)
            ih3 = mod(ibc3,10)
            if(ih1.gt.0.and.ih3.gt.0) then
              ih2 = ih1
            else
              ih2 = 0
            endif
            ibc1 = ibc1/10
            ibc3 = ibc3/10
            if(ibc1.eq.0.or.ibc3.eq.0) then
              ibc2 = 0
            elseif(ibc1.eq.11.and.ibc3.eq.10) then
              ibc2 = ibc3
            elseif(ibc1.eq.11.and.ibc3.eq.1) then
              ibc2 = ibc3
              CalcAlfa = .true.
!              alfa(NPN)=alfa(N1)*(1.-fct) +alfa(N3)*fct
            elseif(ibc1.eq.10.and.ibc3.eq.11) then
              ibc2 = ibc1
            elseif(ibc1.eq.1.and.ibc3.eq.11) then
              ibc2 = ibc1
              CalcAlfa = .true.
!              alfa(NPN)=alfa(N1)*(1.-fct) +alfa(N3)*fct
            elseif(ibc1.eq.11.and.ibc3.eq.11) then
              ibc2 = 11
            elseif(ibc1.eq.10.and.ibc3.eq.10) then
              ibc2 = 10
            elseif(ibc1.eq.ibc3) then
              if(ih2.gt.0) then
                ibc2 = 0
              else
                ibc2 = ibc1
              endif
              if(ibc2.eq.1) then
                CalcAlfa = .true.
!                alfa(NPN)=alfa(N1)*(1.-fct) +alfa(N3)*fct
              endif
            else
              ibc2 = 0
            endif
! *** end ibc loop

            if(CalcAlfa) then
              DY=XYZ(N3,2)-XYZ(N1,2)
              DX=XYZ(N3,1)-XYZ(N1,1)
              IF(ABS(DX).LT.1.e-6) then
                ibc2 = 10
                ALFA(npn)=0.
              else
                ALFA(npn)= atan(DY/DX)
              endif
            endif

            nbc(npn) = 10*ibc2 + ih2

!           if(ih2.gt.0) then
              write(lp,*)' bc at',npn,' code=',nbc(npn),' between',n1,n3
!           endif

            else
              nbc(npn) = 0
            endif
! *** end iboundary loop

            fct = float(kord)/float(iord)
            xyz(NPN,1)=xyz(N1,1)*(1.-fct) +xyz(N3,1)*fct
            xyz(NPN,2)=xyz(N1,2)*(1.-fct) +xyz(N3,2)*fct
            xyz(NPN,3)=xyz(N1,3)*(1.-fct) +xyz(N3,3)*fct
! *** LOOK FOR THE SAME SIDE
210     ICNT=0
        DO 250 J=1,NE
          IF(J.EQ.L) GO TO 250
          DO 200 K=1,NsN,iord
            NN1=IABS(NEN(J,K))
            IF(NN1.NE.N3) GO TO 200
            NN3=K+iord
            IF(NN3.GE.NsN) NN3=1
            NN3=IABS(NEN(J,NN3))
            IF(NN3.NE.N1) GO TO 200
            NN2=IABS(NEN(J,K+iord-kord))
            IF(NN2.NE.0) WRITE(lp,6090) J,K,NN1,NN2,NN3,L,M,N1,NPN,N3
            NEN(J,K+iord-kord)=npn
            ICNT=ICNT+1
            IF(ICNT.GT.1) WRITE(lp,6095) ICNT,L,M,J,K
200       CONTINUE
250     CONTINUE
      enddo
!   *** end kord loop
300   CONTINUE
!   *** add interior nodes
      if(ncn.gt.nsn) then
        nc1=iabs(nen(L,1))
        nc2=iabs(nen(L,1+iord))
        nc3=iabs(nen(L,1+2*iord))
        do j=nsn+1,ncn
          npn=npn+1
          nen(L,j)=-npn
          xyz(npn,1) = xyz(nc1,1)*xnp(j,1) + xyz(nc2,1)*xnp(j,2) + xyz(nc3,1)*xnp(j,3)
          xyz(npn,2) = xyz(nc1,2)*xnp(j,1) + xyz(nc2,2)*xnp(j,2) + xyz(nc3,2)*xnp(j,3)
          xyz(npn,3) = xyz(nc1,3)*xnp(j,1) + xyz(nc2,3)*xnp(j,2) + xyz(nc3,3)*xnp(j,3)
        enddo
      endif

350   CONTINUE

      write(lp,*) ' After expand, np=',npn,' ,ncn=',ncn
! *** write node and triangle file
!     open(83,file='testni.nod',status='unknown')
!     write(83,*) npn
!     do j=1,npn
!       write(83,8301) j,xyz(j,1),xyz(j,2),nbc(j),xyz(j,3)
!     enddo
!     close(83)
8301  format(i5,2(1x,f10.2),1x,i4,1x,f8.2)
!     open(83,file='testni.tri',status='unknown')
!     write(83,*) ne
!     do j=1,ne
!       write(83,8305) (iabs(nen(j,jj)),jj=1,ncn)
!     enddo
!     close(83)
8305  format(20(1x,i5))

!      do j=np+1,npn
!        laste=0
!        lastn=0
!        do nn=1,ne
!          do kk=1,nsn,iord
!            do kkk=1,iord-1
!            n1 = nen(nn,kk+kkk)
!              if(n1.eq.j) then
!                laste = nn
!                lastn = kk+kkk
!              endif
!            enddo
!          enddo
!        enddo
!        if(laste.gt.maxne.or.lastn.gt.maxncn) then
!          write(53,*) ' Error in expand',j,laste,lastn
!          stop
!        endif
!        if(laste.gt.0) then
!          nen(laste,lastn) = -nen(laste,lastn)
!        endif
!      enddo

      np = npn

      return

6090  FORMAT(5X,'SIDE ERROR- ',10I6)
6095  FORMAT(5X,I3,' DUPLICATE SIDES, L,M=',2I5,'  J,K=',2I5)

      END

!***************************************************************
        
