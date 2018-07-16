!***************************************************************

  Module RCMArrays
  
      implicit none

      integer :: ne,np,nsides,npv,npvC,ncn   !nph,nphu,npv  
      integer :: nson,nsed,nsol
      integer :: neqtide=0, neMB=0, nsbc, iOPsol=0
      integer :: noptcount=0, nopt=2, outfileopt=0
      integer :: izcoord=2,izgrid=0,ixycoord
      integer :: jUprofile=0,jCprofile=0, jSprofile=0
      integer, allocatable :: nen(:,:),numsideT(:,:),iside(:,:),IECode(:)
      integer, allocatable :: iends(:,:)
      integer, allocatable :: tsUnodes(:),tsSnodes(:),tsCnodes(:)
      real*8 ::  TET,tet0
      real*8 :: x0off,y0off
      real*8, parameter :: depmin=0.00001
      real*8, parameter :: bigr = 6378136.D0
      real*8, parameter :: pi = 3.141592653589793D0
      real*8, parameter :: rad2deg=180.D0/pi
      real*8, parameter :: deg2rad=pi/180.D0
      real*8, allocatable ::  xp(:),yp(:),zp(:),zdep(:),area(:)
      real*8, allocatable ::  xpt(:),ypt(:),zpt(:),uct(:),vct(:)
      real*8, allocatable ::  eta(:),un(:,:),ut(:,:),wz(:,:),etaMB(:)
      real*8, allocatable ::  uc(:,:),vc(:,:),eqtide(:)
      real*8, allocatable ::  sxy(:,:),sbot(:),sdx(:),sdy(:),sdep(:)
      real*8, allocatable ::  rhv(:),gamma(:),Rn0(:,:) !scratch vectors
      real*8, allocatable ::  qp(:,:),cc(:,:),sigt(:,:),PSU(:,:),TC(:,:)
      real*8, allocatable ::  emax(:),emin(:),spdmax(:),umax(:),vmax(:),t1(:),twet(:)
      character(256) :: OutResFile, OutMaxFile
      character(20) :: outopt
 
  end module RCMArrays

!*********************************************************************

      Program RCM2netCDF
  
      use RCMArrays
  
      implicit none

      integer :: i,j,k,istat,npvgrd,npdim
      integer :: numarg, iargc
      character*256 :: fnamedata='', fnameprofile='', fnamemax=''
      character(18) :: AnalysisTime
      character(len=256) :: arg
      logical :: notOK=.true.
       
!     read time slices from the binary data file
      numarg=iargc()

      if(numarg.eq.0) then
        write(*,*) ' Enter filename for input RiCOM binary file'
        read(*,'(a)') fnamedata
        write(*,*) ' Enter filename for profile points input. CR=none.'
        read(*,'(a)') fnameprofile
      endif
      
      if(numarg.ge.1) then
        i = 1
        call getarg(i,arg)
        fnamedata = arg
      endif
      
      if(numarg.ge.2) then
        i = 2
        call getarg(i,arg)
        fnameprofile = arg
      endif

      open(unit=20,file=fnamedata,status='old',form='unformatted')
      if(fnameprofile.ne.'') open(unit=21,file=fnameprofile,status='old')

      OutResFile = trim(fnamedata)//'.nc'

      read(20) AnalysisTime
      write(*,*) 'AnalysisTime= ',AnalysisTime
      read(20)  !skip windfilename

      read(20) ne,np,nsides,npv,ncn,izcoord,izgrid,ixycoord,x0off,y0off
      read(20) nson,nsed,nsol,neMB,neqtide,iOPsol

      npv = max(npv,1)
      if(npv.gt.1) then 
        nopt= 3
        if(izcoord.lt.2) then
          npvC = npv-1
        else
          npvC = npv
        endif
      else
        nopt= 2
        npvC = 1
      endif

      write(*,*) 'ne,np,nsides,npv=',ne,np,nsides,npv
      
      npvgrd = npvC + 1  !min(1,max(0,izcoord-1))
      npdim = (np+nsides)*npvgrd
        
      ALLOCATE ( xp(np), yp(np), zp(np), nen(ne,ncn), &
          sxy(2,nsides),sbot(nsides),sdx(nsides),sdy(nsides),Rn0(npvgrd,nsides), &
          numsideT(ncn,ne),iside(2,nsides),area(ne),IECode(ne),iends(2,nsides), &
          uc(npvgrd,np),vc(npvgrd,np),zdep(npvgrd),rhv(nsides),gamma(nsides),sdep(nsides), &
          uct(npdim),vct(npdim),xpt(npdim),ypt(npdim),zpt(npdim), &
          eta(ne),un(npvgrd,nsides),ut(npvgrd,nsides),wz(npvc,ne), eqtide(ne), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate ncon main storage arrays',istat
        stop
      endif
        
      if(izcoord.eq.10) then
        ALLOCATE ( etaMB(ne), STAT = istat )
        if(istat.ne.0) then
          write(*,*) 'FATAL ERROR: Cannot allocate etaMB array',istat
          stop
        endif
      endif
        
      if(nson.ne.0) then
        ALLOCATE ( qp(npv,ne), STAT = istat )
        if(istat.ne.0) then
          write(*,*) 'FATAL ERROR: Cannot allocate qp array',istat
          stop
        endif
      endif

      if(nsol.gt.0) then
        ALLOCATE ( sigt(npvc,ne), STAT = istat )
        if(istat.ne.0) then
          write(*,*) 'FATAL ERROR: Cannot allocate sigt array',istat
          stop
        endif
        if(iOPsol.eq.1) then
          ALLOCATE ( PSU(npvc,ne), STAT = istat )
          if(istat.ne.0) then
            write(*,*) 'FATAL ERROR: Cannot allocate PSU array',istat
            stop
          endif
        elseif(iOPsol.eq.2) then
          ALLOCATE ( TC(npvc,ne), STAT = istat )
          if(istat.ne.0) then
            write(*,*) 'FATAL ERROR: Cannot allocate TC array',istat
            stop
          endif
          read(20) (TC(1,j),j=1,ne)
        elseif(iOPsol.eq.3) then
          ALLOCATE ( PSU(npvc,ne),TC(npvc,ne), STAT = istat )
          if(istat.ne.0) then
            write(*,*) 'FATAL ERROR: Cannot allocate PSU,TC arrays',istat
            stop
          endif
        endif
      endif
      
! *** read coordinates
      read(20) (xp(j),j=1,np)
      read(20) (yp(j),j=1,np)
      read(20) (zp(j),j=1,np)
    
      if(npv.gt.1) then
        read(20) (zdep(k),k=1,npvgrd)    !write  z-grid (m)
      endif

! *** read element list
      read(20) ((nen(j,k),k=1,ncn),j=1,ne),(area(j),j=1,ne),(IECode(j),j=1,ne)
      read(20) (sdx(j),j=1,nsides),(sdy(j),j=1,nsides), &
               ((sxy(i,j),i=1,2),j=1,nsides),(sbot(j),j=1,nsides)
      read(20) ((numsideT(i,j),i=1,ncn),j=1,ne),((iside(i,j),i=1,2),j=1,nsides)
      read(20) ((iends(i,j),i=1,2),j=1,nsides)
      
! *** initialize
      tet0 = -1.D0
      un = 0.D0
      ut = 0.D0
      wz = 0.D0
      
      do
        read(20, IOSTAT=istat) TET 
        if(istat.gt.0) then
          write(*,*) 'READ ERROR =',istat,'exiting...'
          exit
        elseif(istat.lt.0) then
          write(*,*) 'END OF FILE - no max min data to process'
          exit
        elseif(tet.lt.tet0) then  !go to writing max/min file
          exit
        endif
        tet0 = tet
        
        read(20, IOSTAT=istat) (eta(j),j=1,ne)
! *** read velocity (un,ut)
        if(istat.eq.0) read(20, IOSTAT=istat) ((un(k,i),k=1,npv),i=1,nsides)
        if(istat.eq.0) read(20, IOSTAT=istat) ((ut(k,i),k=1,npv),i=1,nsides)    
        if(istat.eq.0.and.(npv.gt.1.or.nson.gt.0)) read(20, IOSTAT=istat) ((wz(k,i),k=1,npvc),i=1,ne)

! *** read bottom location if moving bottom
        if(neMB.gt.0.and.noptcount.gt.0) then
          read(20) (zp(j),j=1,np)
        endif
  
! *** read etaMB, uMB if dynamic landslide
        if(izcoord.eq.10) then
          read(20) (etaMB(j),j=1,ne)
          read(20) (un(2,j),j=1,nsides)    !write  unMB (m/s)
          read(20) (ut(2,j),j=1,nsides)    !write  utMB (m/s)
        endif

        if(nson.gt.0.and.istat.eq.0) then
          read(20, IOSTAT=istat) ((qp(k,j),k=1,npv),j=1,ne)
        endif

        if(nsed.gt.0.and.istat.eq.0) then
          read(20, IOSTAT=istat) ((cc(k,j),k=1,npvc+1),j=1,ne)
        endif

        if(nsol.gt.0.and.istat.eq.0) then
          read(20, IOSTAT=istat) ((sigt(k,j),k=1,npvc),j=1,ne)
          if(iOPsol.eq.1) then
            read(20) ((PSU(k,j),k=1,npvc),j=1,ne)
          elseif(iOPsol.eq.2) then
            read(20) ((TC(k,j),k=1,npvc),j=1,ne)
          elseif(iOPsol.eq.3) then
            read(20) ((PSU(k,j),k=1,npvc),j=1,ne)
            read(20) ((TC(k,j),k=1,npvc),j=1,ne)
          endif
        endif

        if(neqtide.gt.0.and.istat.eq.0) then  ! equilibrium tide results
          read(20, IOSTAT=istat) (eqtide(j),j=1,ne)
        endif

        if(istat.ne.0) then
          write(*,*) ' Error reading input file at time=', TET
          exit
        endif
      
! *** write netCDF file
        if(noptcount.eq.0) write(*,*) ' call ncout, nopt=',nopt
        call WritenetCDFData(fname,Quit)

        noptcount = noptcount + 1

      enddo
      
! *** write maxmin results
      if(istat.eq.0) then
        ALLOCATE ( emax(ne),emin(ne),umax(nsides),vmax(nsides),t1(ne),twet(ne), STAT = istat )
        if(istat.ne.0) then
          write(*,*) 'FATAL ERROR: Cannot allocate max arrays',istat
          stop
        endif

        read(20, IOSTAT=istat) (emax(j),j=1,ne),(emin(j),j=1,ne)
        read(20, IOSTAT=istat) (umax(j),j=1,nsides),(vmax(j),j=1,nsides)
        read(20, IOSTAT=istat) (t1(j),j=1,ne),(twet(j),j=1,ne)
        if(istat.ne.0) then
          write(*,*) ' Error reading input max file at time=', TET
          stop
        endif

! *** write netCDF file
        OutMaxFile = trim(fnamedata)//'.max.dat'
        noptcount = 0
        nopt = 4
        write(*,*) ' call ncout, nopt=',nopt
        call WritenetCDFData(fname,Quit)
       
      endif

! *** clean up
      close(20)

      stop
      end

!*********************************************************************  
 
      subroutine WritenetCDFData(fname,Quit)
 
      use RCMArrays
      USE UGrid_netCDFio
     
      implicit none
      
! *** passed varaibles
      character(len=*) :: fname
      logical quit

! *** local variables
      integer :: j,k,nbp
      integer :: np,ne,ncn,ixy,iz,ncount,dcount
      integer :: nindex(itot),nen(3,TotTr)
      integer :: NUMofBDYS, bnode_index(itot), bnode_id(itot), rnodes(2*itot)
      logical :: err
      real*8 :: x0,y0
      character(3) :: uz
      
      quit = .false.
      err = .false.
      
      ixy = ixycoord
      if(scaleY.lt.0.) then
        iz = 0
      else
        iz = 1
      endif
      x0 = x0off
      y0 = y0off
      uz(1:3) = UTMzone(1:3)
      call Create_Grid_netCDF (np,ne,ncn,iXY,iZ,uz,fname,err)
      if(err) then
        quit = .true.
        return
      endif

      if(ncn.eq.3) then
        call Write_Grid_netCDF (np,ne,ncn,x0,y0,xp,yp,zp,code,IEcode,nen,err)
      elseif(ncn.eq.4) then
        call Write_Grid_netCDF (np,ne,ncn,x0,y0,zp,yp,zp,code,IEcode,nen,err)
      endif
      if(err) then
        quit = .true.
        return
      endif
      
      call Boundary_list(np,nl,nbtot,dxray,dyray,code, &
                         nbp,NUMofBDYS,bnode_index, bnode_id,rnodes)

      call Write_Boundary_netCDF (nbp,NUMofBDYS,bnode_index, bnode_id,rnodes,err)

                 
      end subroutine

!************************************************************************
        
        SUBROUTINE Boundary_list(itot,nl,nbtot,dxray,dyray,ncode,&
                                 nbp,nbnd,bnode_index, bnode_id,rnodes)

           IMPLICIT NONE

!       WARNING - Boundary nodes retrieved in RETRONODES are stored here in
!                 local array RNODES( ). The nodes for each boundary are in a
!                 consecutive block, but the order of the nodes within each 
!                 each block is the reverse of the order usual in NODE format
!                 files, i.e. in rnodes( ), the outer boundary is in clockwise 
!                 order and islands are in counterclockwise order. Also, the outer 
!                 boundary is not necessarily the first block & rnodes( ) does not
!                 have the interior nodes in it. If the RETRONODE function in
!                 GRIDGEN is supposed to produce a properly ordered set of nodes,
!                 equivalent to reading in a NODE format file, then 3 further 
!                 steps need to be implemented:
!                  1) Put outer boundary first
!                  2) Reverse order of nodes in each boundary  
!                  3) Add interior nodes to (retro) node array

    
!       Purpose: To form a boundary list for a grid.
!                - Some error detection capability.
!
!       Input variables
!       itot   = number of nodes in grid
!       nbtot   = number adjacent nodes


! *** Passed variables ***

      integer, intent(in) :: itot,nbtot,ncode(itot)
      integer, intent(inout) :: nl(nbtot,itot+1)
      INTEGER, intent(out) :: RNODES(2*itot),nbnd
      integer, intent(out) :: nbp, bnode_index(itot+1), bnode_id(itot)
      real*8, intent(in) :: dxray(itot), dyray(itot)


! *** Local variables ***

      INTEGER NextVer                     ! Integer functions
      integer ListTr(4)
      INTEGER IN, SNB(NBTOT+1)
      INTEGER numverts, ThisJ, NextCcwNbr, NextJ, NextVertex
      INTEGER I,J,K,M
!      LOGICAL ELIVE

!      ListTr(1, ) to ListTr( ,4) are nodes of element in ccw order

      INTEGER NUMN(itot)   ! no. of neighbours of each node

!---  variables added 25/9/98 for RETRONODE operation ------------------------

      integer, parameter :: maxnnb=10000
      integer, save :: NUMperBDY(MAXNNB)
      INTEGER Nfound,NUMofBDYS

!     PASSED variable added 25/9/98 for RETRONODE operation
      LOGICAL :: retrowanted=.true.

!    Local variables added 25/9/98 for RETRONODE operation
      INTEGER firstnode
      LOGICAL followboundary

!     Local variables added 25/9/98 for output of retrieved node file
      INTEGER  sumofnodes, numouter, presum
      REAL*8 minimumx

!---  end of extra variable section ------------------------------------------

! ----BEGIN------------------------------------------------------------

!     number of boundaries so far detected
      NUMofBDYS = 0
!     number of boundary nodes found so far
      Nfound = 0
      NUMperBDY = 0

!     Count sort (ccw from east) the neighbours of each node 
      do I = 1, itot
        numn(I) = 0
        call SORTNBRS(I,NUMN(I),SNB,itot,nbtot,nl,dxray,dyray)
        do J = 1, NUMN(I)
          NL(J,I) = SNB(J)
        enddo
        if(numn(i).lt.nbtot) then
          do j=numn(i)+1,nbtot
            nl(j,i) = 0
          enddo
        endif
      enddo

!      Neighbour lists now sorted ccw, counted and packed left
!      (no need to set unused neighbours spaces at right to zero)

!     Set element counts to zero
      ListTr = 0
 !     LISTERR = 0

! ----------------------------------------------------------------------
!     BEGIN SEARCH FOR ELEMENTS
!     For each node in grid
      DO IN = 1, itot

!       first look after degenerate elements, nodes with 0 or 1 neighbours 
        if (numn(IN).lt.2) then
!           case of node with no or one neighbour.
          cycle
        endif

!       - from here on, node IN assumed to have at least 2 neighbours
!      Begin element check around each node, starting from first neighbour
        j = 1
!        do 60 while (j.le.numn(IN))
        do while (j.le.numn(IN))
          ThisJ = j
!          ELIVE = .true.      ! element still valid

!      begin new element with line from node IN to its Jth neighbour
          ListTr(1) = IN        ! first node of fresh element is node IN
          ListTr(2) = NL(j,IN)
!         kill element if lower-numbered node is encountered      
          if(ListTr(2).lt.IN) then
            j = j + 1
            cycle !go to 60
          endif

          numverts = 2

!        Find  next neighbour of IN in ccw direction for later check
!                                                 of polygon closure

          call CwwNbr(IN,ThisJ,NextCcwNbr,NextJ,numn(in),nbtot,itot,nl)

!         Find consecutive vertices of polygon in ccw direction until
!         polygon closes or terminates for some reason
!         - maybe make NextVer return an error indicator
!         nextvertex = NextVer(IN,NB(IN,j))
          nextvertex = NextVer(IN,NL(j,IN),numn(NL(j,IN)),nl,nbtot,itot)
          ListTr(3) = nextvertex

!         kill element if lower-numbered node is encountered      
          if(ListTr(3).le.IN) then
            j = j + 1
            cycle !go to 60
          endif
          if (nextvertex.eq.NextCcwNbr) then
!         element is triangle
            ListTr(4) = 0
            j = j + 1
            cycle !go to 60
          else
!           ELEMENT IS POLYGON with more than 3 sides or is not closed
            nextvertex = NextVer(NL(j,IN),ListTr(3), &
                                numn(ListTr(3)),nl,nbtot,itot)
!           kill element if lower-numbered node is encountered      
            if(nextvertex.lt.IN) then
              j = j + 1
              cycle !go to 60
            endif
            ListTr(4) = nextvertex

            if (nextvertex.eq.NextCcwNbr) then
!           element is quadrilateral
              j = j + 1
              cycle !go to 60
            else
              if(retrowanted) then
                followboundary = .true.
!          element has 5 or more sides - assume it is a boundary, not an element
                NUMofBDYS = NUMofBDYS + 1
!      store these 4 nodes top down in array RNODE( ) - note that outer boundary
!      will thus be in clockwise order and islands in ccw order. Dirns. will have 
!      to be reversed on output
                do m = 1,4
                  rnodes(Nfound+m) = ListTr(m)
                end do
                Nfound = Nfound + 4
                NUMperBDY(NUMofBDYS) = 4
                firstnode = ListTr(1)

!           look for 5th, 6th, ... nodes on this boundary
                do while(followboundary)
                  nextvertex = NextVer(rnodes(Nfound-1),rnodes(Nfound), &
     &                        numn(rnodes(Nfound)),nl,nbtot,itot)
!             kill element if lower-numbered node is encountered      
                  if(nextvertex.lt.IN) then
                    Nfound = Nfound - NUMperBDY(NUMofBDYS)
                    NUMofBDYS = NUMofBDYS - 1
!                    elcount = elcount - 1
                    j = j + 1
                    cycle !go to 60
                  endif
                  if(nextvertex.ne.firstnode) then
!               boundary is not yet closed, add node to list of boundary nodes
                    Nfound = Nfound + 1
                    rnodes(Nfound) = nextvertex
                    NUMperBDY(NUMofBDYS) = NUMperBDY(NUMofBDYS) + 1

                  else
!             boundary now closed
                    followboundary = .false.               
                  endif
                enddo
              endif
!              elcount = elcount - 1
              j = j + 1
              cycle !go to 60
            endif
          endif

!60      enddo
        enddo

      enddo  !itot loop

!      rnodes( ) holds the various boundaries, each in reverse of required order.
!     
!            First boundary has NUMperBDY(1) nodes, which are in
!              rnodes(1),...,rnodes(NUMperBDY(1))
!            Second boundary has NUMperBDY(2) nodes, which are in
!              rnodes(NUMperBDY(1)+1),....,rnodes(NUMperBDY(1)+(NUMperBDY(2))
!            and so on.

!      First step is to determine which of the NUMofBDYS boundaries is the outer.
!      Any of the extreme values must lie on outer boundary
!        - find minimum x-coordinate over all boundary nodes
        minimumx = 1.E+20
        do k = 1, Nfound
          if(minimumx.lt.dxray(rnodes(k))) minimumx = dxray(rnodes(k))
        end do

!      now find which boundary has the node with minimum x
        sumofnodes = 0
        do k = 1, NUMofBDYS
          do m = sumofnodes+1,sumofnodes+ NUMperBDY(k)
            if(dxray(rnodes(m)).le.minimumx) then
              numouter = k
!           the numouter-th boundary is the outer one
              presum = sumofnodes
!           presum is number of nodes in rnodes( ) prior to start of outer boundaary
!           outer boundary (in ccw order) consists of nodes
!           rnodes(presum+NUMperBNDY(numouter)) to rnodes(presum+1) 
              go to 444
            endif
          end do
          sumofnodes = sumofnodes + NUMperBDY(k)
        end do
444     continue

        sumofnodes = 0
        nbnd = 0
        do k = 1, NUMofBDYS
          nbnd = nbnd + 1
          bnode_index(nbnd) = sumofnodes + 1
          bnode_id(nbnd) = ncode(rnodes(sumofnodes + 1))
          do m = sumofnodes+1,sumofnodes+ NUMperBDY(k)
            if(ncode(rnodes(m)).ne.bnode_id(nbnd)) then !new segment
              nbnd = nbnd + 1
              bnode_index(nbnd) = m
              bnode_id(nbnd) = ncode(rnodes(m))
            endif
          end do
          sumofnodes = sumofnodes + NUMperBDY(k)
        end do
        nbp = sumofnodes
        bnode_index(nbnd+1) = sumofnodes + 1
        

! *** test boundary
      open(41,file='btestout.dat',status='unknown')
      write(41,*) sumofnodes
      write(41,*) NUMofBDYS
      m = 0
      do k=1,NUMofBDYS
        write(41,*) NUMperBDY(k)
        do j=1,NUMperBDY(k)
          m = m + 1
          write(41,*) dxray(rnodes(m)),dyray(rnodes(m)),ncode(rnodes(m))
        enddo
      enddo
      m = 0
      write(41,*) m
      close(41)

      RETURN
      END

!************************************************************************
      
      SUBROUTINE CwwNbr(INODE,ThisJ,NextCcwNbr,NextJ,numin,nbm,nm,nl)

!     Purpose: Given position ThisJ in neighbour list of node INODE
!              routine returns the next neighbour of INODE counter-
!              clockwise, i.e. NextCcwNbr, and its position, NextJ, in
!              neighbour list of INODE, allowing for cyclic nature of 
!              list and assuming neighbours of all nodes already ccw

      IMPLICIT NONE

! *** Passed variables ***
      INTEGER INODE, ThisJ, numin, nbm, nm
      integer nl(nbm,nm+1)

! *** Local variables
      INTEGER NextJ, NextCcwNbr

! ----BEGIN------------------------------------------------------------

      if (ThisJ.lt.numin) then !numn(INODE))then
        NextJ = ThisJ + 1
      else
        NextJ = 1
      endif

      NextCcwNbr = NL(NextJ,INODE)

      RETURN
      END
      
!************************************************************************
      
      INTEGER FUNCTION NextVer(NODE, NBR, numnbr, nl, nbm, nm)

!     Purpose: Given some neighbour, with index NBR, of node NODE -
!              routine returns the next neighbour of NBR clockwise 
!              from NODE (generally not same as next neighbour of NODE 
!              ccw from NBR).
!              Assumes neighbours of all nodes already sorted ccw
! 
!     Use: When finding the vertices of a polygon in ccw direction,
!          next node is NextVer (thisnode,lastnode)

      IMPLICIT NONE

! *** Passed variables ***
      INTEGER NODE, NBR, numnbr, nbm,nm
      integer nl(nbm,nm+1)

! *** Local variables
      INTEGER K

! ----BEGIN------------------------------------------------------------

!     Find where NODE sits in neighbour list of NBR, and take preceding
!                                  neighbour in clockwise direction
      NextVer = 0  ! kill compiler warning
      do k = 1, numnbr  !NUMN(NBR)
        if (NODE.eq.NL(k,NBR)) then
!         thus NODE is kth neighbour of NBR
          if (k.ne.1) then 
            NextVer = NL(k-1,NBR)
          else
            NextVer = NL(numnbr,nbr)  !NUMN(NBR),NBR)
          endif
        endif
      enddo

      RETURN
      END
      
!************************************************************************
      
      SUBROUTINE SORTNBRS(INODE,NBRS,SNB,nm,nbtot,nl,x,y)

!       Purpose: To sort neighbours of current node INODE into
!                counterclockwise order

      IMPLICIT NONE

! *** PASSED VARIABLES

      integer nm,nbtot,nl(nbtot,nm+1)
      integer INODE, NBRS, SNB(NBTOT+1)
      real*8 x(nm),y(nm)
!       INODE - index of current node
!       NBRS  - number of neighbours of current node
!       SNB( ) - sorted neighbours & SNB(NBRS+1) = SNB(1)

! *** LOCAL VARIABLES

      REAL*8 XDIFF,YDIFF,X1,Y1,X2,Y2
      integer I,J,N, TNb1, TNb2
      integer OCT,UNB(NBTOT),NUinOCT(8),NBinOCT(8,nbtot)
!       OCT - octant counter
!       UNB( ) - packed, unsorted neighbours  
!       NUinOCT(I) - No. of neighbours found in octant I
!       NBinOCT(I, ) - Indices of neighbours found in octant I
!       TNb1( ) - temporary first nbr in sort within octant
!       TNb2( ) - temporary second nbr in sort within octant

! ----BEGIN------------------------------------------------------------
!     Find total no. of nbrs of INODE and form packed unsorted list
      NBRS = 0
      unb = 0
      snb = 0
      DO I = 1, NBTOT
        IF(NL(I,INODE).GT.0) THEN
          NBRS = NBRS + 1
          UNB(NBRS) = NL(I,INODE)
        ENDIF
      enddo

!     Set counts of neighbours in octants to zero
!      DO I = 1,8
        NUinOCT = 0
        NBinOCT = 0
!      enddo

!     Sort nbrs into octants, update counts, store nbrs by octant
      DO I = 1, NBRS
        XDIFF = X(UNB(I)) - X(INODE)
        YDIFF = Y(UNB(I)) - Y(INODE)
        CALL FINDOCT(OCT,XDIFF,YDIFF) !INODE,UNB(I))
        NUinOCT(OCT) = NUinOCT(OCT) + 1
        NBinOCT(OCT,NUinOCT(OCT)) = UNB(I)
      enddo

!     If more than one nbr in octant, do ccw sort within octant

      N = 0

!         DO 30 OCT = 1,8
      OCT = 1
      DO 30 WHILE ((OCT .LE. 8) .AND. (N.LT.NBTOT)) 

        IF ( NUinOCT(OCT) .EQ. 1 ) THEN
          N = N + 1
          SNB(N) = NBinOCT(OCT,1)
        ELSEIF (NUinOCT(OCT) .gt. 1 ) THEN
          DO 20 I = 1, NUinOCT(OCT)          
            DO 20 J = 1, NUinOCT(OCT) - I     
              TNb1 = NBinOCT(OCT,J)
              TNb2 = NBinOCT(OCT,J+1)
              X1 = X(TNb1) - X(INODE)
              X2 = X(TNb2) - X(INODE)
              Y1 = Y(TNb1) - Y(INODE)
              Y2 = Y(TNb2) - Y(INODE)
!               Check cw v.ccw order by vector cross-product
              IF(X1*Y2-X2*Y1.LT.0)THEN 
                NBinOCT(OCT,J) = TNb2
                NBinOCT(OCT,J+1) = TNb1
              ENDIF
20        CONTINUE
          DO 25 I = 1, NUinOCT(OCT)
            N = N + 1
            IF (N.LE.NBTOT) THEN
              SNb(N) = NBinOCT(OCT,I)
            ENDIF
25        CONTINUE

        ENDIF

!         Following reset is necessary
        NUinOCT(OCT)=0
        OCT = OCT +1
30    CONTINUE

      SNb(N+1) = SNb(1)
!       
      RETURN
      END

! *********************************************************************
      
      SUBROUTINE FINDOCT(OCT, X, Y ) !P1, P2)

!     Purpose: Finds which octant (45 degree sector) node P2 lies in
!               relative to node P1, by first sorting into quadrants,
!               then rotating axes 45 degrees and sorting into new
!               quadrants. Octants numbered counterclockwise from ENE

      IMPLICIT NONE

! *** PASSED VARIABLES ***
      integer   OCT !, P1, P2
      REAL*8    X, Y

! *** LOCAL VARIABLES ***   
      LOGICAL Up, Rite, UpRite, LeftUP

! ----BEGIN------------------------------------------------------------

      Up = Y .GE. 0
      Rite = X .GE. 0
      UpRite = Y+X .GE. 0
      LeftUp = Y-X .GE. 0

      IF ( Up ) THEN
        IF ( Rite ) THEN
          IF (.NOT. LeftUp ) THEN
            OCT = 1
          ELSE
            OCT = 2
          ENDIF
        ELSE 
          IF ( UpRite ) THEN
            OCT = 3
          ELSE
            OCT = 4
          ENDIF
        ENDIF
      ELSE
        IF (.NOT. Rite ) THEN
          IF ( LeftUp ) THEN
            OCT = 5
          ELSE
            OCT = 6
          ENDIF
        ELSE
          IF (.NOT. UpRite ) THEN
            OCT = 7
          ELSE
            OCT = 8
          ENDIF
        ENDIF
      ENDIF

      end

!********************************************************************************
!************************************************************************
