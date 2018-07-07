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
      integer, allocatable :: connect(:), offset(:)
      integer(1), allocatable :: CellType(:)
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

!***************************************************************

      Program RCM2Tecplot
  
      use RCMArrays
!      use Lib_VTK_IO
  
      implicit none
      
!      include 'tecio.f90'

      integer :: i,j,k,istat,npvgrd,npdim,icount
      integer :: numarg, iargc
      character*256 :: fnamedata='', fnameprofile='', fnamemax=''
      character(18) :: AnalysisTime
      character(len=256) :: arg
      logical :: notOK=.true.
       
!     read time slices from the binary data file
      numarg=iargc()

      if(numarg.eq.0) then
        do while (notOK.eqv..true.)
          write(*,*) ' Enter output option: TECASC (ascii) or TECPLT (plt).'
          read(*,'(a)') outopt
          if(outopt(1:6).eq.'TECASC') then
            outfileopt = 0
            notOK = .false.
          elseif(outopt(1:6).eq.'TECPLT') then
            outfileopt = 1
            notOK = .false.
          else
            write(*,*) ' Unknown option: try again'
          endif
        enddo
        write(*,*) ' Enter filename for input RiCOM binary file'
        read(*,'(a)') fnamedata
        write(*,*) ' Enter filename for profile points input. CR=none.'
        read(*,'(a)') fnameprofile
      endif
      
      if(numarg.ge.1) then
        i = 1
        call getarg(i,arg)
        outopt = arg
        do while (notOK.eqv..true.)
          if(outopt(1:6).eq.'VTKASC') then
            outfileopt = 0
            notOK = .false.
          elseif(outopt(1:6).eq.'VTKBIN') then
            outfileopt = 1
            notOK = .false.
          else
            write(*,*) ' Unknown option: try again'
            write(*,*) ' Enter output option: VTKASC (ascii) or VTKBIN (binary).'
            read(*,'(a)') outopt
          endif
        enddo
      endif
      
      if(numarg.ge.2) then
        i = 2
        call getarg(i,arg)
        fnamedata = arg
      endif
      
      if(numarg.ge.3) then
        i = 3
        call getarg(i,arg)
        fnameprofile = arg
      endif

      open(unit=20,file=fnamedata,status='old',form='unformatted')
      if(fnameprofile.ne.'') open(unit=21,file=fnameprofile,status='old')

!      if(outfileopt.eq.0) then
!        OutResFile = trim(fnamedata)//'.vtu'
!      elseif(outfileopt.eq.1) then
!        OutResFile = trim(fnamedata)//'.vtu'
!      endif
      OutResFile = trim(fnamedata)

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
      npdim = (np+nsides)*(npvC+1)
        
      ALLOCATE ( xp(np), yp(np), zp(np),&
          nen(ne,ncn),connect(2*ne*ncn*npvc),offset(ne*npvc),CellType(ne*npvc), &
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

! *** average (u,v) to centroid
!        do j=1,ne
!          do jj=1,ncn
!            ii = numsideT(jj,j)
!            do kv=1,npv
!              ubar(j,kv) = ubar(j,kv) + un(kv,ii)*sdy(ii) + ut(kv,ii)*sdx(ii)
!              vbar(j,kv) = vbar(j,kv) - un(kv,ii)*sdx(ii) + ut(kv,ii)*sdy(ii)
!            enddo
!          enddo
!        enddo
!        ubar = ubar/float(ncn)
!        vbar = vbar/float(ncn)

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
      
! *** write tecplot file
        if(noptcount.eq.0) write(*,*) ' call vtkout, nopt=',nopt
        call OutputData

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

! *** write tecplot file
        if(outfileopt.eq.0) then
          OutMaxFile = trim(fnamedata)//'.max.vtu'
        elseif(outfileopt.eq.1) then
          OutMaxFile = trim(fnamedata)//'.max.vtu'
        endif
        noptcount = 0
        nopt = 4
        write(*,*) ' call vtkout, nopt=',nopt
        call OutputData
       
      endif

! *** clean up
      close(20)

      stop
      end

!***************************************************************

      subroutine OutputData () !(noptcount, noptwrite)

      USE RCMArrays
      use Lib_VTK_IO

      implicit none
      
!      include 'tecio.f90'
      
      integer, parameter :: ndf0=0
      integer :: n23, npv0, MBout, icount
      integer, save :: npvb=-1
      integer :: i,j,js,k,kv,nn,ncn2,nentmp,mr,npvm,nps,nen1(4,1),nen3(8,1)
      integer, save :: ffmt=0, ftype=0, idbg=1, isdbl=1, i0=0, i1=1,i8=8
      integer :: ZoneType, PVLst(10), VarLoc(10), ShVar(10),ShCon,nptot,netot
      real*8 :: cdep, zz(100),topomin,deptest,zero,zncn
      real*8 :: uu,vv
      real*8 :: bigrI,bigrcI
      character(10) cseq
      character(80) :: vars
      character(4), save :: fseq='0000'


      SELECT CASE (nopt)

      CASE(2)  !eta,C on elements (2d), u,v on edges(3d). Model variable locations.

! *** generate ut for output.
        if(npvb.lt.0) npvb = max(npv,min(2,2*max(0,izcoord-9)))
        Rn0(:,1:nsides) = un(:,1:nsides)
        call  GlobalQInterp ()
        if(allocated(etaMB)) then
          MBout = 1
        else
          MBout = 0
        endif

! *** initialize output files
        if (noptcount.eq.0) then   ! first time through

! *** generate connectivity vector
          icount = 0
          do j=1,ne
            if(ncn.eq.3) then
              do k=1,ncn
                connect(icount+k) = nen(j,k)-1
              enddo
              if(j.eq.1) then
                offset(1) = 3
              else
                offset(j) = offset(j-1) + 3
              endif
              CellType(j) = 5
              icount = icount + 3
            elseif(ncn.eq.4) then
              if(NEN(j,ncn).eq.0) then
                do k=1,ncn
                  connect(icount+k) = nen(j,k)-1
                enddo
                if(j.eq.1) then
                  offset(1) = 3
                else
                  offset(j) = offset(j-1) + 3
                endif
                CellType(j) = 5
                icount = icount + 3
              else
                do k=1,ncn
                  connect(icount+k) = nen(j,k)-1
                enddo
                if(j.eq.1) then
                  offset(1) = 4
                else
                  offset(j) = offset(j-1) + 4
                endif
                CellType(j) = 9
                icount = icount + 4
              endif
            endif
          enddo

!        else   ! *** sequential files
!          return
        endif

        write(fseq,'(I4.4)') noptcount
        write(*,*) noptcount, fseq
        
        if(outfileopt.eq.0) then
          i = VTK_INI_XML( output_format 	= 'ASCII', &
                           filename	=  trim(OutResFile)//fseq//'.vtu', &
                           mesh_topology = 'UnstructuredGrid')
        elseif(outfileopt.eq.1) then
          i = VTK_INI_XML( output_format 	= 'BINARY', &
                           filename	=  trim(OutResFile)//fseq//'.vtu', &
                           mesh_topology = 'UnstructuredGrid')
        endif
          
        nptot = np + nsides
        netot = ne

! *** write the data
        xpt(1:np) = xp(1:np)
        ypt(1:np) = yp(1:np)
        zpt(1:np) = zp(1:np)
        do j=1,nsides
          xpt(np+j) = sxy(1,j)
          ypt(np+j) = sxy(2,j)
          zpt(np+j) = sbot(j)
        enddo

        i = VTK_GEO_XML(NPtot,NE,Xpt,Ypt,Zpt)
        
        i = VTK_CON_XML(NE,CONNECT,OFFSET,CellType)

        zpt = 0.D0
        do j=1,np
          uct(j) = uc(1,j)
          vct(j) = vc(1,j)
        enddo
        do j=1,nsides
          uct(np+j) =  un(1,j)*sdy(j)+ut(1,j)*sdx(j)
          vct(np+j) = -un(1,j)*sdx(j)+ut(1,j)*sdy(j)
        enddo
        i = VTK_DAT_XML('node','OPEN')
        i = VTK_VAR_XML(NPtot,'Vec',uct,vct,zpt)
        if(MBout.gt.0) then
          do j=1,np
            uct(j) = uc(2,j)
            vct(j) = vc(2,j)
          enddo
          do j=1,nsides
            uct(np+j) = un(2,j)*sdy(j)+ut(2,j)*sdx(j)
            vct(np+j) = -un(2,j)*sdx(j)+ut(2,j)*sdy(j)
          enddo
          i = VTK_VAR_XML(NPtot,'VecMB',uct,vct,zpt)
        endif
        i = VTK_DAT_XML('node','CLOSE')

        i = VTK_DAT_XML('cell','OPEN')
        i = VTK_VAR_XML(ne,'eta',eta)
        
        if(nson.gt.0) then
          zpt(1:ne) = wz(1,1:ne)
          i = VTK_VAR_XML(ne,'wz',zpt)
          zpt(1:ne) = qp(1,1:ne)
          i = VTK_VAR_XML(ne,'qp',zpt)
        elseif(nsed.gt.0) then
          zpt(1:ne) = cc(1,1:ne)
          i = VTK_VAR_XML(ne,'cc',zpt)
        elseif(nsol.ne.0) then
          if(iOPsol.eq.0) then
            zpt(1:ne) = sigt(1,1:ne)
            i = VTK_VAR_XML(ne,'sigt',zpt)
          elseif(iOPsol.eq.1) then
            zpt(1:ne) = sigt(1,1:ne)
            i = VTK_VAR_XML(ne,'sigt',zpt)
            zpt(1:ne) = PSU(1,1:ne)
            i = VTK_VAR(ne,'PSU',zpt)
          elseif(iOPsol.eq.2) then
            zpt(1:ne) = sigt(1,1:ne)
            i = VTK_VAR_XML(ne,'sigt',zpt)
            zpt(1:ne) = TC(1,1:ne)
            i = VTK_VAR_XML(ne,'TC',zpt)
          elseif(iOPsol.eq.3) then
            zpt(1:ne) = sigt(1,1:ne)
            i = VTK_VAR_XML(ne,'sigt',zpt)
            zpt(1:ne) = PSU(1,1:ne)
            i = VTK_VAR_XML(ne,'PSU',zpt)
            zpt(1:ne) = TC(1,1:ne)
            i = VTK_VAR_XML(ne,'TC',zpt)
          endif
        elseif(MBout.gt.0) then
          i = VTK_VAR_XML(ne,'etaMB',etaMB)
!          do j=1,np
!            uct(j) = uc(2,j)
!            vct(j) = vc(2,j)
!          enddo
!          do j=1,nsides
!            uct(np+j) = un(2,j)*sdy(j)+ut(2,j)*sdx(j)
!            vct(np+j) = -un(2,j)*sdx(j)+ut(2,j)*sdy(j)
!          enddo
!          zpt = 0.D0
!          i = VTK_DAT(nptot,'node')
!          i = VTK_VAR('vect',nptot,'VecMB',uct,vct,zpt)
        elseif(neqtide.gt.0) then
          i = VTK_VAR_XML(ne,'eqt',eqtide)
        endif
        i = VTK_DAT_XML('cell','CLOSE')
            
        
      CASE(3)  !eta, u, v interpolated to centroid. 2D.

! *** generate ut for output.

        if(npvb.lt.0) npvb = max(npv,min(2,2*max(0,izcoord-9)))
        Rn0(:,1:nsides) = un(:,1:nsides)
        call  GlobalQInterp ()


! *** initialize output files
        if (noptcount.eq.0) then   ! first time through

! *** generate connectivity vector
          icount = 0
          nps = np + nsides
          do j=1,ne
            if(ncn.eq.3) then
              do kv=1,npvC
                do k=1,ncn
                  connect(icount+k) = nen(j,k)-1+kv*nps
                  connect(icount+3+k) = nen(j,k)-1+(kv-1)*nps
                enddo
                if(j.eq.1.and.kv.eq.1) then
                  offset(1) = 6
                else
                  offset(npvc*(j-1)+kv) = offset(npvc*(j-1)+kv-1) + 6
                endif
                CellType(npvc*(j-1)+kv) = 13
                icount = icount + 6
              enddo
            elseif(ncn.eq.4) then
              if(NEN(j,ncn).eq.0) then
                do kv=1,npvC
                  do k=1,ncn
                    connect(icount+k) = nen(j,k)-1+kv*nps
                    connect(icount+3+k) = nen(j,k)-1+(kv-1)*nps
                  enddo
                  if(j.eq.1.and.kv.eq.1) then
                    offset(1) = 6
                  else
                    offset(npvc*(j-1)+kv) = offset(npvc*(j-1)+kv-1) + 6
                  endif
                  CellType(npvc*(j-1)+kv) = 13
                  icount = icount + 6
                enddo
              else
                do kv=1,npvC
                  do k=1,ncn
                    connect(icount+k) = nen(j,k)-1+kv*nps
                    connect(icount+4+k) = nen(j,k)-1+(kv-1)*nps
                  enddo
                  if(j.eq.1.and.kv.eq.1) then
                    offset(1) = 8
                  else
                    offset(npvc*(j-1)+kv) = offset(npvc*(j-1)+kv-1) + 8
                  endif
                  CellType(npvc*(j-1)+kv) = 12
                  icount = icount + 8
                enddo
              endif
            endif
          enddo

        endif

        write(fseq,'(I4.4)') noptcount
        write(*,*) noptcount, fseq
 
        if(outfileopt.eq.0) then
          i = VTK_INI_XML( output_format 	= 'ASCII', &
                           filename	=  trim(OutResFile)//fseq//'.vtu', &
                           mesh_topology = 'UnstructuredGrid')
        elseif(outfileopt.eq.1) then
          i = VTK_INI_XML( output_format 	= 'BINARY', &
                           filename	=  trim(OutResFile)//fseq//'.vtu', &
                           mesh_topology = 'UnstructuredGrid')
        endif
         

        nptot = (np+nsides)*(npvC+1)
        netot = ne*npvc

! *** write the data

        nps = np + nsides
        do k=1,npvc+1
          do j=1,np
              xpt((k-1)*nps+j) = xp(j)
              ypt((k-1)*nps+j) = yp(j)
          enddo
          do j=1,nsides
              xpt(np+(k-1)*nps+j) = sxy(1,j)
              ypt(np+(k-1)*nps+j) = sxy(2,j)
          enddo
        enddo

        if(izcoord.eq.0.or.izcoord.eq.2) then
          do k=1,npvc+1
            do j=1,np
              zpt((k-1)*nps+j) = -zp(j)*zdep(k)
            enddo
            do j=1,nsides
              zpt(np+(k-1)*nps+j) = -sbot(j)*zdep(k)
            enddo
          enddo
        elseif(izcoord.eq.1.or.izcoord.eq.3) then
          do k=1,izgrid-1
            do j=1,np
              zpt((k-1)*nps+j) = -max(zp(j),zdep(izgrid))*zdep(k)
            enddo
            do j=1,nsides
              zpt(np+(k-1)*nps+j) = -max(sbot(j),zdep(izgrid))*zdep(k)
            enddo
          enddo
          do k=izgrid,npvc+1
            do j=1,np
              zpt((k-1)*nps+j) = max(zp(j),zdep(k))
            enddo
            do j=1,nsides
              zpt(np+(k-1)*nps+j) = max(sbot(j),zdep(k))
            enddo
          enddo
        endif

        i = VTK_GEO_XML(NPtot,NEtot,Xpt,Ypt,Zpt)
        
        i = VTK_CON_XML(NEtot,CONNECT,OFFSET,CellType)

!  *** generate velocity vector
        uct = 0.D0
        vct = 0.D0
        do k=1,npv
          do j=1,np
            uct((k-1)*nps+j) = uc(k,j)
            vct((k-1)*nps+j) = vc(k,j)
          enddo
          do j=1,nsides
            uct(np+(k-1)*nps+j) =  un(k,j)*sdy(j)+ut(k,j)*sdx(j)
            vct(np+(k-1)*nps+j) = -un(k,j)*sdx(j)+ut(k,j)*sdy(j)
          enddo
        enddo
! *** average vertical velocity to nodes
        zpt = 0.D0

        i = VTK_DAT_XML('node','OPEN')
        i = VTK_VAR_XML(NPtot,'Vec',uct,vct,zpt)
        i = VTK_DAT_XML('node','CLOSE')

        i = VTK_DAT_XML('cell','OPEN')
!          do j=1,ne
!            do k=1,npvc
!              zpt(npvc*(j-1)+k) = wz(k,j)
!            enddo
!          enddo
!        i = VTK_VAR_XML(ne,'wz',zpt)
        
        if(nson.gt.0) then
          do j=1,ne
            do k=1,npvc
              zpt(npvc*(j-1)+k) = qp(k,j)
            enddo
          enddo
          i = VTK_VAR_XML(netot,'qp',zpt)
        elseif(nsed.gt.0) then
          do j=1,ne
            do k=1,npvc
              zpt(npvc*(j-1)+k) = cc(k,j)
            enddo
          enddo
          i = VTK_VAR_XML(netot,'cc',zpt)
        elseif(nsol.ne.0) then
          if(iOPsol.eq.0) then
            do j=1,ne
              do k=1,npvc
                zpt(npvc*(j-1)+k) = sigt(k,j)
              enddo
            enddo
            i = VTK_VAR_XML(netot,'sigt',zpt)
          elseif(iOPsol.eq.1) then
            do j=1,ne
              do k=1,npvc
                zpt(npvc*(j-1)+k) = sigt(k,j)
              enddo
            enddo
            i = VTK_VAR_XML(netot,'sigt',zpt)
            do j=1,ne
              do k=1,npvc
                zpt(npvc*(j-1)+k) = PSU(k,j)
              enddo
            enddo
            i = VTK_VAR(netot,'PSU',zpt)
          elseif(iOPsol.eq.2) then
            do j=1,ne
              do k=1,npvc
                zpt(npvc*(j-1)+k) = sigt(k,j)
              enddo
            enddo
            i = VTK_VAR_XML(netot,'sigt',zpt)
            do j=1,ne
              do k=1,npvc
                zpt(npvc*(j-1)+k) = TC(k,j)
              enddo
            enddo
            i = VTK_VAR_XML(netot,'TC',zpt)
          elseif(iOPsol.eq.3) then
            do j=1,ne
              do k=1,npvc
                zpt(npvc*(j-1)+k) = sigt(k,j)
              enddo
            enddo
            i = VTK_VAR_XML(netot,'sigt',zpt)
            do j=1,ne
              do k=1,npvc
                zpt(npvc*(j-1)+k) = PSU(k,j)
              enddo
            enddo
            i = VTK_VAR_XML(netot,'PSU',zpt)
            do j=1,ne
              do k=1,npvc
                zpt(npvc*(j-1)+k) = TC(k,j)
              enddo
            enddo
            i = VTK_VAR_XML(netot,'TC',zpt)
          endif
        endif
        i = VTK_DAT_XML('cell','CLOSE')

      CASE(4)  !case 4 is max values

! *** then write it out

          if (noptcount.eq.0) then

! *** generate connectivity vector
            icount = 0
            do j=1,ne
              if(ncn.eq.3) then
                do k=1,ncn
                  connect(icount+k) = nen(j,k)-1
                enddo
                if(j.eq.1) then
                  offset(1) = 3
                else
                  offset(j) = offset(j-1) + 3
                endif
                CellType(j) = 5
                icount = icount + 3
              elseif(ncn.eq.4) then
                if(NEN(j,ncn).eq.0) then
                  do k=1,ncn
                    connect(icount+k) = nen(j,k)-1
                  enddo
                  if(j.eq.1) then
                    offset(1) = 3
                  else
                    offset(j) = offset(j-1) + 3
                  endif
                  CellType(j) = 5
                  icount = icount + 3
                else
                  do k=1,ncn
                    connect(icount+k) = nen(j,k)-1
                  enddo
                  if(j.eq.1) then
                    offset(1) = 4
                  else
                    offset(j) = offset(j-1) + 4
                  endif
                  CellType(j) = 9
                  icount = icount + 4
                endif
              endif
            enddo
          endif

          if(outfileopt.eq.0) then
            i = VTK_INI_XML( output_format 	= 'ASCII', &
                             filename	=  trim(OutMaxFile), &
                             mesh_topology = 'UnstructuredGrid')
          elseif(outfileopt.eq.1) then
            i = VTK_INI_XML( output_format 	= 'BINARY', &
                             filename	=  trim(OutMaxFile), &
                             mesh_topology = 'UnstructuredGrid')
          endif

          !set original coordinates if icoord>0
          if(ixycoord.eq.0) then
!            x00 = -long0 ! 0. !offsets: from input file?
!            y00 = -lat0  ! 0.
            bigrI = rad2deg/bigr
            bigrcI = bigrI/cos(deg2rad*y0off) 
            do j=1,np
              rhv(j) = xp(j)*bigrcI 
              gamma(j) = yp(j)*bigrI
            enddo
            i = VTK_GEO_XML(NP,NE,rhv,gamma,Zp)
          else
            i = VTK_GEO_XML(NP,NE,Xp,Yp,Zp)
          endif

          i = VTK_CON_XML(NE,CONNECT,OFFSET,CellType)

          Rn0 = 0.
          do j=1,ne
            ncn2 = ncn -1 + min0(1,nen(j,ncn))
            zncn = dble(ncn2)
            topomin = min(sbot(numsideT(1,j)),sbot(numsideT(2,j)),&
                          sbot(numsideT(3,j)),sbot(numsideT(ncn2,j)))
            deptest = emax(j) - topomin
            if(deptest.gt.depmin) then
              uu = (umax(numsideT(1,j))+umax(numsideT(2,j))+umax(numsideT(3,j)))/3.
              vv = (vmax(numsideT(1,j))+vmax(numsideT(2,j))+vmax(numsideT(3,j)))/3.
              Rhv(j)= sqrt(uu**2 + vv**2)
            endif
          enddo

          i = VTK_DAT_XML('cell','OPEN')
          i = VTK_VAR_XML(ne,'IE',IEcode)
          i = VTK_VAR_XML(ne,'emax',emax)
          i = VTK_VAR_XML(ne,'emin',emin)
          i = VTK_VAR_XML(ne,'smax',Rhv)
          i = VTK_VAR_XML(ne,'t1',t1)
          i = VTK_VAR_XML(ne,'twet',twet)
          i = VTK_DAT_XML('cell','CLOSE')

      END SELECT
      
      i = VTK_GEO_XML()
      i = VTK_END_XML()

        close(22)

        if(npv.gt.1) then
          ! first time through start new file
          do j=1,jUprofile

            js = numsideT(2,tsUnodes(j))
            if(js.eq.0) cycle

     
            npv0 = npvc+1
            call Getzgrid(js,npv0,zz)

            if(izcoord.gt.1) then
              do k=1,npv0-1
                zz(k) = 0.5D0*(zz(k)+zz(k+1))
              enddo
              npv0 = npv0-1
            endif

            if (noptcount.eq.0) then
              n23 = 23  !+j-1
              write(cseq,"(I2.2)") j
              open(unit=n23,file='profile'//trim(cseq)//'.dat',status='unknown')
              if(nson.eq.0) then
                write(n23,*) 'VARIABLES="z" "u" "v" "w"'
              else
                write(n23,*) 'VARIABLES="z" "u" "v" "w" "q"'
              endif
              write(n23,*) 'ZONE'
              write(n23,*) ' SOLUTIONTIME=',TET

!             ! write the data
              if(nson.eq.0) then
                do k=1,npv0
                  write(n23,'(4(1PE12.4))') zz(k),( un(k,js)*sdy(js)+ut(k,js)*sdx(js)), &
                                 (-un(k,js)*sdx(js)+ut(k,js)*sdy(js)),wz(k,iside(1,js))
                enddo
              else
                do k=1,npv0
                  write(n23,'(4(1PE12.4))') zz(k),( un(k,js)*sdy(js)+ut(k,js)*sdx(js)), &
                       (-un(k,js)*sdx(js)+ut(k,js)*sdy(js)),wz(k,iside(1,js)),qp(k,iside(1,js))
                enddo
              endif
              close(n23)
            else 
          ! append file
              n23 = 23  !+j-1
              write(cseq,"(I2.2)") j
              open(unit=n23,file='profile'//trim(cseq)//'.dat',status='old',position='append')
              write(n23,*) 'ZONE'
              write(n23,*) ' SOLUTIONTIME=',TET

!             ! write the data
              if(nson.eq.0) then
                do k=1,npv0
                  write(n23,'(4(1PE12.4))') zz(k),( un(k,js)*sdy(js)+ut(k,js)*sdx(js)), &
                                 (-un(k,js)*sdx(js)+ut(k,js)*sdy(js)),wz(k,iside(1,js))
                enddo
              else
                do k=1,npv0
                  write(n23,'(4(1PE12.4))') zz(k),( un(k,js)*sdy(js)+ut(k,js)*sdx(js)), &
                       (-un(k,js)*sdx(js)+ut(k,js)*sdy(js)),wz(k,iside(1,js)),qp(k,iside(1,js))
                enddo
              endif
              close(n23)
            endif
          enddo
        endif

        if((nsed.gt.0).and.(npvc.gt.1)) then
        
          js = jCProfile
          call GetzgridE(js,npv0,cdep,zz)

          if (noptcount.eq.0) then
          ! first time through start new file
            open(unit=23,file='profileSed.dat',status='unknown')
            write(23,*) 'VARIABLES="z" "C" '
            write(23,*) 'ZONE'
            write(23,*) ' SOLUTIONTIME=',TET

            do k=1,npvc
!           ! write the data
              write(23,*) zz(k), cc(k,js)
            enddo
            write(23,*) zz(npvc), cc(npvc+1,js)
            if(nsbc.ge.3) then
              write(23,*) zz(npvc), cc(npvc+2,js)
            endif
          else 
          ! append file
            open(unit=23,file='profileSed.dat',status='old',position='append')
            write(23,*) 'ZONE'
            write(23,*) ' SOLUTIONTIME=',TET

            do k=1,npvc
!           ! write the data
              write(23,*) zz(k), cc(k,js)
            enddo
            write(23,*) zz(npvc), cc(npvc+1,js)
            if(nsbc.ge.3) then
              write(23,*) zz(npvc), cc(npvc+2,js)
            endif
          endif
          close(23)
        endif

        if((nsol.ne.0).and.(npvc.gt.1)) then
          
          do j=1,jSprofile

            js = tsSnodes(j)
            call GetzgridE(js,npv0,cdep,zz)
          
            if (noptcount.eq.0) then
          ! first time through start new file
!            do j=1,jSprofile
              n23 = 23  !+j-1
              write(cseq,"(I2.2)") j
              open(unit=n23,file='profileCon'//trim(cseq)//'.dat',status='unknown')
              if(iOPsol.eq.0) then
                write(23,*)'VARIABLES= "Z"  "Sigt" '
              elseif(iOPsol.eq.1) then
                write(23,*)'VARIABLES= "Z"  "Sigt" "PSU" '
              elseif(iOPsol.eq.2) then
                write(23,*)'VARIABLES= "Z"  "Sigt" "T" '
              elseif(iOPsol.eq.3) then
                write(23,*)'VARIABLES= "Z"  "Sigt" "PSU" "T" '
              endif
!              write(23,*) 'VARIABLES="z" "C" '
              write(23,*) 'ZONE'
              write(23,*) ' SOLUTIONTIME=',TET

!              js = tsSnodes(j)
!              cdep=min(xyz(nen(js,1),3),xyz(nen(js,2),3),xyz(nen(js,3),3),xyz(nen(js,ncn),3))
!              ncn2 = ncn -1 + min0(1,nen(js,ncn))
!              cdep = 0. 
!              do jn=1,ncn2
!                cdep = cdep + xyz(nen(js,jn),3)
!              enddo
!              cdep = eta(js)-cdep/dble(ncn2)
!              write(23,*) 'iOP,npvc=',iOPsol,npvc,js,cdep
!             ! write the data
              if(iOPsol.eq.0) then
                do k=1,npv0
                  write(23,'(6(1x,e14.6))') zz(k), sigt(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.1) then
                do k=1,npv0
                  write(23,'(6(1x,e14.6))') zz(k), sigt(k,js), PSU(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.2) then
                do k=1,npv0
                  write(23,'(6(1x,e14.6))') zz(k), sigt(k,js), TC(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.3) then
                do k=1,npv0
                  write(23,'(6(1x,e14.6))') zz(k), sigt(k,js), PSU(k,js), TC(k,js)
                enddo
                close(n23)
              endif
            else 
          ! append file
!            do j=1,jSprofile
              n23 = 23 !+j-1
              write(cseq,"(I2.2)") j
              open(unit=n23,file='profileCon'//trim(cseq)//'.dat',status='old',position='append')
              write(23,*) 'ZONE'
              write(23,*) ' SOLUTIONTIME=',TET

!              js = tsSnodes(j)
!              ncn2 = ncn -1 + min0(1,nen(js,ncn))
!              cdep = 0. 
!              do jn=1,ncn2
!                cdep = cdep + xyz(nen(js,jn),3)
!              enddo
!              cdep = eta(js)-cdep/dble(ncn2)
!             ! write the data
              if(iOPsol.eq.0) then
                do k=1,npv0
                  write(23,'(6(1x,e14.6))') zz(k), sigt(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.1) then
                do k=1,npv0
                  write(23,'(6(1x,e14.6))') zz(k), sigt(k,js), PSU(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.2) then
                do k=1,npv0
                  write(23,'(6(1x,e14.6))') zz(k), sigt(k,js), TC(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.3) then
                do k=1,npv0
                  write(23,'(6(1x,e14.6))') zz(k), sigt(k,js), PSU(k,js), TC(k,js)
                enddo
                close(n23)
              endif
            endif
          enddo
        endif

      RETURN
      END

!*****************************************************************

      SUBROUTINE GlobalQInterp ()

! *** Interpolate velocity to vertices and midsides

      USE RCMArrays
!      USE AdvectArrays

      implicit none

! *** local variables
      integer :: kv,nn,ncn2,j,mr,mr0,sn0,mr2,sn2
      real*8 :: zncn,dnx1,dny1,dnx2,dny2,un1,un2,uu,vv
      real*8 :: detn,areasum

      uc = 0.
      vc = 0.
      ut = 0.
      rhv = 0.

      do nn=1,ne
        ncn2 = ncn -1 + min0(1,nen(nn,ncn))
        zncn = 1./float(ncn2)
        MR0 = numsideT(ncn2,nn)
        sn0=(iside(2,mr0)-2*nn+iside(1,mr0))/(iside(2,mr0)-iside(1,mr0))
        DNx1 =  sn0*sdy(mr0)
        DNy1 = -sn0*sdx(mr0)
        DO J=1,ncn2
          MR = NEN(NN,J)
          rhv(mr) = rhv(mr) + area(nn)*zncn !area weighting
          MR2 = numsideT(j,nn)
          sn2=(iside(2,mr2)-2*nn+iside(1,mr2))/(iside(2,mr2)-iside(1,mr2))
          DNx2 =  sn2*sdy(mr2)
          DNy2 = -sn2*sdx(mr2)
          detn = dnx1*dny2 - dnx2*dny1
          do kv=1,npv
            un1 = sn0*un(kv,mr0)
            un2 = sn2*un(kv,mr2)
            uu =   (dny2*un1-dny1*un2)*area(nn)/detn
            vv = - (dnx2*un1-dnx1*un2)*area(nn)/detn
            uc(kv,mr) = uc(kv,mr) + uu*zncn
            vc(kv,mr) = vc(kv,mr) + vv*zncn
            ut(kv,mr0) = ut(kv,mr0) - 0.5*sn0*(uu*dny1-vv*dnx1)
            ut(kv,mr2) = ut(kv,mr2) - 0.5*sn2*(uu*dny2-vv*dnx2)
          enddo
          mr0 = mr2
          dnx1 = dnx2
          dny1 = dny2
          sn0 = sn2
        enddo
      enddo

      do j=1,np
        if(rhv(j).gt.0.) then !and.nbc(j).eq.0) then
          uc(:,j)  = uc(:,j)/rhv(j)
          vc(:,j) = vc(:,j)/rhv(j)
        else
          uc(:,j)  = 0.
          vc(:,j) = 0.
        endif
      enddo

       do j=1,nsides
         if(iside(2,j).gt.ne) then !openbc
           ut(:,j)  = 0.
         elseif(iside(2,j).eq.0) then !land bcs         
           areasum = area(iside(1,j))
           ut(:,j)  = ut(:,j)/areasum
         else
           areasum = area(iside(1,j)) + area(iside(2,j))
           ut(:,j)  = ut(:,j)/areasum
         endif
       enddo


      return
      end

!*****************************************************************

      SUBROUTINE Getzgrid(j,npv0,zz)

! *** compute node locations in vertical grid

      USE RCMArrays
      
      implicit none

! *** passed variables
      integer, intent(in) :: j
      integer, intent(out) :: npv0
      real*8, intent(out) :: zz(npvC+1)

! *** local variables
      integer :: k
      real*8 :: zbot, depdif, etaside

      zbot = sbot(j)
      etaside = sdep(j) + zbot
      npv0 = npvC+1
      
      if(izgrid.eq.0) then ! sigma, layers
        do k=1,npv0
          zz(k) = etaside + sdep(j)*zdep(k)
        enddo
      else                 ! sigma on z, layers
        depdif = etaside - max(zdep(izgrid),zbot)
        do k=1,izgrid-1
          zz(k) = etaside + depdif*zdep(k)
        enddo
        zz(izgrid) = etaside - depdif
        npv0 = izgrid
        if(zbot.lt.zdep(izgrid)) then
          do k=izgrid+1,npvc+1
            if(zbot.ge.zdep(k)) then
              zz(k:npvC+1) = zbot
              npv0 = k
              exit
            else
              zz(k) = zdep(k)
            endif
          enddo
        endif
      endif

      return
      end

!******************************************************************************
!******************************************************************************

      SUBROUTINE GetzgridE(j,npv0,zbot,zz)

! *** compute node locations in vertical grid for element centers.

      USE RCMArrays
      
      implicit none

! *** passed variables
      integer, intent(in) :: j
      integer, intent(out) :: npv0
      real*8, intent(out) :: zbot, zz(npvC)

! *** local variables
      integer :: k,ncn2
      real*8 :: cdep, wt

      ncn2 = ncn -1 + min0(1,nen(j,ncn))
      wt = 1.D0/dble(ncn2)
      if(izgrid.eq.0) then ! sigma, layers
        cdep = wt*(zp(nen(j,1)) + zp(nen(j,2)) + zp(nen(j,3)))
        if(ncn2.eq.4) then
          cdep = cdep + wt*zp(nen(j,4))
        endif
        zbot = cdep
        cdep = eta(j) - cdep
        do k=1,npvC
          zz(k) = eta(j) + 0.5D0*cdep*(zdep(k)+zdep(k+1))
        enddo
        npv0 = npvC
      else                 ! sigma on z, layers
        zbot = max(zp(nen(j,1)),zp(nen(j,2)),zp(nen(j,3)),zp(nen(j,ncn2)))
        if(zbot.gt.zdep(izgrid)) then
!          wt = 1.D0/dble(ncn2)
          cdep = 0.D0
          do k=1,ncn2
            cdep = cdep + max(zp(nen(j,k)),zdep(izgrid))
          enddo
          zbot = wt*cdep
          cdep = eta(j) - zbot
        else
          cdep = eta(j) - zdep(izgrid)
        endif
        do k=1,izgrid-2
          zz(k) = eta(j) + 0.5D0*cdep*(zdep(k)+zdep(k+1))
        enddo
        npv0 = izgrid-1
        zz(izgrid-1) = eta(j) + 0.5D0*cdep*(zdep(izgrid-1)-1.D0)
        do k=izgrid,npvC
          !set bottom, npv0
          zz(k) = 0.5D0*(zdep(k)+zdep(k+1))
          npv0 = k
          if(zbot.gt.zz(k)) then
            npv0 = k-1
            zz(k:npvC) = zbot
            exit
          endif
        enddo
      endif
      
      return
      end

!******************************************************************************
!*****************************************************************
