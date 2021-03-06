!***************************************************************

  Module RCMArrays
  
      implicit none

      integer :: ne,np,nsides,npv,npvC,ncn   !nph,nphu,npv  
      integer :: nson,nsed,nsol
      integer :: neqtide=0, neMB=0, nsbc, iOPsol=0
      integer :: noptcount=0, nopt=2
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
      real*8, allocatable ::  eta(:),un(:,:),ut(:,:),wz(:,:),etaMB(:)
      real*8, allocatable ::  uc(:,:),vc(:,:),eqtide(:)
      real*8, allocatable ::  sxy(:,:),sbot(:),sdx(:),sdy(:),sdep(:)
      real*8, allocatable ::  rhv(:),gamma(:),Rn0(:,:) !scratch vectors
      real*8, allocatable ::  qp(:,:),cc(:,:),sigt(:,:),PSU(:,:),TC(:,:)
      real*8, allocatable ::  emax(:),emin(:),spdmax(:),umax(:),vmax(:),t1(:),twet(:)
      character*256 OutResFile, OutMaxFile
 
  end module RCMArrays

!***************************************************************

      Program RCM2Tecplot
  
      use RCMArrays
  
      implicit none

      integer :: i,j,k,istat,npvgrd
      integer :: numarg, iargc
      character*256 :: fnamedata='', fnameprofile='', fnamemax=''
      character(18) :: AnalysisTime
      character(len=256) :: arg
       
!     read time slices from the binary data file
      numarg=iargc()

      if(numarg.eq.0) then
        write(*,*) ' Enter filename for input RiCOM binary file'
        read(*,'(a)') fnamedata
        write(*,*) ' Enter filename for profile points input'
        read(*,'(a)') fnameprofile
      elseif(numarg.ge.1) then
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

      OutResFile = trim(fnamedata)//'.dat'

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
      
      npvgrd = npv + min(1,max(0,izcoord-1))
        
      ALLOCATE ( xp(np), yp(np), zp(np), nen(ne,ncn), &
          sxy(2,nsides),sbot(nsides),sdx(nsides),sdy(nsides),Rn0(npvgrd,nsides), &
          numsideT(ncn,ne),iside(2,nsides),area(ne),IECode(ne),iends(2,nsides), &
          uc(npvgrd,np),vc(npvgrd,np),zdep(npvgrd),rhv(nsides),gamma(nsides),sdep(nsides), &
          eta(ne),un(npvgrd,nsides),ut(npvgrd,nsides),wz(npvgrd,ne), eqtide(ne), STAT = istat )
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
        elseif(tet.lt.tet0) then
          exit
        endif
        tet0 = tet
        
        read(20, IOSTAT=istat) (eta(j),j=1,ne)
! *** read velocity (un,ut)
        if(istat.eq.0) read(20, IOSTAT=istat) ((un(k,i),k=1,npv),i=1,nsides)
        if(istat.eq.0) read(20, IOSTAT=istat) ((ut(k,i),k=1,npv),i=1,nsides)    
        if(istat.eq.0.and.(npv.gt.1.or.nson.gt.0)) read(20, IOSTAT=istat) ((wz(k,i),k=1,npv),i=1,ne)

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
        if(noptcount.eq.0) write(*,*) ' call tecout, nopt=',nopt
        call OutputTecplot

        noptcount = noptcount + 1

      enddo

! *** clean up
!      close(20)

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
        OutResFile = trim(fnamedata)//'.max.dat'
        noptcount = 0
        nopt = 4
        write(*,*) ' call tecout, nopt=',nopt
        call OutputTecplot
        
      endif

! *** clean up
      close(20)

      stop
      end

!***************************************************************

      subroutine OutputTecplot () !(noptcount, noptwrite)

      USE RCMArrays

      implicit none
      
      integer, parameter :: ndf0=0
      integer :: n23, npv0, MBout
      integer, save :: npvb=-1
      integer :: j,js,k,kv,nn,ncn2,nentmp,mr,npvm,nps
      real*8 :: cdep, zz(100),topomin,deptest,zero,zncn
      real*8 :: uu,vv
      real*8 :: bigrI,bigrcI
!      real, allocatable, save :: etamax(:),spdmax(:),t1(:),twet(:),twet2(:),twet3(:),etalast(:)
      character(10) cseq

!      integer, parameter :: ndf0=0
!      integer :: noptcount, noptwrite, istat, n23, npv0,MBout
!      integer :: j,js,k,kv,nn,ncn2,nentmp,mr,npvm,nps
!      integer :: js2,jn
!      integer, save :: nopttype, firstcount = 0
!      real*8 :: cdep, zz(100), etaside,speed,topomin,deptest,zero,zncn
!      real*8 :: uu,vv
!      real*8 :: bigrI,bigrcI,x00,y00
!      real*8, allocatable, save :: emax(:),emin(:),spdmax(:),umax(:),vmax(:),t1(:),twet(:)
!      character(10) cseq

! *** Tecplot output- a single Tecplot file with many frames
        if (noptcount.eq.0) then
          ! first time through start new file
            open(unit=22, file=OutResFile,status='unknown')
        else 
          ! append file
            open(unit=22, file=OutResFile,status='old',position='append')
        endif

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

          if (noptcount.eq.0) then
          ! first time through start new file
            if(nson.gt.0) then
                write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "W" "Q" '
            elseif(nsed.gt.0) then
              write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "S" '
            elseif(nsol.ne.0) then
              if(iOPsol.eq.0) then
                write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "Sigt" '
              elseif(iOPsol.eq.1) then
                write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "Sigt" "PSU" '
              elseif(iOPsol.eq.2) then
                write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "Sigt" "T" '
              elseif(iOPsol.eq.3) then
                write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "Sigt" "PSU" "T" '
              endif
            elseif(MBout.gt.0) then
              write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "EMB" "UMB" "VMB" '
            elseif(neqtide.gt.0) then
              write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "EQT" '
            else
              write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" '
            endif
            
            write(22,"('ZONE N=',i7,' E=',i7 )" ) np+nsides,ne

            if(ncn.eq.3) then
              write(22,"(' ZONETYPE=FETRIANGLE DATAPACKING=BLOCK')" )
            elseif(ncn.eq.4) then
              write(22,"(' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK')" )
            endif
            
            if(nsol.ne.0) then
              if(iOPsol.eq.0) then
                write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
              elseif(iOPsol.eq.1) then
                write(22,"(' VARLOCATION=([4,7-8]=CELLCENTERED)')" )
              elseif(iOPsol.eq.2) then
                write(22,"(' VARLOCATION=([4,7-8]=CELLCENTERED)')" )
              elseif(iOPsol.eq.3) then
                write(22,"(' VARLOCATION=([4,7-9]=CELLCENTERED)')" )
              endif
            elseif(nson.ne.0) then
              write(22,"(' VARLOCATION=([4,7,8]=CELLCENTERED)')" )
            elseif(nsed.gt.0) then
              write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
            elseif(neqtide.gt.0.or.MBout.gt.0) then
              write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
            else
              write(22,"(' VARLOCATION=([4]=CELLCENTERED)')" )
            endif
            
            write(22,*) ' SOLUTIONTIME=',TET

            ! write the data
            write(22,'(6(1x,e14.6))') ((xp(j),j=1,np),(sxy(1,j),j=1,nsides),k=1,1) !x
            write(22,'(6(1x,e14.6))') ((yp(j),j=1,np),(sxy(2,j),j=1,nsides),k=1,1) !y
            write(22,'(6(1x,e14.6))') ((zp(j),j=1,np),(sbot(j),j=1,nsides),k=1,1) !z

            write(22,'(6(1x,e14.6))') (eta(j),j=1,ne)

            zero = 0.0
            write(22,'(6(1x,e14.6))') ((uc(1,j),j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,1)
            write(22,'(6(1x,e14.6))') ((vc(1,j),j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,1)

            if(nson.gt.0) then
              write(22,'(6(1x,e14.6))') (wz(1,j),j=1,ne)
              write(22,'(6(1x,e14.6))') (qp(1,j),j=1,ne)
            elseif(nsed.gt.0) then
              write(22,'(6(1x,e14.6))') (cc(1,j),j=1,ne) !(dsg(j),j=1,ne)
            elseif(nsol.ne.0) then
              if(iOPsol.eq.0) then
                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
              elseif(iOPsol.eq.1) then
                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
                write(22,'(6(1x,e14.6))') (PSU(1,j),j=1,ne)
              elseif(iOPsol.eq.2) then
                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
                write(22,'(6(1x,e14.6))') (TC(1,j),j=1,ne)
              elseif(iOPsol.eq.3) then
                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
                write(22,'(6(1x,e14.6))') (PSU(1,j),j=1,ne)
                write(22,'(6(1x,e14.6))') (TC(1,j),j=1,ne)
              endif
            elseif(MBout.gt.0) then
              write(22,'(6(1x,e14.6))') (etaMB(j),j=1,ne)
              write(22,'(6(1x,e14.6))') (uc(2,j),j=1,np),(( un(2,j)*sdy(j)+ut(2,j)*sdx(j)),j=1,nsides)
              write(22,'(6(1x,e14.6))') (vc(2,j),j=1,np),((-un(2,j)*sdx(j)+ut(2,j)*sdy(j)),j=1,nsides)
            elseif(neqtide.gt.0) then
              write(22,'(6(1x,e14.6))') (eqtide(j),j=1,ne)
            endif

            do j=1,ne
            ! and the elements
              if(NEN(j,ncn).eq.0) then
                NENtmp = NEN(j,ncn-1)
              else
                NENtmp = NEN(j,ncn)
              endif
              write(22,*) (NEN(j,k),k=1,ncn-1),NENtmp
            enddo
          else 
          ! append file
            if(neMB.gt.0) then
              write(22,"('ZONE D=(1,2,FECONNECT)')" )
            else
              write(22,"('ZONE D=(1,2,3,FECONNECT)')" )
            endif

            if(ncn.eq.3) then
              write(22,"(' ZONETYPE=FETRIANGLE DATAPACKING=BLOCK')" )
            elseif(ncn.eq.4) then
              write(22,"(' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK')" )
            endif

            if(nsol.ne.0) then
              if(iOPsol.eq.0) then
                write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
              elseif(iOPsol.eq.1) then
                write(22,"(' VARLOCATION=([4,7-8]=CELLCENTERED)')" )
              elseif(iOPsol.eq.2) then
                write(22,"(' VARLOCATION=([4,7-8]=CELLCENTERED)')" )
              elseif(iOPsol.eq.3) then
                write(22,"(' VARLOCATION=([4,7-9]=CELLCENTERED)')" )
              endif
            elseif(nson.ne.0) then
              write(22,"(' VARLOCATION=([4,7,8]=CELLCENTERED)')" )
            elseif(nsed.gt.0) then
              write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
            elseif(neqtide.gt.0.or.neMB.gt.0) then
              write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
            else
              write(22,"(' VARLOCATION=([4]=CELLCENTERED)')" )
            endif

            write(22,*) ' SOLUTIONTIME=',TET

            ! write the data
            if(neMB.gt.0) then
              write(22,'(6(1x,e14.6))') ((zp(j),j=1,np),(sbot(j),j=1,nsides),k=1,1) !z
            endif
            
            write(22,'(6(1x,e14.6))') (eta(j),j=1,ne)

            zero = 0.0
            write(22,'(6(1x,e14.6))') ((uc(1,j),j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,1)
            write(22,'(6(1x,e14.6))') ((vc(1,j),j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,1)

            if(nson.gt.0) then
              write(22,'(6(1x,e14.6))') (wz(1,j),j=1,ne)
              write(22,'(6(1x,e14.6))') (qp(1,j),j=1,ne)
            elseif(nsed.gt.0) then
              write(22,'(6(1x,e14.6))') (cc(1,j),j=1,ne) !(dsg(j),j=1,ne)
            elseif(nsol.ne.0) then
              if(iOPsol.eq.0) then
                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
              elseif(iOPsol.eq.1) then
                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
                write(22,'(6(1x,e14.6))') (PSU(1,j),j=1,ne)
              elseif(iOPsol.eq.2) then
                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
                write(22,'(6(1x,e14.6))') (TC(1,j),j=1,ne)
              elseif(iOPsol.eq.3) then
                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
                write(22,'(6(1x,e14.6))') (PSU(1,j),j=1,ne)
                write(22,'(6(1x,e14.6))') (TC(1,j),j=1,ne)
              endif
            elseif(MBout.gt.0) then
              write(22,'(6(1x,e14.6))') (etaMB(j),j=1,ne)
              write(22,'(6(1x,e14.6))') (uc(2,j),j=1,np),(( un(2,j)*sdy(j)+ut(2,j)*sdx(j)),j=1,nsides)
              write(22,'(6(1x,e14.6))') (vc(2,j),j=1,np),((-un(2,j)*sdx(j)+ut(2,j)*sdy(j)),j=1,nsides)
            elseif(neqtide.gt.0) then
              write(22,'(6(1x,e14.6))') (eqtide(j),j=1,ne)
            endif

          endif
        
        CASE(3)  !eta, u, v interpolated to centroid. 2D.

! *** generate ut for output.
          if(npvb.lt.0) npvb = max(npv,min(2,2*max(0,izcoord-9)))
          Rn0(:,1:nsides) = un(:,1:nsides)
          call  GlobalQInterp ()

!          if(nopttype.eq.1) then !write nodal eta for harmonic analysis
!            write(22) TET,ne,ndf0,ndf0
!            write(22) (eta(j),j=1,ne)

!          elseif(nopttype.eq.2) then !write nodal eta,u,v for harmonic analysis
!            write(22) TET,ne,nsides,npv
!            write(22) (eta(j),j=1,ne)
!            do k=1,npv
!              write(22) (( un(k,j)*sdy(j)+Rn0(k,j)*sdx(j)),j=1,nsides)
!              write(22) ((-un(k,j)*sdx(j)+Rn0(k,j)*sdy(j)),j=1,nsides)
!            enddo

          if (noptcount.eq.0) then
          ! first time through start new file
            if(npv.gt.1) then
              if(nson.gt.0) then
                write(22,*)'VARIABLES="X" "Y" "Z" "U" "V" "W" "Q" '
              elseif(nsed.gt.0) then
                write(22,*)'VARIABLES="X" "Y" "Z" "U" "V" "W" "S" '
              elseif(nsol.ne.0) then
                if(iOPsol.eq.0) then
                  write(22,*)'VARIABLES="X" "Y"  "Z" "U" "V" "W" "Sigt" '
                elseif(iOPsol.eq.1) then
                  write(22,*)'VARIABLES="X" "Y"  "Z" "U" "V" "W" "Sigt" "PSU" '
                elseif(iOPsol.eq.2) then
                  write(22,*)'VARIABLES="X" "Y"  "Z" "U" "V" "W" "Sigt" "T" '
                elseif(iOPsol.eq.3) then
                  write(22,*)'VARIABLES="X" "Y"  "Z" "U" "V" "W" "Sigt" "PSU" "T" '
                endif
              else
                write(22,*)'VARIABLES="X" "Y" "Z" "U" "V" "W" '
              endif
            else
              if(nson.gt.0) then
                write(22,*)'VARIABLES="X" "Y" "Z" "U" "V" "Q" '
              elseif(nsed.gt.0) then
                write(22,*)'VARIABLES="X" "Y" "Z" "U" "V" "S" '
              elseif(nsol.ne.0) then
                if(iOPsol.eq.0) then
                  write(22,*)'VARIABLES="X" "Y"  "Z" "U" "V" "Sigt" '
                elseif(iOPsol.eq.1) then
                  write(22,*)'VARIABLES="X" "Y"  "Z" "U" "V" "Sigt" "PSU" '
                elseif(iOPsol.eq.2) then
                  write(22,*)'VARIABLES="X" "Y"  "Z" "U" "V" "Sigt" "T" '
                elseif(iOPsol.eq.3) then
                  write(22,*)'VARIABLES="X" "Y"  "Z" "U" "V" "Sigt" "PSU" "T" '
                endif
              else
                write(22,*)'VARIABLES="X" "Y" "Z" "U" "V"'
              endif
            endif
      
            if(npv.gt.1) then
              write(22,"('ZONE N=',i7,' E=',i7 )" ) (np+nsides)*(npvC+1),ne*(npvC)
            else
              write(22,"('ZONE N=',i7,' E=',i7 )" ) (np+nsides)*2, ne
            endif

!            if(npv.gt.1) then
              write(22,"(' ZONETYPE=FEBRICK DATAPACKING=BLOCK')" )
!            elseif(ncn.eq.3) then
!              write(22,"(' ZONETYPE=FETRIANGLE DATAPACKING=BLOCK')" )
!            elseif(ncn.eq.4) then
!              write(22,"(' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK')" )
!            endif

            if(npv.gt.1) then
              if(nsol.ne.0) then
                if(iOPsol.eq.0) then
                    write(22,"(' VARLOCATION=([6-7]=CELLCENTERED)')" )
                elseif(iOPsol.eq.1) then
                    write(22,"(' VARLOCATION=([6-8]=CELLCENTERED)')" )
                elseif(iOPsol.eq.2) then
                    write(22,"(' VARLOCATION=([6-8]=CELLCENTERED)')" )
                elseif(iOPsol.eq.3) then
                    write(22,"(' VARLOCATION=([6-9]=CELLCENTERED)')" )
                endif
              else
                write(22,"(' VARLOCATION=([6]=CELLCENTERED)')" )
              endif
            else
              if(nsol.ne.0) then
                if(iOPsol.eq.0) then
                    write(22,"(' VARLOCATION=([6]=CELLCENTERED)')" )
                elseif(iOPsol.eq.1) then
                    write(22,"(' VARLOCATION=([6-7]=CELLCENTERED)')" )
                elseif(iOPsol.eq.2) then
                    write(22,"(' VARLOCATION=([6-7]=CELLCENTERED)')" )
                elseif(iOPsol.eq.3) then
                    write(22,"(' VARLOCATION=([6-8]=CELLCENTERED)')" )
                endif
              elseif(nsed.gt.0.or.nson.ne.0) then
                write(22,"(' VARLOCATION=([6]=CELLCENTERED)')" )
              else
!                write(22,"(' VARLOCATION=([4]=CELLCENTERED)')" )
              endif
            endif

            write(22,*) ' SOLUTIONTIME=',TET

! write x,y
            write(22,'(6(1x,1pe15.7))') ((xp(j),j=1,np),(sxy(1,j),j=1,nsides),k=1,max(2,npvC+1)) !x
            write(22,'(6(1x,1pe15.7))') ((yp(j),j=1,np),(sxy(2,j),j=1,nsides),k=1,max(2,npvC+1)) !y

            rhv = 0.
            do nn=1,ne
              ncn2 = ncn -1 + min0(1,nen(nn,ncn))
              deptest = max(sdep(numsideT(1,nn)),sdep(numsideT(2,nn)),sdep(numsideT(3,nn)),sdep(numsideT(ncn2,nn)))
              if(deptest.le.depmin) cycle
              zncn = 1./dble(ncn2)
              DO J=1,ncn2
                MR = NEN(NN,J)
                gamma(mr) = gamma(mr) + area(nn)*zncn
                rhv(mr) = rhv(mr) + eta(nn)*area(nn)*zncn
              enddo
            enddo
            do j=1,np
!              if(nbc(j).lt.0) then
!                rhv(j) = spec(-nbc(j))
              if(gamma(j).gt.0) then
                rhv(j) = rhv(j)/gamma(j)
              else
                rhv(j) = zp(j)
              endif
            enddo
                       
! *** write z
            if(npv.gt.1) then
              if(izcoord.eq.0.or.izcoord.eq.2) then
                do k=1,npvC+1
                  write(22,'(6(1x,e14.6))') (-zp(j)*zdep(k),j=1,np),(sdep(j)*zdep(k),j=1,nsides) !z
                enddo
              elseif(izcoord.eq.1.or.izcoord.eq.3) then
                do k=1,izgrid-1
                  write(22,'(6(1x,e14.6))') (-max(zp(j),zdep(izgrid))*zdep(k),j=1,np),&
                            (-max(sbot(j),zdep(izgrid))*zdep(k),j=1,nsides) !z
                enddo
                do k=izgrid,npvC+1
                  write(22,'(6(1x,e14.6))') (max(zp(j),zdep(k)),j=1,np),&
                             (max(sbot(j),zdep(k)),j=1,nsides) !z
                enddo
              endif
            elseif(neMB.gt.0) then
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(-sdep(j),j=1,nsides) !surface
              write(22,'(6(1x,e14.6))') (zp(j),j=1,np),(-sdep(j),j=1,nsides) !etaMB
            else
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(-sdep(j),j=1,nsides) !eta
              write(22,'(6(1x,e14.6))') (zp(j),j=1,np),(sbot(j),j=1,nsides) !z
            endif

! *** write u
            zero = 0.0
            write(22,'(6(1x,e14.6))') ((uc(k,j),j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)
            if(izcoord.ge.2) then
              write(22,'(6(1x,e14.6))') (zero,j=1,np),(zero,j=1,nsides)
            endif
            if(npv.eq.1) then
                write(22,'(6(1x,e14.6))') ((uc(k,j),j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)
            endif

! *** write v
            write(22,'(6(1x,e14.6))') ((vc(k,j),j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)
            if(izcoord.ge.2) then
              write(22,'(6(1x,e14.6))') (zero,j=1,np),(zero,j=1,nsides)
            endif
            if(npv.eq.1) then
              write(22,'(6(1x,e14.6))') ((vc(k,j),j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)
            endif

! *** write w
            if(npv.gt.1) then
              write(22,'(6(1x,e14.6))') ((wz(k,j),k=1,npvC),j=1,ne)
            endif

! *** write scalar variables
            npvm = max(1,npvC)
            if(nson.gt.0) then
              write(22,'(6(1x,e14.6))') ((qp(k,j),k=1,npvm),j=1,ne)
            elseif(nsed.gt.0) then
              write(22,'(6(1x,e14.6))') ((cc(k,j),k=1,npvm+1),j=1,ne)
            elseif(nsol.ne.0) then
              if(iOPsol.eq.0) then
                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.1) then
                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
                write(22,'(6(1x,e14.6))') ((PSU(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.2) then
                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
                write(22,'(6(1x,e14.6))') ((TC(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.3) then
                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
                write(22,'(6(1x,e14.6))') ((PSU(k,j),k=1,npvm),j=1,ne)
                write(22,'(6(1x,e14.6))') ((TC(k,j),k=1,npvm),j=1,ne)
              endif
            endif

            ! and the elements

!            write(*,*) ' writing elements, kmin=', kmin

            if(npv.gt.1) then
              nps = np + nsides
              do j=1,ne
              ! and the elements
                if(ncn.eq.4.and.NEN(j,ncn).eq.0) then
                  NENtmp = NEN(j,ncn-1)
                else
                  NENtmp = NEN(j,ncn)
                endif
                do kv=1,npvC
                  write(22,'(8(1x,I8))') (NEN(j,k)+kv*nps,k=1,3),NENtmp+kv*nps,(NEN(j,k)+(kv-1)*nps,k=1,3),NENtmp+(kv-1)*nps
                enddo
              enddo
            else
              nps = np + nsides
              do j=1,ne
                if(NEN(j,ncn).eq.0) then
                  NENtmp = NEN(j,ncn-1)
                else
                  NENtmp = NEN(j,ncn)
                endif
                kv=1
                write(22,'(8(1x,I8))') (NEN(j,k)+kv*nps,k=1,3),NENtmp+kv*nps,(NEN(j,k)+(kv-1)*nps,k=1,3),NENtmp+(kv-1)*nps
              enddo
            endif

          else 
          ! append file
            write(22,"('ZONE D=(1,2,FECONNECT)')" ) 

!            if(npv.gt.1) then
              write(22,"(' ZONETYPE=FEBRICK DATAPACKING=BLOCK')" )
!            elseif(ncn.eq.3) then
!              write(22,"(' ZONETYPE=FETRIANGLE DATAPACKING=BLOCK')" )
!            elseif(ncn.eq.4) then
!              write(22,"(' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK')" )
!            endif

            if(npv.gt.1) then
              if(nsol.ne.0) then
                if(iOPsol.eq.0) then
                    write(22,"(' VARLOCATION=([6-7]=CELLCENTERED)')" )
                elseif(iOPsol.eq.1) then
                    write(22,"(' VARLOCATION=([6-8]=CELLCENTERED)')" )
                elseif(iOPsol.eq.2) then
                    write(22,"(' VARLOCATION=([6-8]=CELLCENTERED)')" )
                elseif(iOPsol.eq.3) then
                    write(22,"(' VARLOCATION=([6-9]=CELLCENTERED)')" )
                endif
              else
                write(22,"(' VARLOCATION=([6]=CELLCENTERED)')" )
              endif
            else
              if(nsol.ne.0) then
                if(iOPsol.eq.0) then
                    write(22,"(' VARLOCATION=([6]=CELLCENTERED)')" )
                elseif(iOPsol.eq.1) then
                    write(22,"(' VARLOCATION=([6-7]=CELLCENTERED)')" )
                elseif(iOPsol.eq.2) then
                    write(22,"(' VARLOCATION=([6-7]=CELLCENTERED)')" )
                elseif(iOPsol.eq.3) then
                    write(22,"(' VARLOCATION=([6-8]=CELLCENTERED)')" )
                endif
              elseif(nsed.gt.0.or.nson.ne.0) then
                write(22,"(' VARLOCATION=([6]=CELLCENTERED)')" )
              else
!                write(22,"(' VARLOCATION=([6]=CELLCENTERED)')" )
              endif
            endif

            write(22,*) ' SOLUTIONTIME=',TET

            ! write the data
!            if(npv.eq.1) write(22,'(6(1x,e14.6))') (eta(j),j=1,ne)

            rhv = 0.
            gamma = 0.
            do nn=1,ne
              ncn2 = ncn -1 + min0(1,nen(nn,ncn))
!              deptest = max(sdep(numsideT(1,nn)),sdep(numsideT(2,nn)),sdep(numsideT(3,nn)),sdep(numsideT(ncn2,nn)))
!              if(deptest.le.depmin) cycle
!              zncn = 1./dble(ncn2)
              DO J=1,ncn2
                MR = NEN(NN,J)
                gamma(mr) = gamma(mr) + area(nn) !*zncn
                rhv(mr) = rhv(mr) + eta(nn)*area(nn) !*zncn
              enddo
            enddo
            do j=1,np
!              if(nbc(j).lt.0) then
!                rhv(j) = spec(-nbc(j))
              if(gamma(j).gt.0) then
                rhv(j) = rhv(j)/gamma(j)
              else
                rhv(j) = zp(j)
              endif
            enddo
                       
! *** write z
            if(npv.gt.1) then
              if(izcoord.eq.0.or.izcoord.eq.2) then
                do k=1,npvC+1
                  write(22,'(6(1x,e14.6))') (-zp(j)*zdep(k),j=1,np),(sdep(j)*zdep(k),j=1,nsides) !z
                enddo
              elseif(izcoord.eq.1.or.izcoord.eq.3) then
                do k=1,izgrid-1
                  write(22,'(6(1x,e14.6))') (-max(zp(j),zdep(izgrid))*zdep(k),j=1,np),&
                            (-max(sbot(j),zdep(izgrid))*zdep(k),j=1,nsides) !z
                enddo
                do k=izgrid,npvC+1
                  write(22,'(6(1x,e14.6))') (max(zp(j),zdep(k)),j=1,np),&
                             (max(sbot(j),zdep(k)),j=1,nsides) !z
                enddo
              endif
            elseif(neMB.gt.0) then
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(-sdep(j),j=1,nsides) !eta
              write(22,'(6(1x,e14.6))') (zp(j),j=1,np),(-sdep(j),j=1,nsides) !etaMB
            else
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(-sdep(j),j=1,nsides) !eta
              write(22,'(6(1x,e14.6))') (zp(j),j=1,np),(sbot(j),j=1,nsides) !z
            endif

! *** write u
            zero = 0.0
            write(22,'(6(1x,e14.6))') ((uc(k,j),j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)
            if(izcoord.ge.2) then
              write(22,'(6(1x,e14.6))') (zero,j=1,np),(zero,j=1,nsides)
            endif
            if(npv.eq.1) then
                write(22,'(6(1x,e14.6))') ((uc(k,j),j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)
            endif

! *** write v
            write(22,'(6(1x,e14.6))') ((vc(k,j),j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)
            if(izcoord.ge.2) then
              write(22,'(6(1x,e14.6))') (zero,j=1,np),(zero,j=1,nsides)
            endif
            if(npv.eq.1) then
              write(22,'(6(1x,e14.6))') ((vc(k,j),j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)
            endif

! *** write w
            if(npv.gt.1) then
              write(22,'(6(1x,e14.6))') ((wz(k,j),k=1,npvC),j=1,ne)
            endif

! *** write scalar variables
            npvm = max(1,npvC)
            if(nson.gt.0) then
              write(22,'(6(1x,e14.6))') ((qp(k,j),k=1,npvm),j=1,ne)
            elseif(nsed.gt.0) then
              write(22,'(6(1x,e14.6))') ((cc(k,j),k=1,npvm+1),j=1,ne)
            elseif(nsol.ne.0) then
              if(iOPsol.eq.0) then
                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.1) then
                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
                write(22,'(6(1x,e14.6))') ((PSU(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.2) then
                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
                write(22,'(6(1x,e14.6))') ((TC(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.3) then
                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
                write(22,'(6(1x,e14.6))') ((PSU(k,j),k=1,npvm),j=1,ne)
                write(22,'(6(1x,e14.6))') ((TC(k,j),k=1,npvm),j=1,ne)
              endif
            endif

          endif
        
        CASE(4)  !case 4 is max values

! *** then write it out

          if (noptcount.eq.0) then

          ! first time through start new file

            write(22,*)'VARIABLES="X" "Y" "Z" "ECode" "EMAX" "EMIN" "SMAX" "T1" "TWET"'

            if(ncn.eq.3) then
              write(22,"('ZONE N=',i7,' E=',i7 )" ) np,ne
              write(22,"(' ZONETYPE=FETRIANGLE DATAPACKING=BLOCK')" )
              write(22,"(' VARLOCATION=([4-9]=CELLCENTERED)')" )
            elseif(ncn.eq.4) then
              write(22,"('ZONE N=',i7,' E=',i7 )" ) np,ne
              write(22,"(' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK')" )
              write(22,"(' VARLOCATION=([4-9]=CELLCENTERED)')" )
            endif

            !set original coordinates if icoord>0
            if(ixycoord.eq.0) then
!              x00 = -long0 ! 0. !offsets: from input file?
!              y00 = -lat0  ! 0.
              bigrI = rad2deg/bigr
              bigrcI = bigrI/cos(deg2rad*y0off) 
              do j=1,np
                rhv(j) = xp(j)*bigrcI 
                gamma(j) = yp(j)*bigrI
              enddo
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np)
              write(22,'(6(1x,e14.6))') (gamma(j),j=1,np)
              write(22,'(6(1x,e14.6))') (zp(j),j=1,np)
            else
!              x00 = 0.
!              y00 = 0.
              write(22,'(6(1x,e14.6))') (xp(j),j=1,np)
              write(22,'(6(1x,e14.6))') (yp(j),j=1,np)
              write(22,'(6(1x,e14.6))') (zp(j),j=1,np)
            endif

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
                Rn0(1,j)= sqrt(uu**2 + vv**2)
              endif
            enddo

            write(22,'(6(1x,i5))') (IECode(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (emax(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (emin(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (Rn0(1,j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (umax(j),j=1,nsides)
!            write(22,'(6(1x,e14.6))') (vmax(j),j=1,nsides)
            write(22,'(6(1x,e14.6))') (t1(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (twet(j),j=1,ne)

            ! and the elements
            do j=1,ne
              if(NEN(j,ncn).eq.0) then
                NENtmp = NEN(j,ncn-1)
              else
                NENtmp = NEN(j,ncn)
              endif
              write(22,*) (NEN(j,k),k=1,ncn-1),NENtmp
            enddo
          
          else 

          ! append file
            if(ncn.eq.3) then
              write(22,"('ZONE D=(1,2,3,FECONNECT)')" ) 
              write(22,"(' ZONETYPE=FETRIANGLE DATAPACKING=BLOCK')" )
              write(22,"(' VARLOCATION=([4-9]=CELLCENTERED)')" )
            elseif(ncn.eq.4) then
              write(22,"('ZONE D=(1,2,3,FECONNECT)')" ) 
              write(22,"(' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK')" )
              write(22,"(' VARLOCATION=([4-9]=CELLCENTERED)')" )
            endif

            ! write the data

            rhv = 0.
            do j=1,ne
              ncn2 = ncn -1 + min0(1,nen(j,ncn))
              topomin = min(sbot(numsideT(1,j)),sbot(numsideT(2,j)),&
                            sbot(numsideT(3,j)),sbot(numsideT(ncn2,j)))
              deptest = emax(j) - topomin
              if(deptest.gt.depmin) then
                uu = (umax(numsideT(1,j))+umax(numsideT(2,j))+umax(numsideT(3,j)))/3.
                vv = (vmax(numsideT(1,j))+vmax(numsideT(2,j))+vmax(numsideT(3,j)))/3.
                rhv(j)=  sqrt(uu**2 + vv**2)
              endif
            enddo

            write(22,'(6(1x,i5))') (IECode(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (emax(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (emin(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (rhv(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (umax(j),j=1,nsides)
!            write(22,'(6(1x,e14.6))') (vmax(j),j=1,nsides)
            write(22,'(6(1x,e14.6))') (t1(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (twet(j),j=1,ne)
          endif
       
        END SELECT

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
