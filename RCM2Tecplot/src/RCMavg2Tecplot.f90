!***************************************************************

  Module RCMArrays
  
      implicit none

      integer :: ne,np,nsides,npv,ncn,npvc=1   !nph,nphu,npv  
      integer :: nson,nsed,nsol,iOPsol
      integer :: neqtide=0, neMB=0, nsbc
      integer :: nopt
      integer :: izcoord=2,izgrid=0
      integer :: jUprofile=0,jCprofile=0, jSprofile=0
      integer, allocatable :: nen(:,:),numsideT(:,:),iside(:,:)
      integer, allocatable :: iends(:,:)
      real ::  TET
      real,parameter :: depmin=0.01
      real, allocatable ::  xp(:),yp(:),zp(:),zdep(:),area(:)
      real, allocatable ::  eta(:),un(:,:),ut(:,:),wz(:,:)
      real, allocatable ::  uc(:,:),vc(:,:),spd(:,:),spdcc(:,:)
      real, allocatable ::  spd2d(:),spd2dcc(:)
      real, allocatable ::  sxy(:,:),sbot(:),sdx(:),sdy(:),sdep(:)
      real, allocatable ::  rhv(:),gamma(:) !scratch vectors
      character*256 OutResFile
 
  end module RCMArrays

!***************************************************************

      Program RCM2Tecplot
  
      use RCMArrays
  
      implicit none

      integer :: i,ii,j,jj,k,kv,it, noptcount, istat,npvgrd,istart,istop
      real ::  fcount,fcountcc,ub,vb,spd1,spd2 
      character*256 fnamedata
      character(18) :: AnalysisTime

!     read time slices from the binary data file

      write(*,*) ' Enter filename for input RiCOM binary file'
      read(*,'(a)') fnamedata
      open(unit=20,file=fnamedata,status='old',form='unformatted')

      write(*,*) ' Enter izcoord'
      read(*,*) izcoord

      write(*,*) ' Enter time slice start index'
      read(*,*) istart

      write(*,*) ' Enter time slice end index'
      read(*,*) istop

      write(*,*) ' Enter output option'
      read(*,*) nopt

      write(*,*) ' Enter filename for output Tecplot ascii file'
      read(*,'(a)') OutResFile
      
      read(20) AnalysisTime
      write(*,*) 'AnalysisTime= ',AnalysisTime
      read(20)  !skip windfilename

      read(20) ne,np,nsides,npv,ncn
      read(20) nson,nsed,nsol,neMB,neqtide

      npv = max(npv,1)

      write(*,*) 'ne,np,nsides,npv=',ne,np,nsides,npv
      
      npvgrd = npv + min(1,max(0,izcoord-1))
        
      ALLOCATE ( xp(np), yp(np), zp(np), nen(ne,ncn), &
          sxy(2,nsides),sbot(nsides),sdx(nsides),sdy(nsides), &
          numsideT(ncn,ne),iside(2,nsides),area(ne),iends(2,nsides), &
          uc(npvgrd,np),vc(npvgrd,np),zdep(npvgrd),rhv(nsides),gamma(nsides),sdep(nsides), &
          eta(ne),un(npv,nsides),ut(npv,nsides),wz(npv,ne), &
          spd(npvgrd,np),spdcc(npvgrd-1,ne),spd2d(np),spd2dcc(ne), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate ncon main storage arrays',istat
        stop
      endif
      
! *** read coordinates
      read(20) (xp(j),j=1,np)
      read(20) (yp(j),j=1,np)
      read(20) (zp(j),j=1,np)
    
      if(npv.gt.1) then
        read(20) (zdep(k),k=1,npvgrd)    !write  z-grid (m)
      endif

! *** read element list
      read(20) ((nen(j,k),k=1,ncn),j=1,ne),(area(j),j=1,ne)
      read(20) (sdx(j),j=1,nsides),(sdy(j),j=1,nsides), &
               ((sxy(i,j),i=1,2),j=1,nsides),(sbot(j),j=1,nsides)
      read(20) ((numsideT(i,j),i=1,ncn),j=1,ne),((iside(i,j),i=1,2),j=1,nsides)
      
!      call GetArea()
      call EdgeSetup()
      
! *** initialize
      un = 0.
      ut = 0.
      if(npv.gt.1) wz = 0.

      do it=1,istart-1
        read(20, IOSTAT=istat) TET
        if(istat.ne.0) then
          exit
        endif
        read(20, IOSTAT=istat) !(eta(j),j=1,ne)
! *** read velocity (un,ut)
        if(istat.eq.0) read(20, IOSTAT=istat) !((un(k,i),k=1,npv),i=1,nsides)
        if(istat.eq.0) read(20, IOSTAT=istat) !((ut(k,i),k=1,npv),i=1,nsides)    
        if(istat.eq.0.and.npv.gt.1) read(20, IOSTAT=istat) !((wz(k,i),k=1,npvgrd-1),i=1,ne)

        if(nson.gt.0.and.istat.eq.0) then
          read(20, IOSTAT=istat)
        endif
        if(nsed.gt.0.and.istat.eq.0) then
          read(20, IOSTAT=istat)
        endif
        if(nsol.gt.0.and.istat.eq.0) then
          read(20, IOSTAT=istat)
        endif
        if(istat.ne.0) then
          write(*,*) ' Error reading input file at time=', TET
          exit
        endif
      enddo

      spd = 0.
      spd2d = 0.
      spd2dcc = 0.
      fcount = 0.
      fcountcc = 0.
      
      do it = 1,istop-istart+1
        read(20, IOSTAT=istat) TET 
        if(istat.ne.0) then
          exit
        endif
        read(20, IOSTAT=istat) (eta(j),j=1,ne)
! *** read velocity (un,ut)
        if(istat.eq.0) read(20, IOSTAT=istat) ((un(k,i),k=1,npv),i=1,nsides)
        if(istat.eq.0) read(20, IOSTAT=istat) ((ut(k,i),k=1,npv),i=1,nsides)    
        if(istat.eq.0.and.npv.gt.1) read(20, IOSTAT=istat) ((wz(k,i),k=1,npvgrd-1),i=1,ne)

! *** average speed
        call  GlobalQInterp ()
        fcount = fcount + 1.
        do j=1,np
          do k=1,npv
            spd1 = sqrt(uc(k,j)**2 + vc(k,j)**2)
            spd(k,j) = spd(k,j) + spd1
            spd2d(j) = spd2d(j) + spd1*(zdep(k)-zdep(k+1))
          enddo
        enddo

! *** average speed to centroid then calculate vertical average
        fcountcc = fcountcc + 1.
        do j=1,ne
          do kv=1,npv
            ub = 0.
            vb = 0.
            do jj=1,ncn
              ii = numsideT(jj,j)
              ub = ub + un(kv,ii)*sdy(ii) + ut(kv,ii)*sdx(ii)
              vb = vb - un(kv,ii)*sdx(ii) + ut(kv,ii)*sdy(ii)
            enddo
            ub = ub/float(ncn)
            vb = vb/float(ncn)
            spd1 = sqrt(ub*ub + vb*vb)
            spdcc(kv,j) = spdcc(kv,j) + spd1
            spd2dcc(j) = spd2dcc(j) + spd1*(zdep(kv)-zdep(kv+1))
          enddo
        enddo

        if(nson.gt.0.and.istat.eq.0) then
          read(20, IOSTAT=istat)
        endif
        if(nsed.gt.0.and.istat.eq.0) then
          read(20, IOSTAT=istat)
        endif
        if(nsol.gt.0.and.istat.eq.0) then
          read(20, IOSTAT=istat)
        endif
        if(istat.ne.0) then
          write(*,*) ' Error reading input file at time=', TET
          exit
        endif
      
      enddo
      
      spd = spd/fcount
      spdcc = spdcc/fcountcc
      spd2d = spd2d/fcount
      spd2dcc = spd2dcc/fcountcc
        
      if(izcoord.eq.2) then        
      !shift from layers to vertices in vertical
        do j=1,np
          spd1 = spd(1,j)
          do k=2,npv
            spd2 = spd(k,j)
            spd(k,j) = 0.5*(spd1 + spd2)
            spd1 = spd2
          enddo
          spd(npv+1,j) = spd1
        enddo
      endif

! *** write tecplot file
        write(*,*) ' call tecout, nopt=',nopt
        call OutputTecplot

! *** clean up
      close(20)
      
      stop
      end

!***************************************************************

      subroutine OutputTecplot()  !(noptcount, noptwrite)

      USE RCMArrays

      implicit none
      
      integer, parameter :: ndf0=0
      integer :: noptcount=0, istat, n23, npv0,npvgrd
      integer :: j,js,k,kk,kv,kmin,n1,n2,nn,ncn2,nentmp,mr,mr0,mr2,npvm,nps
      integer :: nex,ntypex,npx,nprx,ncnx,nenx,iex,js2,jn
      real :: cdep, zz(100), etaside,speed,topomin,deptest,zero,zncn,zlev
      real :: sn0,sn2,DNX1,DNY1,DNX2,DNY2,detn,un1,un2,uu,vv,u,v,w
      real :: bigrI,bigrcI,x00,y00,xc,yc,zc,depth,depdif,zbot
      real, allocatable, save :: etamax(:),spdmax(:),t1(:),twet(:),twet2(:),twet3(:),etalast(:)
      character(10) cseq

! *** Tecplot output- a single Tecplot file with many frames
        if (noptcount.eq.0) then
          ! first time through start new file
            open(unit=22, file=OutResFile,status='unknown')
        else 
          ! append file
            open(unit=22, file=OutResFile,status='old',position='append')
        endif

        SELECT CASE (nopt)

        CASE(1)  !eta,C on elements (2d), u,v on edges(3d). Model variable locations.

! *** generate ut for output.

          call  GlobalQInterp ()

          if (noptcount.eq.0) then
          ! first time through start new file
            if(nsed.gt.0) then
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
            elseif(neqtide.gt.0) then
              write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "EQT" '
            else
              write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "SPD" "SPDCC"'
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
            elseif(nsed.gt.0.or.nson.ne.0) then
              write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
            elseif(neqtide.gt.0) then
              write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
            else
              write(22,"(' VARLOCATION=([4,8]=CELLCENTERED)')" )
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

            write(22,'(6(1x,e14.6))') ((spd2d(j),j=1,np),((zero),j=1,nsides),k=1,1)
            write(22,'(6(1x,e14.6))') (spd2dcc(j),j=1,ne)

            if(nson.gt.0) then
!              write(22,'(6(1x,e14.6))') (qp(j),j=1,ne)
            endif
            if(nsed.gt.0) then
!              write(22,'(6(1x,e14.6))') (cc(1,j),j=1,ne) !(dsg(j),j=1,ne)
            endif
            if(nsol.ne.0) then
              if(iOPsol.eq.0) then
!                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
              elseif(iOPsol.eq.1) then
!                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
!                write(22,'(6(1x,e14.6))') (PSU(1,j),j=1,ne)
              elseif(iOPsol.eq.2) then
!                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
!                write(22,'(6(1x,e14.6))') (TC(1,j),j=1,ne)
              elseif(iOPsol.eq.3) then
!                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
!                write(22,'(6(1x,e14.6))') (PSU(1,j),j=1,ne)
!                write(22,'(6(1x,e14.6))') (TC(1,j),j=1,ne)
              endif
            endif
            if(neqtide.gt.0) then
!                write(22,'(6(1x,e14.6))') (eqtide(j),j=1,ne)
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
            elseif(nsed.gt.0.or.nson.ne.0) then
              write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
            elseif(neqtide.gt.0) then
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
!              write(22,'(6(1x,e14.6))') (qp(j),j=1,ne)
            endif
            if(nsed.gt.0) then
!              write(22,'(6(1x,e14.6))') (cc(1,j),j=1,ne) !(dsg(j),j=1,ne)
            endif
            if(nsol.ne.0) then
              if(iOPsol.eq.0) then
!                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
              elseif(iOPsol.eq.1) then
!                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
!                write(22,'(6(1x,e14.6))') (PSU(1,j),j=1,ne)
              elseif(iOPsol.eq.2) then
!                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
!                write(22,'(6(1x,e14.6))') (TC(1,j),j=1,ne)
              elseif(iOPsol.eq.3) then
!                write(22,'(6(1x,e14.6))') (sigt(1,j),j=1,ne)
!                write(22,'(6(1x,e14.6))') (PSU(1,j),j=1,ne)
!                write(22,'(6(1x,e14.6))') (TC(1,j),j=1,ne)
              endif
            endif
            if(neqtide.gt.0) then
!              write(22,'(6(1x,e14.6))') (eqtide(j),j=1,ne)
            endif

          endif
        
        CASE(2)  !eta, u, v, cc interpolated to vertices. 2D or 3D.
        
        CASE(3)  !eta, u, v interpolated to centroid. 2D only.

! *** generate nodal velocities for output.
          call  GlobalQInterp ()
          call SetDepth ()
          npvgrd = npv + min(1,max(0,izcoord-1))

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
                write(22,*)'VARIABLES="X" "Y" "Z" "U" "V" "spdcc" "spd"'
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


!            write(*,*) 'In output 1, neMB,npMB=',neMB,npMB
      
            if(npv.gt.1) then
              write(22,"('ZONE N=',i7,' E=',i7 )" ) np*npvgrd+nsides*npv,ne*(npvgrd-1)
            else
              write(22,"('ZONE N=',i7,' E=',i7 )" ) np*npvgrd+nsides*npv,ne*(npvgrd-1)
            endif

            write(22,"(' ZONETYPE=FEBRICK DATAPACKING=BLOCK')" )

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
!                write(22,"(' VARLOCATION=([4]=CELLCENTERED)')" ) no cc
              endif
            endif

            write(22,*) ' SOLUTIONTIME=',TET

            ! write the data
            if(npv.gt.1) then
              write(22,'(6(1x,1pe15.7))') ((xp(j),j=1,np),k=1,npvgrd),((sxy(1,j),j=1,nsides),k=1,npv) !x
              write(22,'(6(1x,1pe15.7))') ((yp(j),j=1,np),k=1,npvgrd),((sxy(2,j),j=1,nsides),k=1,npv) !y
            else
              write(22,'(6(1x,1pe15.7))') ((xp(j),j=1,np),k=1,npvgrd),((sxy(1,j),j=1,nsides),k=1,npv) !x
              write(22,'(6(1x,1pe15.7))') ((yp(j),j=1,np),k=1,npvgrd),((sxy(2,j),j=1,nsides),k=1,npv) !y
            endif
            
            rhv = 0.
            gamma = 0.
            do nn=1,ne
              ncn2 = ncn -1 + min0(1,nen(nn,ncn))
              deptest = max(sdep(numsideT(1,nn)),sdep(numsideT(2,nn)),sdep(numsideT(3,nn)),sdep(numsideT(ncn2,nn)))
              if(deptest.le.depmin) cycle
              zncn = 1./float(ncn2)
              DO J=1,ncn2
                MR = NEN(NN,J)
                gamma(mr) = gamma(mr) + area(nn)*zncn
                rhv(mr) = rhv(mr) + eta(nn)*area(nn)*zncn
              enddo
            enddo
            do j=1,np
!              if(nbc(j).lt.0) then
!                rhv(j) = spec(-nbc(j))
!              elseif(gamma(j).gt.0) then
              if(gamma(j).gt.0.) then
                rhv(j) = rhv(j)/gamma(j)
              else
                rhv(j) = zp(j)
              endif
            enddo
                       
            if(npv.gt.1) then
              if(izcoord.eq.0) then
                write(22,'(6(1x,e14.6))') (((rhv(j)+(rhv(j)-zp(j))*zdep(k)),j=1,np),k=1,npvgrd), &
                         ((sdep(j)+sbot(j)+sdep(j)*0.5*(zdep(k)+zdep(k+1)),j=1,nsides),k=1,npv) !z
              elseif(izcoord.eq.1) then
                do k=1,izgrid-1
!                  write(22,'(6(1x,e14.6))') (-amax1(zp(j),zdep(izgrid))*zdep(k),j=1,np),&
!                            (-amax1(sbot(j),zdep(izgrid))*zdep(k),j=1,nsides) !z
                enddo
                do k=izgrid,npv
!                  write(22,'(6(1x,e14.6))') (amax1(zp(j),zdep(k)),j=1,np),&
!                             (amax1(sbot(j),zdep(k)),j=1,nsides) !z
                enddo
              elseif(izcoord.eq.2) then
                write(22,'(6(1x,e14.6))') (((rhv(j)+(rhv(j)-zp(j))*zdep(k)),j=1,np),k=1,npvgrd), &
                         ((sdep(j)+sbot(j)+sdep(j)*0.5*(zdep(k)+zdep(k+1)),j=1,nsides),k=1,npv) !z
              elseif(izcoord.eq.3) then
              elseif(izcoord.eq.4) then
              endif
            elseif(neMB.gt.0) then
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(zp(j),j=1,np), &
                        (sdep(j)+sbot(j),j=1,nsides) !eta,zbot,level
            else
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(zp(j),j=1,np), &
                        (sdep(j)+sbot(j),j=1,nsides) !eta,zbot,level
            endif

            zero = 0.0
            write(22,'(6(1x,e14.6))') ((uc(k,j),j=1,np),k=1,npvgrd), &
                        ((( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)

            write(22,'(6(1x,e14.6))') ((vc(k,j),j=1,np),k=1,npvgrd), &
                        (((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)

            if(npv.gt.1) then
              write(22,'(6(1x,e14.6))') ((spdcc(k,j),k=1,npvgrd-1),j=1,ne)
              write(22,'(6(1x,e14.6))') ((spd(k,j),j=1,np),k=1,npvgrd), &
                          (((zero),j=1,nsides),k=1,npv)
            endif

            npvm = max(1,npv-1)
            if(nson.gt.0) then
!              write(22,'(6(1x,e14.6))') ((qp(j),k=1,npvm),j=1,ne)
            elseif(nsed.gt.0) then
!              write(22,'(6(1x,e14.6))') ((cc(k,j),k=1,npvm),j=1,ne)
            elseif(nsol.ne.0) then
              if(iOPsol.eq.0) then
!                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.1) then
!                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
!                write(22,'(6(1x,e14.6))') ((PSU(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.2) then
!                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
!                write(22,'(6(1x,e14.6))') ((TC(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.3) then
!                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
!                write(22,'(6(1x,e14.6))') ((PSU(k,j),k=1,npvm),j=1,ne)
!                write(22,'(6(1x,e14.6))') ((TC(k,j),k=1,npvm),j=1,ne)
              endif
            endif

            ! and the elements

!            write(*,*) ' writing elements, kmin=', kmin

            if(npv.gt.1) then
              nps = np  ! + nsides
              do j=1,ne
              ! and the elements
                if(ncn.eq.4.and.NEN(j,ncn).eq.0) then
                  NENtmp = NEN(j,ncn-1)
                else
                  NENtmp = NEN(j,ncn)
                endif
                do kv=1,npvgrd-1
                  write(22,'(8(1x,I8))') (NEN(j,k)+kv*nps,k=1,3),NENtmp+kv*nps, &
                                     (NEN(j,k)+(kv-1)*nps,k=1,3),NENtmp+(kv-1)*nps
                enddo
              enddo
            else
              nps = np  ! + nsides
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
              deptest = max(sdep(numsideT(1,nn)),sdep(numsideT(2,nn)),sdep(numsideT(3,nn)),sdep(numsideT(ncn2,nn)))
              if(deptest.le.depmin) cycle
              zncn = 1./float(ncn2)
              DO J=1,ncn2
                MR = NEN(NN,J)
                gamma(mr) = gamma(mr) + area(nn) !*zncn
                rhv(mr) = rhv(mr) + eta(nn)*area(nn) !*zncn
              enddo
            enddo
            do j=1,np
!              if(nbc(j).lt.0) then
!                rhv(j) = spec(-nbc(j))
              if(gamma(j).gt.0.) then
                rhv(j) = rhv(j)/gamma(j)
              else
                rhv(j) = zp(j)
              endif
            enddo
                       
            if(npv.gt.1) then
              if(izcoord.eq.0) then
                write(22,'(6(1x,e14.6))') (((rhv(j)+(rhv(j)-zp(j))*zdep(k)),j=1,np),k=1,npvgrd), &
                         ((sdep(j)+sbot(j)+sdep(j)*0.5*(zdep(k)+zdep(k+1)),j=1,nsides),k=1,npv) !z
              elseif(izcoord.eq.1) then
                do k=1,izgrid-1
!                  write(22,'(6(1x,e14.6))') (-amax1(zp(j),zdep(izgrid))*zdep(k),j=1,np),&
!                            (-amax1(sbot(j),zdep(izgrid))*zdep(k),j=1,nsides) !z
                enddo
                do k=izgrid,npv
!                  write(22,'(6(1x,e14.6))') (amax1(zp(j),zdep(k)),j=1,np),&
!                             (amax1(sbot(j),zdep(k)),j=1,nsides) !z
                enddo
              elseif(izcoord.eq.2) then
                write(22,'(6(1x,e14.6))') (((rhv(j)+(rhv(j)-zp(j))*zdep(k)),j=1,np),k=1,npvgrd), &
                         ((sdep(j)+sbot(j)+sdep(j)*0.5*(zdep(k)+zdep(k+1)),j=1,nsides),k=1,npv) !z
              elseif(izcoord.eq.3) then
              elseif(izcoord.eq.4) then
              endif
            elseif(neMB.gt.0) then
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(zp(j),j=1,np), &
                        (sdep(j)+sbot(j),j=1,nsides) !eta,zbot,level
            else
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(zp(j),j=1,np), &
                        (sdep(j)+sbot(j),j=1,nsides) !eta,zbot,level
            endif

            write(22,'(6(1x,e14.6))') ((uc(k,j),j=1,np),k=1,npvgrd), &
                        ((( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)

            write(22,'(6(1x,e14.6))') ((vc(k,j),j=1,np),k=1,npvgrd), &
                        (((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)

            if(npv.gt.1) then
              write(22,'(6(1x,e14.6))') ((wz(k,j),k=1,npv),j=1,ne)
              write(22,'(6(1x,e14.6))') ((spd(k,j),j=1,np),k=1,npvgrd), &
                          (((zero),j=1,nsides),k=1,npv)
            endif

            npvm = max(1,npv-1)
            if(nson.gt.0) then
!              write(22,'(6(1x,e14.6))') ((qp(j),k=1,npvm),j=1,ne)
            elseif(nsed.gt.0) then
!              write(22,'(6(1x,e14.6))') ((cc(k,j),k=1,npvm),j=1,ne)
            elseif(nsol.ne.0) then
              if(iOPsol.eq.0) then
!                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.1) then
!                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
!                write(22,'(6(1x,e14.6))') ((PSU(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.2) then
!                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
!                write(22,'(6(1x,e14.6))') ((TC(k,j),k=1,npvm),j=1,ne)
              elseif(iOPsol.eq.3) then
!                write(22,'(6(1x,e14.6))') ((sigt(k,j),k=1,npvm),j=1,ne)
!                write(22,'(6(1x,e14.6))') ((PSU(k,j),k=1,npvm),j=1,ne)
!                write(22,'(6(1x,e14.6))') ((TC(k,j),k=1,npvm),j=1,ne)
              endif
            endif

          endif
        
        CASE(4)  !eta, u, v, cc interpolated to vertices. 2D. case 4 is max values

! *** then write it out

          if (noptcount.eq.0) then

          ! first time through start new file

            write(22,*)'VARIABLES="X" "Y" "Z" "ECode" "EMAX" "UMAX" "T1" "TWET" "TWET2" "TWET3" "ELAST"'

            if(npv.gt.1) then
              write(22,"('ZONE N=',i7,' E=',i7 )" ) np*npv,max(ne,ne*(npv-1))
              write(22,"(' ZONETYPE=FEBRICK DATAPACKING=BLOCK')" )
              write(22,"(' VARLOCATION=([4-11]=CELLCENTERED)')" )
            elseif(ncn.eq.3) then
              write(22,"('ZONE N=',i7,' E=',i7 )" ) np,ne
              write(22,"(' ZONETYPE=FETRIANGLE DATAPACKING=BLOCK')" )
              write(22,"(' VARLOCATION=([4-11]=CELLCENTERED)')" )
            elseif(ncn.eq.4) then
              write(22,"('ZONE N=',i7,' E=',i7 )" ) np,ne
              write(22,"(' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK')" )
              write(22,"(' VARLOCATION=([4-11]=CELLCENTERED)')" )
            endif

            !read in original coordinates if icoord>0
              x00 = 0.
              y00 = 0.
              write(22,'(6(1x,e14.6))') (xp(j),j=1,np)
              write(22,'(6(1x,e14.6))') (yp(j),j=1,np)
              write(22,'(6(1x,e14.6))') (zp(j),j=1,np)

!            Rn0 = 0.
            do j=1,ne
              ncn2 = ncn -1 + min0(1,nen(j,ncn))
              topomin = amin1(sbot(numsideT(1,j)),sbot(numsideT(2,j)),&
                            sbot(numsideT(3,j)),sbot(numsideT(ncn2,j)))
              deptest = etamax(j) - topomin
!              if(deptest.gt.depmin) then
!                Rn0(1,j)= (spdmax(numsideT(1,j))+spdmax(numsideT(2,j))+spdmax(numsideT(3,j)))/3.
!              endif
            enddo

!            write(22,'(6(1x,i5))') (IECode(j),j=1,ne)
!!            write(22,'(6(1x,e14.6))') (etamax(j),j=1,ne)
!!            write(22,'(6(1x,e14.6))') (Rn0(1,j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (eta(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (eqtide(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (t1(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (twet(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (twet2(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (twet3(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (etalast(j),j=1,ne)

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
            if(npv.gt.1) then
              write(22,"('ZONE D=(1,2,FECONNECT)')" ) 
              write(22,"(' ZONETYPE=FEBRICK DATAPACKING=BLOCK')" )
              write(22,"(' VARLOCATION=([4-11]=CELLCENTERED)')" )
            elseif(ncn.eq.3) then
              write(22,"('ZONE D=(1,2,3,FECONNECT)')" ) 
              write(22,"(' ZONETYPE=FETRIANGLE DATAPACKING=BLOCK')" )
              write(22,"(' VARLOCATION=([4-11]=CELLCENTERED)')" )
            elseif(ncn.eq.4) then
              write(22,"('ZONE D=(1,2,3,FECONNECT)')" ) 
              write(22,"(' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK')" )
              write(22,"(' VARLOCATION=([4-11]=CELLCENTERED)')" )
            endif

            ! write the data

            rhv = 0.
            do j=1,ne
              ncn2 = ncn -1 + min0(1,nen(j,ncn))
              topomin = amin1(sbot(numsideT(1,j)),sbot(numsideT(2,j)),&
                            sbot(numsideT(3,j)),sbot(numsideT(ncn2,j)))
!              deptest = etamax(j) - topomin
!              if(deptest.gt.depmin) then
!                rhv(j)= (spdmax(numsideT(1,j))+spdmax(numsideT(2,j))+spdmax(numsideT(3,j)))/3.
!              endif
            enddo

!            write(22,'(6(1x,i5))') (IECode(j),j=1,ne)
!!            write(22,'(6(1x,e14.6))') (etamax(j),j=1,ne)
!!            write(22,'(6(1x,e14.6))') (rhv(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (eta(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (eqtide(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (t1(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (twet(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (twet2(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (twet3(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (etalast(j),j=1,ne)
          endif
       
        END SELECT

        close(22)

        if(npv.gt.1) then
          ! first time through start new file
          do j=1,jUprofile

            js = numsideT(1,1)  !tsUnodes(j))
            if(js.eq.0) cycle

            npv0 = npv
            n1 = iside(1,js)
            n2 = iside(2,js)
            if(n2.eq.0) then
              etaside = eta(n1)
            else
              etaside = max(eta(n1),eta(n2))
            endif

            if(izcoord.eq.0) then ! sigma
              do k=1,npv
!                zz(k) = etaside + sdep(js)*zdep(k)
              enddo
            elseif(izcoord.eq.1) then  ! sigma on z
              depdif = etaside - zdep(izgrid)
              do k=1,izgrid-1
                zz(k) = depdif*(1.+zdep(k)) + zdep(izgrid)
              enddo
              zbot = etaside-sdep(js)
              do k=izgrid,npv
                if(zbot.ge.zdep(k)) then
                  zz(k) = zbot
                  npv0 = k
                  exit
                else
                  zz(k) = zdep(k)
                endif
              enddo        
            elseif(izcoord.eq.2) then ! sigma layers
              do k=1,npv
                zz(k) = sbot(js) + sdep(js)*(1.+ 0.5*(zdep(k)+zdep(k+1))) !ref to MSL
              enddo
            endif

            if (noptcount.eq.0) then
              n23 = 23  !+j-1
              write(cseq,"(I2.2)") j
              open(unit=n23,file='profile'//trim(cseq)//'.dat',status='unknown')
              write(n23,*) 'VARIABLES="z" "u" "v" '
              write(n23,*) 'ZONE'
              write(n23,*) ' SOLUTIONTIME=',TET


              do k=1,npv0
!             ! write the data
                write(n23,'(3(1PE12.4))') zz(k),( un(k,js)*sdy(js)+ut(k,js)*sdx(js)), &
                                 (-un(k,js)*sdx(js)+ut(k,js)*sdy(js))
              enddo
              close(n23)
            else 
          ! append file
              n23 = 23  !+j-1
              write(cseq,"(I2.2)") j
              open(unit=n23,file='profile'//trim(cseq)//'.dat',status='old',position='append')
              write(n23,*) 'ZONE'
              write(n23,*) ' SOLUTIONTIME=',TET

              do k=1,npv0
!             ! write the data
                write(n23,'(3(1PE12.4))') zz(k),( un(k,js)*sdy(js)+ut(k,js)*sdx(js)), &
                                 (-un(k,js)*sdx(js)+ut(k,js)*sdy(js))
              enddo
              close(n23)
            endif
          enddo
        endif

        if(nsed.gt.0.and.npvc.gt.1) then
          if (noptcount.eq.0) then
          ! first time through start new file
            open(unit=23,file='profileSed.dat',status='unknown')
            write(23,*) 'VARIABLES="z" "C" '
            write(23,*) 'ZONE'
            write(23,*) ' SOLUTIONTIME=',TET

            js = jCProfile
            js2 = jUProfile
            do j=1,npvc
!           ! write the data
!              write(23,*) sdep(js2)*(1.+zdepC(j)), cc(j,js)
            enddo
!            write(23,*) sdep(js2)*(1.+zdep(npv)), cc(npvc+1,js)
            if(nsbc.ge.3) then
!              write(23,*) sdep(js2)*(1.+zdep(npv)), cc(npvc+2,js)
            endif
          else 
          ! append file
            open(unit=23,file='profileSed.dat',status='old',position='append')
            write(23,*) 'ZONE'
            write(23,*) ' SOLUTIONTIME=',TET

            js = jCProfile
            js2 = jUProfile
            do j=1,npvc
!           ! write the data
!              write(23,*) sdep(js2)*(1.+zdepC(j)), cc(j,js)
            enddo
!            write(23,*) sdep(js2)*(1.+zdep(npv)), cc(npvc+1,js)
            if(nsbc.ge.3) then
!              write(23,*) sdep(js2)*(1.+zdep(npv)), cc(npvc+2,js)
            endif
          endif
          close(23)
        endif

        if((nsol.ne.0).and.(npvc.gt.1)) then
          if (noptcount.eq.0) then
          ! first time through start new file
            do j=1,jSprofile
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

              js = 1  !tsSnodes(j)
!              cdep=amin1(xyz(nen(js,1),3),xyz(nen(js,2),3),xyz(nen(js,3),3),xyz(nen(js,ncn),3))
              ncn2 = ncn -1 + min0(1,nen(js,ncn))
              cdep = 0. 
              do jn=1,ncn2
                cdep = cdep + zp(nen(js,jn))
              enddo
              cdep = eta(js)-cdep/float(ncn2)
!              write(23,*) 'iOP,npvc=',iOPsol,npvc,js,cdep
!             ! write the data
              if(iOPsol.eq.0) then
                do k=1,npvc
!                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.1) then
                do k=1,npvc
!                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), PSU(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.2) then
                do k=1,npvc
!                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), TC(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.3) then
                do k=1,npvc
!                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), PSU(k,js), TC(k,js)
                enddo
                close(n23)
              endif
            enddo
          else 
          ! append file
            do j=1,jSprofile
              n23 = 23 !+j-1
              write(cseq,"(I2.2)") j
              open(unit=n23,file='profileCon'//trim(cseq)//'.dat',status='old',position='append')
!              open(unit=23,file='profileCon.dat',status='old',position='append')
              write(23,*) 'ZONE'
              write(23,*) ' SOLUTIONTIME=',TET

              js = 1  !tsSnodes(j)
              ncn2 = ncn -1 + min0(1,nen(js,ncn))
              cdep = 0. 
              do jn=1,ncn2
                cdep = cdep + zp(nen(js,jn))
              enddo
              cdep = eta(js)-cdep/float(ncn2)
!             ! write the data
              if(iOPsol.eq.0) then
                do k=1,npvc
!                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.1) then
                do k=1,npvc
!                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), PSU(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.2) then
                do k=1,npvc
!                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), TC(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.3) then
                do k=1,npvc
!                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), PSU(k,js), TC(k,js)
                enddo
                close(n23)
              endif
            enddo
          endif
        endif
        
        noptcount = noptcount + 1

! format statements

20    format( 'ZONE N=',i7,' E=',i7,' ET=TRIANGLE F=FEPOINT')

21    format('ZONE D=(1,2,3,FECONNECT) F=FEPOINT ET=TRIANGLE')
22    format('ZONE D=(1,2,FECONNECT) F=FEPOINT ET=TRIANGLE')

30    format( 'ZONE N=',i7,' E=',i7,' ET=QUADRILATERAL F=FEPOINT')

31    format('ZONE D=(1,2,3,FECONNECT) F=FEPOINT ET=QUADRILATERAL')
32    format('ZONE D=(1,2,FECONNECT) F=FEPOINT ET=QUADRILATERAL')

!40    format('ZONE I=',I7,' F=POINT')

!41    format('ZONE D=(1,2,3) F=POINT')

50    format('ZONE N=',i7,' E=',i7,' ET=BRICK F=FEPOINT')

51    format('ZONE D=(1,2,FECONNECT) F=FEPOINT ET=BRICK')

      RETURN
      END

!******************************************************************************

      SUBROUTINE SetDepth

! *** set the side depth

      USE RCMArrays

      implicit none
      
!   local variables
      integer :: j,n1,n2
      real :: depth, etaside, zmax, zmin

!   set sdep

!   Get new ones
      do j=1,nsides
        n1 = iside(1,j)
        n2 = iside(2,j)

        if(n2.le.0.or.n2.gt.ne) then
          etaside = eta(n1)
        else
          etaside = max(eta(n1),eta(n2))
        endif
        
        zmax = max(zp(iends(1,j)),zp(iends(2,j)))
        zmin = min(zp(iends(1,j)),zp(iends(2,j)))

        if(etaside.ge.zmax) then
          depth = etaside - 0.5*(zmax+zmin)  !refdep(j)
        elseif(etaside.gt.zmin) then
          depth = 0.5*(etaside-zmin)**2/(zmax-zmin)
!          depth = etaside - sbot(j)  !force old method
        else
          depth = 0.
        endif
        
        if(depth.gt.depmin) then
          sdep(j) = depth
        else
          sdep(j) = 0.
        endif
      enddo

      end

!*****************************************************************

      SUBROUTINE GetArea ()

! *** compute area of elements

      USE RCMArrays

      implicit none
      
! *** local variables
      integer :: I1,I2,I3,I4,nn,ncne

      write(*,*) 'Finding areas'
        DO nn=1,NE
          I1=Nen(nn,1)
          I2=Nen(nn,2)
          I3=Nen(nn,3)
          Area(nn)=0.5*(xp(I1)*(yp(I2)-yp(I3)) &
     &               +xp(I2)*(yp(I3)-yp(I1)) &
     &               +xp(I3)*(yp(I1)-yp(I2)))
          if(ncn.eq.4) then
            I4=Nen(ncn,nn)
            if(I4.gt.0) then
              Area(nn)=Area(nn) &
     &              +0.5*(xp(I3)*(yp(I4)-yp(I1)) &
     &                   +xp(I4)*(yp(I1)-yp(I3)) &
     &                   +xp(I1)*(yp(I3)-yp(I4)))
            endif
          endif
        enddo

      RETURN
      END

!*****************************************************************

      SUBROUTINE EdgeSetup ()

! *** compute side arrays

      USE RCMArrays

      implicit none
      
! *** local variables
      integer :: L,M,n1,n2,ncn2
      integer :: idxside

      iends = 0
      
! *** Fill side arrays
      DO L=1,NE
        ncn2 = ncn -1 + min0(1,nen(L,ncn))
        DO M=1,ncn2
          idxside = numsideT(M,L)
          if(iends(1,idxside).eq.0) then
!            iside(1,idxside) = L
            N1=NEN(L,M)
            N2=NEN(L,mod(M,ncn2)+1)
            iends(1,idxside) = n1
            iends(2,idxside) = n2
!            sxy(1,idxside)  = 0.5*(xyz(n2,1)+xyz(n1,1))
!            sxy(2,idxside)  = 0.5*(xyz(n2,2)+xyz(n1,2))
!            refdep(idxside) = 0.5*(xyz(n1,3)+xyz(n2,3))
!            slen(idxside)=sqrt((mf1(n2)*(xyz(n1,1)-xyz(n2,1)))**2 +(xyz(n1,2)-xyz(n2,2))**2)
!            sdx(idxside) = (mf1(n2)*(xyz(n2,1)-xyz(n1,1)))/slen(idxside)
!            sdy(idxside) = (xyz(n2,2)-xyz(n1,2))/slen(idxside)
!            if(ncn2.eq.3) then
!              dlinv(idxside) = 2.*area(L)/(3.*slen(idxside))
!            elseif(ncn2.eq.4) then
!              call GetDist2cm(L,M,cmdist)
!              dlinv(idxside) = cmdist  !0.5*slen(idxside)
!            else
!               write(lst6,*) 'Number of element nodes not 3 or 4',L
!              ierr = 72
!              call terminate(ierr, myproc)
              !call exit(72)
!            endif
!          elseif(iside(2,idxside).eq.0) then
!            iside(2,idxside) = L
          endif
        enddo
      enddo
      
      RETURN
      END

!!*****************************************************************

      SUBROUTINE GlobalQInterp ()

! *** Interpolate velocity to vertices and midsides

      USE RCMArrays
!      USE AdvectArrays

      implicit none

! *** local variables
      integer :: kv,nn,ncn2,j,mr,mr0,sn0,mr2,sn2,n1,n2
      real :: zncn,dnx1,dny1,dnx2,dny2,un1,un2,uu,vv,umid,vmid
      real :: detn,areasum

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
!*****************************************************************
