!***************************************************************

  Module RCMArrays
  
      implicit none

      integer :: ne,np,nsides,npv,ncn   !nph,nphu,npv  
      integer :: nson,nsed,nsol
      integer :: nopt
      integer, allocatable :: nen(:,:),numsideT(:,:)
      real, allocatable ::  xp(:),yp(:),zp(:)
      real, allocatable ::  eta(:),un(:,:),ut(:,:)
      real, allocatable ::  sdx(:),sdy(:)
      character*256 OutResFile
 
  end module RCMArrays

!***************************************************************

      Program RCM2Tecplot
  
      use RCMArrays
  
      implicit none

      integer :: noptcount
      character*256 fnamedata
      character(18) :: AnalysisTime

!     read time slices from the binary data file

      write(*,*) ' Enter filename for input RiCOM binary file'
      read(*,'(a)') fnamedata
      open(unit=20,file=fnamedata,status='old',form='unformatted')

      write(*,*) ' Enter output option'
      read(*,'(a)') nopt

      write(*,*) ' Enter filename for output Tecplot ascii file'
      read(*,'(a)') OutResFile
      
      read(20) AnalysisTime
      write(*,*) 'AnalysisTime= ',AnalysisTime
      read(20)  !skip windfilename

      read(20) ne,np,nsides,npv,ncn   !nph,nphu,npv  
      read(20) nson,nsed,nsol

      npv = max(npv,1)

      write(*,*) 'ne,np,nsides,npv=',ne,np,nsides,npv
        
      ALLOCATE ( xp(np), yp(np), zp(np), nen(ne,ncn), &
          sdx(nsides),sdy(nsides),numsideT(ncn,ne), &
          eta(ne),un(npv,np),ut(npv,np), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate ncon main storage arrays',istat
        stop
      endif
      
! *** read coordinates
      read(20) (xp(j),j=1,np)
      read(20) (yp(j),j=1,np)
      read(20) (zp(j),j=1,np)

! *** read element list
      read(20) ((nen(j,k),k=1,ncn),j=1,ne)
      read(20) (sdx(j),j=1,nsides),(sdy(j),j=1,nsides)
      read(20) ((numsideT(i,j),i=1,ncn),j=1,ne)

      do
        read(20, IOSTAT=istat) TET 
        if(istat.ne.0) then
          exit
        endif
        read(20, IOSTAT=istat) (eta(j),j=1,nph)
! *** read velocity (un,ut)
        if(istat.eq.0) read(20, IOSTAT=istat) ((un(kv,i),kv=1,npv),i=1,nphu)
        if(istat.eq.0) read(20, IOSTAT=istat) ((ut(kv,i),kv=1,npv),i=1,nphu)
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
          write(*,*) ' Error reading input file at time='.TET
          exit
        endif
      
! *** write tecplot file
        call OutputTecplot

      enddo

! *** clean up
      close(20)
      
      stop
      end

!***************************************************************

      subroutine OutputTecplot()  !(noptcount, noptwrite)

      USE RCMArrays

      implicit none
      
      integer, parameter :: ndf0=0
      integer :: noptcount, noptwrite, istat, n23, npv0
      integer :: j,js,k,kk,kv,kmin,n1,n2,nn,ncn2,nentmp,mr,mr0,mr2,npvm,nps
      integer :: nex,ntypex,npx,nprx,ncnx,nenx,iex,js2,jn
      integer, save :: nopttype, firstcount = 0
      real :: cdep, zz(100), etaside,speed,topomin,deptest,zero,zncn,zlev
      real :: sn0,sn2,DNX1,DNY1,DNX2,DNY2,detn,un1,un2,uu,vv,u,v,w
      real :: bigrI,bigrcI,x00,y00,xc,yc,zc,depth,depdif,zbot
      real, allocatable, save :: etamax(:),spdmax(:),t1(:),twet(:),twet2(:),twet3(:),etalast(:)
      character(10) cseq

        if(mod(nopt,10).eq.4.and.noptwrite.eq.0) then ! compute max values

          if(firstcount.eq.0) then ! allocate and initialize storagearrays
            ALLOCATE (etamax(ne),spdmax(nsides),t1(ne),twet(ne),twet2(ne),twet3(ne),etalast(ne), STAT = istat )
            if(istat.ne.0) then
              write(*,*) 'FATAL ERROR: Cannot allocate etamax vector arrays'
              stop 71
            endif
            firstcount = 1
            etamax=-99999.
            spdmax=0.
            t1 = -1.
            twet = 0.
            twet2 = 0.
            twet3 = 0.
            etalast = -99999.
          endif

          do j=1,nsides
            if(sdep(j).gt.depmin) then
              speed = sqrt(Un(1,j)**2 + Ut(1,j)**2)
              spdmax(j) = amax1(spdmax(j),speed)
            endif
          enddo

          do j=1,ne
            ncn2 = ncn -1 + min0(1,nen(j,ncn))
!            topomin = amin1(refdep(numsideT(1,j)),refdep(numsideT(2,j)),&
!                            refdep(numsideT(3,j)),refdep(numsideT(ncn2,j)))
            topomin = amin1(xyz(nen(j,1),3),xyz(nen(j,2),3),xyz(nen(j,3),3),&
                        xyz(nen(j,ncn2),3))
            deptest = eta(j) - topomin
            if(deptest.gt.depmin) then
              etamax(j) = amax1(etamax(j),eta(j))
              if(t1(j).lt.0.) t1(j) = idx*dt0
              twet(j) = twet(j) + dt0
              if(deptest.gt.0.5) twet2(j) = twet2(j) + dt0
              if(deptest.gt.1.5) twet3(j) = twet3(j) + dt0
!            elseif(twet(j).le.0..and.mod(idx,100).eq.0) then
!              eta(j) = topomin
            endif
            etalast(j) = eta(j)
          enddo

          return

        endif

! *** Tecplot output- a single Tecplot file with many frames
        if (noptcount.eq.0) then
          nopttype = nopt/10
          nopt = mod(nopt,10)
          ! first time through start new file
          if(nopttype.gt.0) then
            open(unit=22, file=OutResFile,status='unknown',form='unformatted')
          else
            open(unit=22, file=OutResFile,status='unknown')
          endif
        else 
          ! append file
          if(nopttype.gt.0) then
            open(unit=22, file=OutResFile,status='old',position='append',form='unformatted')
          else
            open(unit=22, file=OutResFile,status='old',position='append')
          endif
        endif

        SELECT CASE (nopt)

        CASE(1)  !eta,C on elements (2d), u,v on edges(3d). Model variable locations.

! *** generate ut for output.
          Rn0(:,1:nsides) = un(:,1:nsides)
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
            elseif(nsed.gt.0.or.nson.ne.0) then
              write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
            elseif(neqtide.gt.0) then
              write(22,"(' VARLOCATION=([4,7]=CELLCENTERED)')" )
            else
              write(22,"(' VARLOCATION=([4]=CELLCENTERED)')" )
            endif
            
            write(22,*) ' SOLUTIONTIME=',TET

            ! write the data
            write(22,'(6(1x,e14.6))') ((xyz(j,1),j=1,np),(sxy(1,j),j=1,nsides),k=1,1) !x
            write(22,'(6(1x,e14.6))') ((xyz(j,2),j=1,np),(sxy(2,j),j=1,nsides),k=1,1) !y
            write(22,'(6(1x,e14.6))') ((xyz(j,3),j=1,np),(refdep(j),j=1,nsides),k=1,1) !z

            write(22,'(6(1x,e14.6))') (eta(j),j=1,ne)

            zero = 0.0
            write(22,'(6(1x,e14.6))') ((Rn0Gx(j,1),j=1,np),(( un(k,j)*sdy(j)+Rt0G(j,k)*sdx(j)),j=1,nsides),k=1,1)
            write(22,'(6(1x,e14.6))') ((Rn0Gy(j,1),j=1,np),((-un(k,j)*sdx(j)+Rt0G(j,k)*sdy(j)),j=1,nsides),k=1,1)

            if(nson.gt.0) then
              write(22,'(6(1x,e14.6))') (qp(j),j=1,ne)
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
              write(22,'(6(1x,e14.6))') ((xyz(j,3),j=1,np),(refdep(j),j=1,nsides),k=1,1) !z
            endif
            
            write(22,'(6(1x,e14.6))') (eta(j),j=1,ne)

            zero = 0.0
            write(22,'(6(1x,e14.6))') ((Rn0Gx(j,1),j=1,np),(( un(k,j)*sdy(j)+Rt0G(j,k)*sdx(j)),j=1,nsides),k=1,1)
            write(22,'(6(1x,e14.6))') ((Rn0Gy(j,1),j=1,np),((-un(k,j)*sdx(j)+Rt0G(j,k)*sdy(j)),j=1,nsides),k=1,1)

            if(nson.gt.0) then
              write(22,'(6(1x,e14.6))') (qp(j),j=1,ne)
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
            elseif(neqtide.gt.0) then
                write(22,'(6(1x,e14.6))') (eqtide(j),j=1,ne)
            endif

          endif
        
        CASE(2)  !eta, u, v, cc interpolated to vertices. 2D or 3D.

          rhv = 0.
          gamma = 0.
          Rn0 = 0.
          dz = 0.
          if(idx.eq.0) Au=sdep
          do nn=1,ne
            ncn2 = ncn -1 + min0(1,nen(nn,ncn))
            deptest = amax1(sdep(numsideT(1,nn)),sdep(numsideT(2,nn)),sdep(numsideT(3,nn)),sdep(numsideT(ncn2,nn)))
            if(deptest.le.depmin) cycle
            zncn = 1./float(ncn2)
            DO J=1,ncn2
              MR = NEN(NN,J)
              gamma(mr) = gamma(mr) + area(nn)*zncn
              rhv(mr) = rhv(mr) + eta(nn)*area(nn)*zncn
              do kk=1,npv
                Rn0(kk,mr) = Rn0(kk,mr) + Aiz(kk,nn)*area(nn)*zncn  !w
              enddo
            enddo
          enddo

          Aiz(:,1:nsides) = Rn0  !w
          Rn0 = 0.
          do nn=1,ne
            ncn2 = ncn -1 + min0(1,nen(nn,ncn))
            zncn = 1./float(ncn2)
            MR0 = numsideT(ncn2,nn)
            sn0=(iside(2,mr0)-2*nn+iside(1,mr0))/(iside(2,mr0)-iside(1,mr0))
            DNx1 =  sn0*sdy(mr0)
            DNy1 = -sn0*sdx(mr0)
            DO J=1,ncn2
              MR = NEN(NN,J)
!              gamma(mr) = gamma(mr) + area(nn)*zncn
!              rhv(mr) = rhv(mr) + eta(nn)*area(nn)*zncn
              MR2 = numsideT(j,nn)
              sn2=(iside(2,mr2)-2*nn+iside(1,mr2))/(iside(2,mr2)-iside(1,mr2))
              DNx2 =  sn2*sdy(mr2)
              DNy2 = -sn2*sdx(mr2)
              detn = dnx1*dny2 - dnx2*dny1
              do kk=1,npv
                un1 = sn0*un(kk,mr0)
                un2 = sn2*un(kk,mr2)
                uu =   (dny2*un1-dny1*un2)*area(nn)*zncn/detn
                vv = - (dnx2*un1-dnx1*un2)*area(nn)*zncn/detn
                Rn0(kk,mr) = Rn0(kk,mr) + uu
                dz(kk,mr)  = dz(kk,mr)  + vv
              enddo
              mr0 = mr2
              dnx1 = dnx2
              dny1 = dny2
              sn0 = sn2
            enddo
          enddo

          do j=1,np
            if(nbc(j).lt.0) then
              rhv(j) = spec(-nbc(j))
            elseif(gamma(j).gt.0) then
              rhv(j) = rhv(j)/gamma(j)
            else
              rhv(j) = xyz(j,3)
            endif
            do kk=1,npv
              if(nbc(j).eq.1) then
                Rn0(kk,j) = 0.
                dz(kk,j)  = 0.
                Aiz(kk,j) = 0.
              elseif(nbc(j).gt.1.and.nbc(j).lt.7000) then
                js = nbc(j)
                Rn0(kk,j) = spec(js)
                dz(kk,j)  = spec(js+1)
                Aiz(kk,j) = 0.
              elseif(gamma(j).gt.0.) then
                Rn0(kk,j) = Rn0(kk,j)/gamma(j)
                dz(kk,j)  = dz(kk,j)/gamma(j)
                Aiz(kk,j) = Aiz(kk,j)/gamma(j)
              else
                Rn0(kk,j) = 0.
                dz(kk,j)  = 0.
                Aiz(kk,j) = 0.
              endif
            enddo
          enddo

          if(nopttype.eq.1) then !write nodal eta for harmonic analysis
            write(22) TET,np,ndf0,ndf0
            write(22) (rhv(j),j=1,np)

          elseif(nopttype.eq.2) then !write nodal eta,u,v for harmonic analysis
            write(22) TET,np,np,npv
            write(22) (rhv(j),j=1,np)
            do k=1,npv
              write(22) (Rn0(k,j),j=1,np)
              write(22) (dz(k,j),j=1,np)
            enddo

          elseif (noptcount.eq.0) then
          ! first time through start new file

            if(npv.gt.1) then
              write(22,*)'VARIABLES="X" "Y" "Z" "U" "V" "W" '
            elseif(iwind.ge.2) then
              write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" "Ps" "Tx" "Ty" '
            else
              write(22,*)'VARIABLES="X" "Y" "Z" "ETA" "U" "V" '
            endif
            if(npv.gt.1) then
              write(22,50) np*npv,max(ne,ne*(npv-1))
            elseif(ncn.eq.3) then
              write(22,20) np,ne
            elseif(ncn.eq.4) then
              write(22,30) np,ne
            endif
            write(22,*) ' SOLUTIONTIME=',TET

            ! write the data
            if(npv.gt.1) then
              do k=1,npv
                do j=1,np
                  zlev = rhv(j)+amax1(0.,(rhv(j)-xyz(j,3)))*zdep(k)
                  u = Rn0(k,j)
                  v = dz(k,j)
                  w = Aiz(k,j)
                  write(22,'(3(1x,e14.6),3(1x,f10.5))') xyz(j,1),xyz(j,2),zlev,u,v,w
                enddo
              enddo

              do j=1,ne
              ! and the elements
                if(ncn.eq.4.and.NEN(j,ncn).eq.0) then
                  NENtmp = NEN(j,ncn-1)
                else
                  NENtmp = NEN(j,ncn)
                endif
                do kv=1,npv-1
                  write(22,'(8(1x,I7))') (NEN(j,k)+kv*np,k=1,3),NENtmp+kv*np,(NEN(j,k)+(kv-1)*np,k=1,3),NENtmp+(kv-1)*np
                enddo
              enddo
            
            else
              if(iwind.ge.2) then
                do j=1,np
                  u = Rn0(1,j)
                  v = dz(1,j)
                  write(22,'(3(1x,e14.6),4(1x,f10.5),2(1pe12.4))') xyz(j,1),xyz(j,2),xyz(j,3),rhv(j),u,v,Ps(j),tsx(j),tsy(j)
                enddo
              else
                do j=1,np
                  u = Rn0(1,j)
                  v = dz(1,j)
                  write(22,'(3(1x,e14.6),3(1x,f10.5))') xyz(j,1),xyz(j,2),xyz(j,3),rhv(j),u,v
                enddo
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
            endif
          else 
          ! append file
            if(npv.gt.1) then
              write(22,51) 
            elseif(ncn.eq.3) then
              if(neMB.gt.0.or.(icase.ge.1.and.icase.le.4)) then
                write(22,22)
              else 
                write(22,21)
              endif
            elseif(ncn.eq.4) then
              if(neMB.gt.0.or.(icase.ge.1.and.icase.le.4)) then
                write(22,32)
              else 
                write(22,31)
              endif
            endif
            write(22,*) ' SOLUTIONTIME=',TET

            ! write the data
            if(npv.gt.1) then
              do k=1,npv
                do j=1,np
                  zlev = rhv(j)+amax1(0.,(rhv(j)-xyz(j,3)))*zdep(k)
                  u = Rn0(k,j)
                  v = dz(k,j)
                  w = Aiz(k,j)
                  write(22,'(e14.6,3(1x,f10.5))') zlev,u,v,w
                enddo
              enddo
            else
              if(neMB.gt.0.or.(icase.ge.1.and.icase.le.4)) then
                do j=1,np
                  u = Rn0(1,j)
                  v = dz(1,j)
                  write(22,'(e14.6,3(1x,f10.5))') xyz(j,3),rhv(j),u,v
                enddo
              else
                if(iwind.ge.2) then
                  do j=1,np
                    u = Rn0(1,j)
                    v = dz(1,j)
                    write(22,'(4(1x,f10.5),2(1pe12.4))') rhv(j),u,v,ps(j),tsx(j),tsy(j)
                  enddo
                else
                  do j=1,np
                    u = Rn0(1,j)
                    v = dz(1,j)
                    write(22,'(3(1x,f10.5))') rhv(j),u,v
                  enddo
                endif
              endif
            endif
          endif
        
        CASE(3)  !eta, u, v interpolated to centroid. 2D only.

! *** generate ut for output.
          Rn0(:,1:nsides) = un(:,1:nsides)
          call  GlobalQInterp ()

          if(nopttype.eq.1) then !write nodal eta for harmonic analysis
            write(22) TET,ne,ndf0,ndf0
            write(22) (eta(j),j=1,ne)

          elseif(nopttype.eq.2) then !write nodal eta,u,v for harmonic analysis
            write(22) TET,ne,nsides,npv
            write(22) (eta(j),j=1,ne)
            do k=1,npv
              write(22) (( un(k,j)*sdy(j)+Rn0(k,j)*sdx(j)),j=1,nsides)
              write(22) ((-un(k,j)*sdx(j)+Rn0(k,j)*sdy(j)),j=1,nsides)
            enddo

          elseif (noptcount.eq.0) then
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


!            write(*,*) 'In output 1, neMB,npMB=',neMB,npMB
      
            if(npv.gt.1) then
              write(22,"('ZONE N=',i7,' E=',i7 )" ) (np+nsides)*npv,ne*(npv-1)
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

            ! write the data
            write(22,'(6(1x,1pe15.7))') ((xyz(j,1),j=1,np),(sxy(1,j),j=1,nsides),k=1,max(2,npv)) !x
            write(22,'(6(1x,1pe15.7))') ((xyz(j,2),j=1,np),(sxy(2,j),j=1,nsides),k=1,max(2,npv)) !y

            rhv = 0.
            do nn=1,ne
              ncn2 = ncn -1 + min0(1,nen(nn,ncn))
              deptest = amax1(sdep(numsideT(1,nn)),sdep(numsideT(2,nn)),sdep(numsideT(3,nn)),sdep(numsideT(ncn2,nn)))
              if(deptest.le.depmin) cycle
              zncn = 1./float(ncn2)
              DO J=1,ncn2
                MR = NEN(NN,J)
                gamma(mr) = gamma(mr) + area(nn)*zncn
                rhv(mr) = rhv(mr) + eta(nn)*area(nn)*zncn
              enddo
            enddo
            do j=1,np
              if(nbc(j).lt.0) then
                rhv(j) = spec(-nbc(j))
              elseif(gamma(j).gt.0) then
                rhv(j) = rhv(j)/gamma(j)
              else
                rhv(j) = xyz(j,3)
              endif
            enddo
                       
            if(npv.gt.1) then
              if(izcoord.eq.0.or.izcoord.eq.2) then
                do k=1,npv
                  write(22,'(6(1x,e14.6))') (-xyz(j,3)*zdep(k),j=1,np),(sdep(j)*zdep(k),j=1,nsides) !z
                enddo
              elseif(izcoord.eq.1) then
                do k=1,izgrid-1
                  write(22,'(6(1x,e14.6))') (-amax1(xyz(j,3),zdep(izgrid))*zdep(k),j=1,np),&
                            (-amax1(refdep(j),zdep(izgrid))*zdep(k),j=1,nsides) !z
                enddo
                do k=izgrid,npv
                  write(22,'(6(1x,e14.6))') (amax1(xyz(j,3),zdep(k)),j=1,np),&
                             (amax1(refdep(j),zdep(k)),j=1,nsides) !z
                enddo
              endif
            elseif(npMB.gt.0) then
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(-sdep(j),j=1,nsides) !surface
              write(22,'(6(1x,e14.6))') (xyz(j,3),j=1,np),(-sdep(j),j=1,nsides) !etaMB
            else
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(-sdep(j),j=1,nsides) !eta
              write(22,'(6(1x,e14.6))') (xyz(j,3),j=1,np),(refdep(j),j=1,nsides) !z
            endif

            zero = 0.0
            write(22,'(6(1x,e14.6))') ((Rn0Gx(j,k),j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)
            if(npv.eq.1) then
              if(neMB.gt.0) then
                write(22,'(6(1x,e14.6))') ((Rn0Gx(j,k),j=1,np),( Qsn(j),j=1,nsides),k=1,npv)
              else
                write(22,'(6(1x,e14.6))') ((Rn0Gx(j,k),j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)
              endif
            endif
            write(22,'(6(1x,e14.6))') ((Rn0Gy(j,k),j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)
            if(npv.eq.1) then
              write(22,'(6(1x,e14.6))') ((Rn0Gy(j,k),j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)
            endif

            if(npv.gt.1) then
              write(22,'(6(1x,e14.6))') ((Aiz(k,j),k=1,npv-1),j=1,ne)
            endif

            npvm = max(1,npv-1)
            if(nson.gt.0) then
              write(22,'(6(1x,e14.6))') ((qp(j),k=1,npvm),j=1,ne)
            elseif(nsed.gt.0) then
              write(22,'(6(1x,e14.6))') ((cc(k,j),k=1,npvm),j=1,ne)
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
                do kv=1,npv-1
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
!              deptest = amax1(sdep(numsideT(1,nn)),sdep(numsideT(2,nn)),sdep(numsideT(3,nn)),sdep(numsideT(ncn2,nn)))
!              if(deptest.le.depmin) cycle
!              zncn = 1./float(ncn2)
              DO J=1,ncn2
                MR = NEN(NN,J)
                gamma(mr) = gamma(mr) + area(nn) !*zncn
                rhv(mr) = rhv(mr) + eta(nn)*area(nn) !*zncn
              enddo
            enddo
            do j=1,np
              if(nbc(j).lt.0) then
                rhv(j) = spec(-nbc(j))
              elseif(gamma(j).gt.0) then
                rhv(j) = rhv(j)/gamma(j)
              else
                rhv(j) = xyz(j,3)
              endif
            enddo
                       
            if(npv.gt.1) then
              if(izcoord.eq.0.or.izcoord.eq.2) then
                do k=1,npv
                  write(22,'(6(1x,e14.6))') (-xyz(j,3)*zdep(k),j=1,np),(sdep(j)*zdep(k),j=1,nsides) !z
                enddo
              elseif(izcoord.eq.1) then
                do k=1,izgrid-1
                  write(22,'(6(1x,e14.6))') (-amax1(xyz(j,3),zdep(izgrid))*zdep(k),j=1,np),&
                            (-amax1(refdep(j),zdep(izgrid))*zdep(k),j=1,nsides) !z
                enddo
                do k=izgrid,npv
                  write(22,'(6(1x,e14.6))') (amax1(xyz(j,3),zdep(k)),j=1,np),&
                             (amax1(refdep(j),zdep(k)),j=1,nsides) !z
                enddo
              endif
            elseif(npMB.gt.0) then
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(-sdep(j),j=1,nsides) !eta
              write(22,'(6(1x,e14.6))') (xyz(j,3),j=1,np),(-sdep(j),j=1,nsides) !etaMB
            else
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np),(-sdep(j),j=1,nsides) !eta
              write(22,'(6(1x,e14.6))') (xyz(j,3),j=1,np),(refdep(j),j=1,nsides) !z
            endif

            zero = 0.0
            write(22,'(6(1x,e14.6))') ((Rn0Gx(j,k),j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)
            if(npv.eq.1) then
              if(neMB.gt.0) then
                write(22,'(6(1x,e14.6))') ((Rn0Gx(j,k),j=1,np),( Qsn(j),j=1,nsides),k=1,npv)
              else
                write(22,'(6(1x,e14.6))') ((Rn0Gx(j,k),j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)
              endif
            endif
            write(22,'(6(1x,e14.6))') ((Rn0Gy(j,k),j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)
            if(npv.eq.1) then
              write(22,'(6(1x,e14.6))') ((Rn0Gy(j,k),j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)
            endif
!            write(22,'(6(1x,e14.6))') ((zero,j=1,np),(( un(k,j)*sdy(j)+ut(k,j)*sdx(j)),j=1,nsides),k=1,npv)
!            write(22,'(6(1x,e14.6))') ((zero,j=1,np),((-un(k,j)*sdx(j)+ut(k,j)*sdy(j)),j=1,nsides),k=1,npv)

            if(npv.gt.1) then
              write(22,'(6(1x,e14.6))') ((Aiz(k,j),k=1,npv-1),j=1,ne)
            endif

            npvm = max(1,npv-1)
            if(nson.gt.0) then
              write(22,'(6(1x,e14.6))') ((qp(j),k=1,npvm),j=1,ne)
            elseif(nsed.gt.0) then
              write(22,'(6(1x,e14.6))') ((cc(k,j),k=1,npvm),j=1,ne)
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
            if(icoord.gt.0) then
!              OPEN(UNIT=19,file=trim(OutputPath)//GridFileName,status='OLD',FORM='UNFORMATTED')      
!              READ (19) NEx,NTYPEx,NPx,NPRx,NCNx
!              READ (19) ((NENx,K=1,NCNx),J=1,NEx),(IEx,J=1,NEx),(rhv(J),J=1,NPx),(gamma(J),J=1,NPx) 
!              close (19)

              x00 = -long0 ! 0. !offsets: from input file?
              y00 = -lat0  ! 0.
              bigrI = rad2deg/bigr
              bigrcI = bigrI/cos(deg2rad*lat0) 
              do j=1,np
                rhv(j) = xyz(j,1)*bigrcI 
                gamma(j) = xyz(j,2)*bigrI
              enddo
              write(22,'(6(1x,e14.6))') (rhv(j),j=1,np)
              write(22,'(6(1x,e14.6))') (gamma(j),j=1,np)
              write(22,'(6(1x,e14.6))') (xyz(j,3),j=1,np)
            else
              x00 = 0.
              y00 = 0.
              write(22,'(6(1x,e14.6))') (xyz(j,1),j=1,np)
              write(22,'(6(1x,e14.6))') (xyz(j,2),j=1,np)
              write(22,'(6(1x,e14.6))') (xyz(j,3),j=1,np)
            endif

            Rn0 = 0.
            do j=1,ne
              ncn2 = ncn -1 + min0(1,nen(j,ncn))
              topomin = amin1(refdep(numsideT(1,j)),refdep(numsideT(2,j)),&
                            refdep(numsideT(3,j)),refdep(numsideT(ncn2,j)))
              deptest = etamax(j) - topomin
              if(deptest.gt.depmin) then
                Rn0(1,j)= (spdmax(numsideT(1,j))+spdmax(numsideT(2,j))+spdmax(numsideT(3,j)))/3.
              endif
            enddo

            write(22,'(6(1x,i5))') (IECode(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (etamax(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (Rn0(1,j),j=1,ne)
            write(22,'(6(1x,e14.6))') (eta(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (eqtide(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (t1(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (twet(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (twet2(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (twet3(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (etalast(j),j=1,ne)

            ! and the elements
            do j=1,ne
              if(NEN(j,ncn).eq.0) then
                NENtmp = NEN(j,ncn-1)
              else
                NENtmp = NEN(j,ncn)
              endif
              write(22,*) (NEN(j,k),k=1,ncn-1),NENtmp
            enddo
          
            ! write the data at the centroid
            open(unit=23,file=trim(OutputPath)//'MaxCentroid.dat',status='unknown')

            do nn=1,ne
              ncn2 = ncn -1 + min0(1,nen(nn,ncn))
              zncn = 1./float(ncn2)
              xc = 0.
              yc = 0.
              zc = 0.
              speed = 0.
              if(icoord.gt.0) then
                DO k=1,ncn2
                  j=nen(nn,k)
                  xc = xc + rhv(j)*zncn
                  yc = yc + gamma(j)*zncn
                  zc = zc + xyz(j,3)*zncn
                  speed = speed + spdmax(numsideT(k,nn))*zncn
                enddo
              else
                DO k=1,ncn2
                  j=nen(nn,k)
                  xc = xc + xyz(j,1)*zncn
                  yc = yc + xyz(j,2)*zncn
                  zc = zc + xyz(j,3)*zncn
                  speed = speed + spdmax(numsideT(k,nn))*zncn
                enddo
              endif
              depth = etamax(nn)-zc
              if(depth.le.0.) then
                speed = 0.
                depth = 0.
              endif
              write(23,'(2(1x,1Pe14.7),10(1x,1Pe12.5))') xc+x00,yc+y00,zc,depth,&
                            etamax(nn),speed,t1(nn),twet(nn),twet2(nn),twet3(nn),etalast(nn)
            enddo
            close(23)


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
              topomin = amin1(refdep(numsideT(1,j)),refdep(numsideT(2,j)),&
                            refdep(numsideT(3,j)),refdep(numsideT(ncn2,j)))
              deptest = etamax(j) - topomin
              if(deptest.gt.depmin) then
                rhv(j)= (spdmax(numsideT(1,j))+spdmax(numsideT(2,j))+spdmax(numsideT(3,j)))/3.
              endif
            enddo

            write(22,'(6(1x,i5))') (IECode(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (etamax(j),j=1,ne)
!            write(22,'(6(1x,e14.6))') (rhv(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (eta(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (eqtide(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (t1(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (twet(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (twet2(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (twet3(j),j=1,ne)
            write(22,'(6(1x,e14.6))') (etalast(j),j=1,ne)
          endif
       
        END SELECT

        close(22)

        if(npv.gt.1) then
          ! first time through start new file
          do j=1,jUprofile

            js = numsideT(1,tsUnodes(j))
            if(js.eq.0) cycle

            npv0 = npv
            n1 = iside(1,js)
            n2 = iside(2,js)
            if(n2.eq.0) then
              etaside = eta(n1)
            else
              etaside = amax1(eta(n1),eta(n2))
            endif

            if(izcoord.eq.0) then ! sigma
              do k=1,npv
                zz(k) = etaside + sdep(js)*zdep(k)
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
                zz(k) = refdep(js) + sdep(js)*(1.+ 0.5*(zdep(k)+zdep(k+1))) !ref to MSL
              enddo
            endif

            if (noptcount.eq.0) then
              n23 = 23  !+j-1
              write(cseq,"(I2.2)") j
              open(unit=n23,file=trim(OutputPath)//'profile'//trim(cseq)//'.dat',status='unknown')
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
              open(unit=n23,file=trim(OutputPath)//'profile'//trim(cseq)//'.dat',status='old',position='append')
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

        if((nsed.gt.0).and.(npvc.gt.1)) then
          if (noptcount.eq.0) then
          ! first time through start new file
            open(unit=23,file=trim(OutputPath)//'profileSed.dat',status='unknown')
            write(23,*) 'VARIABLES="z" "C" '
            write(23,*) 'ZONE'
            write(23,*) ' SOLUTIONTIME=',TET

            js = jCProfile
            js2 = jUProfile
            do j=1,npvc
!           ! write the data
              write(23,*) sdep(js2)*(1.+zdepC(j)), cc(j,js)
            enddo
            write(23,*) sdep(js2)*(1.+zdep(npv)), cc(npvc+1,js)
            if(nsbc.ge.3) then
              write(23,*) sdep(js2)*(1.+zdep(npv)), cc(npvc+2,js)
            endif
          else 
          ! append file
            open(unit=23,file=trim(OutputPath)//'profileSed.dat',status='old',position='append')
            write(23,*) 'ZONE'
            write(23,*) ' SOLUTIONTIME=',TET

            js = jCProfile
            js2 = jUProfile
            do j=1,npvc
!           ! write the data
              write(23,*) sdep(js2)*(1.+zdepC(j)), cc(j,js)
            enddo
            write(23,*) sdep(js2)*(1.+zdep(npv)), cc(npvc+1,js)
            if(nsbc.ge.3) then
              write(23,*) sdep(js2)*(1.+zdep(npv)), cc(npvc+2,js)
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
              open(unit=n23,file=trim(OutputPath)//'profileCon'//trim(cseq)//'.dat',status='unknown')
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

              js = tsSnodes(j)
!              cdep=amin1(xyz(nen(js,1),3),xyz(nen(js,2),3),xyz(nen(js,3),3),xyz(nen(js,ncn),3))
              ncn2 = ncn -1 + min0(1,nen(js,ncn))
              cdep = 0. 
              do jn=1,ncn2
                cdep = cdep + xyz(nen(js,jn),3)
              enddo
              cdep = eta(js)-cdep/float(ncn2)
!              write(23,*) 'iOP,npvc=',iOPsol,npvc,js,cdep
!             ! write the data
              if(iOPsol.eq.0) then
                do k=1,npvc
                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.1) then
                do k=1,npvc
                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), PSU(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.2) then
                do k=1,npvc
                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), TC(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.3) then
                do k=1,npvc
                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), PSU(k,js), TC(k,js)
                enddo
                close(n23)
              endif
            enddo
          else 
          ! append file
            do j=1,jSprofile
              n23 = 23 !+j-1
              write(cseq,"(I2.2)") j
              open(unit=n23,file=trim(OutputPath)//'profileCon'//trim(cseq)//'.dat',status='old',position='append')
!              open(unit=23,file='profileCon.dat',status='old',position='append')
              write(23,*) 'ZONE'
              write(23,*) ' SOLUTIONTIME=',TET

              js = tsSnodes(j)
              ncn2 = ncn -1 + min0(1,nen(js,ncn))
              cdep = 0. 
              do jn=1,ncn2
                cdep = cdep + xyz(nen(js,jn),3)
              enddo
              cdep = eta(js)-cdep/float(ncn2)
!             ! write the data
              if(iOPsol.eq.0) then
                do k=1,npvc
                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.1) then
                do k=1,npvc
                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), PSU(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.2) then
                do k=1,npvc
                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), TC(k,js)
                enddo
                close(n23)
              elseif(iOPsol.eq.3) then
                do k=1,npvc
                  write(23,'(6(1x,e14.6))') cdep*zdepC(k), sigt(k,js), PSU(k,js), TC(k,js)
                enddo
                close(n23)
              endif
            enddo
          endif
        endif

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
!*****************************************************************
