
      Program RiCOMHarmonicAnalysis

!     This version of the tidal heights analysis program has
!     been adapted for the analysis of hourly output from
!     RiCOM model runs.

!     all the er* arrays contain residual mean square errors
!     or standard errors of estimate of y

      implicit none

      real*8, parameter :: PI=3.1415926536, TWOPI=2.*PI, fac=pi/180.
      integer :: i,ii,ia,ib,ik,j,jj,jk,jkm,jm,jmkm,k,kk,kjm,kv,m,mm,m2,m2m2,m1m1,nn
      integer :: ntime,nskip,nskip1,nc,ninf,ihr0,nprof,np,npv, nson,nsed,nsol
      integer :: nlev, npvx, nfile, ncon, ncon2, nph, nphu,ndf,kkount,infz,nii
      integer :: iex,nex,ntypex,npx,nprx,node,nentmp
      integer :: ne,nnodes,nsides,ncn,izcoord,npvgrd,neMB,neqtide
      real*8 tet, arg
      real a(400)  !,rhs(ncon2),rhu(ncon2),rhv(ncon2)
      real :: x,y,ab,ba1,wj,wk,wd,ws,sind,sins,td,ts,gg,errz0,aa,pp
      real :: ersf,eruf,ervf,ain,phin,dep
      integer, allocatable ::  npt(:),nen(:,:),numsideT(:,:)
      real, allocatable ::  xp(:),yp(:),zp(:),zdep(:)
      real, allocatable ::  un(:,:),ut(:,:)
      real, allocatable ::  sdx(:),sdy(:)
      real, allocatable ::  zeta(:),ubar(:,:),vbar(:,:),zs(:)
      real, allocatable ::  erz(:),erub(:),ervb(:)
      real, allocatable ::  ampz(:,:),phz(:,:),ampub(:,:)
      real, allocatable ::  phub(:,:),ampvb(:,:),phvb(:,:)
      real, allocatable ::  freq(:),f(:),vpu(:)
      real, allocatable ::  rhs(:),rhu(:),rhv(:)
      real, allocatable ::  rhs2(:,:),rhu2(:,:),rhv2(:,:)
      real, allocatable ::  amaj(:),amin(:),g(:),ainc(:)
      real, allocatable ::  freqref(:),freqinf(:),finf(:),vpuinf(:),ratio(:),phdif(:)
      integer :: istat,iarray
      character*4, allocatable :: name(:),namref(:),naminf(:) !  ncon)

      integer :: numarg, iargc, nconx=20
      real*8 :: fnx(20),nux(20),freqx(20)
!      character*4 :: namex(20)
      integer :: isp(50)
      character*4 icn
      character*256 fname,fnamedata
      character(18) :: AnalysisTime
      logical(4) resOK, openfile

      write(*,*) ' Start analysis'
!      PI=3.1415926536
!      TWOPI=2.*PI
!      fac=pi/180.

!     open control file

      numarg=iargc()
      
      if(numarg.ge.1) then
        i = 1
        call getarg(i,fname)
        open(unit=3,file=fname,status='old',iostat=istat)
        if(istat.ne.0) then
          write(*,*), 'Unable to open input control file'
          stop
        endif
      else
        ResOK = OpenFile(3,'Open Control File',fname,"Control file(*.dat),*.dat;All files(*.*),*.*;")

        if(.not.resOK) then
          stop
        endif
      endif

!     read control file

      read(3,*) ntime,nskip,nc,ninf,ihr0,nprof,np
!     ntime = total number of hourly time steps
!     nskip = number to skip initially
!     nc = number of constituents
!     ihr0 = reference hour for first record
!      nlev = level number to be processed
!     ninf = number of inferred constituents
!     nprof = 0 analyse horzontal layer
!           = 1 analyse profiles
!           = 2 analyse points and only output PEST file.
!     np = number of profiles to analyze (or vertical level if nprof=0)
      write(*,*) ' number of hourly time steps =', ntime
      write(*,*) ' number of steps to skip initially =', nskip
      write(*,*) ' number of constituents =', nc
      write(*,*) ' Index for hour 0 of analysis =', ihr0
      if(nprof.eq.0) then
        nlev = np
        write(*,*) ' Index for vertical level =', nlev
      else
        write(*,*) ' number of profiles =', np
        ALLOCATE (npt(np), STAT = istat ) 
        if(istat.ne.0) then
          write(*,*) 'FATAL ERROR: Cannot allocate npt storage array'
          stop
        endif
        read(3,*) (npt(i),i=1,np)
!        if(nprof.eq.2) then  !read sigma data
!                ResOK = OpenFile(5,'Open sigma File',fname,"data file(*.dat),*.dat;All files(*.*),*.*;")
!          read(5,*) npvx
!          ALLOCATE (zs(npvx), STAT = istat ) 
!          if(istat.ne.0) then
!            write(*,*) 'FATAL ERROR: Cannot allocate zs storage array'
!            stop
!          endif
!          do j=1,npvx
!            read(5,*) zs(j)
!          enddo
!        endif
      endif

      nfile=ntime
      nskip1=nskip+1
      ncon = nc
      ncon2 = 2*ncon

      ALLOCATE (name(ncon),freq(ncon),f(ncon),vpu(ncon),&
             namref(ncon),freqref(ncon),finf(ncon),vpuinf(ncon),&
             naminf(ncon),freqinf(ncon),ratio(ncon),phdif(ncon),STAT = istat ) 
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate freq storage arrays',istat
        stop
      endif

      do i=1,nc
        read(3,'(a4,6x,3f10.7)') name(i),freq(i),f(i),vpu(i)
      enddo
!      assume the frequency is in cycles/hr and vpu is in degrees

!    read in inference info if requested
      if(ninf.gt.0) then
        do i=1,ninf
          read(3,798) namref(i),freqref(i),naminf(i),freqinf(i),&
                      ratio(i),phdif(i),finf(i),vpuinf(i)
        end do !i
      endif
798   format(a4,6x,f10.7,a4,6x,f10.7,4f7.3)

!     read time slices from the binary data file

      read(3,'(a)') fnamedata
      open(unit=20,file=fnamedata,status='old',form='unformatted')
      
      read(20) AnalysisTime
      write(*,*) 'AnalysisTime= ',AnalysisTime
      
! *** convert phases for the start of the time series for known constituents
      call ConvertPhases(AnalysisTime,ncon,f,vpu,name,freq)
      
! *** convert phases for the start of the time series for known inf constituents
      if(ninf.gt.0) then
        call ConvertPhases(AnalysisTime,ninf,finf,vpuinf,naminf,freqinf)
      endif
      
      read(20)  !skip windfilename

      read(20) ne,nnodes,nsides,npv,ncn,izcoord   !nph,nphu,npv  
      read(20) nson,nsed,nsol,neMB,neqtide
      nph = ne
      nphu = nsides
      npv = max(npv,1)
      if(nprof.eq.0) then
        np = nph !**
        nn = nph*npv
      else
        nn = np*npv
      endif
! just keep going      rewind (20)

      write(*,*) 'nph,nphu,nnodes,npv=',nph,nphu,nnodes,npv
      write(*,*) 'ncon,np*npv=',ncon,nn
        
      ALLOCATE ( xp(nnodes), yp(nnodes), zp(nnodes), nen(ne,ncn), &
          sdx(nsides),sdy(nsides),numsideT(ncn,ne), &
          un(npv,nphu),ut(npv,nphu), &
          zeta(nph),ubar(nph,npv),vbar(nph,npv), &
          ampz(nn,ncon),phz(nn,ncon),ampub(nn,ncon), &
          phub(nn,ncon),ampvb(nn,ncon),phvb(nn,ncon), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate ncon main storage arrays',istat
        stop
      endif
      
! *** read coordinates
      read(20) (xp(j),j=1,nnodes)
      read(20) (yp(j),j=1,nnodes)
      read(20) (zp(j),j=1,nnodes)
    
      if(npv.gt.1) then
        npvgrd = npv + min(1,max(0,izcoord-1))
        ALLOCATE ( zdep(npvgrd), STAT = istat )
        if(istat.ne.0) then
          write(*,*) 'FATAL ERROR: Cannot allocate zdep',istat
          stop
        endif
        read(20) (zdep(k),k=1,npvgrd)    !write  z-grid (m)
      endif

! *** read element list
      read(20) ((nen(j,k),k=1,ncn),j=1,ne)
      read(20) (sdx(j),j=1,nsides),(sdy(j),j=1,nsides)
      read(20) ((numsideT(i,j),i=1,ncn),j=1,ne)

      write(*,*) 'ncon2=',ncon2

      ALLOCATE (rhs(ncon2),rhu(ncon2),rhv(ncon2),  & 
          rhs2(nn,ncon2),rhu2(nn,ncon2),rhv2(nn,ncon2), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate ncon2 main storage arrays',istat
        stop
      endif

      ALLOCATE (erz(nn),erub(nn),ervb(nn), &
          amaj(nn),amin(nn),g(nn),ainc(nn), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate nn main storage arrays',istat
        stop
      endif

      zeta = 0.
      ubar = 0.
      vbar = 0.
      rhs2=0.
      rhu2=0.
      rhv2=0.
      ampz=0.
      phz=0.
      ampub=0.
      phub=0.
      ampvb=0.
      phvb=0.
      amaj=0.
      amin=0.
      g=0.
      ainc=0.
      
      if(nskip.gt.0) then
        do k=1,nskip
          read(20) TET  !,nph,nphu,npv
          read(20)
          read(20)
          read(20)    
          if(npv.gt.1) then
            read(20)    
          endif
          if(k.gt.1.and.neMB.gt.0) then
            read(20)    
          endif
          if(nson.gt.0) then
            read(20)
          endif
          if(nsed.gt.0) then
            read(20)
          endif
          if(nsol.gt.0) then
            read(20)
          endif
          if(neqtide.gt.0) then
            read(20)
          endif
       enddo
      endif

      do k=nskip1,ntime
!        l=1
        read(20) TET 
        read(20) (zeta(j),j=1,nph)
! *** read velocity (un,ut)
        read(20) ((un(kv,i),kv=1,npv),i=1,nphu)
        read(20) ((ut(kv,i),kv=1,npv),i=1,nphu)
        if(npv.gt.1) then
          read(20)    
        endif
! *** average (u,v) to centroid
        do j=1,ne
          do jj=1,ncn
            ii = numsideT(jj,j)
            do kv=1,npv
              ubar(j,kv) = ubar(j,kv) + un(kv,ii)*sdy(ii) + ut(kv,ii)*sdx(ii)
              vbar(j,kv) = vbar(j,kv) - un(kv,ii)*sdx(ii) + ut(kv,ii)*sdy(ii)
            enddo
          enddo
        enddo
        ubar = ubar/float(ncn)
        vbar = vbar/float(ncn)
        if(k.gt.1.and.neMB.gt.0) then
          read(20)    
        endif
        if(nson.gt.0) then
          read(20)
        endif
        if(nsed.gt.0) then
          read(20)
        endif
         if(nsol.gt.0) then
          read(20)
        endif
        if(neqtide.gt.0) then
          read(20)
        endif
        
        do j=1,nc
          jj=nc+j
          arg=freq(j)*twopi*(ihr0+k)
          if(nprof.eq.0) then
            do kv=1,npv  
              do i=1,nph
                ik = (kv-1)*nph + i
                rhu2(ik,j)=rhu2(ik,j)+ubar(i,kv)*cos(arg)
                rhu2(ik,jj)=rhu2(ik,jj)+ubar(i,kv)*sin(arg)
                rhv2(ik,j)=rhv2(ik,j)+vbar(i,kv)*cos(arg)
                rhv2(ik,jj)=rhv2(ik,jj)+vbar(i,kv)*sin(arg)
                rhs2(ik,j)=rhs2(ik,j)+zeta(i)*cos(arg)
                rhs2(ik,jj)=rhs2(ik,jj)+zeta(i)*sin(arg)
              enddo
            enddo
          else
            do i=1,np 
              do kv=1,npv  
                ik = (i-1)*npv + kv
                rhs2(ik,j)=rhs2(ik,j)+zeta(npt(i))*cos(arg)
                rhs2(ik,jj)=rhs2(ik,jj)+zeta(npt(i))*sin(arg)
                rhu2(ik,j)=rhu2(ik,j)+ubar(npt(i),kv)*cos(arg)
                rhu2(ik,jj)=rhu2(ik,jj)+ubar(npt(i),kv)*sin(arg)
                rhv2(ik,j)=rhv2(ik,j)+vbar(npt(i),kv)*cos(arg)
                rhv2(ik,jj)=rhv2(ik,jj)+vbar(npt(i),kv)*sin(arg)
              enddo
            enddo
          endif
        enddo
      enddo

      close(unit=20)
      write(6,*) ' time series read'

!    set up the least square matrix
!***********************************************************************
!*  GENERATE THE MATRIX BY SUMMING THE CONTRIBUTIONS FROM EACH OF THE
!*  CONTINUOUS BLOCKS OF DATA.
!
      m=nc
      M2=2*M
      M2M2=(M2*M2+M2)/2
      M1M1=(3*M*M+M)/2+1
      MM=(M*M+M)/2
      DO 160 JK=1,M2M2
160   A(JK)=0.0

!***********************************************************************
!*  IA AND IB ARE RESPECTIVELY THE FIRST AND LAST READINGS IN A
!*  CONTINUOUS BLOCK OF DATA.

      IA=ihr0+nskip+1
      ib=ia+nfile-nskip-1
      AB=IA+IB
      BA1=IB-IA+1

!***********************************************************************
!*  TRIGONOMETRIC IDENTITIES AVOID THE SUMMATIONS REQUIRED IN
!*  CALCULATING THE MATRIX ELEMENTS.
!*  ONLY THE LOWER TRIANGLE (INCLUDING DIAGONAL) OF THE MATRIX NEED BE
!*  GENERATED.

      DO 260 J=1,M
        WJ=PI*Freq(J)
        SIND=BA1
        SINS=BA1
        IF(J.GT.1) SINS=SIN(2.*BA1*WJ)/SIN(2.*WJ)
240     DO 255 K=J,M
          JK=M2*(J-1)+K-(J*J-J)/2
          JKM=JK+M
          JM=J+M
          KJM=M2*(K-1)+JM-(K*K-K)/2
          JMKM=M2*(JM-1)+K+M-(JM*JM-JM)/2
          WK=PI*Freq(K)
          WD=WJ-WK
          WS=WJ+WK
          IF(J.EQ.K) GO TO 250
          SIND=SIN(BA1*WD)/SIN(WD)
          SINS=SIN(BA1*WS)/SIN(WS)
250       TD=COS(AB*WD)*SIND
          TS=COS(AB*WS)*SINS
          A(JK)=A(JK)+.5*(TD+TS)
          A(JMKM)=A(JMKM)+.5*(TD-TS)
          TD=SIN(AB*WD)*SIND
          TS=SIN(AB*WS)*SINS
          IF(J.EQ.K) GO TO 255
          A(JKM)=A(JKM)+.5*(-TD+TS)
255     A(KJM)=A(KJM)+.5*(TD+TS)
260   CONTINUE
!***********************************************************************
!*  S(1) IS FORCED TO ZERO BY SETTING THIS MATRIX ELEMENT TO LARGE VALUE

  300 A(M1M1)=Nfile-nskip


!    check out the matrix
      write(6,*) ' lower triangle of matrix: m2m2=',m2m2
      write(6,251) (a(jk),jk=1,m2m2)
251      format(10e13.6)
      CALL TRIANG(A,RHS,gg,M2)
      write(6,*) ' matrix condition=',gg
      IF(gg.GT.0.) GO TO 371
      WRITE(6,372) gg
  372 FORMAT(' THE LEAST SQUARES MATRIX IS ILL-CONDITIONED AND COULD NOT BE TRIANGULARIZED: g=',f10.5)
      STOP
  371 CONTINUE
!    calculate the inverse matrix diagonal element corresponding
!    to the constant term. This will enable us to calculate an estimate
!    of the standard error for this element. See pg 40 in my old
!    tidal height analysis manual.
      do 426 j=1,nc
        jj=nc+j
        rhs(j)=0.
        rhs(jj)=0.
426      continue
      rhs(1)=1.
      GG=1.005
      CALL SOLUTN(A,RHs,GG,M2)
      errz0=sqrt(rhs(1))
      write(6,*) ' errz0=',errz0

!    now create the right hand sides and solve for the elevation
!    and vertically averaged velocities
      ndf=nfile-nskip-2*nc
      nskip1=nskip+1
      do 400 i=1,nn
        erz(i)=0.
        erub(i)=0.
        ervb(i)=0.
        do 401 j=1,nc
          jj=nc+j
!          rhs2(i,j)=0.
!          rhs2(i,jj)=0.
!          rhu2(i,j)=0.
!          rhu2(i,jj)=0.
!          rhv2(i,j)=0.
!          rhv2(i,jj)=0.
!        do 402 k=nskip1,nfile
!          arg=freq(j)*twopi*(ihr0+k)
!          rhu2(i,j)=rhu2(i,j)+ubar(i,k)*cos(arg)
!          rhu2(i,jj)=rhu2(i,jj)+ubar(i,k)*sin(arg)
!          rhv2(i,j)=rhv2(i,j)+vbar(i,k)*cos(arg)
!          rhv2(i,jj)=rhv2(i,jj)+vbar(i,k)*sin(arg)
!          rhs2(i,j)=rhs2(i,j)+zeta(i,k)*cos(arg)
!402        rhs2(i,jj)=rhs2(i,jj)+zeta(i,k)*sin(arg)
          rhs(j)=rhs2(i,j)
          rhs(jj)=rhs2(i,jj)
          rhu(j)=rhu2(i,j)
          rhu(jj)=rhu2(i,jj)
          rhv(j)=rhv2(i,j)
          rhv(jj)=rhv2(i,jj)
401     continue
        CALL SOLUTN(A,RHU,GG,M2)
        CALL SOLUTN(A,RHV,GG,M2)
        CALL SOLUTN(A,RHS,GG,M2)
!    write(6,*) ' calculated vertically integrated solutions'

        kkount = 0
        do 403 j=1,nc
          jj=nc+j
!    save results for error estimate
          rhs2(i,j)=rhs(j)
          rhs2(i,jj)=rhs(jj)
          rhu2(i,j)=rhu(j)
          rhu2(i,jj)=rhu(jj)
          rhv2(i,j)=rhv(j)
          rhv2(i,jj)=rhv(jj)

!    calculate amp/phase
          x=rhs(j)
          y=rhs(jj)
          ampz(i,j)=sqrt(x*x+y*y)/f(j)
          if(ampz(i,j).lt.1.e-10) then
            phz(i,j)=0.
          else
            phz(i,j)=vpu(j)+atan2(y,x)*360./twopi
            if(phz(i,j).lt.0.) phz(i,j)=phz(i,j)+360.
            if(phz(i,j).gt.360.) phz(i,j)=phz(i,j)-360.
          end if

!    check for inference
          infz=0
          if(ninf.gt.0) then
            do k=1,ninf
              if(name(j).eq.namref(k)) then
                infz = 1
                freqref(k) = freq(j)
              endif
            end do
          endif

          if(infz.eq.1) then
            kk=k
            infz=1
            kkount=kkount+1
!    1st undo the astronomical argument correction as subroutine infer
!    assumes it hasn't been done yet
            aa=ampz(i,j)*f(j)
            pp=phz(i,j)-vpu(j)
            call infer(nfile,namref(kk),freqref(kk),aa,pp,f(j),vpu(j),&
              naminf(kk),freqinf(kk),ain,phin,finf(kk),vpuinf(kk),ratio(kk),phdif(kk))
            ain=ain/finf(kk)
            phin=phin+vpuinf(kk)
            if(phin.lt.0.) phin=phin+360.
            if(phin.gt.360.) phin=phin-360.
            ampz(i,nc+kkount)=ain
            phz(i,nc+kkount)=phin
            ampz(i,j)=aa/f(j)
            phz(i,j)=pp+vpu(j)
            if(phz(i,j).lt.0.) phz(i,j)=phz(i,j)+360.
            if(phz(i,j).gt.360.) phz(i,j)=phz(i,j)-360.
          end if !infz

          x=rhu(j)
          y=rhu(jj)
          ampub(i,j)=sqrt(x*x+y*y)/f(j)
          if(ampub(i,j).gt.1.e-8) then
            phub(i,j)=vpu(j)+atan2(y,x)*360./twopi
            if(phub(i,j).lt.0.) phub(i,j)=phub(i,j)+360.
            if(phub(i,j).gt.360.) phub(i,j)=phub(i,j)-360.
          else
            phub(i,j)=0.
          end if

!    inference
          if(ninf.gt.0.and.infz.eq.1) then
!     we already know kk and kkount from the elevations
!    1st undo the astronomical argument correction as subroutine infer
!    assumes it hasn't been done yet
            aa=ampub(i,j)*f(j)
            pp=phub(i,j)-vpu(j)
            call infer(nfile,namref(kk),freqref(kk),aa,pp,f(j),vpu(j),&
               naminf(kk),freqinf(kk),ain,phin,finf(kk),vpuinf(kk),ratio(kk),phdif(kk))
            ain=ain/finf(kk)
            phin=phin+vpuinf(kk)
            if(phin.lt.0.) phin=phin+360.
            if(phin.gt.360.) phin=phin-360.
            ampub(i,nc+kkount)=ain
            phub(i,nc+kkount)=phin
            ampub(i,j)=aa/f(j)
            phub(i,j)=pp+vpu(j)
            if(phub(i,j).lt.0.) phub(i,j)=phub(i,j)+360.
            if(phub(i,j).gt.360.) phub(i,j)=phub(i,j)-360.
          end if !ninf

          x=rhv(j)
          y=rhv(jj)
          ampvb(i,j)=sqrt(x*x+y*y)/f(j)
          if(ampvb(i,j).gt.1.e-8) then
            phvb(i,j)=vpu(j)+atan2(y,x)*360./twopi
            if(phvb(i,j).lt.0.) phvb(i,j)=phvb(i,j)+360.
            if(phvb(i,j).gt.360.) phvb(i,j)=phvb(i,j)-360.
          else
            phvb(i,j)=0.
          end if

!    inference
          if(ninf.gt.0.and.infz.eq.1) then
!     we already know kk and kkount from the elevations
!    1st undo the astronomical argument correction as subroutine infer
!    assumes it hasn't been done yet
            aa=ampvb(i,j)*f(j)
            pp=phvb(i,j)-vpu(j)
            ain=ain/finf(kk)
            phin=phin+vpuinf(kk)
            if(phin.lt.0.) phin=phin+360.
            if(phin.gt.360.) phin=phin-360.
            ampvb(i,nc+kkount)=ain
            phvb(i,nc+kkount)=phin
            ampvb(i,j)=aa/f(j)
            phvb(i,j)=pp+vpu(j)
            if(phvb(i,j).lt.0.) phvb(i,j)=phvb(i,j)+360.
            if(phvb(i,j).gt.360.) phvb(i,j)=phvb(i,j)-360.
          end if !ninf


403     continue
400   continue

!    calculate the residual mean square error
      open(unit=20,file=fnamedata,status='old',form='unformatted')
      read(20) !analysistime
      read(20) !windfilename
      read(20) !dim
      read(20) !dim
      read(20) !x
      read(20) !y
      read(20) !z
      if(npv.gt.1) then
        read(20) !zdep
      endif
      read(20) !list  
      read(20) !sdx,sdy
      read(20) !numsideT
      if(nskip.gt.0) then
        do k=1,nskip
          read(20) TET !,npx,npx,npvx
          read(20)  !eta
          read(20)  !u
          read(20)  !v
          if(k.gt.1.and.neMB.gt.0) then
            read(20)  !z moving    
          endif
          if(npv.gt.1) then
            read(20) !w
          endif
          if(nson.gt.0) then
            read(20)
          endif
          if(nsed.gt.0) then
            read(20)
          endif
          if(nsol.gt.0) then
            read(20)
          endif
          if(neqtide.gt.0) then
            read(20)
          endif
        enddo
      endif

      ubar = 0.
      vbar = 0.
      do 420 k=nskip1,nfile
        read(20) TET 
        read(20) (zeta(j),j=1,nph)
! *** read velocity (un,ut)
        read(20) ((un(kv,i),kv=1,npv),i=1,nphu)
        read(20) ((ut(kv,i),kv=1,npv),i=1,nphu)
        if(npv.gt.1) then
          read(20) !w
        endif
! *** average (u,v) to centroid
        do j=1,ne
          do jj=1,ncn
            ii = numsideT(jj,j)
            do kv=1,npv
              ubar(j,kv) = ubar(j,kv) + un(kv,ii)*sdy(ii) + ut(kv,ii)*sdx(ii)
              vbar(j,kv) = vbar(j,kv) - un(kv,ii)*sdx(ii) + ut(kv,ii)*sdy(ii)
            enddo
          enddo
        enddo
        ubar = ubar/float(ncn)
        vbar = vbar/float(ncn)
        if(k.gt.1.and.neMB.gt.0) then
          read(20)  !z moving    
        endif
        if(nson.gt.0) then
          read(20)
        endif
        if(nsed.gt.0) then
          read(20)
        endif
        if(nsol.gt.0) then
          read(20)
        endif
        if(neqtide.gt.0) then
          read(20)
        endif

        do i=1,nn
          ersf=0.
          eruf=0.
          ervf=0.
          do 421 j=1,nc
            arg=freq(j)*twopi*(ihr0+k)
            jj=nc+j
            ersf=ersf+rhs2(i,j)*cos(arg)
            ersf=ersf+rhs2(i,jj)*sin(arg)
            eruf=eruf+rhu2(i,j)*cos(arg)
            eruf=eruf+rhu2(i,jj)*sin(arg)
            ervf=ervf+rhv2(i,j)*cos(arg)
            ervf=ervf+rhv2(i,jj)*sin(arg)
421          continue
          if(nprof.eq.0) then
            nii = mod(i-1,nph)+1
            kv = (i-1)/nph + 1
          else
            ii = 1+(i-1)/npv
            kv = mod(i-1,npv)+1
            nii = npt(ii)
          endif
          erz(i)=erz(i)+(zeta(nii)-ersf)**2
          erub(i)=erub(i)+(ubar(nii,kv)-eruf)**2
          ervb(i)=ervb(i)+(vbar(nii,kv)-ervf)**2
        enddo
420   continue

      close(unit=20)

!400    continue
      do 422 i=1,nn
        erz(i)=sqrt(erz(i)/ndf)
        erub(i)=sqrt(erub(i)/ndf)
        ervb(i)=sqrt(ervb(i)/ndf)
422   continue
      write(6,*) ' completed zeta,ubar,vbar analysis'
      write(6,*) ' erz(1),eru(1),erv(1)=',erz(1),erub(1),ervb(1)
!404    continue



!    real lel file for coordinates
!      read(3,'(a)') fname
!      open(unit=25,file=fname,status='old',form='unformatted')
!      READ (25) NEx,NTYPEx,NPx,NPRx,NCNx
! *** allocate arrays
!      ALLOCATE (nen(NEx,NCNx), STAT = istat )
!      if(istat.ne.0) then
!        write(*,*) 'FATAL ERROR: Cannot allocate element storage arrays'
!        stop
!      endif

!      nen = 0
!      READ (25) ((nen(j,k),K=1,NCNx),J=1,NEx),(IEx,J=1,NEx),(ubar(J,1),J=1,NPx),&
!                 (vbar(J,1),J=1,NPx),(zeta(J),J=1,NPx)
!      close(unit=25)

!    write out the results in Tecplot format
      read(3,'(a)') fname
      if(nprof.eq.0.or.nprof.eq.1) then
        do j=1,nc
          ia = LEN_TRIM(name(j))
          open(unit=25,file=name(j)(1:ia)//fname,status='unknown',form='formatted')
!          fac=3.1415926/180.
          if(nprof.eq.0) then
            write(25,*)  'VARIABLES= "x" "y" "z" "e" "ep" "u" "up" "v" "vp" "umaj" "umin" "ph" "dir"'
            write(25,*) ' ZONE T="',name(j),'" '
            write(25,"(' N=',i7,' E=',i7 )" )  nnodes,ne  !nn
!            write(25,*)  '  F=FEPOINT ET=TRIANGLE'
            if(ncn.eq.3) then
              write(25,"(' ZONETYPE=FETRIANGLE DATAPACKING=BLOCK')" )
            elseif(ncn.eq.4) then
              write(25,"(' ZONETYPE=FEQUADRILATERAL DATAPACKING=BLOCK')" )
            endif
            write(25,"(' VARLOCATION=([4-13]=CELLCENTERED)')" )
          else
            write(25,*) 'VARIABLES= "s" "e" "ep" "u" "up" "v" "vp" "umaj" "umin" "ph" "dir"'
          endif
! ** calculate ellipse parameters
          if(nphu.gt.0) then
            call ellipse(nn,ampub(1,j),phub(1,j),ampvb(1,j),phvb(1,j),amaj,amin,g,ainc)
          endif
!          do i=1,nn
          if(nprof.eq.0) then
!              if(npv.gt.1) then
!                ii = mod(i-1,nph) + 1
!                write(25,510) ampub(i,j)*cos(fac*phub(i,j)),ampvb(i,j)*cos(fac*phvb(i,j)),&
!                                 ampz(i,j)*cos(fac*phz(i,j))
!              else
            write(25,'(6(1x,e14.6))') (xp(i),i=1,nnodes)
            write(25,'(6(1x,e14.6))') (yp(i),i=1,nnodes)
            write(25,'(6(1x,e14.6))') (zp(i),i=1,nnodes)
            write(25,'(6(1x,e14.6))') (ampz(i,j),i=1,ne)  !nn)
            write(25,'(6(1x,e14.6))') (phz(i,j),i=1,ne)  !nn)
            write(25,'(6(1x,e14.6))') (ampub(i,j),i=1,ne)  !nn)
            write(25,'(6(1x,e14.6))') (phub(i,j),i=1,ne)  !nn)
            write(25,'(6(1x,e14.6))') (ampvb(i,j),i=1,ne)  !nn)
            write(25,'(6(1x,e14.6))') (phvb(i,j),i=1,ne)  !nn)
            write(25,'(6(1x,e14.6))') (amaj(i),i=1,ne)  !nn)
            write(25,'(6(1x,e14.6))') (amin(i),i=1,ne)  !nn)
            write(25,'(6(1x,e14.6))') (g(i),i=1,ne)  !nn)
            write(25,'(6(1x,e14.6))') (Ainc(i),i=1,ne)  !nn)
!              endif
          elseif(nprof.eq.1) then  !write profiles
            do i=1,nn
              if(mod(i,npv).eq.1) then
                node = (i-1)/npv + 1
                write(25,*)  ' ZONE T="',npt(node),'"'
              endif
!              if(nprof.eq.1) then
                dep = -float(mod(i-1,npv))
!              else
!              dep = 1. !zdep(mod(i-1,npv)+1)
!              endif
              write(25,500) dep,ampz(i,j),phz(i,j),ampub(i,j),phub(i,j),ampvb(i,j),phvb(i,j),&
                          amaj(i),amin(i),g(i),Ainc(i)
            enddo
          endif
      
          if(nprof.eq.0) then  !.and.npv.le.1) then
            do jj=1,ne
              if(NEN(jj,ncn).eq.0) then
                NENtmp = NEN(jj,ncn-1)
              else
                NENtmp = NEN(jj,ncn)
              endif
              write(25,*) (NEN(jj,k),k=1,ncn) !,NENtmp
            enddo
          endif
          
          close(unit=25)
        enddo
      elseif(nprof.eq.2) then  !write PEST file
! *** write amp,phase for all constituents at each point
!        open(26,file='outHAAP.dat',status='unknown',form='formatted')
        open(26,file=fname,status='unknown',form='formatted')
        do i=1,nn
          if(mod(i-1,npv).eq.0) then
            node = (i-1)/npv + 1
            write(26,520) npt(node),((ampz(i,j),phz(i,j)),j=1,nc)
          endif
        enddo
        close(26)
! *** write real, imag amp for all constituents at each point
        open(26,file='outHARI.dat',status='unknown',form='formatted')
        do i=1,nn
          if(mod(i-1,npv).eq.0) then
            node = (i-1)/npv + 1
            write(26,520) npt(node),((ampz(i,j),phz(i,j)),j=1,nc)
          endif
        enddo
        close(26)
      endif
500 format(e14.7,1x,E14.7,1x,e12.5,1x,15f8.3)
510 format(15(1x,1PE12.4))
520 format(i8,15(1x,F7.4,1x,F7.2))

    stop
    end

!***********************************************************************
      SUBROUTINE CHLSKY(A,F,G,N)

!***********************************************************************
!*    ROUTINE TO SOLVE SYMMETRIC POSITIVE DEFINITE MATRICES FOUND IN
!*    LINEAR LEAST SQUARES PROBLEMS.  THE CHOLESKY SQUARE METHOD (ALSO
!*    ATTRIBUTED TO BANACHIEWICZ) IS USED.  SEE LEAST SQUARES PROGRAM-
!*    MING FOR BAND MATRICES BY B. EVANS. COMPUTATIONAL METHODS FOR
!*    LINEAR ALGEBRA BY FADDEEV AND FADDEEVA PAGE 145. ADVANCES IN
!*    COMPUTERS VOL 2 PAGE 57.  COMPUTER SOLUTION OF LINEAR ALGEBRAIC
!*    SYSTEMS BY FORSYTHE AND MOLER PAGE 114.
!*
!*    DIRECTIONS...
!*    1... CALL TRIANG WITH THE UPPER TRIANGLE OF THE MATRIX STORED
!*    BY ROWS IN THE ONE DIMENSIONAL ARRAR A.  THE PRICE OF THIS
!*    STORAGE SAVING IS EXTRA PROCESSING TIME DUE TO A COMPLEX
!*    SUBSCRIPTING. IN TRIANG, A IS FACTORIZED AND THE NEW MATRIX
!*    REPLACES A.  THE MATRIX CONDITION IS RETURNED IN G.  IF (G.LE.0),
!*    THEN G IS SET TO -I (ROW NUMBER) AND AN IMMEDIATE RETURN IS MADE.
!*    2...CALL SOLUTN WITH THE RIGHT HAND SIDE ARRAY IN F. THE SOLUTION
!*    IS RETURNED IN F.

!*    DEC/76 ... INDEXING SIMPLIFIED TO IMPROVE SPEED.
!*    IF F(I)=0 FOR ALL I EXCEPT F(IENTRY)=1.0, THEN SETTING G TO
!*    IENTRY CAUSES THE PROGRAM TO CALCULATE ONLY THE PARTIAL SOLUTION
!*    F(IENTRY)...F(N). THIS IS USEFUL FOR GETTING THE
!*    DIAGONAL ELEMENTS OF THE INVERSE MATRIX WHICH ARE RELATED TO THE
!*    EXPECTED ERRORS IN THE REGRESSION COEFFICIENTS AND ALLOWS A SPEED
!*    IMPROVEMENT OF ABOUT 3X. FOR OTHER TYPES OF RHS THE FULL SOLUTION
!*    IS FOUND. TO CANCEL THIS FEATURE SET IENTRY=1 FOR ALL TIME.

!*    TIMING FOR...TRIANGULARIZE + SOLVE ONE RHS + GET INVERSE MATRIX
!*    DIAGONAL ELEMENTS (N=124, CYBER/74)
!*       OLD TIME = 3.0+.09+12.0=15.9 SEC.
!*       NEW TIME = 2.0+.06+ 3.5= 5.56 SEC.=37 PERCENT OF OLD.

      DIMENSION A(400),F(50)

!***********************************************************************
!*    CHOLESKY FACTORIZATION OF THE MATRIX A.
      ENTRY TRIANG(A,F,G,N)
      NP1=N+1
10    IF(A(1))20,20,30
20    G=-1.
      RETURN
30    A(1)= SQRT(A(1))
      DO 40 J=2,N
40    A(J)=A(J)/A(1)
      G=1.
      II=1
50    DO 100 I=2,N
      IP1=I+1
      IM1=I-1
      II=II+N+2-I
      AII=A(II)
      LI=I-N
      DO 60 L=1,IM1
      LI=LI+NP1-L
60    AII=AII-A(LI)*A(LI)
      IF(AII)65,65,70
65    G=-I
      RETURN
70    G=G*AII/A(II)
      A(II)= SQRT(AII)
      IF(I-N)75,100,75
75    DO 90 J=IP1,N
      IJ=II+J-I
      AIJ=A(IJ)
      JMI=J-I
      LI=I-N
      DO 80 L=1,IM1
      LI=LI+NP1-L
      LJ=LI+JMI
80    AIJ=AIJ-A(LI)*A(LJ)
90    A(IJ)=AIJ/A(II)
100   CONTINUE
      G= SQRT(G)
      RETURN

!***********************************************************************
!*    SOLUTION OF THE TRIANGULARIZED MATRIX EQUATION
!*    G NOW DENOTES WHETHER OR NOT AN EXPECTED ERROR CALCULATION IS
!*    DESIRED; THAT IS, IF THE RIGHT HAND SIDE OF THE MATRIX EQUATION
!*    IS A VECTOR WITH VALUE 1 IN ONLY ONE LOCATION AND ZEROES
!*    EVERYWHERE ELSE.  IN THESE CASES, THE FORWARD AND BACKWARD
!*    SUBSTITUTION METHODS CAN BE ABBREVIATED AND RESULT IN AN APPROX.
!*    3 FACTOR TIME SAVING.  A G VALUE EQUAL TO K INDICATES SUCH A
!*    CASE WITH THE 1 IN THE K-TH POSITION.  A G VALUE EQUAL TO 1
!*    INDICATES THAT A STANDARD MATRIX SOLUTION METHOD IS REQUIRED.
      ENTRY SOLUTN(A,F,G,N)

!***********************************************************************
!     FORWARD SOLUTION OF RIGHT HAND SIDE

      IENTRY=INT(G)
      NP1=N+1
      DO 220 I=IENTRY,N
      IM1=I-1
      FI=F(I)
      LI=((N+N-IENTRY)*(IENTRY-1))/2+I
      IF(I.EQ.IENTRY) GO TO 220
      DO 210 L=IENTRY,IM1
      FI=FI-A(LI)*F(L)
  210 LI=LI+N-L
220   F(I)=FI/A(LI)

!***********************************************************************
!*    BACKWARD SUBSTITUTION FOR FINAL SOLUTION

230   NN=(N*(N+1))/2
      F(N)=F(N)/A(NN)
      IF(IENTRY.EQ.N) RETURN
      II=NN
      IENTB=N+1-IENTRY
      DO 260 IB=2,IENTB
      I=N+1-IB
      II=II+I-NP1
      IP1=I+1
      FI=F(I)
      III=II-I
      DO 250 L=IP1,N
      IL=III+L
250   FI=FI-A(IL)*F(L)
260   F(I)=FI/A(II)
      RETURN
      END

!***********************************************************************

      SUBROUTINE INFER(N,KONA,SIGA,AAN,EPSAN,fan,vuan,KONI,SIGI,AIN,&
               EPSIN,fin,vuin,R,ZETA)

!***********************************************************************
!*  THIS SUBROUTINE INFERS THE AMPLITUDE & PHASE OF A CONSTITUENT NOT
!*  INCLUDED IN THE ANALYSIS FROM ONE THAT WAS. THE AMPLITUDE & PHASE OF
!*  THE CONSTITUENT USED FOR THE INFERENCE IS APPROPRIATELY ADJUSTED.
!
!*     KONAN & SIGAN ARE THE NAME AND FREQUENCY OF A CONSTITUENT TO BE
!*                   USED FOR INFERENCE.
!*     KONIN & SIGIN ARE THE NAME AND FREQUNCY OF A CONSTITUENT TO BE
!*                   INFERRED FROM KONAN.
!*     R IS THE AMPLITUDE RATIO OF KONIN TO KONAN.
!*     ZETA IS THE GREENWICH PHASE LAG OF THE MAIN CONSTITUENT -
!*          GREENWICH PHASE LAG OF THE INFERRED CONSTITUENT.
!*  KONAN,SIGAN,KONIN,SIGIN,R,ZETA ARE READ FROM DEVICE MTD WITH THE
!*  FORMAT (2(4X,A5,E16.10),2F10.3). THEIR END IS DENOTED BY A BLANK
!*  CARD.

      character*4 KONI,KONA
      real R,ZETA
      DATA KBLANK/'     '/

!***********************************************************************
!*  PERFORM INFERENCE

      PI=3.1415926536
      TWOPI=2.*PI
!    convert all input angles to cycles
      ZETAC=ZETA/360.
      vuanc=vuan/360.
      vuinc=vuin/360.
      epsan=epsan/360.
!      CALL VUF (KH,KONAN(L),VUAN,FAN,XLAT)
!      CALL VUF (KH,KONIN(L),VUIN,FIN,XLAT)
      SIGDIF=SIGA-SIGI
      DSIGIA=SIGDIF*N*PI
      IBANG=int(SIGDIF*N*.5D0+.0001)
      BANG=(SIGDIF*N*0.50-IBANG)*TWOPI
      V=SIN(BANG)/DSIGIA
      FACTOR=R*V*FIN/FAN
      ETA=(VUINc-VUANc+ZETAC)*TWOPI
      STERM=FACTOR*SIN(ETA)
      CTERM=1.+FACTOR*COS(ETA)
!    write(25,*) ' siga,v,factor,eta,sterm,cterm=',siga,v,factor,eta,sterm,cterm
      if(abs(sterm).lt.1.e-6.and.abs(cterm).lt.1.e-6) then
      div=0.
      dg=0.
      else
        DIV=SQRT(CTERM*CTERM+STERM*STERM)
        DG=ATAN2(STERM,CTERM)/TWOPI
      end if
      EPSAN=EPSAN+DG
      AAN=AAN/DIV
      EPSIN=EPSAN+(VUANc-ZETAC-VUINc)
      AIN=(AAN/FAN)*R*FIN
!    convert all angles back to degrees
      epsan=epsan*360.
      epsin=epsin*360.
    
      RETURN
      END
      
!***********************************************************************
