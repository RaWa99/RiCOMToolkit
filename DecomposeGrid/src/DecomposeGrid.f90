! *********************************************************************

      PROGRAM DecomposeGrid

! *********************************************************************
! ***  DecomposeGrid
! ***  Decomposes grid for parallel execution of RiCOM, based on
! ***             partition table generated by Metis 4.
! ***  Reads preprocessed binary grid file for RiCOM (prepared by PreMod5)
! ***             Then divides grid in ncpu parts in separate directories.
! ***  Reference: Jacek
! ***  Written by: Roy A. Walters
! ***  
! ***  IMPORTANT NOTE:
! ***  This program is provided for personal use and may not be distributed
!      to anyone else without the author's explicit permission.
! *********************************************************************

! *********************************************************************
! *** Exit Codes
! *** 70 unable to open input file
! *** 71 problem allocating arrays
! *** 72 Unknown parameter value/invalid input
! *********************************************************************

      USE MainData
      
      implicit none

! *** local variables
      integer ivar(20)
      real rvar(5)
      real ctime0, stime0, ctime1, stime1, etime

      include 'RCMversion.h'

! *** initialization
      call ticktock(ctime0,stime0)
      CALL DataInput(ivar,rvar)
      call PartitionGrid()
      call decompose()
      call decompose_restart()
      call OutputGrids
      call ticktock(ctime1,stime1)
      etime = ctime1 - ctime0
      write(6,*) ' Elapsed proc time in input=', etime
      etime = stime1 - stime0
      write(6,*) ' Elapsed sys time in input=', etime


      STOP
      END

!*****************************************************************
!*****************************************************************

      SUBROUTINE DataInput(ivar, rvar)

! *** Read grid data created by PROGRAM PREMOD, read boundary data,
! *** and initialize model.

      USE MainData

      implicit none

! *** passed variables
      integer ivar(*)
      real rvar(*)

! *** local variables
      integer :: i,j,k,istat,ierr,nn,ncn2,n1,n2,js
      integer nnp,nndf,npvx,ndfv,iopt,ndim
      logical ResOK, openbinfile
      CHARACTER FNAME*256, ans*1

! *** START
      write(*,*) ' A grid decomposition program for the model'
      write(*,*) ' RiCOM Version ',trim(version),' (River and Coastal Ocean Model)'
      write(*,*) ' created by Roy Walters,'
      write(*,*) ' Internet:  rawalters@gmail.com'
      write(*,*) ' Program size: Allocatable arays. A message will appear'
      write(*,*) '               if there is insufficient memory'
      write(*,*) ' '


! *** Initial values
      TET = 0.

      ResOK = OpenBinFile(22,'Open pre-processed Grid File',fname, &
     &    "Input file(*.rcm),*.rcm;All files(*.*),*.*;")

      if(.not.resOK) then
        write(*,*), 'Unable to open input file'
        call exit(71)
      endif

!      OPEN(UNIT=22,file=GridFileName,status='OLD',FORM='UNFORMATTED')      
      READ (22) NE,NTYPE,NP,NPR,NCN,nsides,nnbr,izup,ifront

!  Main arrays for grid
      ALLOCATE (nen(ne,ncn),xyz(np,3),area(ne),nbc(np),alfa(np), &
        ieadj(5,ne),numsideT(4,ne),IECode(ne), STAT = istat )

      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate main storage arrays'
        call exit(71)
      endif

      READ (22) ((NEN(J,K),K=1,NCN),J=1,NE),(IECode(j),J=1,NE), &
     &         ((XYZ(J,K),J=1,NP),K=1,3),(ALFA(j),J=1,NP),         &
     &         (NBC(J),J=1,NP)
      
      if(nsides.gt.0) then
!  Side-based arrays
        ALLOCATE (sdep(0:nsides),slen(nsides),refdep(nsides),sdx(nsides),sdy(nsides), &
          sxy(2,nsides),dlinv(nsides),iside(3,nsides),iends(2,nsides), STAT = istat )
        if(istat.ne.0) then
          write(*,*) 'FATAL ERROR: Cannot allocate side-based storage arrays'
          call exit(71)
        endif

        READ(22) (area(j),j=1,ne),((ieadj(k,j),k=1,5),j=1,ne),((numsideT(k,j),j=1,ne),k=1,4) 
        READ(22) ((iside(k,j),k=1,2),j=1,nsides),((sxy(k,j),k=1,2),j=1,nsides)
        READ(22) (refdep(j),j=1,nsides),(slen(j),j=1,nsides)
        READ(22) (sdx(j),j=1,nsides),(sdy(j),j=1,nsides),(dlinv(j),j=1,nsides)
        READ(22,IOSTAT=istat) ((iends(k,j),k=1,2),j=1,nsides)
! *** temporarily fill iends until put in PreMod
        if(istat.ne.0) then
          do nn=1,ne
            ncn2 = ncn -1 + min0(1,nen(nn,ncn))
            do j=1,ncn2
              N1=NEN(nn,j)
              N2=NEN(nn,mod(j,ncn2)+1)
              js = numsideT(j,nn)
              iends(1,js ) = n1
              iends(2,js) = n2
            enddo
          enddo
        endif
      else
        nsides = 0
      endif
      CLOSE(unit=22)
 
      write(6,*) ' ne,ncn,ntype=',ne,ncn,ntype
      write(6,*) ' np,nsides,nnbr=',np,nsides,nnbr

! *** restart file
      write(*,*) ' Decompose restart file? (Y/N)'
      READ(*,'(a)') ans
      
      if(ans(1:1).eq.'Y'.or.ans(1:1).eq.'y') then

        ResOK = OpenBinFile(2,'Open restart File',fname, &
          "Input file(*.bin),*.bin;All files(*.*),*.*;")

        if(.not.resOK) then
          write(*,*), 'Unable to open restart file'
          call exit(71)
        endif
        
        irst = 1

! *** If iwind = 3 (i.e. NZLAM wind/pressure input) read in Analysis time for restart file
        IF ((iwind .eq. 2) .or. (iwind .eq. 3)) THEN
          READ(2) TET,nnp,nndf,npvx,ndfv,iopt,iAnalysisTime
		  write(*,*), 'Analysis time from restart file', iAnalysisTime
        ELSE
          READ(2) TET,nnp,nndf,npvx,ndfv,iopt
        END IF
 
        ndim = max(nsides,ne)
        npv = npvx
        allocate ( eta(0:ndim),un(npv,nsides),ut(npv,nsides),ut1(npv,nsides),&
                   ut2(npv,nsides), STAT = istat )
        if(istat.ne.0) then
          write(*,*) 'FATAL ERROR: Cannot allocate dependent variables'
          call exit(71)
        endif
        
        READ(2) (eta(I),I=1,nnp),(un(1,i),i=1,nnp),(ut(1,i),i=1,nnp)
!        eta = amax1(eta-0.1,elev)
        if(npvx.gt.1) then
          do k=2,npvx
            read(2,end=188) (un(k,j),j=1,nnp),(ut(k,j),j=1,nnp)
          enddo
        endif
        ut1 = ut
        ut2 = ut
        read(2,end=188) ((ut1(k,j),k=1,npvx),j=1,nnp)
        read(2,end=188) ((ut2(k,j),k=1,npvx),j=1,nnp)
188     CLOSE (unit=2)
      else
        irst = 0
      endif

      RETURN

      END

!*****************************************************************
!*****************************************************************

      Subroutine PartitionGrid()

      USE MainData
      
      implicit none

! *** local variables
      integer j, k, istat
      integer etype, numflag, edgecut
      integer, allocatable ::  eptr(:), eind(:), vwgt(:), vsize(:) 
!      integer, allocatable :: epart(:), npart(:)
      real, allocatable :: tpwgts(:)
      logical ResOK, openfile
      CHARACTER FNAME*256, ans*1

      write(*,*) ' Partition grid (Y/N)?'
      READ(*,'(a)') ans
      
      if(ans(1:1).eq.'Y'.or.ans(1:1).eq.'y') then
      
        write(*,*) ' Enter the number of partitions'
        read(*,*) nproc

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
            eind(3*(j-1)+k) = nen(j,k)
          enddo
        enddo

        call metis_partmeshdual(ne,np,eind,etype,numflag,nproc,&
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
        open(21,file='partitiontable.dat', status='unknown')
        do j=1,ne
          write(21,*) epart(j)
        enddo
        close(21)

      else 
        ResOK = OpenFile(22,'Open Partition Table File',fname, &
     &    "Input file(*.lel),*.lel;All files(*.*),*.*;")

        if(.not.resOK) then
          write(*,*), 'Unable to open partition file'
          call exit(71)
        endif

        ALLOCATE (epart(ne), STAT = istat )
        if(istat.ne.0) then
          write(*,*) 'FATAL ERROR: Cannot allocate partition array'
          call exit(71)
        endif
        
        do j=1,ne
          read(22,*) epart(j)
        enddo
        close(22)
      
      endif


      RETURN
      END

!*****************************************************************
!*****************************************************************

      Subroutine decompose()

      USE MainData
      
      implicit none

! *** local variables
      integer j,j2,jp,k,kk,istat,ncn2,nn,ep,nepmax,neh
      integer nspmax,mr,mr2,ns,nppmax,eG
      logical found

! *** determine size of partition arrays and form local to global map

      ALLOCATE (nep(nproc),nehalo(nproc),npp(nproc),nphalo(nproc),&
                nsp(nproc),nshalo(nproc),elemapG2L(2,ne), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate decompose arrays'
        call exit(71)
      endif
      
      !debug
      open(99,file='basedata.dat',status='unknown')
      write(99,*) ' global element data=ne,npart,nen(1:4),numside(1:4)'
      do j=1,ne
        write(99,'(10i5)') j,epart(j),(nen(j,k),k=1,ncn),(numsideT(k,j),k=1,ncn)
      enddo
      
      nep = 0
      nehalo = 0
      elemapG2L = 0
      
      do j=1,ne
        nn = nep(epart(j)) + 1
        elemapG2L(1,j) = epart(j)
        elemapG2L(2,j) = nn
!        elemapL2G(nn,epart(j)) = j
        nep(epart(j)) = nn
        !  count halo elements
        ncn2 = ncn -1 + min0(1,nen(j,ncn))
        do k=1,ncn2
!          mr = numsideT(k,j)
          j2 = ieadj(k+1,j)
          if(j2.gt.0) then
!            write(*,*) ' j,k,mr,j2,p1,p2',j,k,mr,j2,epart(j),epart(j2)
            if(epart(j).ne.epart(j2)) then
              nehalo(epart(j)) = nehalo(epart(j)) + 1
            endif
          endif
        enddo
      enddo
      
      write(*,*) ' nep=',(nep(j),j=1,nproc)
      write(*,*) ' nehalo=',(nehalo(j),j=1,nproc)
      
! *** form local to global map for elements

      nepmax = maxval(nep) + maxval(nehalo)
      ALLOCATE (elemapL2G(nepmax,nproc), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate L2G map array'
        call exit(71)
      endif
      
      elemapL2G = 0
      nehalo = 0
      
      do j=1,ne
        ep = elemapG2L(1,j)
        nn = elemapG2L(2,j)
        elemapL2G(nn,ep) = j
        ! add halo elements
        nn = nep(ep)
        ncn2 = ncn -1 + min0(1,nen(j,ncn))
        do k=1,ncn2
!          mr = numsideT(k,j)
          j2 = ieadj(k+1,j)
          if(j2.gt.0) then
!            write(*,*) ' j,k,mr,j2,p1,p2',j,k,mr,j2,epart(j),epart(j2)
            if(epart(j).ne.epart(j2)) then
              neh = nehalo(ep) + 1
              elemapL2G(nn+neh,ep) = j2
              nehalo(ep) = neh
            endif
          endif
        enddo
      enddo
      
      write(*,*) ' nep=',(nep(j),j=1,nproc)
      write(*,*) ' nehalo=',(nehalo(j),j=1,nproc)

      !debug
      write(99,*) ' local element data= npart,j,mapL2G'
      do jp=1,nproc
        write(99,*) ' jp, nep, nehalo=',nep(jp),nehalo(jp)
        do j=1,nep(jp)+nehalo(jp)
          write(99,'(10i5)') jp,j,elemapL2G(j,jp)
        enddo
      enddo

      
! *** form local to global map for sides

      nspmax = 4*(maxval(nep) + maxval(nehalo))
      ALLOCATE (sidemapL2G(nspmax,nproc), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate L2G map array'
        call exit(71)
      endif
      
      nsp = 0
      nshalo = 0
      sidemapL2G = 0
      
      do j=1,ne
        ep = epart(j)
        ns = nsp(ep)
        ncn2 = ncn -1 + min0(1,nen(j,ncn))
        do k=1,ncn2
          mr = numsideT(k,j)
          sidemapL2G(ns+1,ep) = mr
          ns = ns + 1
!          write(*,*) ' ns,ep,sidemapL2G=',ns,ep,sidemapL2G(ns,ep) 
        enddo
        nsp(ep) = ns
      enddo
      
      write(*,*) ' nsp=',(nsp(j),j=1,nproc)

      !debug
      write(99,*) ' local side data= npart,j,sidemapL2G'
           
      ! eliminate duplicates
      do jp=1,nproc
        ns = 1
        do k=2,nsp(jp)
          mr = sidemapL2G(k,jp)
          found = .false.
          do kk=1,k-1
            mr2 = sidemapL2G(kk,jp)
            if(mr2.eq.mr) then
              found = .true.
              exit
            endif
          enddo
          if(.not.found) then
            ns = ns + 1
            sidemapL2G(ns,jp) = mr
 !           write(*,*) ' jp,k,kk,mr,mr2,ns=',jp,k,kk,mr,mr2,ns
          endif
        enddo
        nsp(jp) = ns
        ! add sides for halo
        do j=nep(jp)+1,nep(jp)+nehalo(jp)
          eG = elemapL2G(j,jp)
          ncn2 = ncn -1 + min0(1,nen(eG,ncn))
          do k=1,ncn2
            mr = numsideT(k,eG)
            sidemapL2G(ns+1,jp) = mr
            ns = ns + 1
!            write(*,*) ' add ns,ep,sidemapL2G=',ns,jp,sidemapL2G(ns,jp) 
          enddo
        enddo
        nshalo(jp) = ns - nsp(jp)
        ! eliminate duplicates
        ns = nsp(jp)
        do k=nsp(jp)+1,nsp(jp)+nshalo(jp)
          mr = sidemapL2G(k,jp)
          found = .false.
          do kk=1,k-1
            mr2 = sidemapL2G(kk,jp)
            if(mr2.eq.mr) then
              found = .true.
              exit
            endif
          enddo
          if(.not.found) then
            ns = ns + 1
            sidemapL2G(ns,jp) = mr
 !           write(*,*) ' jp,k,kk,mr,mr2,ns=',jp,k,kk,mr,mr2,ns
          endif
        enddo
        nshalo(jp) = ns - nsp(jp)
        write(99,*) ' jp, nsp, nshalo=',nsp(jp),nshalo(jp)
        do j=1,nsp(jp)+nshalo(jp)
          write(99,'(10i5)') jp,j,sidemapL2G(j,jp)
        enddo
      enddo
     
      write(*,*) ' nsp=',(nsp(j),j=1,nproc)
      write(*,*) ' nshalo=',(nshalo(j),j=1,nproc)
     
! *** form local to global map for nodes

      nppmax = 4*(maxval(nep) + maxval(nehalo))
      ALLOCATE (nodemapL2G(nppmax,nproc), STAT = istat )
      if(istat.ne.0) then
        write(*,*) 'FATAL ERROR: Cannot allocate L2G map array'
        call exit(71)
      endif
      
      npp = 0
      nphalo = 0
      nodemapL2G = 0
      
      do j=1,ne
        ep = elemapG2L(1,j)
        nn = elemapG2L(2,j)
        ns = npp(ep)
        ncn2 = ncn -1 + min0(1,nen(j,ncn))
        do k=1,ncn2
          mr = nen(j,k)
          nodemapL2G(ns+1,ep) = mr
          ns = ns + 1
        enddo
        npp(ep) = ns
      enddo
      
      write(*,*) ' npp=',(npp(j),j=1,nproc)

      !debug
      write(99,*) ' local node data= npart,j,nodemapL2G'
      
      ! eliminate duplicates
      do jp=1,nproc
        ns = 1
        do k=2,npp(jp)
          mr = nodemapL2G(k,jp)
          found = .false.
          do kk=1,k-1
            mr2 = nodemapL2G(kk,jp)
            if(mr2.eq.mr) then
              found = .true.
              exit
            endif
          enddo
          if(.not.found) then
            ns = ns + 1
            nodemapL2G(ns,jp) = mr
          endif
        enddo
        npp(jp) = ns
        ! add nodes for halo
        do j=nep(jp)+1,nep(jp)+nehalo(jp)
          eG = elemapL2G(j,jp)
          ncn2 = ncn -1 + min0(1,nen(eG,ncn))
          do k=1,ncn2
            mr = nen(eG,k)
            nodemapL2G(ns+1,jp) = mr
            ns = ns + 1
          enddo
        enddo
        nphalo(jp) = ns - npp(jp)
        ! eliminate duplicates
        ns = npp(jp)
        do k=npp(jp)+1,npp(jp)+nphalo(jp)
          mr = nodemapL2G(k,jp)
          found = .false.
          do kk=1,k-1
            mr2 = nodemapL2G(kk,jp)
            if(mr2.eq.mr) then
              found = .true.
              exit
            endif
          enddo
          if(.not.found) then
            ns = ns + 1
            nodemapL2G(ns,jp) = mr
          endif
        enddo
        nphalo(jp) = ns - npp(jp)
        write(99,*) ' jp, npp, nphalo=',npp(jp),nphalo(jp)
        do j=1,npp(jp)+nphalo(jp)
          write(99,'(10i5)') jp,j,nodemapL2G(j,jp)
        enddo
      enddo
      
      write(*,*) ' npp=',(npp(j),j=1,nproc)
      write(*,*) ' nphalo=',(nphalo(j),j=1,nproc)


      RETURN
      END

!*****************************************************************
!*****************************************************************

      Subroutine decompose_restart()

      USE MainData
      
      implicit none

! *** local variables

      RETURN
      END

!*****************************************************************
!*****************************************************************

      Subroutine OutputGrids

      USE MainData
      
      implicit none

! *** local variables
      integer j,jp,k,istat
      integer nG,eG,sG,ncn2
      integer nepmax,nppmax,nspmax,netot,nptot,nstot
      integer, allocatable :: nenp(:,:),ieadjp(:,:),numsideTp(:,:),IECodep(:),nbcp(:)
      integer, allocatable :: isidep(:,:),iendsp(:,:)
      real, allocatable :: xyzp(:,:),alfap(:)
      real, allocatable ::  Areap(:),sdepp(:),slenp(:),refdepp(:),sdxp(:),sdyp(:)
      real, allocatable ::  sxyp(:,:), dlinvp(:)
      character*6 :: numpe='PE0000'

!  allocate arrays for grid partition
        nepmax = maxval(nep) + maxval(nehalo)
        nppmax = maxval(npp) + maxval(nphalo)
        ALLOCATE (elemapG2Lp(ne),nodemapG2L(np),sidemapG2L(nsides),&
          nenp(nepmax,4),xyzp(nppmax,3),areap(nepmax),nbcp(nppmax),alfap(nppmax), &
          ieadjp(5,nepmax),numsideTp(4,nepmax),IECodep(nepmax), STAT = istat )

        if(istat.ne.0) then
          write(*,*) 'FATAL ERROR: Cannot allocate main output arrays'
          call exit(71)
        endif
        
        if(nsides.gt.0) then
          nspmax = maxval(nsp) + maxval(nshalo)
          ALLOCATE (sdepp(0:nspmax),slenp(nspmax),refdepp(nspmax),sdxp(nspmax),sdyp(nspmax), &
            sxyp(2,nspmax),dlinvp(nspmax),isidep(3,nspmax),iendsp(2,nspmax), STAT = istat )
          if(istat.ne.0) then
            write(*,*) 'FATAL ERROR: Cannot allocate side-based storage arrays',istat
            call exit(71)
          endif
        endif

! *** form local domain arrays and write out      

      do jp=1,nproc
      
        if(jp.eq.2) numpe(6:6)='1'
                
        netot = nep(jp) + nehalo(jp)
        nptot = npp(jp) + nphalo(jp)
        nstot = nsp(jp) + nshalo(jp)

        ! form local mappings
        elemapG2Lp = 0
        do j=1,netot
          nG = elemapL2G(j,jp)
          elemapG2Lp(nG) = j
        enddo
        nodemapG2L = 0
        do j=1,nptot
          nG = nodemapL2G(j,jp)
          nodemapG2L(nG) = j
        enddo
        sidemapG2L = 0        
        do j=1,nstot
          nG = sidemapL2G(j,jp)
          sidemapG2L(nG) = j
        enddo
        
        ! form local arrays for output
        nenp = 0
        iecodep = 1
        ieadjp = 0
        numsideTp = 0
        do j=1,netot
          eG = elemapL2G(j,jp)
          areap(j) = area(eG)
          iecodep(j) = iecode(eG)
          ieadjp(1,j) = j
          ncn2 = ncn -1 + min0(1,nen(eG,ncn))
          do k=1,ncn2
            nG = nen(eG,k)
            nenp(j,k) = nodemapG2L(nG)
            nG = ieadj(k+1,eG)
            if(nG.gt.0) ieadjp(k+1,j) = elemapG2Lp(nG)
            nG = numsideT(k,eG)
            numsideTp(k,j) = sidemapG2L(nG)
          enddo
        enddo
        
        xyzp = 0
        alfap = 0
        nbcp = 0
        do j=1,nptot
          nG = nodemapL2G(j,jp)
          xyzp(j,1:3) = xyz(nG,1:3)
          alfap(j) = alfa(nG)
          nbcp(j) = nbc(nG)
        enddo
        
! *** write it out
        OPEN(UNIT=22,file='GridFile'//numpe//'.xye',status='unknown')      
        
        write(22,'(a)') '#XYE'
        write(22,*) nptot,netot
        do j=1,nptot
          write(22,'(3e12.4,i5)') (xyzp(j,k),k=1,3),nbcp(j)
        enddo
        do j=1,netot
          write(22,'(4i7,i5)') (nenp(j,k),k=1,4),iecode(j)
        enddo
        close(22)
        
        OPEN(UNIT=22,file='gridfile'//numpe//'.bin',status='unknown',FORM='UNFORMATTED')      
        write(22) netot,nehalo(jp),nptot,nphalo(jp),NCN,nstot,nshalo(jp),izup,ifront
       
        write(22) ((NENp(J,K),K=1,NCN),J=1,netot),(IECodep(j),J=1,netot), &
              ((XYZp(J,K),J=1,NPtot),K=1,3),(ALFAp(j),J=1,NPtot),(NBCp(J),J=1,NPtot)
      
        if(nsides.gt.0) then
!  Side-based arrays
                         
          ! form local arrays
          isidep = 0
          sdepp(0)=sdep(0)
          do j=1,nstot  !p(jp)
            sG = sidemapL2G(j,jp)
            sdepp(j)=sdep(sG)
            slenp(j)=slen(sG)
            refdepp(j)=refdep(sG)
            sdxp(j)=sdx(sG)
            sdyp(j)=sdy(sG)
            dlinvp(j)=dlinv(sG)
            sxyp(1:2,j)=sxy(1:2,sG)
            do k=1,2
              nG = iside(k,sG)
              if(nG.gt.0) isidep(k,j)=elemapG2Lp(nG)
              nG = iends(k,sG)
              iendsp(k,j)=nodemapG2L(nG)
            enddo
            isidep(3,j)=iside(3,sG)
          enddo
          ! resort isides for halo
          do j=nsp(jp)+1,nstot
            if(isidep(1,j).eq.0) then
              isidep(1,j) = isidep(2,j)
              isidep(2,j) = 0
            endif
          enddo
        
          ! write it out
          write(22) (areap(j),j=1,netot),((ieadjp(k,j),k=1,5),j=1,netot),((numsideTp(k,j),j=1,netot),k=1,4) 
          write(22) ((isidep(k,j),k=1,2),j=1,nstot),((sxyp(k,j),k=1,2),j=1,nstot)
          write(22) (refdepp(j),j=1,nstot),(slenp(j),j=1,nstot)
          write(22) (sdxp(j),j=1,nstot),(sdyp(j),j=1,nstot),(dlinvp(j),j=1,nstot)
          write(22,IOSTAT=istat) ((iendsp(k,j),k=1,2),j=1,nstot)
        endif
! *** write map information
        write(22) nproc,jp
        write(22) (elemapG2Lp(j),j=1,ne),(elemapL2G(j,jp),j=1,netot)
        write(22) (nodemapG2L(j),j=1,ne),(nodemapL2G(j,jp),j=1,nptot)
        write(22) (sidemapG2L(j),j=1,ne),(sidemapL2G(j,jp),j=1,nstot)

        CLOSE(unit=22)

! *** write base output
        OPEN(UNIT=22,file='baseout'//numpe//'.dat',status='unknown')      
        write(22,*) 'netot,nehalo(jp),nptot,nphalo(jp),NCN,nstot,nshalo(jp),izup,ifront'
        write(22,*) netot,nehalo(jp),nptot,nphalo(jp),NCN,nstot,nshalo(jp),izup,ifront
        write(22,*) 'elements j,nen(j,1:ncn),iecode(j)'
        do j=1,netot
          write(22,'(10i6)') j,(nenp(j,k),k=1,4),iecodep(j)
        enddo
        write(22,*) 'elements j,area,ieadj(1:5,j),numside(1:4,j)'
        do j=1,netot
          write(22,'(i6,e12.4,10i6)') j,areap(j),(ieadjp(k,j),k=1,5),(numsideTp(k,j),k=1,4)
        enddo
        write(22,*) 'nodes j,xyz(j,1:3),alfa(j)'
        do j=1,nptot
          write(22,'(i6,4e12.4)') j,(xyzp(j,k),k=1,3),alfap(j)
        enddo
        write(22,*) 'sides j,iside(1:3,j),iends(1:2,j)'
        do j=1,nstot
          write(22,'(10i6)') j,(isidep(k,j),k=1,3),(iendsp(k,j),k=1,2)
        enddo
        write(22,*) 'sides j,iside(1:3,j),iends(1:2,j)'
        do j=1,nstot
          write(22,'(i6,8e12.4)') j,(sxyp(k,j),k=1,2),refdepp(j),slenp(j),sdxp(j),sdyp(j),dlinvp(j)
        enddo
        close(22)
      enddo

      RETURN
      END

!***************************************************************
!*****************************************************************

        subroutine ticktock(time1,time2)

          integer(8) ic8,crate8
          real time1,time2

! *** processor time
          call cpu_time(time1)
! *** system time
          call system_clock(count=ic8,count_rate=crate8)
          time2 = float(ic8)/float(crate8)

        return
        end

!***************************************************************
!***************************************************************

