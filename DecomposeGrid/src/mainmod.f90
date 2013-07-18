!******************************************************************************
!******************************************************************************

      MODULE MainData

      SAVE

! Parameters

! Scalar data
      integer nproc, numcomm
      integer  NE,NP,NPR,NCN,ntype,nsides,npv
      integer  nnbr, numhalo
      integer  irst
      integer iwind, iwindx,iAnalysisTime
      integer izup,ifront
      real TET
! 1D data
      integer, allocatable ::  nep(:), nehalo(:), npp(:), nphalo(:), nsp(:), nshalo(:)
      integer, allocatable ::  nbc(:), IECode(:), inehalo(:), ncomm(:)
      integer, allocatable ::  epart(:), npart(:),halomap(:,:)
      real, allocatable ::  Area(:), alfa(:)
      real, allocatable ::  sdep(:),slen(:),refdep(:),sdx(:),sdy(:)
      real, allocatable ::  dlinv(:)
      real, allocatable ::  eta(:)

! 2D data
      integer, allocatable ::  nen(:,:), ieadj(:,:)
      integer, allocatable ::  numsideT(:,:), iside(:,:),iends(:,:) !numside(:,:),
      integer, allocatable ::  elemapG2L(:,:), nodemapG2L(:,:), sidemapG2L(:,:)
      integer, allocatable ::  elemapL2G(:,:), nodemapL2G(:,:),sidemapL2G(:,:)
      integer, allocatable ::  icomm(:,:), ecomm(:,:), commid(:,:)
      integer, allocatable ::  ine(:,:), ene(:,:),ins(:,:), ens(:,:)
      real, allocatable ::  sxy(:,:), xyz(:,:)
      real, allocatable ::  un(:,:), ut(:,:), wn(:,:), ut1(:,:), ut2(:,:) !AB3 Coriolis

! 3D data
      integer, allocatable ::  ipoly(:,:,:), epoly(:,:,:)
      integer, allocatable ::  iedge(:,:,:), eedge(:,:,:)


! Character data
      character*80 OutResFile
      character*(128) GridFileName, WindFileName
      character(120) :: WindSourceFileName
      character(18) :: AnalysisTime
      character(30) :: rstFile
      character(2) :: nrstchar
      character(15) :: version

      END

!******************************************************************************
!******************************************************************************
