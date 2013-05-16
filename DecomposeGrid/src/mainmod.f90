!******************************************************************************
!******************************************************************************

      MODULE MainData

      SAVE

! Parameters

! Scalar data
      integer nproc
      integer  NE,NP,NPR,NDF,NCN,ntype,ntypes,ntypec,ndfe2,nsides,nsidesbc
      integer  ICVG,IDX,NPRT,itest,nitn,ncon,nsed,jsed,nsol,nnbr,nbx,npMB,neMB
      integer  npv, npvc,izcoord,izgrid
      integer  NOPT,NEQS,idpdn,isolve,NSon,iqp,nsbc,irst,irstout,nrstout,icase,nitsol,NMit
      integer iwind, iwindx,iAnalysisTime,neqtide
      integer icoord,izup,ifront,maxrow,maxsto,ielmSR
      real TET
! 1D data
      integer, allocatable ::  nep(:), nehalo(:), npp(:), nphalo(:), nsp(:), nshalo(:)
      integer, allocatable ::  nbc(:),IndexQ(:), IECode(:),numsideeq(:)
      integer, allocatable ::  epart(:), npart(:)
      real, allocatable ::  Area(:), alfa(:)
      real, allocatable ::  sdep(:),slen(:),refdep(:),sdx(:),sdy(:)
      real, allocatable ::  dlinv(:)
      real, allocatable ::  z(:),zdep(:),zdepC(:)
      real, allocatable ::  eta(:), qp(:) 

! 2D data
      integer, allocatable ::  nen(:,:), ieadj(:,:)
      integer, allocatable ::  numsideT(:,:), iside(:,:), isidebc(:,:),iends(:,:) !numside(:,:),
      integer, allocatable ::  elemapG2L(:,:), nodemapG2L(:), sidemapG2L(:)
      integer, allocatable ::  elemapL2G(:,:), nodemapL2G(:,:), sidemapL2G(:,:)
      real, allocatable ::  ampl(:,:),phase(:,:)
      real, allocatable ::  sxy(:,:), xyz(:,:), FRC(:,:)
      real, allocatable ::  un(:,:), ut(:,:), wn(:,:), ut1(:,:), ut2(:,:) !AB3 Coriolis

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
