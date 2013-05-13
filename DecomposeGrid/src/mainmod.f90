!******************************************************************************
!******************************************************************************

      MODULE SISL3DMOD

      SAVE

! Parameters
      real, parameter :: GRAV = 9.80621, cpheat = 4187.
      real,parameter :: Nuw=0.0000011
      real,parameter :: kappa=0.4
      real,parameter :: rhow=1020., rhos=2650.
      real,parameter :: bigr = 6378136.

! Scalar data
      integer  NE,NP,NPR,NDF,NCN,ntype,ntypes,ntypec,ndfe2,nsides,nsidesbc
      integer  ICVG,IDX,NPRT,itest,nitn,ncon,nsed,jsed,nsol,nnbr,nbx,npMB,neMB
      integer  npv, npvc,izcoord,izgrid
      integer  NOPT,NEQS,idpdn,isolve,NSon,iqp,nsbc,irst,irstout,nrstout,icase,nitsol,NMit
      integer ifr, itn, iwn, ivfr, ivSf, ihfr, iFDrag, irampa, irampq
      integer iwind, iwindx,iAnalysisTime,neqtide
      integer ntsdata,ntsskip,ntsUdata,ntsUskip,ntsCdata,ntsCskip,ntsSdata,ntsSskip
      integer noptskip,noptstart
      integer jPProfile, jUProfile, jCProfile, jSProfile
      integer icoord,izup,ifront,maxrow,maxsto,ielmSR
      real TMAX,DELT,dt0,tscale,CT,TET,depmin,gamma0,SXR,SYR,theta,theta0,thetaq, zminq
      real Tload, drhodx, drhody
      real elev  !, fcor
!      real tsx,tsy,ps  ! data for expansion of met input
      real rad2deg, deg2rad
      real*8 lat0, long0
! 1D data
      integer, allocatable ::  nbc(:),IndexQ(:), IndexQ1(:), IndexQ2(:), IECode(:),numsideeq(:)
      integer, allocatable ::  iepart(:)
      real, allocatable ::  Area(:),wetarea(:),mf1(:),fcor(:)
      real, allocatable ::  sdep(:),slen(:),refdep(:),sdx(:),sdy(:)
      real, allocatable ::  dlinv(:),gamma(:),gammaFD(:) !,sv(:)
      real, allocatable ::  z(:),zdep(:),zdepC(:)
      real, allocatable ::  eta(:), qp(:) 
      real, allocatable ::  SPEC(:),SPEC0(:),TPeriod(:)  !ampl(:),phase(:)
      real, allocatable ::  tsx(:),tsy(:),ps(:)  ! data for expansion of met input
      real, allocatable ::  Qnext(:),Qpast(:),Qstart(:),Qend(:)
      real, allocatable ::  QSedStart(:),QSedEnd(:),SedSPEC(:)
      real, allocatable ::  QConStart(:),QConEnd(:)!,ConSPEC(:)
      real, allocatable ::  eqtide(:)

! 2D data
      integer, allocatable ::  nen(:,:), ieadj(:,:)
      integer, allocatable ::  numsideT(:,:), iside(:,:), isidebc(:,:),iends(:,:) !numside(:,:),
      real, allocatable ::  ampl(:,:),phase(:,:)
      real, allocatable ::  sxy(:,:), xyz(:,:), FRC(:,:)
      real, allocatable ::  un(:,:), ut(:,:), wn(:,:), ut1(:,:), ut2(:,:) !AB3 Coriolis
      real, allocatable ::  Rn0(:,:), dz(:,:), AIZ(:,:)
      real, allocatable ::  CulSto(:,:), CSed(:,:)
      real, allocatable ::  cc(:,:), cS(:,:), Dhs(:,:), Dhc(:,:)
      real, allocatable ::  ConSPEC(:,:)

! Solution arrays
      integer, allocatable ::  NKE(:)
      real, allocatable ::  CME(:,:), FE(:)
      real, allocatable ::  rhv(:)

! PCG Solver
      real, allocatable :: ae(:),r(:),p(:),Mp(:),au(:)

! Character data
      character*80 OutResFile
      character*(128) GridFileName, WindFileName
      character(120) :: WindSourceFileName
      character(18) :: AnalysisTime
      character(30) :: rstFile
      character(2) :: nrstchar
      character(15) :: version

! Time series data
!      integer, parameter :: maxnpts=10010, maxcolts=15
      real, allocatable :: tsdat(:,:), tsUdat(:,:), tsCdat(:,:), tsSdat(:,:)
      integer, allocatable :: tsnodes(:), tsUnodes(:), tsCnodes(:), tsSnodes(:)
      integer, allocatable :: tsUlevel(:), tsClevel(:), tsSlevel(:)


      END

!******************************************************************************
!******************************************************************************

      MODULE ProfileArrays

      SAVE

! 1D data
      real, allocatable :: zz(:),A(:),B(:),C(:),D1(:),D2(:),Avp(:),FDrag(:),Kvp(:)

      END

!******************************************************************************
!******************************************************************************
