
    Module MainArrays

!     Variable names

!  np       - number of vertices (nodes) in grid
!  ne       - number of elements in grid
!  nsides   - number of edges in grid
!  X(I) - x-coordordinate of vertex I 
!  Y(I) - y-coordordinate of vertex I
!  Z(I) - bottom elevation of vertex I
!  NCODE(I)  - computational code for nodes
!  nen(i,j) - element list 
!  ECODE(I)  - computational code for elements
      
!  scalar arrays
      INTEGER np, ne, nsides, nnbr  
      integer igridtype
      real*8 x0off, y0off, scaleX, scaleY

!  1D arrays
      integer, allocatable ::  ECode(:),NCode(:)
      real, allocatable ::  x(:),y(:),z(:)

! 2D arrays
      integer, allocatable ::  nen(:,:), NL(:,:)

    End Module
