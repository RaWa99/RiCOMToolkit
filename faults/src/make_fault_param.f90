 PROGRAM make_fault_param

INTEGER NF,NX,NY,i,j,unitin,unitout,code,Ntot
real*8 dx,dy,latend,lonend,depthend,slipcm,alpha,mu,lambda,loncen,latcen,depthcen
character(256) :: infile,outfile,g1,g2,g3,g4,g5,g6
! 1D datatxt 
real*8 :: edgelon(5),edgelat(5),edgedepth(5)
! 2D data
real*8, allocatable :: latedge(:,:),lonedge(:,:),depthedge(:,:)
real*8, allocatable :: slip(:,:),rake(:,:),strike(:,:),dip(:,:)


write(*,'(a)') 'Name of input file?: '
read(*,'(a)') infile
!write(*,'(a)') 'Name of output file?: '
!read(*,'(a)') outfile
outfile='fault.param'


lambda=2.d10
mu=2.d10

unitin=13
unitout=14
open(unit=unitin,file=trim(infile),status='OLD',form='FORMATTED')

open(unit=unitout,file=trim(outfile),status='REPLACE',form='FORMATTED')

read(unitin,6011) NF
Ntot=0
do i=1,NF
   read(unitin,'(a16,i4,a18,i4,a4,f6.2,a15,i4,a4,f6.2,a2)',IOSTAT=code) g1,j,g2,NX,g3,dx,g4,NY,g5,dy,g6
   if (i.ne.j) then
      write(*,*) 'Problem, mismatch with fault number'
   endif
   read(unitin,'(a26,i6)',IOSTAT=code) g1,j  ! Check 
   if (i.ne.j) then
      write(*,*) 'Problem, mismatch with fault number'
   endif
   read(unitin,'(a)') g1
   Ntot=NX*NY
   allocate(latedge(Ny+1,Nx+1),lonedge(Ny+1,Nx+1),depthedge(Ny+1,Nx+1))
   allocate(slip(Ny,Nx),rake(Ny,Nx),strike(Ny,Nx),dip(Ny,Nx))
   write(unitout,'(i6)') Ntot
   do j=1,5
      read(unitin,'(3(f15.5))') edgelon(j),edgelat(j),edgedepth(j)
   enddo
   if (edgelon(1).ne.edgelon(5)) write(*,*) 'mismatch with edge longitude'
   if (edgelat(1).ne.edgelat(5)) write(*,*) 'mismatch with edge longitude'
   if (edgedepth(1).ne.edgedepth(5)) write(*,*) 'mismatch with edge longitude'
   read(unitin,'(a)') g1
   do j=1,Ny
      do k=1,Nx
         read(unitin,*) latedge(j,k), lonedge(j,k),depthedge(j,k),slipcm,rake(j,k),strike(j,k),dip(j,k)
         slip(j,k)=slipcm/1.d2
      enddo
      k=Nx+1
      alpha=dble(Ny+1-j)/dble(Ny)
      latedge(j,k)=edgelat(2)*alpha+edgelat(3)*(1.d0-alpha)
      lonedge(j,k)=edgelon(2)*alpha+edgelon(3)*(1.d0-alpha)
      depthedge(j,k)=edgedepth(2)*alpha+edgedepth(3)*(1.d0-alpha)
   enddo
   j=Ny+1
   do k=1,Nx
      alpha=dble(Nx+1-k)/dble(Nx)
      latedge(j,k)=edgelat(4)*alpha+edgelat(3)*(1.d0-alpha)
      lonedge(j,k)=edgelon(4)*alpha+edgelon(3)*(1.d0-alpha)
      depthedge(j,k)=edgedepth(4)*alpha+edgedepth(3)*(1.d0-alpha)
   enddo
   k=Nx+1
   latedge(j,k)=edgelat(3)
   lonedge(j,k)=edgelon(3)
   depthedge(j,k)=edgedepth(3)
   do j=1,Ny
      do k=1,Nx
         loncen=0.25d0*(lonedge(j,k)+lonedge(j+1,k)+lonedge(j,k+1)+lonedge(j+1,k+1))
         if (loncen .lt. 0.d0) then
            loncen=loncen+3.6d2
         endif
         latcen=0.25d0*(latedge(j,k)+latedge(j+1,k)+latedge(j,k+1)+latedge(j+1,k+1))
         depthcen=0.25d0*(depthedge(j,k)+depthedge(j+1,k)+depthedge(j,k+1)+depthedge(j+1,k+1))
         write(unitout,6012) latcen, loncen,dx*1.d3,dy*1.d3,depthcen*1.d3,dip(j,k),strike(j,k),lambda,mu,slip(j,k),rake(j,k)
      enddo
   enddo
   deallocate(latedge,lonedge,depthedge,slip,rake,strike,dip)
enddo
close(unitin)
close(unitout)
write(*,'(a,i6)') 'Ntot = ',Ntot
if (NF .gt. 1) then
   write(*,'(a,i6,a)') 'There are ',NF,' faults.  The fault param file will need to be adjusted'
endif

6011 FORMAT('#Total number of fault_segments=    ',i)
6012 FORMAT(f12.5,f12.5,f10.2,f10.2,f10.2,f8.2,f8.2,es10.2,es10.2,f8.2,f8.2)
END PROGRAM make_fault_param



