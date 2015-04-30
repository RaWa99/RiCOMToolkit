!! bin2netcdf
!!
!! FORTRAN90 program that converts a binary grid file to a NetCDF
!! file.
!!
!! Mark Cheeseman, NIWA 
!! Feburary, 2012
!!=====================================================================

  PROGRAM bin2netcdf

     use netcdf
     implicit none

     character(len=40) :: GridFileName
     character(len=43) :: filename
  
     integer :: NE, NP, NTYPE, nsides, NCN, NPR, nnbr, izup, ifront, &
                maxrow, maxsto, j, k, ierr, ncid
     integer, dimension(8)                :: dimID
     integer, dimension(14)               :: varID
     integer, dimension(:),   allocatable :: IEcode, nbc
     integer, dimension(:,:), allocatable :: nen, ieadj, numsideT, iside

     real(kind=4)                              :: alfa
     real(kind=4), dimension(:),   allocatable :: area, refdep, slen, &
                                                  dlinv, sdx, sdy
     real(kind=4), dimension(:,:), allocatable :: xyz, sxy

  !
  ! Set the name of the binary gridfile to be read
  !--------------------------------------------------------------------
      GridFileName = 'nzlam061123f_2m_up.lel'
      write(*,*) "converting grid file: ", GridFileName

  !
  ! Open the binary grid file and read in important sizing and 
  ! configuration variables
  !--------------------------------------------------------------------
      open( unit=22, file=GridFileName, status='OLD', &
            form='UNFORMATTED')      
      read(22) NE, NTYPE, NP, NPR, NCN, nsides, nnbr, izup, ifront, &
               maxrow, maxsto

  !
  ! Allocate memory for the working arrays and variables
  !--------------------------------------------------------------------
      allocate( nen(NE,NCN) )
      allocate( IEcode(NE) )
      allocate( xyz(NP,3) )
      allocate( nbc(NP) )
      allocate( area(ne) )
      allocate( ieadj(5,ne) )
      allocate( numsideT(4,ne) )
      allocate( iside(2,nsides) )
      allocate( sxy(2,nsides) ) 
      allocate( refdep(nsides), slen(nsides) )
      allocate( sdx(nsides),sdy(nsides),dlinv(nsides) )

  !
  ! Read in some important working arrays
  !--------------------------------------------------------------------
      read(22) ((nen(j,k),k=1,NCN),j=1,NE), &
               (IEcode(j),j=1,NE),          &
               ((xyz(j,k),j=1,NP),k=1,3),   &
               (alfa,j=1,NP),               &
               (nbc(j),j=1,NP)

      if ( nsides>0 ) then
         read(22) (area(j),j=1,ne),            &
                  ((ieadj(k,j),k=1,5),j=1,ne), &
                  ((numsideT(k,j),j=1,ne),k=1,4)
         read(22) ((iside(k,j),k=1,2),j=1,nsides),&
                  ((sxy(k,j),k=1,2),j=1,nsides)
         read(22) (refdep(j),j=1,nsides), &
                  (slen(j),j=1,nsides)
         read(22) (sdx(j),j=1,nsides), &
                  (sdy(j),j=1,nsides), &
                  (dlinv(j),j=1,nsides)
      endif

      close(22)
 
  !
  ! Create the NetCDF grid definition file
  !--------------------------------------------------------------------
      filename = trim(GridFileName)//'.nc'
      ierr = nf90_create( filename, nf90_clobber, ncid )

  !
  ! Add some global attributes 
  !--------------------------------------------------------------------
      ierr = nf90_put_att( ncid, nf90_global, "title", &
                           "RiCOM Grid File" )
      ierr = nf90_put_att( ncid, nf90_global, "num_vertices", NP)
      ierr = nf90_put_att( ncid, nf90_global, "num_elements", NE)
      if ( ncn==3 ) then
         ierr = nf90_put_att( ncid, nf90_global, "element_shape", &
                              "triangle" )
      elseif ( ncn==4 ) then
         ierr = nf90_put_att( ncid, nf90_global, "element_shape", &
                              "quadrangles" )
      endif
      ierr = nf90_put_att( ncid, nf90_global, &
                           "max_num_neighbors_per_vertex", &
                           nnbr )
      if ( izup==0 ) then
         ierr = nf90_put_att( ncid, nf90_global, "z_positive_dir", &
                             "downwards" )
      elseif ( izup==1 ) then
         ierr = nf90_put_att( ncid, nf90_global, "z_positive_dir", &
                             "upwards" )
      endif
      ierr = nf90_put_att( ncid, nf90_global, "IFRONT", ifront )
      ierr = nf90_put_att( ncid, nf90_global, "NTYPE",  ntype )

  !
  ! Create the required dimensions
  !--------------------------------------------------------------------
      ierr = nf90_def_dim( ncid, 'NE',         ne, dimID(1) )
      ierr = nf90_def_dim( ncid, 'NCN',       ncn, dimID(2) )
      ierr = nf90_def_dim( ncid, 'NSIDES', nsides, dimID(3) )
      ierr = nf90_def_dim( ncid, 'NP',         np, dimID(4) )
      ierr = nf90_def_dim( ncid, 'dim2',        2, dimID(5) )
      ierr = nf90_def_dim( ncid, 'dim3',        3, dimID(6) )
      ierr = nf90_def_dim( ncid, 'dim4',        4, dimID(7) )
      ierr = nf90_def_dim( ncid, 'dim5',        5, dimID(8) )

  !
  ! Define the necessary variable datasets
  !--------------------------------------------------------------------
      ierr = nf90_def_var( ncid, "IEcode",  nf90_int, (/dimID(1)/), &
                           varID(1) )
      ierr = nf90_put_att( ncid, varID(1), "long_name", &
                           "element code number" )

      ierr = nf90_def_var( ncid, "NBC",     nf90_int, (/dimID(4)/), &
                           varID(2) )
      ierr = nf90_put_att( ncid, varID(2), "long_name", &
                           "boundary codes for each element" )

      ierr = nf90_def_var( ncid, "Area",   nf90_real, (/dimID(1)/), &
                           varID(3) )
      ierr = nf90_put_att( ncid, varID(3), "long_name", &
                           "area of each element" )

      ierr = nf90_def_var( ncid, "refdep", nf90_real, (/dimID(3)/), &
                           varID(4) )
      ierr = nf90_put_att( ncid, varID(4), "long_name", &
                           "reference depth for each side of every element" )

      ierr = nf90_def_var( ncid, "slen",   nf90_real, (/dimID(3)/), &
                           varID(5) )
      ierr = nf90_put_att( ncid, varID(5), "long_name", &
                           "length of each element side" )

      ierr = nf90_def_var( ncid, "sdx",    nf90_real, (/dimID(3)/), &
                           varID(6) )
      ierr = nf90_put_att( ncid, varID(6), "long_name", &
               "X components of the unit vector from each element side" )

      ierr = nf90_def_var( ncid, "sdy",    nf90_real, (/dimID(3)/), &
                           varID(7) )
      ierr = nf90_put_att( ncid, varID(7), "long_name", &
               "Y components of the unit vector from each element side" )

      ierr = nf90_def_var( ncid, "dlinv",  nf90_real, (/dimID(3)/), &
                           varID(8) )
      ierr = nf90_put_att( ncid, varID(8), "long_name", &
               "inverse normal distance between 2 neighboring elements" )

      ierr = nf90_def_var( ncid, "NEN",    nf90_int,  (/dimID(1:2)/),&
                           varID(9) )
      ierr = nf90_put_att( ncid, varID(9), "long_name", &
               "vertex numbering" )

      ierr = nf90_def_var( ncid, "XYZ",    nf90_real, &
                           (/dimID(4),dimID(6)/), varID(10) )
      ierr = nf90_put_att( ncid, varID(10), "long_name", &
               "X,Y,Z position of each vertex" )

      ierr = nf90_def_var( ncid, "IEADJ",   nf90_int, &
                           (/dimID(8),dimID(1)/), varID(11) )
      ierr = nf90_put_att( ncid, varID(11), "long_name", &
               "element adjacency list" )

      ierr = nf90_def_var( ncid, "numsideT",nf90_int, &
                           (/dimID(7),dimID(1)/), varID(12) )
      ierr = nf90_put_att( ncid, varID(12), "long_name", &
               "side numbering for each element" )

      ierr = nf90_def_var( ncid, "iside",   nf90_int, &
                           (/dimID(5),dimID(3)/), varID(13) )
      ierr = nf90_put_att( ncid, varID(13), "long_name", &
               "element numbering per side" )

      ierr = nf90_def_var( ncid, "SXY",    nf90_real, &
                           (/dimID(5),dimID(3)/), varID(14) )
      ierr = nf90_put_att( ncid, varID(14), "long_name", &
               "X,Y coordinates of midpoint of each side" )
      ierr = nf90_enddef( ncid )

  !
  ! Write the necessary data to hard disk
  !--------------------------------------------------------------------
      ierr = nf90_put_var( ncid, varID(1), IEcode )
      ierr = nf90_put_var( ncid, varID(2), nbc )
      ierr = nf90_put_var( ncid, varID(3), Area )
      ierr = nf90_put_var( ncid, varID(4), refdep )
      ierr = nf90_put_var( ncid, varID(5), slen )
      ierr = nf90_put_var( ncid, varID(6), sdx )
      ierr = nf90_put_var( ncid, varID(7), sdy )
      ierr = nf90_put_var( ncid, varID(8), dlinv )
      ierr = nf90_put_var( ncid, varID(9), nen )
      ierr = nf90_put_var( ncid, varID(10), xyz )
      ierr = nf90_put_var( ncid, varID(11), ieadj )
      ierr = nf90_put_var( ncid, varID(12), numsideT )
      ierr = nf90_put_var( ncid, varID(13), iside )
      ierr = nf90_put_var( ncid, varID(14), sxy )

      ierr = nf90_close( ncid )

  !
  ! Free allocated memory
  !--------------------------------------------------------------------
      deallocate( nen,IEcode,xyz,nbc )
      deallocate( area,ieadj,numsideT )
      deallocate( iside,sxy )
      deallocate( refdep,slen,sdx,sdy,dlinv )

  END PROGRAM bin2netcdf
