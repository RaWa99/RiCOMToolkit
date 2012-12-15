      program reorder
c
c
      parameter ( maxnp=1000000, maxne=2*maxnp )
      integer itime, e2, snode, nc, oldpro, newpro, kdum
      integer xadj(maxnp+1), adj(16*maxnp), mask(maxnp)
      integer ls(maxnp), xls(maxnp+1), hlevel(maxnp)
      integer nnn(maxnp), iw(3*maxnp+1)
      integer npn(3*maxne), xnpn(maxne+1)
      integer nen(maxne), ien(maxne), npe(maxne)
      integer in(maxne,4), mtrl(maxne)
      real et
      character*80 filnam, timestr*11, ans*1
      logical exis, all, feswms
	logical(4) resOK,openfile
c
c   get filenames for incidence list and reordered lists
c
      write(*,*) ' This program reorders nodes and/or elements'
      write(*,*) '  for minimum bandwidth and rms bandwidth'
      write(*,*) ' max nodes=',maxnp,' max elements=',maxne
c      itime = mclock()
c      et = float(itime)/100.
   13 continue
      feswms = .false.
      np = 0
      ne = 0
c
c *** read input data
340   continue
      write(*,*)'Read element list '
c      read(5,99,err=13,end=13) filnam

c      INQUIRE (FILE=FilNam, EXIST=Exis)
c      IF (Exis) THEN
c      else
c        PRINT *,'FILE DOES NOT EXIST - QUIT(Y/N)?'
c        READ(5,'(A)') ANS
c        IF (ANS.EQ.'Y'.or.ans.eq.'y') then
c          stop
c        else
c          GOTO 340
c        endif
c      ENDIF
c
c      open(unit=12,file=filnam,status='old')

      ResOK = OpenFile(12,'Open Input Element File',filnam,
     &  "Element file(*.el*),*.el*;All files(*.*),*.*;")
c
440   continue
      write(*,*)'Browse for re-ordered element filename '
c      read(5,99,err=13,end=13) filnam
c      INQUIRE (FILE=FilNam, EXIST=Exis)
c      IF (Exis) THEN
c        PRINT *,'FILE ALREADY EXISTS - RENAME(Y/N)?'
c        READ(5,'(A)') ANS
c        IF (ANS.EQ.'Y'.or.ans.eq.'y') then
c          GOTO 440
c        endif
c      ENDIF
c
c      open(unit=14,file=filnam,status='new')

      ResOK = OpenFile(14,'Open NEW ReOrdered Element File',filnam,
     &  "Element file(*.elr),*.elr;All files(*.*),*.*;")
c
      read(12,99,err=31,end=31) filnam
C      rewind 12
      write(6,*) 'Your file has the following first line - ',
     1  'there should be either 3 or 5 integers'
      write(6,99) filnam
      go to 33
   31 continue
      write(*,*) ' Error reading file'
      go to 13
   33 continue
      rewind 12
c
      if(feswms) then
        inum = 0
        ntype = 3
        imat = 2
      else
C
C       This bit changed by Derek Goring Nov 12 1996
C        write(6,5)
C        read(5,*) inum,ntype,imat
C
C
        NINT=ICHECK(FILNAM)
	  write(*,*) ' There are',nint,
     *        ' integers in the above line'
        IF(NINT .EQ. 3)THEN
           IMAT=0
        ELSEIF(NINT .EQ. 5)THEN
           IMAT=2
        ELSE
           WRITE(6,*)'The triangle file is incorrect'
           STOP
        ENDIF
C
C	Set the default options for the other paras
C
	  INUM=0
	  NTYPE=3	  
	
      endif
c
      do j=1,maxne
        if(inum.eq.0 .and. imat.eq.0)then
          read(12,*,end=14,err=14) (in(j,i),i=1,ntype)
        elseif(inum.eq.0 .and. imat.eq.1)then
          read(12,*,end=14,err=14) (in(j,i),i=1,ntype),mtrl(j)
        elseif(inum.eq.0 .and. imat.eq.2)then
          read(12,*,end=14,err=14) (in(j,i),i=1,ntype+1),mtrl(j)
        elseif(inum.eq.1 .and. imat.eq.0)then
          read(12,*,end=14,err=14) k,(in(k,i),i=1,ntype)
        elseif(inum.eq.1 .and. imat.eq.1)then
          read(12,*,end=14,err=14) k,(in(k,i),i=1,ntype),mtrl(k)
        endif
        ne = ne + 1
      enddo
      go to 14
16    continue
      write(*,*) ' Error reading file'
      go to 13
   14 write(6,*)' Element file currently contains: ',ne,'entries'
      close(12)
      xnpn(1) = 1
      iadj = 0
	do i=1,ne
	  if(in(i,4).gt.0) then
	    ntype = 4
	  else
	    ntype = 3
	  endif
        npe(i)=ntype
        xnpn(i+1) = xnpn(i) + ntype
c	  xnpn(i) = ntype*(i-1) + 1
        do j=1,ntype
c          npn(ntype*(i-1)+j) = in(i,j)
          npn(xnpn(i)+j-1) = in(i,j)
        enddo 
        i4 = in(i,4)
        if(ntype.eq.3) then
          np = max(np,in(i,1),in(i,2),in(i,3))
	  else
          np = max(np,in(i,1),in(i,2),in(i,3),in(i,4))
	  endif
	  iadj = iadj + ntype*(ntype-1)
      enddo
      write(*,*) ' Maximum node number:', np
      nnp = np
      nne = ne
c      xnpn(ne+1) = xnpn(ne) + ntype
      inpn = xnpn(ne+1)-1
c      iadj = ne*ntype*(ntype-1)
      all = .true.
	WRITE(6,*)'Reordering nodes'
      call graph(nnp,nne,inpn,npn,xnpn,iadj,adj,xadj)
c
      e2 = xadj(np+1) -1
      call diamtr(nnp,e2,adj,xadj,mask,ls,xls,hlevel,snode,nc)
c
      call label(nnp,e2,adj,xadj,nnn,iw,oldpro,newpro)
c
      write(*,*) ' Old profile width =', oldpro
      write(*,*) ' New profile width =', newpro
c
      call reseq2(nnn,npn,xnpn,nen,ien,npe,ntype,nne,nnp,ntype,all)
c
      kdum=0
	ntype = 3
      do j=1,ne
        jj = ien(j)
        if(inum.eq.0 .and. imat.eq.0)then
c          jj = ien(j)
          write(14,140)(in(jj,i),i=1,ntype)
        elseif(inum.eq.0 .and. imat.eq.1)then
          write(14,140) (in(jj,i),i=1,ntype), mtrl(jj)
        elseif(inum.eq.0 .and. imat.eq.2)then
          write(14,140) (in(jj,i),i=1,ntype+1), mtrl(jj)
        elseif(inum.eq.1 .and. imat.eq.0)then
          write(14,140) j,(in(jj,i),i=1,ntype)
        elseif(inum.eq.1 .and. imat.eq.1)then
          write(14,140) j,(in(jj,i),i=1,ntype),mtrl(jj)
        endif
      enddo
      close(14)
c
c      itime = mclock()
c      et = float(itime)/100. - et
c      write(*,*) ' elapsed time= ', et, ' sec'
c
    5 format(1x,' enter desired element list format of i, j, k where'
     *,/,12x,'i = 1 for element number included, 0 otherwise.'
     *,/,12x,'j = maximum number of nodes per element and'
     *,/,12x,'k = 2 for FESWMS format file,'
     *,/,12x,'k = 1 for material parameter included, 0 otherwise.')
    7 format(2(2x,i4))
    8 format(a8)
    9 format(1x,'your datafile has an unacceptable format. select anothe
     *r.')
   99 format(a80)
  140 format(10(1x,i6))
  141 format(3(2x,1pe13.6))
  142 format(2x,i5,3(2x,1pe13.6))
       end
************************************************************************
      subroutine graph(n,ne,inpn,npn,xnpn,iadj,adj,xadj)
************************************************************************
*
*     PURPOSE:
*     --------
*
*     Form adjacency list for a graph corresponding to a finite element
*     mesh
*
*     INPUT:
*     -----
*
*     N      - Number of nodes in graph (finite element mesh)
*     NE     - Number of elements in finite element mesh
*     INPN   - Length of NPN = XNPN(NE+1)-1
*     NPN    - List of node numbers for each element
*     XNPN   - Index vector for NPN
*            - nodes for element I are found in NPN(J), where
*              J = XNPN(I), XNPN(I)+1, ..., XNPN(I+1)-1
*     IADJ   - Length of vector ADJ
*            - Set IADJ = NE*NEN*(NEN-1) for a mesh of a single type of
*              element with NEN nodes
*            - IADJ = (NEN(1)*(NEN(1)-1)+,.....,+NEN(NE)*(NEN(NE)-1))
*              for a mesh of elements with varying numbers of nodes
*     ADJ    - Undefined
*     XADJ   - Undefined
*
*     OUTPUT:
*     -------
*
*     N      - Unchanged
*     NE     - Unchanged
*     INPN   - Unchanged
*     NPN    - Unchanged
*     XNPN   - Unchanged
*     IADJ   - Unchanged
*     ADJ    - Adjacency list for all nodes in graph
*            - List of length 2E where E is the number of edges in 
*              the graph (note that 2E = XADJ(N+1)-1 )
*     XADJ   - Index vector for ADJ
*            - Nodes adjacent to node I are found in ADJ(J), where
*              J = XADJ(I), XADJ(I)+1, ..., XADJ(I+1)-1
*            - Degree of node I given by XADJ(I+1)-XADJ(I)
*
*     NOTES:
*     ------
*
*     This routine typically requires about 25 percent elbow room for
*     assembling the ADJ list (i.e. IADJ/2E is typically around 1.25).
*     In some cases, the elbow room may be larger (IADJ/2E is slightly 
*     less than 2 for the 3-noded triangle) and in other cases it may be
*     zero (IADJ/2E = 1 for bar elements)
*
*     PROGRAMMER:             Scott Sloan
*     -----------
*
*     LAST MODIFIED:          1 March 1991        Scott Sloan
*     --------------
*
*     COPYRIGHT 1989:         Scott Sloan
*     ---------------         Department of Civil Engineering
*                             University of Newcastle
*                             NSW 2308
*
************************************************************************
      integer i,j,k,l,m,n
      integer ne
      integer iadj,inpn,nen1
      integer jstop,jstrt,lstop,lstrt,mstop,mstrt,nodej,nodek
      integer adj(iadj),npn(inpn)
      integer xadj(n+1),xnpn(ne+1)
*
*     Initialise the adjacency list and its index vector
*
c      write(*,*) ' Enter graph'
      do 8511 i = 1,iadj
         adj(i) = 0
 8511 continue
      do 8513 i = 1,n
        xadj(i) = 0
 8513 continue
*
*     Estimate the degree of each node (always an overestimate)
*
      do 8515 i = 1,ne
         jstrt = xnpn(i)
         jstop = xnpn(i+1)-1
         nen1 = jstop-jstrt
         do 8517 j = jstrt,jstop
            nodej = npn(j)
            xadj(nodej) = xadj(nodej)+nen1
 8517    continue
 8515 continue
*
*     Reconstruct XADJ to point to start of each set of neighbours
*
      l = 1
      do 8519 i = 1,n
         l = l+xadj(i)
         xadj(i) = l-xadj(i)
 8519 continue
      xadj(n+1) = l
*
*     Form adjacency list (which may contain zeros)
*
      do 8521 i = 1,ne
         jstrt = xnpn(i)
         jstop = xnpn(i+1)-1
         do 8523 j = jstrt,jstop-1
            nodej = npn(j)
            lstrt = xadj(nodej)
            lstop = xadj(nodej+1)-1
            do 8525 k = j+1,jstop
               nodek = npn(k)
               do 8527 l = lstrt,lstop
                  if(adj(l).eq.nodek)goto 70
                  if(adj(l).eq.0)goto 55
 8527          continue
               write(6,1000)
               stop
   55          continue
               adj(l) = nodek
               mstrt = xadj(nodek)
               mstop = xadj(nodek+1)-1
               do 8529 m = mstrt,mstop
                  if(adj(m).eq.0)goto 65
 8529          continue
               write(6,1000)
               stop
   65          continue
               adj(m) = nodej
   70          continue
 8525       continue
 8523    continue
 8521 continue
*
*     Strip any zeros from adjacency list
*
      k = 0
      jstrt = 1
      do 8531 i = 1,n
         jstop = xadj(i+1)-1
         do 8533 j = jstrt,jstop
            if(adj(j).eq.0) go to 8534
            k = k+1
            adj(k) = adj(j)
 8533    continue
 8534    continue
         xadj(i+1) = k+1
         jstrt = jstop+1
 8531 continue
*
*     Error message
*
 1000 format(//,1x,'***error in graph***', //,1x,'cannot assemble node a
     &djacency list', //,1x,'check npn and xnpn arrays')
      end
************************************************************************
      subroutine diamtr(n,e2,adj,xadj,mask,ls,xls,hlevel,snode,nc)
************************************************************************
*
*     PURPOSE:
*     --------
*
*     Find nodes which define a psuedo-diameter of a graph and store
*     distances from end node
*
*     INPUT:
*     ------
*
*     N      - The total number of nodes in the graph
*     E2     - Twice the number of edges in the graph  = XADJ(N+1)-1
*     ADJ    - Adjacency list for all nodes in the graph
*            - List of length 2E where E is the number of edges in 
*              the graph and 2E = XADJ(N+1)-1
*     XADJ   - Index vector for ADJ
*            - Nodes adjacent to node I are found in ADJ(J), where
*              J = XADJ(I), XADJ(I)+1, ...,XADJ(I+1)-1
*            - Degree of node I given by XADJ(I+1)-XADJ(I)
*     MASK   - Masking vector for graph
*            - Visible nodes have MASK = 0, node invisible otherwise
*     LS     - Undefined
*     XLS    - Undefined
*     HLEVEL - Undefined
*     SNODE  - Undefined
*     NC     - Undefined
*
*     OUTPUT:
*     ------
*
*     N      - Unchanged
*     E2     - Unchanged
*     ADJ    - Unchanged
*     XADJ   - Unchanged
*     MASK   - List of distances of nodes from the end node
*     LS     - List of nodes in the component
*     XLS    - Not used
*     HLEVEL - Not used
*     SNODE  - Starting node for numbering
*     NC     - The number of nodes in this component of graph
*
*     SUBROUTINES CALLED:  ROOTLS, ISORTI
*     -------------------
*
*     NOTE:      SNODE and ENODE define a pseudo-diameter
*     -----
*
*     PROGRAMMER:             Scott Sloan
*     -----------
*
*     LAST MODIFIED:          1 March 1991      Scott Sloan
*     --------------
*
*     COPYRIGHT 1989:         Scott Sloan
*     ---------------         Department of Civil Engineering
*                             University of Newcastle
*                             NSW 2308
*
************************************************************************
      integer i,j,n
      integer e2,nc
      integer node
      integer depth,enode,hsize,istop,istrt,jstop,jstrt,snode,width
      integer degree,ewidth,mindeg,sdepth
      integer ls(n)
      integer adj(e2),xls(n+1)
      integer mask(n),xadj(n+1)
      integer hlevel(n)
c
C      write(*,*) ' Enter diamtr'
*
*     Choose first guess for starting node by min degree
*     Ignore nodes that are invisible (MASK ne 0)
*
      mindeg = n
      do 8511 i = 1,n
         if(mask(i) .eq. 0)then
            degree = xadj(i+1)-xadj(i)
            if(degree.lt.mindeg)then
               snode = i
               mindeg = degree
            end if
         end if
 8511 continue
*
*     Generate level structure for node with min degree
*
      call rootls(n,snode,n+1,e2,adj,xadj,mask,ls,xls,sdepth,width)
*
*     Store number of nodes in this component
*
      nc = xls(sdepth+1)-1
*
*     Iterate to find start and end nodes
*     
   15 continue
*
*     Store list of nodes that are at max distance from starting node
*     Store their degrees in XLS
*
      hsize = 0
      istrt = xls(sdepth)
      istop = xls(sdepth+1)-1
      do 8513 i = istrt,istop
         node = ls(i)
         hsize = hsize+1
         hlevel(hsize) = node
         xls(node) = xadj(node+1)-xadj(node)
 8513 continue
*
*     Sort list of nodes in ascending sequence of their degree
*     Use insertion sort algorithm
*
      if(hsize.gt.1)call isorti(hsize,hlevel,n,xls)
*
*     Remove nodes with duplicate degrees
*
      istop = hsize
      hsize = 1
      degree = xls(hlevel(1))
      do 8515 i = 2,istop
         node = hlevel(i)
         if(xls(node) .ne. degree)then
            degree = xls(node)
            hsize = hsize+1
            hlevel(hsize) = node
         endif
 8515 continue
*
*     Loop over nodes in shrunken level
*
      ewidth = nc+1
      do 8517 i = 1,hsize
         node = hlevel(i)
*        
*        Form rooted level structures for each node in shrunken level
*        
         call rootls(n,node,ewidth,e2,adj,xadj,mask,ls,xls,depth,width)
         if(width.lt.ewidth)then
*           
*           Level structure was not aborted during assembly
*           
            if(depth.gt.sdepth)then
*              
*              Level structure of greater depth found
*              Store new starting node, new max depth, and begin 
*              a new iteration
*
               snode = node
               sdepth = depth
               goto 15
            endif
*           
*           Level structure width for this end node is smallest so far
*           store end node and new min width
*           
            enode = node
            ewidth = width
         end if
 8517 continue
*
*     Generate level structure rooted at end node if necessary
*
      if(node .ne. enode)then
         call rootls(n,enode,nc+1,e2,adj,xadj,mask,ls,xls,depth,width)
      endif
*
*     Store distances of each node from end node
*
      do 8519 i = 1,depth
         jstrt = xls(i)
         jstop = xls(i+1)-1
         do 8521 j = jstrt,jstop
            mask(ls(j)) = i-1
 8521    continue
 8519 continue
      end
************************************************************************
      subroutine isorti(nl,list,nk,key)
************************************************************************
*
*     PURPOSE:
*     --------
*
*     Order a list of integers in ascending sequence of their keys 
*     using insertion sort
*
*     INPUT:
*     ------
*
*     NL   - Length of LIST
*     LIST - A list of integers
*     NK   - Length of KEY (NK must be ge NL)
*     KEY  - A list of integer keys
*
*     OUTPUT:
*     -------
*
*     NL   - Unchanged
*     LIST - A list of integers sorted in ascending sequence of KEY
*     NK   - Unchanged
*     KEY  - Unchanged
*
*     NOTE:    Efficient for short lists only (NL lt 20)
*     -----
*
*     PROGRAMMER:             Scott Sloan
*     -----------
*
*     LAST MODIFIED:          1 March 1991     Scott Sloan
*     --------------
*
*     COPYRIGHT 1989:         Scott Sloan
*     ---------------         Department of Civil Engineering
*                             University of Newcastle
*                             NSW 2308
*
************************************************************************
      integer i,j,t
      integer nl,nk
      integer value
      integer key(nk)
      integer list(nl)
*
C      write(*,*) ' Enter sorti'
      do 8511 i = 2,nl
         t = list(i)
         value = key(t)
         do 8513 j = i-1,1,-1
            if(value .ge. key(list(j)))then
               list(j+1) = t
               goto 20
            endif
            list(j+1) = list(j)
 8513    continue
         list(1) = t
   20    continue
 8511 continue
      end
************************************************************************
      subroutine label(n,e2,adj,xadj,nnn,iw,oldpro,newpro)
************************************************************************
*
*     PURPOSE:
*     --------
*
*     Label a graph for small profile and rms wavefront
*
*     INPUT:
*     ------
*
*     N      - Total number of nodes in graph
*     E2     - Twice the number of edges in the graph = XADJ(N+1)-1
*     ADJ    - Adjacency list for all nodes in graph
*            - List of length 2E where E is the number of edges in
*              the graph and 2E = XADJ(N+1)-1
*     XADJ   - Index vector for ADJ
*            - Nodes adjacent to node I are found in ADJ(J), where
*              J = XADJ(I), XADJ(I)+1, ..., XADJ(I+1)-1
*            - Degree of node I given by XADJ(I+1)-XADJ(I)
*     NNN    - Undefined
*     IW     - Undefined
*     OLDPRO - Undefined
*     NEWPRO - Undefined
*
*     OUTPUT:
*     -------
*
*     N      - Unchanged
*     E2     - Unchanged
*     ADJ    - Unchanged
*     XADJ   - Unchanged
*     NNN    - List of new node numbers
*            - New number for node I given by NNN(I)
*            - If original node numbers give a smaller profile then
*              NNN is set so that NNN(I) = I for I = 1,N
*     IW     - Not used
*     OLDPRO - Profile using original node numbering
*     NEWPRO - Profile for new node numbering
*            - If original profile is smaller than new profile, then
*              original node numbers are used and NEWPRO = OLDPRO
*
*     SUBROUTINES CALLED:  DIAMTR, NUMBER, PROFIL
*     -------------------
*
*     PROGRAMMER:             Scott Sloan
*     -----------
*
*     LAST MODIFIED:          1 March 1991        Scott Sloan
*     --------------
*
*     COPYRIGHT 1989:         Scott Sloan
*     ---------------         Department of Civil Engineering
*                             University of Newcastle
*                             NSW 2308
*
***********************************************************************
      integer i,n
      integer e2,i1,i2,i3,nc
      integer snode
      integer lstnum,newpro,oldpro
      integer iw(3*n+1)
      integer adj(e2),nnn(n)
      integer xadj(n+1)
*
*     Set all new node numbers = 0
*     This is used to denote all visible nodes
*
C      write(*,*) ' Enter label'
      do 8511 i = 1,n
         nnn(i) = 0
 8511 continue
*
*     Define offsets
*
      i1 = 1
      i2 = i1+n
      i3 = i2+n+1
*
*     Loop while some nodes remain unnumbered
*
      lstnum = 0
 8513 if (lstnum.lt.n) then
*        
*        Find end points of p-diameter for nodes in this component
*        Compute distances of nodes from end node
*        
         call diamtr(n,e2,adj,xadj,nnn,iw(i1),iw(i2),iw(i3),snode,nc)
*        
*        Number nodes in this component
*
         call number(n,nc,snode,lstnum,e2,adj,xadj,nnn,iw(i1),iw(i2))
      go to 8513 
      endif
*
*     Compute profiles for old and new node numbers
*
      call profil(n,nnn,e2,adj,xadj,oldpro,newpro)
*
*     Use original numbering if it gives a smaller profile
*
      if(oldpro.lt.newpro)then
         do 8515 i = 1,n
            nnn(i) = i
 8515    continue
         newpro = oldpro
      end if
      end
************************************************************************
      subroutine number(n,nc,snode,lstnum,e2,adj,xadj,s,q,p)
************************************************************************
*
*     PURPOSE:
*     --------
*
*     Number nodes in component of graph for small profile and rms
*     wavefront
*
*     INPUT:
*     ------
*
*     N      - Number of nodes in graph
*     NC     - Number of nodes in component of graph
*     SNODE  - Node at which numbering starts
*     LSTNUM - Count of nodes which have already been numbered
*     E2     - Twice tne number of edges in the graph = XADJ(N+1)-1
*     ADJ    - Adjacency list for all nodes in graph
*            - List of length 2E where E is the number of edges in
*              the graph and 2E = XADJ(N+1)-1
*     XADJ   - Index vector for ADJ
*            - Nodes adjacent to node I are found in ADJ(J), where
*              J = XADJ(I), XADJ(I)+1, ..... , XADJ(I+1)-1
*     S      - List giving the distance of each node in this
*              component from the end node
*     Q      - List of nodes which are in this component
*            - Also used to store queue of active or preactive nodes
*     P      - Undefined
*
*     OUTPUT:
*     -------
*
*     N      - Unchanged
*     NC     - Unchanged
*     SNODE  - Unchanged
*     LSTNUM - Count of numbered nodes (input value incremented by NC)
*     E2     - Unchanged
*     ADJ    - Unchanged
*     XADJ   - Unchanged
*     S      - List of new node numbers 
*            - New number for node I is S(I)
*     Q      - Not used
*     P      - Not used
*
*     NOTES:
*     ------
*
*     S also serves as a list giving the status of the nodes
*     during the numbering process:
*     S(I) gt 0 indicates node I is postactive
*     S(I) =  0 indicates node I is active 
*     S(I) = -1 indicates node I is preactive
*     S(I) = -2 indicates node I is inactive
*     P is used to hold the priorities for each node
*
*     PROGRAMMER:             Scott Sloan
*     -----------
*
*     LAST MODIFIED:          1 March 1991    Scott Sloan
*     --------------
*
*     COPYRIGHT 1989:         Scott Sloan
*     ---------------         Department of Civil Engineering
*                             University of Newcastle
*                             NSW 2308
*
************************************************************************
      integer i,j,n
      integer e2,nc,nn,w1,w2
      integer nbr
      integer nxt,node,prty
      integer jstop,jstrt,istop,istrt,nabor,snode
      integer addres,lstnum,maxprt
      integer p(n),q(nc),s(n)
      integer adj(e2)
      integer xadj(n+1)
*
      parameter (w1 = 1, w2 = 2)
*
C      write(*,*) ' Enter number'
*     Initialise priorities and status for each node in this component
*     Initial priority = W1*DIST - W2*DEGREE     where:
*     W1     = a positive weight
*     W2     = a positive weight
*     DEGREE = initial current degree for node
*     DIST   = distance of node from end node
*
      do 8511 i = 1,nc
         node = q(i)
         p(node) = w1*s(node)-w2*(xadj(node+1)-xadj(node)+1)
         s(node) = -2
 8511 continue
*
*     Insert starting node in queue and assign it a preactive status
*     NN is the size of queue
*
      nn = 1
      q(nn) = snode
      s(snode) = -1
*
*     Loop while queue is not empty
*
 8513 if (nn .gt. 0) then
*        
*        Scan queue for node with max priority
*        
         addres = 1
         maxprt = p(q(1))
         do 8515 i = 2,nn
            prty = p(q(i))
            if(prty.gt.maxprt)then
               addres = i
               maxprt = prty
            end if
 8515    continue
*
*        NEXT is the node to be numbered nxt
*
         nxt = q(addres)
*        
*        Delete node NEXT from queue
*       
         q(addres) = q(nn)
         nn = nn-1
         istrt = xadj(nxt)
         istop = xadj(nxt+1)-1
         if(s(nxt).eq.-1)then
*           
*           Node NEXT is preactive, examine its neighbours
*           
            do 8517 i = istrt,istop
*              
*              Decrease current degree of neighbour by -1
*
               nbr = adj(i)
               p(nbr) = p(nbr)+w2
*              
*              Add neighbour to queue if it is inactive
*              assign it a preactive status
*
               if(s(nbr).eq.-2)then
                  nn = nn+1
                  q(nn) = nbr
                  s(nbr) = -1
               end if
 8517       continue
         end if
*        
*        Store new node number for node NEXT
*        Status for node NEXT is now postactive
*        
         lstnum = lstnum+1
         s(nxt) = lstnum
*        
*        Search for preactive neighbours of node NEXT
*        
         do 8519 i = istrt,istop
            nbr = adj(i)
            if(s(nbr).eq.-1)then
*              
*              Decrease current degree of preactive neighbour by -1
*              assign neighbour an active status
*              
               p(nbr) = p(nbr)+w2
               s(nbr) = 0
*              
*              Loop over nodes adjacent to preactive neighbour
*              
               jstrt = xadj(nbr)
               jstop = xadj(nbr+1)-1
               do 8521 j = jstrt,jstop
                  nabor = adj(j)
*                 
*                 Decrease current degree of adjacent node by -1
*                 
                  p(nabor) = p(nabor)+w2
                  if(s(nabor) .eq. -2)then
*                    
*                    Insert inactive node in queue with a preactive status
*                    
                     nn = nn+1
                     q(nn) = nabor
                     s(nabor) = -1
                  end if
 8521          continue
            end if
 8519    continue
      go to 8513 
      endif
      end
************************************************************************
      subroutine profil(n,nnn,e2,adj,xadj,oldpro,newpro)
************************************************************************
*
*     PURPOSE:
*     --------
*
*     Compute the profiles using both original and new node numbers
*
*     INPUT:
*     ------
*
*     N      - Total number of nodes in graph
*     NNN    - List of new node numbers for graph
*            - New node number for node I is given by NNN(I)
*     E2     - Twice the number of edges in the graph = XADJ(N+1)-1
*     ADJ    - Adjacency list for all nodes in graph
*            - List of length 2E where E is the number of edges in 
*              the graph and 2E = XADJ(N+1)-1
*     XADJ   - Index vector for ADJ
*            - Nodes adjacent to node I are found in ADJ(J), where
*              J = XADJ(I), XADJ(I)+1, ..., XADJ(I+1)-1
*            - Degree of node I given by XADJ(I+1)-XADJ(I)
*     OLDPRO - Undefined
*     NEWPRO - Undefined
*
*     OUTPUT:
*     -------
*
*     N      - Unchanged
*     NNN    - Unchanged
*     E2     - Unchanged
*     ADJ    - Unchanged
*     XADJ   - Unchanged
*     OLDPRO - Profile with original node numbering
*     NEWPRO - Profile with new node numbering
*
*     NOTE:      Profiles include diagonal terms
*     -----
*
*     PROGRAMMER:             Scott Sloan
*     -----------
*
*     LAST MODIFIED:          13 August 1991     Scott Sloan
*     --------------
*
*     COPYRIGHT 1989:         Scott Sloan
*     ---------------         Department of Civil Engineering
*                             University of Newcastle
*                             NSW 2308
*
***********************************************************************
      integer i,j,n
      integer e2
      integer jstop,jstrt
      integer newmin,newpro,oldmin,oldpro
      integer adj(e2),nnn(n)
      integer xadj(n+1)
*
C      write(*,*) ' Enter profil'
*     Set profiles and loop over each node in graph
*
      oldpro = 0
      newpro = 0
      do 8511 i = 1,n
         jstrt = xadj(i)
         jstop = xadj(i+1)-1
         oldmin = i
         newmin = nnn(i)
*        
*        Find lowest numbered neighbour of node I
*        (using both old and new node numbers)
*        
         do 8513 j = jstrt,jstop
            oldmin = min(oldmin,adj(j))
            newmin = min(newmin,nnn(adj(j)))
 8513    continue
*        
*        Update profiles
*        
         oldpro = oldpro+idim(i,oldmin)
         newpro = newpro+idim(nnn(i),newmin)
 8511 continue
*
*     Add diagonal terms to profiles
*
      oldpro = oldpro+n
      newpro = newpro+n
      end
************************************************************************
      subroutine rootls(n,root,maxwid,e2,adj,xadj,mask,ls,xls,depth,widt
     &   h)
************************************************************************
*
*     PURPOSE:
*     --------
*
*     Generate rooted level structure using a FORTRAN 77 implementation
*     of the algorithm given by George and Liu
*
*     INPUT:
*     ------
*
*     N      - Number of nodes
*     ROOT   - Root node for level structure
*     MAXWID - Max permissible width of rooted level structure
*            - Abort assembly of level structure if width is ge MAXWID
*            - Assembly ensured by setting MAXWID = N+1
*     E2     - Twice the number of edges in the graph = XADJ(N+1)-1
*     ADJ    - Adjacency list for all nodes in graph
*            - List of length 2E where E is the number of edges in 
*              the graph and 2E = XADJ(N+1)-1
*     XADJ   - Index vector for ADJ
*            - Nodes adjacent to node I are found in ADJ(J), where
*              J = XADJ(I), XADJ(I)+1, ..., XADJ(I+1)-1
*            - Degree of node I is XADJ(I+1)-XADJ(I)
*     MASK   - Masking vector for graph
*            - Visible nodes have MASK = 0
*     LS     - Undefined
*     XLS    - Undefined
*     DEPTH  - Undefined
*     WIDTH  - Undefined
*
*     OUTPUT:
*     -------
*
*     N      - Unchanged
*     ROOT   - Unchanged
*     MAXWID - unchanged
*     E2     - Unchanged
*     ADJ    - Unchanged
*     XADJ   - Unchanged
*     MASK   - Unchanged
*     LS     - List containing a rooted level structure
*            - List of length NC
*     XLS    - Index vector for LS
*            - Nodes in level I are found in LS(J), where
*              J = XLS(I), XLS(I)+1, ..., XLS(I+1)-1
*            - List of max length NC+1
*     DEPTH  - Number of levels in rooted level structure
*     WIDTH  - Width of rooted level structure
*
*     NOTE:  If WIDTH ge MAXWID then assembly has been aborted
*     -----
*
*     PROGRAMMER:             Scott Sloan
*     -----------
*
*     LAST MODIFIED:          1 March 1991      Scott Sloan
*     --------------
*
*     COPYRIGHT 1989:         Scott Sloan
*     ---------------         Department of Civil Engineering
*                             University of Newcastle
*                             NSW 2308
*
************************************************************************
      integer i,j,n
      integer e2,nc
      integer nbr
      integer node,root
      integer depth,jstop,jstrt,lstop,lstrt,lwdth,width
      integer maxwid
      integer ls(n)
      integer adj(e2),xls(n+1)
      integer mask(n),xadj(n+1)
*
C      write(*,*) ' Enter rootls'
*     Initialisation
*
      mask(root) = 1
      ls(1) = root
      nc = 1
      width = 1
      depth = 0
      lstop = 0
      lwdth = 1
 8511 if (lwdth.gt.0) then
*        
*        LWDTH is the width of the current level
*        LSTRT points to start of current level
*        LSTOP points to end of current level
*        NC counts the nodes in component
*        
         lstrt = lstop+1
         lstop = nc
         depth = depth+1
         xls(depth) = lstrt
*        
*        Generate next level by finding all visible neighbours
*        of nodes in current level
*        
         do 8513 i = lstrt,lstop
            node = ls(i)
            jstrt = xadj(node)
            jstop = xadj(node+1)-1
            do 8515 j = jstrt,jstop
               nbr = adj(j)
               if(mask(nbr) .eq. 0)then
                  nc = nc+1
                  ls(nc) = nbr
                  mask(nbr) = 1
               end if
 8515       continue
 8513    continue
*        
*        Compute width of level just assembled and the width of the
*        level structure so far
*        
         lwdth = nc-lstop
         width = max(lwdth,width)
*        
*        Abort assembly if level structure is too wide
*        
         if(width .ge. maxwid) goto 35
      go to 8511 
      endif
      xls(depth+1) = lstop+1
*
*     Reset MASK = 0 for nodes in the level structure
*
   35 continue
      do 8517 i = 1,nc
         mask(ls(i)) = 0
 8517 continue
      end
************************************************************************
      subroutine reseq2(newnum,npn,xnpn,nen,ien,npe,maxnod,
     *                  net,nodes,nvn,all)
c
      integer newnum(nodes),npn(3*net),nen(net),npe(net),ien(net)
      integer xnpn(net+1)
      logical all
c*****************************************************************
c    subroutine reseq2 - resequences element nubers to minimise
c                        the frontwidth. reorders the elements in
c                        an ascending sequence of their lowest 
c                        numbered nodes
c*****************************************************************
c   
c      nen   - array containing new element numbers. the address in
c              this array indicates the old element number, e.g.
c              nen(1)=6 means that the new number for old element one
c              is 6. dimension equal to net.
c
c      ien   - array containing initial element numbers. the address
c              in this array indicates the new element numbers, e.g.
c              ien(6)=1 means that the old number for new element six
c              is one. dimension equal to net.
c 
C      write(*,*) ' Enter reseq2'
	WRITE(*,*)'Re-ordering elements -- please wait'
       do 10 i=1,net
10     nen(i)=0
       kount=0
c
c  loop over each new node number
c  loop only over corner nodes if all=.false.
c
       do 40 i=1,nodes
c 
c  loop over each element
c  skip to next element if already numbered
c
       do 30 j=1,net
       if(nen(j).gt.0) go to 30
       nn=npe(j)
       if(.not.all)nn=nvn
c       i1=(j-1)*maxnod
       i1=xnpn(j)-1
c
c  loop over each node in element
c  use only corner nodes if all=.false.
c  assume that corner nodes are listed first in nodal definition
c  vectors if all=.false.
c
       do 20 k=1,nn
       n=npn(i1+k)
       n=newnum(n)
       if(n.ne.i) go to 20
       kount=kount+1
       nen(j)=kount
       ien(kount)=j
       if(kount.eq.net) go to 50
       go to 30
20     continue
30     continue
40     continue
50     return
       end
C
        INTEGER FUNCTION ICHECK(STRING)
C
C       Routine to check the number of integers in a string
C       and stop execution if they are not all integers
C
C       Uses the subroutine UNRAVEL
C
C       Written by Derek Goring Nov 12 1996
C
        character*80 string
        integer numtype,nint,nreal
        REAL X(100)
        INTEGER IX(100)
        CALL UNRAVEL(STRING,NUMTYPE,NINT,NREAL,IX,X)
        IF(NUMTYPE .NE. 1)THEN
                WRITE(6,*)'The triangle file contains non-integers'
                STOP
        ENDIF

        icheck = nint

        RETURN
        END
C-------------------------------------------------------------------
      SUBROUTINE UNRAVEL(teststr,numtype,nints,nreals,ix,x)
C     Written: 14 Oct 93 by Falconer Henry, based on routines by
C              Luc Cuypers and Todd Warnes.
C     Purpose: Extracts integers OR reals from character string 'teststr' 
C              - Assumes that string is all REAL*4 in F or E format  
C              or all INTEGER*4 , in free format. Fields can be
C              separated by a blank, a slash (/) or comma.
C              - Decides whether record contains reals or integers
C              by looking for a period (.) anywhere in the record.
C              If even one period is found on line, all numbers 
C              in record will be read as reals, any integer being
C              being converted to a real.
C              - Will also detect a blank record. 
C              - Returns :
C                 numtype = 0 if blank record
C                         = 1 if integers found
C                         = 2 if reals found
C                 nints   = number of integers found
C                 nreals  = number of reals found  
C                 ix(1),... ix(nints) = integers found
C                 x(1),... x(nreals)  = reals found

      real x(*)                      ! reals found in record
      integer ix(*)                  ! integers found in record
      integer nlr                    ! length of record
      integer realscan, intscan      ! integer functions
      integer nreals                 ! number of reals found in record
      integer nints                  ! number of integers found in record
      integer numtype                ! data type, returned by s/r unravel
      character*(*) teststr          ! current record as character string
      character*1 perchar            ! period character
      logical index                  ! Fortran function
      external blank1
      logical*1 blank1                 ! function
      external len_trim
      integer len_trim               ! function
      external isinsetc
      logical isinsetc               ! function
c
      perchar = '.'
c     check for blank record - if found, warn and quit
      numtype = -1  
      nints = 0
      nreals = 0
      if (blank1(teststr)) then
        numtype = 0
      else
c     check if any periods (.) in record;
        nlr = len_trim(teststr)
        if (isinsetc(perchar,.true.,teststr,nlr)) then        
c         there is a period in the string - assume string all reals
          nreals = realscan (teststr, x)
          numtype = 2
        else
c         there is no period in string - assume all integers
          nints = intscan (teststr,ix)
          numtype = 1
        endif
      endif
      return 
      end
C------------------------------------------------------------------------
      logical*1 function BLANK1 (string)
*
* Luc E. Cuypers, April 24, 1990.
*
* BLANK1 becomes true if STRING is all blanks.
*
* User must declare function as follows :
*
*    external   BLANK1
*    logical*1  BLANK1
*
* STRING may be any length.
*
      implicit      NONE
      character*(*) string
      character     BLANKCHAR, NULLCHAR, a
      parameter    (BLANKCHAR = CHAR (32), NULLCHAR = CHAR (0))
c      integer*2     MaxLen, i
      integer*4     MaxLen, i
*
* Determine length of input STRING
*
      MaxLen = LEN (string)
*
      BLANK1 = .false.
      do i = 1, MaxLen
        a = string(i:i)
        if (a .ne. BLANKCHAR  .and.  a. ne. NULLCHAR) return
      end do
      BLANK1 = .true.
      return
      END
C-------------------------------------------------------------------------
      integer function INTSCAN (string,  ivector)
*     input/output :      o        i        o
*     type          :     int      char   int array
*
* Scans a string for integer numbers. The number of integers found is returned
* via INTSCAN. The integers themselves are stored in VECTOR. Make sure it's
* declared large enough!
*
* Within STRING, integers must be separated by blanks, slash (/), or
* comma (,). No other characters (except a minus sign) may appear within 
* STRING!
*
      implicit      NONE
      character*(*) string
      integer       ivector(*)
*
      external      LEN_TRIM
      integer       LEN_TRIM
      external      ISINSETC, INTCHR, BLANK1
      logical       ISINSETC
      logical*1     BLANK1
      integer       INTCHR, slen, zero, count, i, j, l
      character     delimit(3)
*
      data  delimit /  ' ',',','/'  /
      data  zero    / 0 /
*
* i/j are the beginning & end positions of the numbers found in STRING
*
      slen  = LEN_TRIM (string)
      if (slen .eq. zero) then
          INTSCAN = zero
          return
      end if
      i     = zero
      j     = zero
      count = zero
      do l = 1, slen
         if (ISINSETC(string(l:l),.false.,delimit,3)) then
*            character is a delimiter.
             if (j.ne.zero) then
*                the 1st blank following a 'minus' sign is ignored...
                 if (BLANK1(string(l:l)).and.string(l-1:l-1).eq.'-') 
     &               go to 888
                 count = count + 1
                 ivector(count) = INTCHR(string(i:j))
                 i = zero
                 j = zero
             end if
             go to 888
           else
             if (i.eq.zero) then
                 i = l
                 j = l
               else
                 j = l
             end if
         end if
888      continue
      end do
*
* Last integer...
*
      if (j.ne.zero) then
          count = count + 1
          ivector(count) = INTCHR(string(i:j))
      end if
*
      INTSCAN = count
      return
      end
C-------------------------------------------------------------------------
*
*  Extracts numbers from string.
*  The numbers are stored in a real array (vector). The number of numbers
*  found is returned via the function value (integer).
*
*  Adapted from subroutine 'SCAN' (A. Blaskovich  06-jul-87) by L. Cuypers
*
      integer function REALSCAN (string, vector)
*     input/output :      o         i      o
*     type         :     int       char   real
*
      implicit      none
      external      upperc
      external      LEN_TRIM 		!added by RFH
      character*(*) string
      character     c,lastc,nextc,upperc,fmt*10
      logical       sos
      integer       digits,ptr(30,2),slen,i,j,k,zero /0/
      integer       LEN_TRIM		! added by RFH
      real          vector(*)
*
*  initialization
      slen      = LEN_TRIM(string)
      realscan  = zero
      if (slen .eq. zero) return
      sos       = .false.
      digits    = zero
      i         = zero
      c         = ' '
*
*  scan contents of string
10    lastc = c
      i     = i + 1
      if (i.le.slen + 1) then
         c = upperc(string(i:i))
         nextc = upperc(string(i+1:i+1))
         if (c.eq.'+'.or.c.eq.'-') then
            if (sos) then
               if (lastc.eq.'E'.or.lastc.eq.'D') goto 10
               goto 30
            endif
            goto 20
         else if (c.eq.' '.or.c.eq.',') then
            if (sos) goto 30
            goto 10
         else if (c.eq.'E'.or.c.eq.'D') then
            if (sos.and.(nextc.eq.'-'.or.nextc.eq.'+'
     &              .or.(nextc.ge.'0'.and.'9'.ge.nextc))) goto 10
            if (sos) goto 30
            goto 10
         else if ((c.ge.'0'.and.'9'.ge.c).or.c.eq.'.') then
            if (c.ne.'.') digits = digits + 1
            if (sos.and.c.eq.'.'.and.lastc.eq.'.') goto 30
            if (sos) goto 10
            goto 20
         else
            if (sos) goto 30
            goto 10
         endif
*
*  beginning of number detected
20       sos = .true.
         ptr(realscan+1,1) = i
         goto 10
*
*  end of number detected
30       if (digits.gt.zero) then
            digits    = zero
            sos       = .false.
            realscan  = realscan + 1
            ptr(realscan,2) = i - 1
            endif
         if ((c.eq.'+'.or.c.eq.'-'.or.c.eq.'.').and.i.ne.slen) goto 20
         goto 10
      endif
*
*  assign values to vector
      if (realscan.gt.zero) then
         do 60 i = 1, realscan
            j = ptr(i,1)
            k = ptr(i,2)
            write (fmt,50) k - j + 1
50          format('(f',i2,'.0)')
            read  (string(j:k),fmt) vector(i)
60       continue
      endif
      return
      end
C--------------------------------------------------------------------
      logical function ISINSETC (char,ignorecase,charset,n)
*     input/output :       o      i       i         i    i
*     type         :      log    char    log       char int
*                                                  array
*
* Is CHAR in the character set CHARSET ? N is the number of elements in
* CHARSET. Upper/lower case may be ignored (set IGNORECASE true).
*
      implicit  logical (a - z)
      external  UPPERC
      character char, upperchar, charset(*), UPPERC
      integer   n, i
      logical   ignorecase
*
      ISINSETC = .FALSE.
      if (ignorecase) then
          upperchar = UPPERC(char)
          do 10 i = 1, n
            if (upperchar.eq.UPPERC(charset(i))) then
                ISINSETC = .true.
                return
            end if
   10     continue
        else
          do 20 i = 1, n
            if (char.eq.charset(i)) then
                ISINSETC = .true.
                return
            end if
   20     continue
      end if
      return
      end
C------------------------------------------------------------------------
      integer function LEN_TRIM(AU)
*
* Purpose       : returns the position of the last non-blank character in
*                 a string.
*
* Written by    : T. Warnes  December 09, 1986
* Revised       :
*
      CHARACTER*(*) AU
      INTEGER I
*
      I = LEN(AU)
      if (I .eq. 0) then
        I = 1
        AU = ' '
        goto 20
      endif
10    if (AU(I:I) .ne. ' ' .OR. I .eq. 1) goto 20
        I = I - 1
        goto 10
20    continue
      LEN_TRIM = I
      return
      end
C-------------------------------------------------------------------
      CHARACTER*(*) FUNCTION UPPERC(STRING)
*
* Converts lowercase characters in a string to uppercase
*
* Calling pgm must declare : EXTERNAL        UPPERC
*                            CHARACTER*(len) UPPERC
*
      implicit      NONE
      character*(*) string
c      integer*2     Flen, Slen, i, ascii, offset, a, z
      integer*4     Flen, Slen, i, ascii, offset, a, z
      external len_trim				! added by RFH
      integer len_trim
      parameter    (a = 97, z = 122, offset = 32)
*
      Flen = LEN (upperc)
      Slen = LEN_TRIM (string)
      if (Slen.gt.Flen) stop
     &  'function UPPERC : argument string is longer than function'
*
      UPPERC = ' '
      do i = 1, Slen
         ascii = ICHAR (string(i:i))
         if ((ascii.ge.a).and.(ascii.le.z)) then
             UPPERC(i:i) = CHAR (ascii - offset)
           else
             UPPERC(I:I) = string(i:i)
         end if
      end do
      return
      end
C------------------------------------------------------------------------
      INTEGER FUNCTION INTCHR(STRING)
*     input/output :     o       i
*     type         :    int     char
*
* LUC E. CUYPERS, MAY 23, 1984.
*
* CONVERTS A (SUB)STRING FROM CHARACTER TO INTEGER.
* THE INPUT STRING MAY CONTAIN A MINUS OR PLUS SIGN.
* LEADING, TRAILING OR EMBEDDED BLANKS IN THE INPUT STRING ARE IGNORED.
* IN ADDITION, ANY CHARACTERS PRECEDING THE '+' OR '-' ARE IGNORED.
* CHARACTER DATA (OTHER THAN DIGITS 0 TO 9) FOLLOWING THE '+' OR '-'
* RESULT IN ERROR CONDITION.
*
* ERROR CONDITION RESULTS IN AN ERROR MESSAGE AND INTCHR = 99999
*
* CALLS FUNCTION COMPRS (REMOVES BLANKS) 
*
      CHARACTER*(*)  STRING
      EXTERNAL       COMPRS
      CHARACTER*20   COMPRS
      CHARACTER*20   WORK
      CHARACTER      A
      LOGICAL        MINUS
*
      INTEGER    ASCII0, ASCII9
      CHARACTER  BLANK
      PARAMETER (ASCII0 = 48,
     &           ASCII9 = 57,
     &           BLANK  = ' ',
     &           IFLAG  = -99999,
     &           LNP    = 6)
*
   50 FORMAT(' FUNCT. INTCHR : INPUT STRING CONTAINS CHARACTER DATA')
*
      MINUS = .FALSE.
      INTCHR = 0
      WORK = COMPRS(STRING,1,NCHAR)
*
      K    = 0
      K = INDEX(WORK,'-')
      IF (K.NE.0) THEN
          MINUS = .TRUE.
        ELSE
          K = INDEX(WORK,'+')
*         MINUS IS ALREADY .FALSE.
      END IF
*
* Check for garbage ahead of the '_' or '+' character
*
      IF (K.GT.1) THEN
          DO 5 I = 1,K-1
            IF (WORK(I:I).NE.BLANK) THEN
                WRITE (LNP,50)
                INTCHR = IFLAG
                RETURN
            END IF
    5     CONTINUE
      END IF
*
      ISTART = 1 + K
      DO 10 I = ISTART,NCHAR
          A   = WORK(I:I)
          IA  = ICHAR(A)
          IF ((IA.LT.ASCII0).OR.(IA.GT.ASCII9)) THEN
              WRITE (LNP,50)
              INTCHR = IFLAG
              RETURN
          END IF
          INTCHR = INTCHR*10 + (IA - ASCII0)
   10 CONTINUE
*
      IF (MINUS) INTCHR = - INTCHR
*
      RETURN
      END
C-----------------------------------------------------------------
      CHARACTER*(*) FUNCTION COMPRS(CHARAY,NIN,NOUT)
*     input/output :           o       i    i    o
*     type         :          char   char  int  int
*                                    array
*
* LUC E. CUYPERS, MAY 26,1984.
*
* FUNCTION COMPRESSES (REMOVES BLANKS FROM) INPUT ARRAY 'CHARAY'.
* 'NIN' ELEMENTS OF 'CHARAY' ARE COMPRESSED INTO 'COMPRS'.
*     CHARAY : (INPUT) CHARACTER ARRAY, ARRAY ELEMENT, OR (SUB)STRING.
*     NIN    : (INPUT) NUMBER OF ELEMENTS IN 'CHARAY' TO BE PROCESSED.
*              SPECIFY 1 TO PROCESS A SINGLE (SUB)STRING.
*     NOUT   : (OUTPUT) NR. OF CHARACTERS STORED IN 'COMPRS'.
*              0 IF INPUT STRING CONTAINED NO NON-BLANK CHARACTERS.
*
* THE SUBROUTINE IS ALSO USEFUL TO DETERMINE THE NUMBER OF NON-BLANK
* CHARACTERS IN A STRING OR CHARACTER ARRAY.
*
* CALLING PGM MUST DECLARE AS FOLLOWS :
*     EXTERNAL COMPRS
*     CHARACTER*(NR. OF CHAR. NEEDED) COMPRS
*
* CALLS FUNCTION CHRINT 
*
* IF AN ERROR CONDITION ARISES, COMPRS IS RETURNED AS A BLANK STRING.
*
      external      CHRINT
      character*12  CHRINT
      character*64  error1
      character*(*) charay(*)
      character*4   insert
      character     a,blank
      integer       nIn, nOut, MaxChar, idum
      integer      lnp, i, j, InLen, Loc, one
      parameter    (blank = ' ', lnp = 6, one = 1)
*
      MaxChar = LEN (COMPRS)
*
      COMPRS = blank
*
      loc    = one
      nout   = 0
      do i = one, nIn
          inlen = LEN_TRIM (charay(i))
          do j = one, inlen
              a = charay(i)(j:j)
              if (a.ne.blank) then
                  COMPRS(loc:loc) = a
                  loc = loc + one
              end if
          end do
      end do
      nOut = loc - one
*
* Write error message & blank output string if number of characters
* exceeds declared function length.
*
      if (nout.gt.MaxChar) then
          insert = CHRINT (MaxChar,-1,idum)
          error1 = '('' function COMPRS : input array has over '
     &             //insert//' char.'')'
          write (lnp,error1)
          COMPRS = blank
      end if
      return
      end
C------------------------------------------------------------------------
      CHARACTER*(*) FUNCTION CHRINT(NUMBER,IDIGIT,NCHAR)
*     input/output :           o      i      i      o
*     type         :          char   int    int    int
*
* LUC E. CUYPERS, MAY 26, 1984
*
* CONVERTS AN INTEGER NUMBER TO CHARACTER.
*
*      NUMBER : (INPUT) POS. OR NEG. INTEGER NUMBER
*      IDIGIT : (OPTIONAL INPUT) NUMBER OF DIGITS IN 'NUMBER', NEGATIVE
*               OR ZERO IF UNKNOWN. IF IDIGIT IS LESS THAN THE ACTUAL
*               NUMBER OF DIGITS  ERROR CONDITION RESULTS & CHRINT IS RE-
*               TURNED BLANK. IF IDIGIT IS LARGER THAN THE ACTUAL NUMBER
*               OF DIGITS, THE CHARACTER STRING WILL HAVE ZEROS IN THE
*               LEFTMOST POSITIONS AS IN EXAMPLE 1.
*               IF SPECIFIED, IDIGIT SHOULD BE LESS THAN OR EQUAL TO THE
*               DECLARED FUNCTION LENGTH (MINUS 1 TO LEAVE ROOM FOR A
*               MINUS SIGN, IF THERE IS ANY).
*      NCHAR  : (OUTPUT) NR. OF MEANINGFUL CHARACTERS IN 'CHRINT'.
*               IF IDIGIT IS SPECIFIED, NCHAR = IDIGIT (IDIGIT+1 FOR
*               NEGATIVE NUMBERS).
*               0 FOR ERROR CONDITION
*
* EXAMPLE : WITH NUMBER = 3      AND  IDIGIT = 3,      CHRINT = 003
*           WITH NUMBER = 3      AND  IDIGIT = 1,      CHRINT = 3
*           WITH NUMBER = 3      AND  IDIGIT =-1 OR 0, CHRINT = 3
*
* CALLING PGM ***MUST*** DECLARE FUNCTION AS FOLLOWS :
*     EXTERNAL CHRINT
*     CHARACTER*(NUMBER OF CHARACTERS NEEDED) CHRINT
*
*
* ERROR CONDITION RESULTS IN  CHRINT SET TO BLANKS & NCHAR = 0
*
      CHARACTER    BLANK
      CHARACTER*20 BUFFER
      INTEGER      ASCII0
      PARAMETER   (ASCII0 = 48, BLANK = ' ', LNP = 6)
      LOGICAL      MINUS
*
   50 FORMAT(' FUNCTION CHRINT : INPUT NUMBER ABOVE 10**10!')
   75 FORMAT(' FUNCTION CHRINT : INPUT PARAMETER ''IDIGIT'' TOO SMALL')
  100 FORMAT(' FUNCTION CHRINT : FUNCTION TOO SHORT TO STORE',
     &       ' OUTPUT STRING')
*
      MINUS  = .FALSE.
      CHRINT = BLANK
      NCHAR  = 0
      MAXLEN = LEN(CHRINT)
*
      NUM = IABS(NUMBER)
      IF (NUMBER.LT.0) MINUS = .TRUE.
*
      IF (IDIGIT.GT.MAXLEN) THEN
          WRITE(LNP,100)
          RETURN
      END IF
*
      IF (IDIGIT.LE.0) THEN
          ITEST = 1
          DO 10 I = 1,10
            ITEST = ITEST*10
            IF (NUM.LT.ITEST) THEN
                NCHAR = I
                GO TO 15
            END IF
   10     CONTINUE
        ELSE
          NCHAR = IDIGIT
      END IF
      IF (NCHAR.EQ.0) THEN
          WRITE(LNP,50)
          RETURN
      END IF
*
   15 CONTINUE
*
      IF (NCHAR.GT.MAXLEN) THEN
          WRITE (LNP,100)
          NCHAR  = 0
          RETURN
      END IF
*
* Generate string starting at rightmost position
*
      NEWNUM = NUM
      DO 20 I = NCHAR,1,-1
          NUMD10 = NEWNUM/10
          ITEMP  = NEWNUM - NUMD10*10
          CHRINT(I:I) = CHAR((ITEMP + ASCII0))
          NEWNUM = NUMD10
   20 CONTINUE
*
* Check if IDIGIT was specified large enough
*
      IF (NEWNUM.NE.0) THEN
          WRITE (LNP,75)
          CHRINT = BLANK
          NCHAR  = 0
          RETURN
      END IF
*
      IF (MINUS) THEN
          NCHAR  = NCHAR + 1
          IF (NCHAR.LE.MAXLEN) THEN
              BUFFER(1:NCHAR) = CHRINT
c(1:NCHAR)
              CHRINT          = '-'//BUFFER(1:NCHAR)
            ELSE
              WRITE (LNP,100)
              CHRINT = BLANK
              NCHAR  = 0
              RETURN
          END IF
      END IF
*
      RETURN
      END
