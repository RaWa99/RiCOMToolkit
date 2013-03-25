!*********************************************************************************

      logical function OpenFile(nunit, prompt, fname, template)

!+
!+      logical function OpenFile(nunit, prompt, fname, template)
!+ Call Sequence:
!+              open(unit=lun, file=name, status='unknown',...)
!+ Purpose:  To select the name of an existing or new file, so the caller can
!+           then open the file for writing.
!+ Givens :  character*(*)      prompt - contains a message displayed with
!+                              the request for a file name.
!+           character*(*)      template - contains file name template
!+                              including wild cards etc. 
!+                              May require system dependent implementation.
!+                              May be ignored completely.
!+ Returns : logical            PigGetOpenFileName 
!+                                      is .TRUE. if the user selected a 
!+                                      valid filename, which is returned
!+                                      in name. The file must exist.
!+                                      is .FALSE. if the user cancelled
!+                                      the operation.
!+           character*(*)      name - contains a valid filename of an
!+                              existing file.
!+ Effects:  Gets the name of an existing file from the user.
!+
      integer nunit
      character(*) fname, prompt, template
      logical(4) resOK

      write(*,*) prompt
      read(*,'(a)') FNAME

      OPEN(UNIT=nunit,file=FNAME,status='unknown')

      OpenFile = .true. !resOK

      return
      end

!*********************************************************************************

      logical function OpenBinFile(nunit, prompt, fname, template)

!+
!+      logical function OpenFile(nunit, prompt, fname, template)
!+ Call Sequence:
!+              open(unit=lun, file=name, status='unknown',...)
!+ Purpose:  To select the name of an existing or new file, so the caller can
!+           then open the file for writing.
!+ Givens :  character*(*)      prompt - contains a message displayed with
!+                              the request for a file name.
!+           character*(*)      template - contains file name template
!+                              including wild cards etc. 
!+                              May require system dependent implementation.
!+                              May be ignored completely.
!+ Returns : logical            PigGetOpenFileName 
!+                                      is .TRUE. if the user selected a 
!+                                      valid filename, which is returned
!+                                      in name. The file must exist.
!+                                      is .FALSE. if the user cancelled
!+                                      the operation.
!+           character*(*)      name - contains a valid filename of an
!+                              existing file.
!+ Effects:  Gets the name of an existing file from the user.
!+
      integer nunit
      character(*) fname, prompt, template
      logical(4) resOK

      write(*,*) prompt
      read(*,'(a)') FNAME

      OPEN(UNIT=nunit,file=FNAME,form='unformatted',status='unknown')

      OpenBinFile = .true. !resOK

      return
      end

!*********************************************************************************

      logical function OpenFUFile(nunit, prompt, fname, template)

!   USE DFLIB
!+
!+      logical function OpenFile(nunit, prompt, fname, template)
!+ Call Sequence:
!+              open(unit=lun, file=name, status='unknown',...)
!+ Purpose:  To select the name of an existing or new file, so the caller can
!+           then open the file for writing.
!+ Givens :  character*(*)      prompt - contains a message displayed with
!+                              the request for a file name.
!+           character*(*)      template - contains file name template
!+                              including wild cards etc. 
!+                              May require system dependent implementation.
!+                              May be ignored completely.
!+ Returns : logical            PigGetOpenFileName 
!+                                      is .TRUE. if the user selected a 
!+                                      valid filename, which is returned
!+                                      in name. The file must exist.
!+                                      is .FALSE. if the user cancelled
!+                                      the operation.
!+           character*(*)      name - contains a valid filename of an
!+                              existing file.
!+ Effects:  Gets the name of an existing file from the user.
!+
      integer nunit
      character(*) fname, prompt, template
      logical(4) resOK
      character bin

      write(*,*) prompt
      read(*,'(a)') FNAME

      write(*,*) 'Is file formatted or unformatted (F/U)?'
      read(*,'(a)') bin
!   inquire(name=fName,FORM=bin)

      if(bin.eq.'U'.or.bin.eq.'u') then
        open(nunit, file=fName, form='unformatted',status='unknown')
        OpenFUFile = .true. !resOK
      elseif(bin.eq.'F'.or.bin.eq.'f') then
        OPEN(UNIT=nunit,file=FNAME,status='unknown')
        OpenFUFile = .true. !resOK
      else
        OpenFUFile = .false. !resOK
      endif

      return
      end

!*********************************************************************************
