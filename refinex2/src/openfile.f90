!*********************************************************************************

      logical function OpenFile(nunit, prompt, fname, template)

	USE DFLIB
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
	character *(*) fname, prompt, template
!     Following 7 lines for Windows browsing
	integer(4) len, len2
	logical(4) resOK
	character(256) pathbuf
	character(2) fdrive
	character(256) fdir
	character(8) fname8
	character(4) fext

	call setmessageqq( template, QWIN$MSG_FILEOPENDLG )
	open(nunit, file=' ')
	inquire(unit=nunit,name=fName)
	len = fullpathqq(fname, pathbuf)
	len2 = splitpathqq(pathbuf,fdrive,fdir,fname8,fext)
    fnlen = len_trim(fname8)
	fname = fname8(:fnlen)//fext
	resOK = changedirqq(fdrive//fdir(:len2))
	OpenFile = resOK

!     Following leaves window open (needs USE DFLIB above)
    len = SETEXITQQ(QWIN$EXITPERSIST)

	return

	end

!*********************************************************************************
      
      logical function OpenBinFile(nunit, prompt, fname, template)

	USE DFLIB
!+
!+      logical function OpenBinFile(nunit, prompt, fname, template)
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
	character *(*) fname, prompt, template
!     Following 7 lines for Windows browsing
	integer(4) len, len2
	logical(4) resOK
	character(256) pathbuf
	character(2) fdrive
	character(256) fdir
	character(8) fname8
	character(4) fext

	call setmessageqq( template, QWIN$MSG_FILEOPENDLG )
	open(nunit, file=' ', form='unformatted')
	inquire(unit=nunit,name=fName)
	len = fullpathqq(fname, pathbuf)
	len2 = splitpathqq(pathbuf,fdrive,fdir,fname8,fext)
    fnlen = len_trim(fname8)
	fname = fname8(:fnlen)//fext
	resOK = changedirqq(fdrive//fdir(:len2))
	OpenBinFile = resOK

!     Following leaves window open (needs USE DFLIB above)
    len = SETEXITQQ(QWIN$EXITPERSIST)

	return

	end
      
!*********************************************************************************

