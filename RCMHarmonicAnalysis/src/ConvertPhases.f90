! *********************************************************************

subroutine ConvertPhases(AnalysisTime,ncon,fn,nu,name,freq)

! *********************************************************************
! ***  subroutiine to adjust nodal correction factor fn and phase nu=(V0+u)+w*t
! ***  for the actual forecast period. All times are referenced to GMT.
! ***  V0 is referenced to the time 000hr 1/1/1900
! ***  t is elapsed time in hours from 000hr 1/1/1900 to the analysis time
! *********************************************************************

implicit none

! *** passed variables
integer, intent(in)    :: ncon 
character*4, intent(in) :: name(ncon)
real*8, intent(inout)       :: freq(ncon)
real, intent(inout)      :: fn(ncon),nu(ncon)
character(18), intent(in) ::  AnalysisTime


! *** Local variables
integer ::  j, k
integer, parameter :: ntidesmx = 13
integer ::  ifu(ntidesmx)=(/6,10,6,7,3,4,10,4,6,6,6,8,10/)
real    ::   fm(ntidesmx), mu(ntidesmx)
real    ::   v0(ntidesmx)=(/6.328,0.,63.686,200.382,10.190,356.138,349.810,53.496,121.044,12.656,315.289,128.970,1.031/)
real*8  ::    w(ntidesmx)=(/28.9841042373,30.0000000000,28.4397295415,30.0821372786,15.0410686393,13.9430366980, &
                    14.9589313607,13.3986609022,27.89533548458,27.9682084746,28.5125831704,29.5284789331,29.9589333224/)
real*8  ::  thours, dph, diffdph, v1(ntidesmx)
character*4 ::  tide(ntidesmx)=(/'M2  ','S2  ','N2  ','K2  ',&
                                 'K1  ','O1  ','P1  ','Q1  ',&
                                 '2N2 ','MU2 ','NU2 ','L2  ','T2  '/)

! *** calculate hours from 1/1/1900 00:00
      thours = 0.
      call GetTheBloodyElapsedHours(AnalysisTime, thours)

! *** calculate f and u, then add V0+time offset
      v1 = v0 + w*(thours)
      v1 = mod(v1+360.D0,360.D0)
      call fmu(thours,ntidesmx,fm,mu,ifu)

! *** load the constituents according to their frequency

      do j=1,ncon
        do k=1,ntidesmx
          if(name(j).eq.tide(k)) then !found one
            fn(j)=fm(k)
            nu(j) = v1(k) + mu(k)
            freq(j) = w(k)/360.D0
          endif
        enddo
      enddo

      nu = mod(nu+360.,360.)

!      do j=1,ntidesmx
!        write(lst,*) fm(j),mu(j),v0(j),v1(j),tide(j)
!      enddo

      do j=1,ncon
        write(*,'(4x,a4,2x,f10.8,f10.7,f10.4)') name(j),freq(j),fn(j),nu(j)
      enddo

end subroutine


! *********************************************************************

subroutine fmu(thours,ncon,fn,nu,ifu)

!	f and mu corrections for thours from 1/1/1900

implicit none

! *** passed variables
integer, intent(in)    :: ncon,ifu(ncon)
real, intent(out)      :: fn(ncon),nu(ncon)
real*8, intent(in)     :: thours

! *** Local variables
real :: f(24), u(24)
real an,p,xx,yy,pi
real*8 dan,dp

pi = 2.*asin(1.)
dan=259.156d0-0.00220641d0*thours
dan=dmod(dan,360.d0)
an=dan*pi/180.

dp=334.384d0+0.0046418249d0*thours
dp=dmod(dp,360.d0)
p=dp*pi/180.

      f(1)=1.000-0.13*cos(an)+0.0013*cos(2.*an)
      u(1)=0.
      f(2)=1.0429+0.4135*cos(an)-0.004*cos(2.*an)
      u(2)=-23.74*sin(an)+2.68*sin(2.*an)-0.38*sin(3.*an)
      f(3)=1.006+0.115*cos(an)-0.0088*cos(2.*an)+0.0006*cos(3.*an)
      u(3)=-8.86*sin(an)+0.68*sin(2.*an)-0.07*sin(3.*an)
      f(4)=1.0089+0.1871*cos(an)-0.0147*cos(2.*an)+0.0014*cos(3.*an)
      u(4)=10.80*sin(an)-1.35*sin(2.*an)+0.19*sin(3.*an)
      f(5)=1.0129+0.1676*cos(an)-0.0170*cos(2.*an)+0.0016*cos(3.*an)
      u(5)=-12.94*sin(an)+1.34*sin(2.*an)-0.19*sin(3.*an)
      f(6)=1.0004-0.0373*cos(an)+0.0002*cos(2.*an)
      u(6)=-2.14*sin(an)
      f(7)=1.0241+0.2863*cos(an)+0.0083*cos(2.*an)-0.0015*cos(3.*an)
      u(7)=-17.74*sin(an)+0.68*sin(2.*an)-0.04*sin(3.*an)
      xx=1.0000-0.2505*cos(2.*p)-0.1102*cos(2.*p-an)-0.0156*cos(2.*p-2.*an)-0.0370*cos(an)
      yy=-0.2505*sin(2.*p)-0.1102*sin(2.*p-an)-0.0156*sin(2.*p-2.*an)- 0.0370*sin(an)
      f(8)=sqrt(xx*xx+yy*yy)
      u(8)=atan(yy/xx)*180./pi
      xx=2.*cos(p)+0.4*cos(p-an)
      yy=sin(p)+0.2*sin(p-an)
      f(9)=sqrt(xx*xx+yy*yy)
      u(9)=atan(yy/xx)*180./pi
      f(10)=1.
      u(10)=0.
      f(11)=1.1027+0.6504*cos(an)+0.0317*cos(2.*an)-0.0014*cos(3.*an)
      u(11)=-36.68*sin(an)+4.02*sin(2.*an)-0.57*sin(3.*an)
      f(12)=f(4)*f(4)
      u(12)=2.*u(4)
      f(13)=f(6)*f(6)
      u(13)=2.*u(6)
      f(14)=f(6)*f(7)
      u(14)=u(6)+u(7)
      f(15)=f(3)*f(3)
      u(15)=2.*u(3)
      f(16)=f(4)
      u(16)=-u(4)
      f(17)=f(6)
      u(17)=-u(6)
      f(18)=f(6)*f(4)
      u(18)=u(6)+u(4)
      f(19)=sqrt(f(6)**3.)
      u(19)=1.5*u(6)
      f(20)=f(6)*f(3)
      u(20)=u(6)+u(3)
      f(21)=f(6)**3.
      u(21)=3.*u(6)
      f(22)=f(14)*f(7)
      u(22)=u(14)+u(7)
      f(23)=f(13)
      u(23)=0.
      f(24)=f(13)*f(13)
      u(24)=2.*u(13)

      fn = f(ifu)
      nu = u(ifu)

      end subroutine

! *********************************************************************
      subroutine GetTheBloodyElapsedHours(AnalysisTime, thours)
      
      implicit none

! *** passed variables
      character(*) AnalysisTime
      real*8 thours

! *** local variables
      integer iyy,imm,idd,icc,kd,kd0,ihr,imn

! *** get reference day 000h 1/1/1900
      iyy = 0
      imm = 1
      idd = 1
      icc = 19
      call GDAY(idd,imm,iyy,icc,kd0)
      kd0 = (kd0-1)*24

! *** parse analysistime
      read(AnalysisTime(1:2),'(I2)')   icc
      read(AnalysisTime(3:4),'(I2)')   iyy
      read(AnalysisTime(6:7),'(I2)')   imm
      read(AnalysisTime(9:10),'(I2)')  idd
      read(AnalysisTime(12:13),'(I2)') ihr
      read(AnalysisTime(15:16),'(I2)') imn

      call GDAY(idd,imm,iyy,icc,kd)
      kd = (kd-1)*24 + ihr -kd0

      thours = float(kd) + float(imn)/60.

      return
      end

! *********************************************************************

      SUBROUTINE GDAY(IDD,IMM,IYY,ICC,KD)

!     GIVEN DAY,MONTH,YEAR AND CENTURY(EACH 2 DIGITS), GDAY RETURNS
!     THE DAY#, KD BASED ON THE GREGORIAN CALENDAR.
!     THE GREGORIAN CALENDAR, CURRENTLY 'UNIVERSALLY' IN USE WAS
!     INITIATED IN EUROPE IN THE SIXTEENTH CENTURY. NOTE THAT GDAY
!     IS VALID ONLY FOR GREGORIAN CALENDAR DATES.

!     KD=1 CORRESPONDS TO JANUARY 1, 0000
	
! 	Note that the Gregorian reform of the Julian calendar 
!	omitted 10 days in 1582 in order to restore the date
!	of the vernal equinox to March 21 (the day after
!	Oct 4, 1582 became Oct 15, 1582), and revised the leap 
!	year rule so that centurial years not divisible by 400
!	were not leap years.

!C   THIS ROUTINE WAS WRITTEN BY EUGENE NEUFELD, AT IOS, IN JUNE 1990.

      implicit none

! *** passed variables
      integer idd,imm,iyy,icc,kd
      integer, save :: lst = 6

! *** local variables
      INTEGER NDP(13)
      INTEGER NDM(12)
      DATA NDP/0,31,59,90,120,151,181,212,243,273,304,334,365/
      DATA NDM/31,28,31,30,31,30,31,31,30,31,30,31/

!  TEST FOR INVALID INPUT:
      IF(ICC.LT.0)THEN
        WRITE(lst,*) ' INPUT ERROR. ICC = ',ICC
        STOP 72
      ENDIF
      IF(IYY.LT.0.OR.IYY.GT.99)THEN
        WRITE(lst,*) ' INPUT ERROR. IYY = ',IYY
        STOP 72
      ENDIF
      IF(IMM.LE.0.OR.IMM.GT.12)THEN
        WRITE(lst,*) ' INPUT ERROR. IMM = ',IMM
        STOP 72
      ENDIF
      IF(IDD.LE.0)THEN
        WRITE(lst,*) ' INPUT ERROR. IDD = ',IDD
        STOP 72
      ENDIF
      IF(IMM.NE.2.AND.IDD.GT.NDM(IMM))THEN
        WRITE(lst,*) ' INPUT ERROR. IDD = ',IDD
        STOP 72
      ENDIF
      IF(IMM.EQ.2.AND.IDD.GT.29)THEN
        WRITE(lst,*) ' INPUT ERROR. IDD = ',IDD
        STOP 72
      ENDIF
      IF(IMM.EQ.2.AND.IDD.GT.28.AND.((IYY/4)*4-IYY.NE.0.OR.(IYY.EQ.0.AND.(ICC/4)*4-ICC.NE.0)))THEN
        WRITE(lst,*) ,IDD
        STOP 72
      ENDIF
!5000  FORMAT(' INPUT ERROR. ICC = ',I7)
!5010  FORMAT(' INPUT ERROR. IYY = ',I7)
!5020  FORMAT(' INPUT ERROR. IMM = ',I7)
!5030  FORMAT(' INPUT ERROR. IDD = ',I7)

! *** CALCULATE DAY# OF LAST DAY OF LAST CENTURY:
      KD = ICC*36524 + (ICC+3)/4

! *** CALCULATE DAY# OF LAST DAY OF LAST YEAR:
      KD = KD + IYY*365 + (IYY+3)/4

! *** ADJUST FOR CENTURY RULE:
! *** (VIZ. NO LEAP-YEARS ON CENTURYS EXCEPT WHEN THE 2-DIGIT
! *** CENTURY IS DIVISIBLE BY 4.)
      IF(IYY.GT.0.AND.(ICC-(ICC/4)*4).NE.0) KD=KD-1
! *** KD NOW TRULY REPRESENTS THE DAY# OF THE LAST DAY OF LAST YEAR.

! *** CALCULATE DAY# OF LAST DAY OF LAST MONTH:
      KD = KD + NDP(IMM)

! *** ADJUST FOR LEAP YEARS:
      IF(IMM.GT.2.AND.((IYY/4)*4-IYY).EQ.0.AND.((IYY.NE.0).OR.(((ICC/4)*4-ICC).EQ.0)))   KD=KD+1
! *** KD NOW TRULY REPRESENTS THE DAY# OF THE LAST DAY OF THE LAST MONTH.

! *** CALCULATE THE CURRENT DAY#:
      KD = KD + IDD

      RETURN
      end

! *********************************************************************
