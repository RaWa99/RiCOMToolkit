
    subroutine ellipse(np,u,up,v,vp,amaj,amin,g,ainc)

    integer np
    real u(np),up(np),v(np),vp(np)
    real amaj(np),amin(np),g(np),ainc(np)

!	calculate the ellipse parameters
	pi=3.1415926536
	fac=pi/180.

	do j=1,np
!	  cx = vel(j,3)
!	  sx = vel(j,4)
!	  cy = vel(j,5)
!	  sy = vel(j,6)
	  cx = u(j)*cos(up(j)*fac)
	  sx = u(j)*sin(up(j)*fac)
	  cy = v(j)*cos(vp(j)*fac)
	  sy = v(j)*sin(vp(j)*fac)
	  cxpsy = 0.5*(cx+sy)
	  cymsx = 0.5*(cy-sx)
	  cxmsy = 0.5*(cx-sy)
	  cypsx = 0.5*(cy+sx)
	  apl = sqrt(cxpsy**2+cymsx**2)
	  amn = sqrt(cxmsy**2+cypsx**2)
	  if(apl.lt.1.e-5) then
	    eps=0.
	  else
	    epl=atan2(cymsx,cxpsy)/fac
	  end if
	  if(amn.lt.1.e-5) then
	    emin=0.
	  else
	    emin=atan2(cypsx,cxmsy)/fac
	  end if

	  amaj(j)=apl+amn
	  amin(j)=apl-amn

	  gpl = -epl
	  gmin = emin
	  if(gmin-gpl.lt.0.) gmin=gmin+360.
	  g(j)=0.5*(gpl+gmin)
	  ainc(j)=0.5*(gmin-gpl)
	  if( g(j).gt.360.) g(j)=g(j)-360.
	  if( g(j).lt.0.) g(j)=g(j)+360.

	enddo

    end
