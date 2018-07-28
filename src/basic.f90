!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Invert u from dq/dy
! Not currently used; but perhaps useful as a *reference*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program basic
  integer, parameter :: jmax = 801
  common /array/  dqdy1(jmax),dqdy2(jmax)
  common /brray/  u1(jmax),u2(jmax),ub(jmax),uc(jmax)
  common /crray/  al(jmax),gam(jmax),zeta(jmax),eta(jmax)
  common /drray/  pb(jmax),pc(jmax)
  eld = 800000.  !(m)
  beta = 1.7e-11 
  width = 80000000.   !(m)
  dy = width/float(jmax-1)
  bwidth = 1000000. !(m)
  pi = acos(-1.)
  const = 1.5e-4/(bwidth*sqrt(2.*pi))
  dd = 0.
  ! dd = 500000. !(m)
  do j =1,jmax
    y = dy*float(j-1)
    y1 = y-0.5*width-dd
    y2 = y-0.5*width+dd
    o1 = y1*y1/(bwidth*bwidth)
    o2 = y2*y2/(bwidth*bwidth)
    dqdy1(j) = beta + const*exp(-o1)
    dqdy2(j) = beta - const*exp(-o2)
    pb(j) = -(dqdy1(j)+dqdy2(j)-2.*beta)
    pc(j) = -(dqdy1(j)-dqdy2(j))
    write(*,*) j,y,dqdy1(j),dqdy2(j),dqdy1(j)-dqdy2(j)
  enddo
  al(jmax-1) = 0.
  gam(jmax-1) = 0.
  zeta(jmax-1) = 0.
  eta(jmax-1) = 0.
  do j = jmax-1,2,-1
    al(j-1) = -1./(al(j)-2.)
    gam(j-1) = (pb(j)*dy*dy-gam(j))/(al(j)-2.)
    zeta(j-1) = -1./(zeta(j)-2.-2.*dy*dy/(eld*eld))
    eta(j-1) = (pc(j)*dy*dy-eta(j))/  &
      (zeta(j)-2.-2.*dy*dy/(eld*eld))
  enddo
  ub(1) = 0.
  uc(1) = 0.
  do j = 1,jmax-1
    ub(j+1) = al(j)*ub(j)+gam(j)
    uc(j+1) = zeta(j)*uc(j)+eta(j)
  enddo
  do j = 1,jmax
    u1(j) = 0.5*(ub(j)+uc(j)) 
    u2(j) = 0.5*(ub(j)-uc(j)) 
    write(*,*) j,u1(j),u2(j)
  enddo
  stop
end
