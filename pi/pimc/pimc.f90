program pimc
!
!use prec
implicit none
!

real(kind=2) :: ran3
external ran3
integer(kind=3),parameter :: N=10**10
integer(kind=3) :: i,d,s1,s2
real(kind=2) :: x,y,r,p,c
real(kind=2),parameter :: pi=3.141592653589798d00
!
d=0
!
s1=2
s2=3
do i=1,N
  x=Ran3(s1)
  y=Ran3(s2)
  r=x*x+y*y
  if (r<=1) then
    d=d+1
  end if
end do
!
p=4.d00*dble(d)/dble(N)
c=p/pi
!
write(*,*) 'pi=',p
write(*,*) 'pi/pi=',c
!
stop
end program pimc