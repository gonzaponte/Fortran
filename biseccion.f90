subroutine biseccion(f,min,max,p)
!
integer,parameter :: int = SELECTED_INT_KIND(9)       !Define un enteiro
integer,parameter ::  dp = SELECTED_REAL_KIND(15,307) !Define un real
!
character (len=100),intent(in) :: f
real(kind=dp),intent(in) :: min,max,p
!
integer(kind=int) :: i,size
real(kind=dp) :: pinv,delta
real(kind=dp),allocatable :: x(:),y(:)
!
pinv=1/p
delta=max-min
size=nint(delta*pinv)
!
allocate (x(size),y(size))
!
x=min-1
!
do i=1:size
  x(i)=x(i)+dble(i)*p
enddo
!
y=f
return
end program