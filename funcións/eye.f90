program eye2(a)

real (kind=2) :: n
real (kind=2)intent(inout) :: a

a=0

n=sqrt(size(a))
do i=1,n
  a(i,i)=1
enddo

return
end subroutine eye2