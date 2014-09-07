program eye3(a)

real (kind=2) :: n
real (kind=2)intent(inout) :: a

a=0

n=(size(a))^(1.0/3.0)
do i=1,n
  a(i,i)=1
enddo

return
end subroutine eye3
