subroutine escalar(x,y,z)

real (kind=2), intent(in) :: x,y
real (kind=2), intent(our) :: z
real (kind=2) :: s

if (size(x)/=size(y)) then
  write (*,*) 'Erro: os vectores deben ter a mesma dimensión0'
  return
endif

z=0
do i=1,size(x)
  z=z+x(i)*y(i)
enddo

return
end subroutine escalar