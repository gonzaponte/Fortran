subroutine zeros(v)
real (kind=2),intent(inout) :: v
integer (kind=3) :: i
do i=1,size(v)
  	v(i)=0
enddo
return
endsubroutine zeros