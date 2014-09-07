subroutine big(matrix,a,b,maior)
!
!A subrutina big busca o maior elemento en valor absoluto dunha matriz bidimensional de dimensi�ns
!a x b. Os par�metros de entrada son a matriz e as dimensi�ns, mentres que a sa�da � o maior ele-
!mento en valor absoluto.
!
use prec
implicit none
!
!Declaro as variables de entrada:
!-a,b: 	  dimensi�ns da matriz a x b
!-matrix: matriz � que lle buscamos o maior elemento
!
integer (kind=int),intent(in) :: a,b
real (kind=dp),intent(in) :: matrix(a,b)
!
!Agora defino a variable de sa�da:
!-maior: � o elemento de maior valor absoluto
!
real (kind=dp),intent(out) :: maior
!
!Finalmente declaro variables auxiliares:
!-x: 	variable auxiliar para asignar momentaneamente a cada elemento da matriz
!-i,j:  variables para correr en bucle e acceder a cada elemento da matriz
!
real (kind=dp) :: x
integer (kind=int) :: i,j
!
!Po�o a cero o elemento maior, de xeito que calquera elemento non nulo lle cambiar� ese valor. Os 
!bucles recorren t�dalas columnas para cada fila.
!
maior=0
do i=1,a
  do j=1,b
!
!	  Agisno a compo�ente (i,j) da matriz � variable x
!
      x=matrix(i,j)
!
!	  Se o elemento en cuesti�n � maior que o que te�o gardado como maior, substit�eo, se non non
!	  se fai nada.
!
      if (abs(x)>maior) then
         maior=abs(x)
      endif
!
   enddo
enddo
!
!Fin.
!
return
end subroutine big