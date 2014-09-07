subroutine big(matrix,a,b,maior)
!
!A subrutina big busca o maior elemento en valor absoluto dunha matriz bidimensional de dimensións
!a x b. Os parámetros de entrada son a matriz e as dimensións, mentres que a saída é o maior ele-
!mento en valor absoluto.
!
use prec
implicit none
!
!Declaro as variables de entrada:
!-a,b: 	  dimensións da matriz a x b
!-matrix: matriz á que lle buscamos o maior elemento
!
integer (kind=int),intent(in) :: a,b
real (kind=dp),intent(in) :: matrix(a,b)
!
!Agora defino a variable de saída:
!-maior: é o elemento de maior valor absoluto
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
!Poño a cero o elemento maior, de xeito que calquera elemento non nulo lle cambiará ese valor. Os 
!bucles recorren tódalas columnas para cada fila.
!
maior=0
do i=1,a
  do j=1,b
!
!	  Agisno a compoñente (i,j) da matriz á variable x
!
      x=matrix(i,j)
!
!	  Se o elemento en cuestión é maior que o que teño gardado como maior, substitúeo, se non non
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