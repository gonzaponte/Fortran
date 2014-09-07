subroutine potencialLJ(N,ldp,ldpi,L33,Rc2,Rc6,rx,ry,rz,cor,U,phiv,phivv,fx,fy,fz)
!
!A subrutina potencialLJ calcula a enerxía potencial dun sistema, a súa primeira e segunda derivada
!e a aceleración resultante a partir do número de partículas, o lado da caixa o raio de corte e as
!posicións iniciais.
!
!Uso o módulo prec para definir un enteiro e un real en dobre precisión. Afirmo que non hai variables
!implícitas.
!
use prec
implicit none
!
!Declaro as variables de entrada:
!-N: 		número de partículas do sistema
!-ldp:  	lado do cubo de simulación en dobre precisión
!-ldpi: 	inversa de ldp
!-L33:  	un tercio do anterior
!-Rc2:  	raio de corte do potencial ao cadrado
!-Rc6:  	raio de corte elevado á sexta potencia
!-cor:  	factor para as correccións á enerxía
!-rx,ry,rz: compoñentes do vector posición de tódalas partículas.
!
integer (kind=int),intent(in) :: N
real (kind=dp),intent(in) :: ldp,ldpi,L33,Rc2,Rc6,cor
real (kind=dp),dimension(N),intent(in) :: rx,ry,rz
!
!Agora as variables de saída:
!-U: 		enerxía potencial do sistema
!-phiv: 	primeira derivada da enerxía potencial
!-phivv: 	segunda derivada da enerxía potencial
!-fx,fy,fz: compoñentes da aceleración de tódalas partículas
!
real (kind=dp),intent(out) :: U,phiv,phivv
real (kind=dp),dimension(N),intent(out) :: fx,fy,fz
!
!E finalmente variables a usar dentro da subrutina:
!-i,j: 		  variables auxiliares para correr nos bucles
!-dx,dy,dz:   compoñentes do vector distancia entre dúas partículas
!-D2: 		  módulo ao cadrado dese vector
!-D2i: 		  inverso do anterior
!-D6: 		  inversa da sexta potencia da distancia
!-D12:  	  inversa da duodécima potencia da distancia.
!-fmod: 	  parte común escalar dos vectores de aceleración.
!-DeltaU: 	  corrección á enerxía potencial
!-Deltaphiv:  corrección á suma da primeira derivada da enerxía potencial
!-Deltaphivv: corrección á suma da segunda derivada da enerxía potencial
!
integer (kind=int) :: i,j
real (kind=dp) :: dx,dy,dz,D2,D2i,D6,D12,fmod,DeltaU,Deltaphiv,Deltaphivv
!
!A continuación,calculo a enerxía potencial, as súas derivadas e a aceleración. Poño a cero tódolos
!acumuladores necesarios para calcular esas cantidades.
!
U=0.d00
phiv=0.d00
phivv=0.d00
fx=0.d00
fy=0.d00
fz=0.d00
!
!Corro o bucle en i dende un ata N-1 e en j de i+1 ata N para non calcular
!dúas veces o mesmo.
!
do i=1,N-1
  do j=i+1,N
!
!	Resto os vectores de cada par de partículas.
!    
    dx=rx(i)-rx(j)
    dy=ry(i)-ry(j)
    dz=rz(i)-rz(j)
!
!	As seguintes liñas redefinen a distancia cunha partícula como a da súa imaxe máis próxima.
!    
    dx=dx-ldp*dnint(dx*ldpi)
	dy=dy-ldp*dnint(dy*ldpi)
	dz=dz-ldp*dnint(dz*ldpi)
!
!	Calculo o módulo ao cadrado e as inversas do módulo á sexta e a duodécima potencia.
!    
    D2=dx*dx+dy*dy+dz*dz
!
!	A seguinte condición impón que se a distancia entre dúas partículas é menor  que o raio de 
!	corte, a teñamos en conta para a enerxía potencial, derivadas e aceleracións. Se cumpre a 
!	condición, sumará a enerxía potencial correspondente ao potencial de Lennard - Jones e fará
!	o mesmo coas derivadas e fmod. Calcúlanse as aceleracións do seguinte xeito:
!	Para a partícula i (a que estamos considerando) súmaselle a aceleración producida pola par-
!	tícula j, que é fmod*distancia. Como r_ij = - r_ji e fmod é o mesmo, para a partícula j rés-
!	taselle fmod*distancia porque esta é a aceleración producida pola partícula i. Ao facer o 
!	cálculo deste xeito e acumular o resultado non é necesario facer que os dous bucles vaian 
!	dende 1 ata N.
!    
    if (D2<=Rc2) then
		D2i=1/D2                		!Inversa do módulo ao cadrado
      	D6=D2i*D2i*D2i          		!Inversa do módulo á sexta potencia
      	D12=D6*D6               		!Inversa do módulo á duodécima potencia
      	U=U+D12-D6              		!Enerxía potencial
      	phiv=phiv-2.d00*D12+D6      	!Primeira derivada
     	phivv=phivv+26.d00*D12-7.d00*D6 !Segunda derivada
      	fmod=24.d00*(2.d00*D12-D6)*D2i  !Parte común da forza
      	fx(i)=fx(i)+fmod*dx     		!Compoñente x da forza sobre i debida a j
      	fx(j)=fx(j)-fmod*dx     		!Compoñente x da forza sobre j debida a i
      	fy(i)=fy(i)+fmod*dy     		!Compoñente y da forza sobre i debida a j
      	fy(j)=fy(j)-fmod*dy     		!Compoñente y da forza sobre j debida a i
      	fz(i)=fz(i)+fmod*dz     		!Compoñente z da forza sobre i debida a j
      	fz(j)=fz(j)-fmod*dz     		!Compoñente z da forza sobre j debida a i
    endif
!    
  enddo
enddo
!
!Calculo agora as correccións á enerxía potencial e as súas derivadas e introduzo as correccións
!nos seus valores. Ter en conta que Rc3 e Rc6 son as inversas da 3ª e 6ª potencia.
!
DeltaU=cor*(Rc6/3.d00-1.d00)/6.d00
U=4.d00*U+DeltaU
!
Deltaphiv=cor*(1.d00-2.d00*Rc6/3.d00)
phiv=(24.d00*phiv+Deltaphiv)*L33
!
Deltaphivv=cor*(26.d00*Rc6/3.d00-7.d00)
phivv=((24.d00*phivv+Deltaphivv)*L33-2.d00*phiv)*L33
!
!Fin.
!
return
end subroutine potencialLJ