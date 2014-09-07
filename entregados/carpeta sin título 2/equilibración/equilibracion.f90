program equilibracion
!
!Introd�cese o programa equilibraci�n que leva o estado inicial creado polo programa redefcc ao
!equilibrio. Ten como par�metros de entrada os de sa�da do programa redefcc:
!
!Defino o que � un enteiro e o que � un real en doble precisi�n mediante o m�dulo prec e afirmo que
!non hai variables impl�citas.
!
use prec
implicit none
!
!Declaro variables:
!
!-N:    	n�mero de part�culas do sistema
!-k: 		n�mero de celdas por eixo
!-L: 		lado da cubo de simulaci�n
!-pasos: 	n�mero de pasos da simulaci�n
!-gravar: 	define cada cando se gravan as enerx�as
!-i: 		variable para correr no bucle
!-ndp: 		n�mero de part�culas en dobre precisi�n
!-ldp: 		lado da cubo en dobre precisi�n
!-ldpi: 	inverso de ldp
!-a: 		lado de cada celda
!-rho: 		densidade de part�culas
!-L3: 		volume da caixa de simulaci�n
!-L33: 		constante definida como 1/3 da inversa do volume
!-Rc2: 		raio de corte do potencial ao cadrado
!-Rc3: 		raio de corte do potencial ao cubo
!-Rc6:  	raio de corte do potencial � sexta
!-cor:      constante para as correcci�ns � enerx�a potencial e derivadas
!-U:        enerx�a potencial
!-phiv:     primeira derivada da enerx�a potencial
!-phivv:    segunda derivada da enerx�a potencial
!-T:        enerx�a cin�tica
!-E:        enerx�a total
!-tempo:    variable para saber en que paso da simulaci�n me atopo
!-dt:       paso temporal da simulaci�n
!-Px,Py,Pz: momento total do sistema
!-t1,t2: 	variables para controlar o tempo de execuci�n do programa
!-rx,ry,rz: posici�ns das N part�culas
!-vx,vy,vz: velocidades das N part�culas
!-fx,fy,fz: aceleraci�ns das N part�culas
!
integer(kind=int) :: N,k,L,Rc,pasos,gravar,i
real(kind=dp) :: ndp,ldp,ldpi,a,rho,L3,L33,Rc2,Rc3,Rc6,cor,U,phiv,phivv,T,E,tempo,dt,Px,Py,Pz,t1,t2
real (kind=dp),parameter :: pi=3.14159265458979
real(kind=dp),dimension(500) :: rx,ry,rz,vx,vy,vz,fx,fy,fz
!
!Defino os formatos que vou usar
!
8000	format(1pe19.12,8(1x,e19.12))
8001	format(1pe19.12,4(1x,e19.12))
8002	format(i8,i8,i8,i8,i8,i8)
8003	format(1pe19.12)
8004	format(1pe19.12,1x,e19.12,1x,e19.12,1x,e19.12)
!
call cpu_time(t1)
!
!Collo os datos que solta o programa redefcc e as�gnoos �s variables do programa.
!
open(10,file='N�meros.dat')
	read (10,8002) N,k,L,Rc,pasos,gravar              !Enteiros
    read (10,8003) dt                                 !Paso temporal
    read (10,8000) ndp,ldp,ldpi,a,rho,L3,Rc2,Rc3,Rc6  !Reais 1
	read (10,8001) U,phiv,phivv,T,E                   !Reais 2
close(10)
!
open(20,file='Matrixbin.dat',form='unformatted')
	read(20) rx,ry,rz,vx,vy,vz,fx,fy,fz
close(20)
!
!Calculo unhas constantes que usarei m�is adiante
!
cor=16.d00*pi*rho*ndp*Rc3
L33=1/(3*L3)
!
!Agora ven o programa en si. O seguinte bucle avanza no tempo:
!
do i=1,pasos
!
!	Chamo a subrutina verletv que utiliza o algoritmo verlet velocity que me da as novas posici�ns,
!	velocidades e aceleraci�ns, a enerx�a cin�tica, potencial e as derivadas desta �ltima a partir
!	do n�mero de part�culas, o lado da caixa e o raio de corte.
!
	call verletv(N,ldp,ldpi,L33,Rc2,Rc6,cor,rx,ry,rz,vx,vy,vz,fx,fy,fz,dt,U,phiv,phivv)
!
!	A seguinte condici�n di que se i � un m�ltiplo enteiro de veces a variable gravar me escriba nun
!	ficheiro o tempo de simulaci�n e as enerx�as total, cin�tica e potencial.
!
	if (mod(i,gravar)==0) then
!
!		Calculo a enerx�a cin�tica e total e o paso de tempo no que me atopo.
!
		T=sum(vx*vx+vy*vy+vz*vz)*0.5d00
    	E=T+U
		tempo=dt*dble(i)
!
!		Gravo as enerx�as e o tempo de simulaci�n
!
    	open (30,file='Resultados.dat',access='append')
    		write (30,8004) tempo,E,T,U
    	close (30)
!
!		Gravo as velocidades en binario a partir da metade da simulaci�n.
!
		if (i>=(pasos/2)) then
        	open (40,file='uves.dat',form='unformatted',position='append')
           		write (40) vx,vy,vz
        	close (40)
    	endif
!
  	endif
!
!	Aqu� remata o bucle.
!
enddo
!
!Calculo o momento total e am�soo en pantalla, ao ser un sistema illado dar� case cero.
!
Px=sum(vx)
Py=sum(vy)
Pz=sum(vz)
!
write (*,*) 'P=',Px,Py,Pz
!
!Finalmente gravo as posici�ns, velocidades e aceleraci�ns no equilibrio noutro ficheiro.
!
open (50,file='RVA.dat')
	do i=1,N
    	write (50,8000) rx(i),ry(i),rz(i),vx(i),vy(i),vz(i),fx(i),fy(i),fz(i)
    enddo
close(50)
!
open (60,file='RVAbin.dat',form='unformatted')
    	write (60) rx,ry,rz,vx,vy,vz,fx,fy,fz
close(60)
!
call cpu_time(t2)
!
write (*,*) 'O programa tardou:',t2-t1
!
!Fin.
!
stop
end program