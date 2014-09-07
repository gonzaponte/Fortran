program dinmol
!
!Introd�cese o programa dinmol que usa o estado en equilibrio creado polo programa equilibraci�n
!para calcular propiedades do sistema como presi�n, capacidade calor�fica, etc. Os par�metros de
!entrada lense dos datos de sa�da do programa de equilibraci�n.
!
!Defino o que � un enteiro e o que � un real en doble precisi�n mediante o m�dulo prec e afirmo que
!non hai variables impl�citas.
!
use prec
implicit none
!
!Declaro variables:
!
!-N:    	  n�mero de part�culas do sistema
!-k: 		  n�mero de celdas por eixo
!-L: 		  lado da cubo de simulaci�n
!-pasos: 	  n�mero de pasos da simulaci�n
!-gravar: 	  define cada cando se gravan as enerx�as
!-i: 		  variable para correr no bucle
!-ndp: 		  n�mero de part�culas en dobre precisi�n
!-ldp: 		  lado da cubo en dobre precisi�n
!-ldpi: 	  inverso de ldp
!-a: 		  lado de cada celda
!-rho: 		  densidade de part�culas
!-L3: 		  volume da caixa de simulaci�n
!-L3inv:	  inversa do volume
!-L33: 		  constante definida como 1/3 da inversa do volume
!-Rc2: 		  raio de corte do potencial ao cadrado
!-Rc3: 		  inversa do raio de corte do potencial ao cubo
!-Rc6:  	  inversa do raio de corte do potencial � sexta
!-cor:        constante para as correcci�ns � enerx�a potencial e derivadas
!-U:          enerx�a potencial
!-phiv:       primeira derivada da enerx�a potencial
!-phivv:      segunda derivada da enerx�a potencial
!-T:          enerx�a cin�tica
!-E:          enerx�a total
!-tempo:      variable para saber en que paso da simulaci�n me atopo
!-dt:         paso temporal da simulaci�n
!-Px,Py,Pz:   momento total do sistema
!-Tinv:		  inversa da enerx�a cin�tica
!-Tinvphiv:	  cociente entre a primeira derivada da enerx�a potencial e a enerx�a cin�tica
!-Tinvphiv2:  cociente entre o cadrado da primeira derivada da enerx�a potencial e a enerx�a cin�tica
!-###media:	  media da variable ###
!-Temp:		  temperatura
!-Cv:		  capacidade calor�fica a volume constante
!-Presion:	  presi�n
!-alphaE:	  coeficiente de dilataci�n a enerx�a constante
!-gamma:	  gamma de Gruneisen
!-kS:		  coeficiente de compresibilidade adiab�tico
!-t1,t2: 	  variables para controlar o tempo de execuci�n do programa
!-pasosi:	  inversa do n�mero de pasos total
!-pi:		  n�mero pi
!-f:		  graos de liberdade
!-rx,ry,rz:   posici�ns das N part�culas
!-vx,vy,vz:   velocidades das N part�culas
!-fx,fy,fz:   aceleraci�ns das N part�culas
!
integer(kind=int) :: N,k,L,Rc,pasos,gravar,i
real(kind=dp) :: ndp,ldp,ldpi,a,rho,L3,L3inv,L33,Rc2,Rc3,Rc6,cor,U,phiv,phivv,T,E,tempo,dt,Px,Py,Pz
real(kind=dp) :: Tinv,Tinvphiv,Tinvphiv2
real(kind=dp) :: Tmedia,Tinvmedia,Umedia,phivmedia,phivvmedia,Tinvphivmedia,Tinvphiv2media
real(kind=dp) :: Temp,Cv,Presion,alphaE,gamma,kS,t1,t2,pasosi
real (kind=dp),parameter :: pi=3.14159265458979d00,f=3.d00*499.d00
real(kind=dp),dimension(500) :: rx,ry,rz,vx,vy,vz,fx,fy,fz
!
!Defino os formatos que vou usar
!
8000	format(1pe19.12,8(1x,e19.12))
8002	format(i8,i8,i8,i8,i8,i8)
8003	format(1pe19.12)
8004	format(1pe19.12,9(1x,e19.12))
!
call cpu_time(t1)
!
!Collo os datos que solta o programa redefcc e as�gnoos �s variables do programa.
!
open(10,file='N�meros.dat')
	read (10,8002) N,k,L,Rc,pasos,gravar              !Enteiros
    read (10,8003) dt                                 !Paso temporal
    read (10,8000) ndp,ldp,ldpi,a,rho,L3,Rc2,Rc3,Rc6  !Reais
close(10)
!
!Collo o �ltimo arguivo de posici�ns velocidades e aceleraci�ns para o programa.
!
open (20,file='RVAbin.dat',form='unformatted')
    	read (20) rx,ry,rz,vx,vy,vz,fx,fy,fz
close(20)
!
!Calculo unhas constantes que me far�n falta.
!
cor=16.d00*pi*rho*ndp*Rc3
L3inv=1.d00/L3
L33=L3inv/3.d00
pasosi=1.d00/dble(pasos)
!
!Agora ven o programa en si. Anulo acumuladores de medias.
!
Tmedia=0
Umedia=0
phivmedia=0
phivvmedia=0
Tinvmedia=0
Tinvphiv2media=0
Tinvphivmedia=0
!
!Avanzo no tempo:
!
do i=1,pasos
!
!	Chamo a subrutina verletv que utiliza o algoritmo verlet velocity que me da as novas posici�ns,
!	velocidades e aceleraci�ns, a enerx�a cin�tica, potencial e as derivadas desta �ltima a partir
!	do n�mero de part�culas, o lado da caixa e o raio de corte.
!
	call verletv(N,ldp,ldpi,L33,Rc2,Rc6,cor,rx,ry,rz,vx,vy,vz,fx,fy,fz,dt,U,phiv,phivv)
!
!	Calculo a enerx�a cin�tica, enerx�a total, inversa da cin�tica, primeira derivada da potencial
!	entre a enerx�a cin�tica, o cadrado da primeira derivada da potencial entre a cin�tica.
!	Estas son as variables que me fan falta para o c�lculo das propiedades do sistema.
!
	T=sum(vx*vx+vy*vy+vz*vz)*0.5d00
    E=T+U
    Tinv=1.d00/T
    Tinvphiv=Tinv*phiv
    Tinvphiv2=Tinvphiv*phiv
!
!	Sumo os valores �s s�as respectivas medias.
!
    Tmedia=Tmedia+T
	Umedia=Umedia+U
	phivmedia=phivmedia+phiv
	phivvmedia=phivvmedia+phivv
	Tinvmedia=Tinvmedia+Tinv
   	Tinvphivmedia=Tinvphivmedia+Tinvphiv
	Tinvphiv2media=Tinvphiv2media+Tinvphiv2
!  
!	A seguinte condici�n di que se o �ndice i � un m�ltiplo enteiro de veces a variable gravar
!	faga unha serie de instrucci�ns.
!
	if (mod(i,gravar)==0) then
!		
!		Calculo o tempo de simulaci�n
!
		tempo=dt*dble(i)
!
!		Gravo as enerx�as e o tempo de simulaci�n para saber se todo vai ben
!
    	open (30,file='Resultados_n.dat',access='append')
    		write (30,8004) tempo,E,T,U
    	close (30)
!
!		Gravo as posici�ns, velocidades e aceleraci�ns para o seguinte programa, a�nda que estas
!		�ltimas non fan falta.
!		
        open (40,file='RVA_bin_n.dat',form='unformatted',position='append')
    		write (40) rx,ry,rz,vx,vy,vz,fx,fy,fz
    	close (40)
!
  	endif
!
enddo
!
!Divido as variables de medias entre o n�mero de veces que sumei para obter as medias propiamente
!ditas.
!
Tmedia			= Tmedia		 * pasosi
Umedia			= Umedia		 * pasosi
phivmedia		= phivmedia		 * pasosi
phivvmedia		= phivvmedia	 * pasosi
Tinvmedia		= Tinvmedia		 * pasosi
Tinvphivmedia	= Tinvphivmedia  * pasosi
Tinvphiv2media	= Tinvphiv2media * pasosi
!
!Calculo as constantes caracter�sticas do sistema coas ecuaci�ns en variables reducidas.
!
Temp	= 2.d00*Tmedia/f
Cv		= 1.d00+(2.d00/f-1.d00)*Tmedia*Tinvmedia;Cv=1.d00/Cv
Presion	= ndp*Temp*L3inv-phivmedia
alphaE	= ((1.d00-2.d00/f)*Tmedia*Tinvphivmedia-phivmedia)*L3;alphaE=1.d00/alphaE
gamma	= ndp/Cv+L3*(0.5d00*f-1.d00)*(phivmedia*Tinvmedia-Tinvphivmedia)
kS		= ndp*Temp*L3inv*(1.d00+2*gamma-ndp/Cv)+L3*phivvmedia
kS		= kS-L3*(0.5d00*f-1.d00)*(Tinvphiv2media-2.d00*phivmedia*Tinvphivmedia+Tinvmedia*phivmedia*phivmedia)
kS		= 1.d00/kS
!
!Calculo o momento total e am�soo en pantalla
!
Px=sum(vx)
Py=sum(vy)
Pz=sum(vz)
!
write (*,*) 'P=',Px,Py,Pz
!
!Finalmente gravo as constantes obtidas sempre no mesmo arquivo.
!
open (50,file='Constantes.dat',access='append')
	write (50,8004) Temp,Presion,alphaE,Cv,gamma,kS
close(50)
!
!Gravo as posici�ns, velocidades e aceleraci�ns para a seguinte simulaci�n. As que estaban
!gravadas nese arquivo xa non me interesan.
!
open (20,file='RVAbin.dat',form='unformatted')
    write (20) rx,ry,rz,vx,vy,vz,fx,fy,fz
close(20)
!
call cpu_time(t2)
!
write (*,*) 'O programa tardou:',t2-t1
!
!Fin.
!
stop
end program