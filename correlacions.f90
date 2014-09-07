program correlacions
!
!Introd�cese o programa correlaci�ns que calcula a desviaci�n cuadr�tica media, as correlaci�ns das
!velocidades e a funci�n de distribuci�n radial.
!
!Defino precisi�n e digo que non hai variables impl�citas.
!
use prec
implicit none
!
!Declaro variables:
!
!-N:    	  	n�mero de part�culas do sistema
!-k: 		  	n�mero de celdas por eixo
!-L: 		  	lado da cubo de simulaci�n
!-Rc:		  	raio de corte do potencial
!-i,j,m:	  	variables para correr bucles
!-r:		  	variable para asignar a posici�n na funci�n de distribuci�n radial
!-L2pri:		lonxitude dos vectores da fdr
!-ndp: 		  	n�mero de part�culas en dobre precisi�n
!-ndpi:		  	inversa ddo n�mero de part�culas
!-ldp: 		  	lado da cubo en dobre precisi�n
!-ldpi: 	  	inversa de ldp
!-a: 		  	lado de cada celda
!-rho: 		  	densidade de part�culas
!-L3: 		  	volume da caixa de simulaci�n
!-dt:         	paso temporal da simulaci�n
!-t1,t2: 	  	variables para controlar o tempo de execuci�n do programa
!-Ntotal:		n�mero de pasos temporais totais
!-ttotal:		n�mero de pasos temporais avanzados dende o de referencia
!-paso:			paso temporal para coller tempos de referencia
!-pr:		 	precisi�n � hora de asignar posici�ns na funci�n de distribuci�n radial
!-pri:		 	inversa de pr
!-pi:		  	n�mero pi
!-rx,ry,rz:   	posici�ns das N part�culas
!-vx,vy,vz:   	velocidades das N part�culas
!-fx,fy,fz:   	aceleraci�ns das N part�culas
!-drx,dry,drz: 	variables para almacenar a distancia en cada eixo entre d�as part�culas
!-dv:			variable para usar nas correlaci�ns de velocidades
!-modr,modv:	promedio a t�dalas part�culas a un tempo dado
!-corrr,corrv:	vectores de desprazamento cuadr�tico medio e correlaci�ns de velocidades
!-D:			variable para a obtenci�n do coeficiente de difusi�n
!-x,y,z:		variables auxiliares para o c�lculo de distancias entre d�as part�culas
!-rimx,y,z:		distancia en cada eixo entre d�as part�culas
!-dim2:			m�dulo ao cadrado da distancia
!-gasid:		factor dependente de r de gas ideal
!-fgas:			factor constante de gas ideal
!-corte:		factor de corte de distancia que se considera para a fdr
!-fdr:			funci�n de distribuci�n radial non normalizada
!-fdrn:			funci�n de distribuci�n radial normalizada
!
integer(kind=int) :: N,k,L,Rc,i,j,m,r,L2pri
real(kind=dp) :: ndp,ndpi,ldp,ldpi,a,rho,L3,dt,t1,t2
integer(kind=int),parameter :: Ntotal=5000,ttotal=300,paso=10
real (kind=dp),parameter :: pr=0.001,pri=1/pr,pi=3.1415926535898d00
real(kind=dp),dimension(Ntotal,500) :: rx,ry,rz,vx,vy,vz,fx,fy,fz
real(kind=dp),dimension(500) :: drx,dry,drz,dv
real(kind=dp),dimension(ttotal) :: modr,modv,tempo,corrr,corrv
real(kind=dp) :: D,x,y,z,rimx,rimy,rimz,dim2,gasid,fgas,corte
integer(kind=int),allocatable,dimension(:) :: fdr
real(kind=dp),allocatable,dimension(:) :: fdrn

!
!Defino os formatos que vou usar
!
8000	format(1pe19.12,5(1x,e19.12))
8001	format(1pe19.12,1x,e19.12,1x,e19.12)
8002	format(i8,i8,i8,i8)
8003	format(1pe19.12)
8004	format(i8,1pe19.12)
!
call cpu_time(t1)
!
!Collo variables interesantes para a simulaci�n, algunhas non son necesarias
!
open(10,file='N�meros.dat')
	read (10,8002) N,k,L,Rc 						  !Enteiros
    read (10,8003) dt                                 !Paso temporal
    read (10,8000) ndp,ldp,ldpi,a,rho,L3			  !Reais
close(10)
!
!Leo as posici�ns, velocidades e aceleraci�ns dos arquivos de sa�da po programa de din�mica molecular
!As aceleraci�ns non son necesarias no programa, pero como o arquivo as cont�n l�oas para poder ler
!ben as outras variables.
!
open(20,file='RVA_bin.dat',form='unformatted')
	do i=1,5000
		read(20) rx(i,:),ry(i,:),rz(i,:),vx(i,:),vy(i,:),vz(i,:),fx(i,:),fy(i,:),fz(i,:)
	enddo
close(20)
!
!Calculo variables que far�n falta m�is adiante
!
ndpi=1.d00/ndp
L2pri=nint(0.5d00*ldp*pri)
fgas=4.d00/3.d00*pi*rho
!
!Anulo os vectores acumuladores para comezar o c�lculo das correlaci�ns.
!
corrr=0
corrv=0
!
!Escollo un tempo de referencia e calculo o desprazamento cuadr�tico medio desa part�cula para os
!seguintes "ttotal" instantes de tempo e tam�n as correlaci�ns das velocidades. Promedio sobre t�-
!dalas part�culas e acumulo o vector nunha variable.
!
!
do i=0,Ntotal-ttotal,paso					!Escollo o tempo de refencia cada "paso"
   do j=i,i+ttotal-1						!Para cada instante de tempo posterior fago as seguintes contas:
     m=j-i+1								!Dime o instante de tempo no que me atopo
     drx=rx(j+1,:)-rx(i+1,:);drx=drx*drx	!Calculo a distancia en x e el�voa ao cadrado
     dry=ry(j+1,:)-ry(i+1,:);dry=dry*dry	!   "    "     "     "  y "   "    "    "
     drz=rz(j+1,:)-rz(i+1,:);drz=drz*drz	!   "    "     "     "  z "   "    "    "
     modr(m)=sum(drx+dry+drz)				!Calculo o m�dulo cadrado da distancia promediado a t�dalas part�culas
     dv=vx(j+1,:)*vx(i+1,:)					!Multiplico a compo�ente x da velocidade do tempo de referencia pola de outro tempo
     dv=dv+vy(j+1,:)*vy(i+1,:)				!Repito para y e acumulo
     dv=dv+vz(j+1,:)*vz(i+1,:)				!Repito para z e acumulo
     modv(m)=sum(dv)						!Promedio a t�dalas part�culas
   enddo
   corrr=corrr+modr							!Acumulo o vector nunha variable para facer a gr�fica
   corrv=corrv+modv							!  "     "   "      "      "      "     "   "    "
enddo
!
!Introduzo os factores divisores:
!
!-ndpi:						inversa do n�mero de part�culas porque promediei sobre estas
!-(Ntotal-ttotal)/paso:		n�mero de instantes temporais tomados de referencia e, polo tanto acumulados
!
corrr=corrr*ndpi*dble(paso)/dble(Ntotal-ttotal)
corrv=corrv*ndpi*dble(paso)/dble(Ntotal-ttotal)
!
!Monto o vector de instantes de tempo para representar os vectores anteriores fronte a este
!
do i=0,ttotal-1
  tempo(i+1)=dble(i)*dt
enddo
!
!Agora calculo o coefifiente de difusi�n do sistema que � 1/6 da pendente cando t-> inf do desprazamento
!cuadr�tico medio fronte ao tempo. Para iso collo os valores da pendente a partir da metade e promedio.
!
D=0							!Anulo o acumulador
!
do i=ttotal/2,ttotal		!Faise unha divisi�n enteira, pero dame igual coller unha m�is ou unha menos
  D=D+corrr(i)/tempo(i) 	!Acumulo a pendente
enddo
!
D=D/dble(3*ttotal)			!Divido por o n�mero de termos (ttotal/2) e por 6, en total divido por 3�ttotal
!
!Gardo nun arquivo en ASCII o tempo, o desprazamento cuadr�tico medio e a correlaci�n das velocidades
!
open(30,file='correlacions.dat')
	do i=1,ttotal
    	write(30,8001) tempo(i),corrr(i),corrv(i)
    enddo
close(30)
!
!Agora paso ao c�lculo da funci�n de distribuci�n radial.
!
!Uso unha variable que define o cadrado do l�mite a ter en conta na distancias entre part�culas
!
corte=0.5d00*(ldp-pr);corte=corte*corte
!
!Asigno memoria para os vectores que me servir�n para acumular a funci�n de distribuci�n radial
!
allocate(fdr(0:L2pri-1),fdrn(0:L2pr)-1)
!
!Anulo o acumulador da funci�n de distrbuci�n total
!
fdrn=0
!
!Comezo o bucle.Para cada instante de tempo terei unha funci�n de ditribuci�n radial que acumularei
!e ao final promedio sobre todas.
!
do i=1,Ntotal
  fdr=0			!Anulo o acumulador da funci�n de distribuci�n dun tempo dado
  do j=1,N-1	!Corro sobre t�dalas part�culas (menos a �ltima) e t�moa de referencia
!
    x=rx(i,j)	!Collo a posici�n en x desta
    y=ry(i,j)	!  "   "    "     "  y   "
    z=rz(i,j)	!  "   "    "     "  z   "
!
    do m=j+1,N	!Corro por t�dalas part�culas restantes (as que non deixei atr�s)
!
      rimx=rx(i,m)-x;rimx=rimx-ldp*dnint(rimx*ldpi)		!Calculo a distancia en x da imaxe m�is pr�xima coa part�cula de referencia
      rimy=ry(i,m)-y;rimy=rimy-ldp*dnint(rimy*ldpi) 	!   "    "     "     "  y  "   "    "      "     "		"	  "		 "
      rimz=rz(i,m)-z;rimz=rimz-ldp*dnint(rimz*ldpi) 	!   "    "     "     "  z  "   "    "      "     "		"	  "	 	 "
!
      dim2=rimx*rimx+rimy*rimy+rimz*rimz	!Calculo o m�dulo ao cadrado da distancia
!
      if(dim2<=corte)then					!Se esta distancia � menor que o corte que establec�n ao principio:
        	r=nint(dsqrt(dim2)*pri)			!Asigno esa distancia a un valor da fdr
            fdr(r)=fdr(r)+2					!Acumulo 2 part�culas nesa distancia (a part�cula j e a k)
      endif
!
    enddo
  enddo
  fdrn=fdrn+dble(fdr)						!Acumulo a funci�n de distribuci�n do tempo i na total
enddo
!
!Agora divido polo factor dependente de r de gas ideal.
!
do i=1,L2pri
  gasid=((dble(i)+0.5d00)*pr)**3.d00-((dble(i)-0.5d00)*pr)**3.d00
  fdrn(i)=fdrn(i)/gasid
enddo
!
!Normalizo a funci�n de distribuci�n:
!-fgas:		divido polo factor constante de gas ideal
!-ndpi:		multiplico pola inversa do n�mero de part�culas
!-Ntotal:	divido polo n�mero de instantes temporais acumulados
!
fdrn=fdrn*ndpi/(fgas*dble(Ntotal))
!
!Gardo a distancia e os valores da funci�n de distribuci�n radial nun arquivo en ASCII.
!
open(40,file='fdr.dat')
	do i=1,L2pri
    	write(40,8004) i,fdrn(i)
    enddo
close(40)
!
deallocate(fdr,fdrn)
!
!Gardo o coeficiente de difusi�n.
!
open(50,file='Coef_difusion.dat')
	write(50,8003) D
close(50)
!
call cpu_time(t2)
write(*,*) t2-t1
!
!Fin.
!
stop
end program