program correlacions
!
!Introdúcese o programa correlacións que calcula a desviación cuadrática media, as correlacións das
!velocidades e a función de distribución radial.
!
!Defino precisión e digo que non hai variables implícitas.
!
use prec
implicit none
!
!Declaro variables:
!
!-N:    	  	número de partículas do sistema
!-k: 		  	número de celdas por eixo
!-L: 		  	lado da cubo de simulación
!-Rc:		  	raio de corte do potencial
!-i,j,m:	  	variables para correr bucles
!-r:		  	variable para asignar a posición na función de distribución radial
!-L2pri:		lonxitude dos vectores da fdr
!-ndp: 		  	número de partículas en dobre precisión
!-ndpi:		  	inversa ddo número de partículas
!-ldp: 		  	lado da cubo en dobre precisión
!-ldpi: 	  	inversa de ldp
!-a: 		  	lado de cada celda
!-rho: 		  	densidade de partículas
!-L3: 		  	volume da caixa de simulación
!-dt:         	paso temporal da simulación
!-t1,t2: 	  	variables para controlar o tempo de execución do programa
!-Ntotal:		número de pasos temporais totais
!-ttotal:		número de pasos temporais avanzados dende o de referencia
!-paso:			paso temporal para coller tempos de referencia
!-pr:		 	precisión á hora de asignar posicións na función de distribución radial
!-pri:		 	inversa de pr
!-pi:		  	número pi
!-rx,ry,rz:   	posicións das N partículas
!-vx,vy,vz:   	velocidades das N partículas
!-fx,fy,fz:   	aceleracións das N partículas
!-drx,dry,drz: 	variables para almacenar a distancia en cada eixo entre dúas partículas
!-dv:			variable para usar nas correlacións de velocidades
!-modr,modv:	promedio a tódalas partículas a un tempo dado
!-corrr,corrv:	vectores de desprazamento cuadrático medio e correlacións de velocidades
!-D:			variable para a obtención do coeficiente de difusión
!-x,y,z:		variables auxiliares para o cálculo de distancias entre dúas partículas
!-rimx,y,z:		distancia en cada eixo entre dúas partículas
!-dim2:			módulo ao cadrado da distancia
!-gasid:		factor dependente de r de gas ideal
!-fgas:			factor constante de gas ideal
!-corte:		factor de corte de distancia que se considera para a fdr
!-fdr:			función de distribución radial non normalizada
!-fdrn:			función de distribución radial normalizada
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
!Collo variables interesantes para a simulación, algunhas non son necesarias
!
open(10,file='Números.dat')
	read (10,8002) N,k,L,Rc 						  !Enteiros
    read (10,8003) dt                                 !Paso temporal
    read (10,8000) ndp,ldp,ldpi,a,rho,L3			  !Reais
close(10)
!
!Leo as posicións, velocidades e aceleracións dos arquivos de saída po programa de dinámica molecular
!As aceleracións non son necesarias no programa, pero como o arquivo as contén léoas para poder ler
!ben as outras variables.
!
open(20,file='RVA_bin.dat',form='unformatted')
	do i=1,5000
		read(20) rx(i,:),ry(i,:),rz(i,:),vx(i,:),vy(i,:),vz(i,:),fx(i,:),fy(i,:),fz(i,:)
	enddo
close(20)
!
!Calculo variables que farán falta máis adiante
!
ndpi=1.d00/ndp
L2pri=nint(0.5d00*ldp*pri)
fgas=4.d00/3.d00*pi*rho
!
!Anulo os vectores acumuladores para comezar o cálculo das correlacións.
!
corrr=0
corrv=0
!
!Escollo un tempo de referencia e calculo o desprazamento cuadrático medio desa partícula para os
!seguintes "ttotal" instantes de tempo e tamén as correlacións das velocidades. Promedio sobre tó-
!dalas partículas e acumulo o vector nunha variable.
!
!
do i=0,Ntotal-ttotal,paso					!Escollo o tempo de refencia cada "paso"
   do j=i,i+ttotal-1						!Para cada instante de tempo posterior fago as seguintes contas:
     m=j-i+1								!Dime o instante de tempo no que me atopo
     drx=rx(j+1,:)-rx(i+1,:);drx=drx*drx	!Calculo a distancia en x e elévoa ao cadrado
     dry=ry(j+1,:)-ry(i+1,:);dry=dry*dry	!   "    "     "     "  y "   "    "    "
     drz=rz(j+1,:)-rz(i+1,:);drz=drz*drz	!   "    "     "     "  z "   "    "    "
     modr(m)=sum(drx+dry+drz)				!Calculo o módulo cadrado da distancia promediado a tódalas partículas
     dv=vx(j+1,:)*vx(i+1,:)					!Multiplico a compoñente x da velocidade do tempo de referencia pola de outro tempo
     dv=dv+vy(j+1,:)*vy(i+1,:)				!Repito para y e acumulo
     dv=dv+vz(j+1,:)*vz(i+1,:)				!Repito para z e acumulo
     modv(m)=sum(dv)						!Promedio a tódalas partículas
   enddo
   corrr=corrr+modr							!Acumulo o vector nunha variable para facer a gráfica
   corrv=corrv+modv							!  "     "   "      "      "      "     "   "    "
enddo
!
!Introduzo os factores divisores:
!
!-ndpi:						inversa do número de partículas porque promediei sobre estas
!-(Ntotal-ttotal)/paso:		número de instantes temporais tomados de referencia e, polo tanto acumulados
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
!Agora calculo o coefifiente de difusión do sistema que é 1/6 da pendente cando t-> inf do desprazamento
!cuadrático medio fronte ao tempo. Para iso collo os valores da pendente a partir da metade e promedio.
!
D=0							!Anulo o acumulador
!
do i=ttotal/2,ttotal		!Faise unha división enteira, pero dame igual coller unha máis ou unha menos
  D=D+corrr(i)/tempo(i) 	!Acumulo a pendente
enddo
!
D=D/dble(3*ttotal)			!Divido por o número de termos (ttotal/2) e por 6, en total divido por 3·ttotal
!
!Gardo nun arquivo en ASCII o tempo, o desprazamento cuadrático medio e a correlación das velocidades
!
open(30,file='correlacions.dat')
	do i=1,ttotal
    	write(30,8001) tempo(i),corrr(i),corrv(i)
    enddo
close(30)
!
!Agora paso ao cálculo da función de distribución radial.
!
!Uso unha variable que define o cadrado do límite a ter en conta na distancias entre partículas
!
corte=0.5d00*(ldp-pr);corte=corte*corte
!
!Asigno memoria para os vectores que me servirán para acumular a función de distribución radial
!
allocate(fdr(0:L2pri-1),fdrn(0:L2pr)-1)
!
!Anulo o acumulador da función de distrbución total
!
fdrn=0
!
!Comezo o bucle.Para cada instante de tempo terei unha función de ditribución radial que acumularei
!e ao final promedio sobre todas.
!
do i=1,Ntotal
  fdr=0			!Anulo o acumulador da función de distribución dun tempo dado
  do j=1,N-1	!Corro sobre tódalas partículas (menos a última) e tómoa de referencia
!
    x=rx(i,j)	!Collo a posición en x desta
    y=ry(i,j)	!  "   "    "     "  y   "
    z=rz(i,j)	!  "   "    "     "  z   "
!
    do m=j+1,N	!Corro por tódalas partículas restantes (as que non deixei atrás)
!
      rimx=rx(i,m)-x;rimx=rimx-ldp*dnint(rimx*ldpi)		!Calculo a distancia en x da imaxe máis próxima coa partícula de referencia
      rimy=ry(i,m)-y;rimy=rimy-ldp*dnint(rimy*ldpi) 	!   "    "     "     "  y  "   "    "      "     "		"	  "		 "
      rimz=rz(i,m)-z;rimz=rimz-ldp*dnint(rimz*ldpi) 	!   "    "     "     "  z  "   "    "      "     "		"	  "	 	 "
!
      dim2=rimx*rimx+rimy*rimy+rimz*rimz	!Calculo o módulo ao cadrado da distancia
!
      if(dim2<=corte)then					!Se esta distancia é menor que o corte que establecín ao principio:
        	r=nint(dsqrt(dim2)*pri)			!Asigno esa distancia a un valor da fdr
            fdr(r)=fdr(r)+2					!Acumulo 2 partículas nesa distancia (a partícula j e a k)
      endif
!
    enddo
  enddo
  fdrn=fdrn+dble(fdr)						!Acumulo a función de distribución do tempo i na total
enddo
!
!Agora divido polo factor dependente de r de gas ideal.
!
do i=1,L2pri
  gasid=((dble(i)+0.5d00)*pr)**3.d00-((dble(i)-0.5d00)*pr)**3.d00
  fdrn(i)=fdrn(i)/gasid
enddo
!
!Normalizo a función de distribución:
!-fgas:		divido polo factor constante de gas ideal
!-ndpi:		multiplico pola inversa do número de partículas
!-Ntotal:	divido polo número de instantes temporais acumulados
!
fdrn=fdrn*ndpi/(fgas*dble(Ntotal))
!
!Gardo a distancia e os valores da función de distribución radial nun arquivo en ASCII.
!
open(40,file='fdr.dat')
	do i=1,L2pri
    	write(40,8004) i,fdrn(i)
    enddo
close(40)
!
deallocate(fdr,fdrn)
!
!Gardo o coeficiente de difusión.
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