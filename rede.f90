program rede
 
!Introdúcese o programa fcc que prepara as condicións iniciais colocando partículas nun cubo en rede
!fcc. Os parámetros de entrada serán os seguintes:
!
!-N: número de partículas do sistema
!-L: lado do cubo de simulación
!-E: enerxía total do sistema
!-Rc: raio de corte do potencial

!Defino o que é un enteiro e o que é un real en doble precisión mediante o módulo prec e afirmo que
!non hai variables implícitas.

use prec
implicit none


!Declaro variables escalares a usar:

!-s,i,x,y,z: variables auxiliares para bucles
!-Parámetro pi=3.14159265458979
!-k: número de celdas por eixo
!-ndp: número de partículas en dobre precisión
!-ldp: lado da caixa de simulación en dobre precisión
!-ldpi: inversa de ldp
!-L3: volume do cubo de simulación
!-L33: inversa do triple do volume
!-rho: densidade de partículas
!-kdp: número de celdas por eixo en dobre precisión
!-a: lado de cada celda
!-U: enerxía potencial total
!-phiv: 1ª derivada da enerxía potencial
!-phivv: 2ª derivada da enerxía potencial
!-T0: enerxía cinética "aleatoria"
!-T: enerxía cinética real
!-E: enerxía total
!-factor: factor que usarei para reescalar a enerxía cinética
!-Px,Py,Pz: compoñentes do momento total do sistema
!-Ax,Ay,Az: compoñentes da aceleración total do sistema
!-Rc2: raio de corte ao cadrado
!-Rc3: inversa do raio de corte ao cubo
!-Rc6: inversa do raio de corte elevado a seis
!-cor: factor a usar nas correccións da enerxía potencial e as súas derivadas

integer (kind=int),parameter :: N=500,L=10,Rc=5
real (kind=dp),parameter :: pi=3.14159265458979
integer (kind=int) :: s,i,k,x,y,z
real (kind=dp) :: ndp,ldp,ldpi,L3,L33,rho,kdp,a,U,phiv,phivv,T0,T,E,factor,Px,Py,Pz,Ax,Ay,Az
real (kind=dp) :: Rc2,Rc3,Rc6,cor

!Declaro variables vectoriais a usar:

!-rx,ry,rz: compoñentes dos vectores posición de cada partícula
!-vx,vy,vz: compoñentes dos vectores velocidade de cada partícula
!-fx,fy,fz: compoñentes dos vectores aceleración de cada partícula
!-difx,dify,difz: compoñentes dos vectores de desprazamentos infinitesimais das partículas

real (kind=dp), dimension(N) :: rx,ry,rz,vx,vy,vz,fx,fy,fz,difx,dify,difz

!Defino variables que usarei máis adiante.

ndp=dble(N)               !Número de partículas en dobre precisión
ldp=dble(L)               !Lado da caixa en dobre precisión
ldpi=1.d00/ldp            !Inversa do lado da caixa
Rc2=dble(Rc*Rc)           !Raio de corte ao cadrado
Rc3=1/(Rc2*dble(Rc))      !Inversa do raio de corte ao cubo
Rc6=Rc3*Rc3               !Inversa da sexta potencia do raio de corte

!Calculo as variables propias da caixa de simulación que usarei no programa:

L3=ldp*ldp*ldp                   !Volume
rho=ndp*ldpi*ldpi*ldpi           !Densidade
kdp=(0.25d00*ndp)**(1.d00/3.d00) !Nº de celdas por eixo en dobre precisión
a=ldp/kdp                        !Lado da celda
k=nint(kdp)                      !Nº de celdas por eixo en formato enteiro

!Calculo outros factores para usar máis adiante

L33=1.d00/(3*L3)                 !=1/3V
cor=16.d00*pi*rho*ndp*Rc3        !=16·pi·rho·N/Rc^3

!Comprobo que os parámetros sexan adecuados.

if ((rho>1.d00) .or. (rho<0.1d00)) then
	write(*,*)'Erro: a densidade toma un valor extremo.'
    pause
    stop
endif

if ((kdp-dble(k))>0.001d00) then
 	 write (*,*) 'Erro: N debe ser da forma N=4k^3.)'
 	 pause
 	 stop
endif

if (dble(Rc)>0.5d00*ldp) then
	write(*,*)'Erro: o raio de corte non pode ser maior que a metade da caixa de simulación.'
    pause
    stop
endif

!Para a composición das coordenadas, collo unha celda e poño partículas en (0,0,0),(a/2,a/2,0),
!(a/2,0,a/2) e (0,a/2,a/2). Repito a operación desprazándome por tódolos x's, y's e z's. Polo tanto
!as coordenadas das partículas dunha celda calquera son as representadas abaixo. O desprazamento de
!calquera coordenada é, obviamente, o tamaño da celda: "a". Os bucles corren no número de celda co-
!mezando en 0 para facilitar a asignación de coordenadas. En consecuencia a última celda debe ser a
!k-1. As coordenadas das celdas deben ir de 0 ata L-a (última celda) porque a celda situada en x=L,
!y=L ou z=L colocaría partículas fóra da caixa.

!Defino o índice s que recorrerá tódolos elementos dos vectores posición. Póñoo a cero para
!que comece no elemento 1.

s=0
do z=0,k-1
    do y=0,k-1
        do x=0,k-1
            s=s+1
            rx(s)=x*a
            ry(s)=y*a
            rz(s)=z*a
            s=s+1
            rx(s)=(x+0.5d00)*a
            ry(s)=(y+0.5d00)*a
            rz(s)=z*a
            s=s+1
            rx(s)=(x+0.5d00)*a
            ry(s)=y*a
            rz(s)=(z+0.5d00)*a
            s=s+1
            rx(s)=x*a
            ry(s)=(y+0.5d00)*a
            rz(s)=(z+0.5d00)*a
        enddo
    enddo
enddo

!Introduzo unha condición de erro por se non se colocaron tódalas partículas.

if (s/=N) then
  write (*,*) 'Non se colocaron tódalas partículas'
  write (*,*) '   s=',s,'e   N=',N
endif

!Agora pretendo mover tódalas partículas unha distancia "infinitesimal" da súa posición inicial de
!rede fcc. Para iso creo uns vectores de números pseudo-aleatorios entre 0 e 1.

call random_number(difx)
call random_number(dify)
call random_number(difz)

!Agora reescalo estes vectores para que vaian entre -0.001 e 0.001 e desprazo as coordenadas iniciais
!nesa cantidade. Canto máis grande sexa a caixa menos significativo será.

rx=rx+0.001d00*(2*difx-1)
ry=ry+0.001d00*(2*dify-1)
rz=rz+0.001d00*(2*difz-1)

!A continuación,calculo a enerxía potencial. Para iso chamo á subrutina potencialLJ que ten como
!variables de entrada N,ldp,ldpi,L33,Rc2,Rc6,rx,ry,rz,cor e devolve U,phiv,phivv,fx,fy,fz. A ex-
!plicación inclúese na subrutina.

call potencialLJ(N,ldp,ldpi,L33,Rc2,Rc6,rx,ry,rz,cor,U,phiv,phivv,fx,fy,fz)

!A continuación realizo a asignación de velocidades aleatorias a cada partícula. Para iso xenero
!tres vectores de números pseudo-aleatorios entre 0 e 1 que son as compoñentes das velocidades para
!cada partícula.

call random_number(vx)
call random_number(vy)
call random_number(vz)

!Coas seguintes ordes reescalo as velocidades para que vaian de -1 a 1.

vx=2*vx-1
vy=2*vy-1
vz=2*vz-1

!Calculo o momento total que como son variables reducidas non é máis que as sumas de compoñentes da
!velocidade.

Px=sum(vx)
Py=sum(vy)
Pz=sum(vz)

!Como non será cero (e debe selo), desprazamos as compoñentes das velocidades nunha cantidade igual
!á súa media (coincide por ser variables reducidas) para que o momento total sexa nulo. Debe selo
!ao ser a colectividade microcanónica e, polo tanto, un sistema illado.

vx=vx-Px/N
vy=vy-Py/N
vz=vz-Pz/N

!Calculo a enerxía cinética "aleatoria" resultante

T0=0.5*sum(vx*vx+vy*vy+vz*vz)

!Amoso en pantalla os valores obtidos para a enerxía cinética e potencial e pido que se introduza
!a enerxía total.

write (*,*) '   U=',U
write (*,*) '   T=',T0
write (*,*) 'Insertar enerxía total'
read (*,*) E

!Calculo a nova enerxía cinética e se a enerxía introducida non é adecuada dará un erro.

T=E-U

if(T<0) then
  write (*,*) 'Erro: enerxía cinética negativa'
  stop
endif

!Reescalo as velocidades mediante as seguintes liñas. Obviamente o momento total segue sendo nulo.

factor=sqrt(T/T0)

vx=factor*vx
vy=factor*vy
vz=factor*vz

!Calculo o momento (de novo) e a aceleración total:

Px=sum(vx)
Py=sum(vy)
Pz=sum(vz)

Ax=sum(fx)
Ay=sum(fy)
Az=sum(fz)

!Amoso en pantalla tódalas variables interesantes da simulación:

write (*,*) ' N=',N           !Número de partículas
write (*,*) ' a=',a           !Lado de cada celda
write (*,*) ' k=',k           !Número de celdas por eixo
write (*,*) ' rho=',rho       !Densidade de partículas
write (*,*) ' L=',L           !Lado da caixa de simulación
write (*,*) ' V=',L3          !Volume da caixa de simulación
write (*,*) ' U=',U           !Enerxía potencial total
write (*,*) ' dU/dV=',phiv    !Primeira derivada da enerxía potencial
write (*,*) ' d2U/dV2=',phivv !Segunda derivada da enerxía potencial
write (*,*) ' T=',T           !Enerxía cinética
write (*,*) ' E=',E           !Enerxía total
write (*,*)
write (*,*) ' P=',Px,Py,Pz    !Momento total
write (*,*) ' F=',Ax,Ay,Az    !Aceleración total

!Defino os formatos a usar no gravado de arquivos.

9000	format(1pe13.6,8(1x,e13.6))
8000	format(1pe19.12,8(1x,e19.12))
8001	format(1pe19.12,4(1x,e19.12))
8002	format(4(i8))

!Gravo en ASCII os vectores de posicións, velocidades e aceleracións
!As tres primeiras columnas son x,y,z, as tres seguintes vx,vy,vz e as tres últimas ax,ay,az.

open(10,file='Matrix.dat')
	do i=1,N
	  write(10,9000) rx(i),ry(i),rz(i),vx(i),vy(i),vz(i),fx(i),fy(i),fz(i)
    enddo
close(10)

!Gravo en binario os vectores de posicións, velocidades e aceleracións.

open(20,file='Matrixbin.dat',form='unformatted')
	write(20) rx,ry,rz,vx,vy,vz,fx,fy,fz
close(20)

!Gravo en ASCII as variables escalares de interés:

open(30,file='Números.dat')
	write (30,8002) N,k,L,Rc                           !Enteiros
    write (30,8000) ndp,ldp,ldpi,a,rho,L3,Rc2,Rc3,Rc6  !Reais 1
	write (30,8001) U,phiv,phivv,T,E                   !Reais 2
close(30)

!Fin.

stop

end program rede