program rede
 
!Introd�cese o programa fcc que prepara as condici�ns iniciais colocando part�culas nun cubo en rede
!fcc. Os par�metros de entrada ser�n os seguintes:
!
!-N: n�mero de part�culas do sistema
!-L: lado do cubo de simulaci�n
!-E: enerx�a total do sistema
!-Rc: raio de corte do potencial

!Defino o que � un enteiro e o que � un real en doble precisi�n mediante o m�dulo prec e afirmo que
!non hai variables impl�citas.

use prec
implicit none


!Declaro variables escalares a usar:

!-s,i,x,y,z: variables auxiliares para bucles
!-Par�metro pi=3.14159265458979
!-k: n�mero de celdas por eixo
!-ndp: n�mero de part�culas en dobre precisi�n
!-ldp: lado da caixa de simulaci�n en dobre precisi�n
!-ldpi: inversa de ldp
!-L3: volume do cubo de simulaci�n
!-L33: inversa do triple do volume
!-rho: densidade de part�culas
!-kdp: n�mero de celdas por eixo en dobre precisi�n
!-a: lado de cada celda
!-U: enerx�a potencial total
!-phiv: 1� derivada da enerx�a potencial
!-phivv: 2� derivada da enerx�a potencial
!-T0: enerx�a cin�tica "aleatoria"
!-T: enerx�a cin�tica real
!-E: enerx�a total
!-factor: factor que usarei para reescalar a enerx�a cin�tica
!-Px,Py,Pz: compo�entes do momento total do sistema
!-Ax,Ay,Az: compo�entes da aceleraci�n total do sistema
!-Rc2: raio de corte ao cadrado
!-Rc3: inversa do raio de corte ao cubo
!-Rc6: inversa do raio de corte elevado a seis
!-cor: factor a usar nas correcci�ns da enerx�a potencial e as s�as derivadas

integer (kind=int),parameter :: N=500,L=10,Rc=5
real (kind=dp),parameter :: pi=3.14159265458979
integer (kind=int) :: s,i,k,x,y,z
real (kind=dp) :: ndp,ldp,ldpi,L3,L33,rho,kdp,a,U,phiv,phivv,T0,T,E,factor,Px,Py,Pz,Ax,Ay,Az
real (kind=dp) :: Rc2,Rc3,Rc6,cor

!Declaro variables vectoriais a usar:

!-rx,ry,rz: compo�entes dos vectores posici�n de cada part�cula
!-vx,vy,vz: compo�entes dos vectores velocidade de cada part�cula
!-fx,fy,fz: compo�entes dos vectores aceleraci�n de cada part�cula
!-difx,dify,difz: compo�entes dos vectores de desprazamentos infinitesimais das part�culas

real (kind=dp), dimension(N) :: rx,ry,rz,vx,vy,vz,fx,fy,fz,difx,dify,difz

!Defino variables que usarei m�is adiante.

ndp=dble(N)               !N�mero de part�culas en dobre precisi�n
ldp=dble(L)               !Lado da caixa en dobre precisi�n
ldpi=1.d00/ldp            !Inversa do lado da caixa
Rc2=dble(Rc*Rc)           !Raio de corte ao cadrado
Rc3=1/(Rc2*dble(Rc))      !Inversa do raio de corte ao cubo
Rc6=Rc3*Rc3               !Inversa da sexta potencia do raio de corte

!Calculo as variables propias da caixa de simulaci�n que usarei no programa:

L3=ldp*ldp*ldp                   !Volume
rho=ndp*ldpi*ldpi*ldpi           !Densidade
kdp=(0.25d00*ndp)**(1.d00/3.d00) !N� de celdas por eixo en dobre precisi�n
a=ldp/kdp                        !Lado da celda
k=nint(kdp)                      !N� de celdas por eixo en formato enteiro

!Calculo outros factores para usar m�is adiante

L33=1.d00/(3*L3)                 !=1/3V
cor=16.d00*pi*rho*ndp*Rc3        !=16�pi�rho�N/Rc^3

!Comprobo que os par�metros sexan adecuados.

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
	write(*,*)'Erro: o raio de corte non pode ser maior que a metade da caixa de simulaci�n.'
    pause
    stop
endif

!Para a composici�n das coordenadas, collo unha celda e po�o part�culas en (0,0,0),(a/2,a/2,0),
!(a/2,0,a/2) e (0,a/2,a/2). Repito a operaci�n despraz�ndome por t�dolos x's, y's e z's. Polo tanto
!as coordenadas das part�culas dunha celda calquera son as representadas abaixo. O desprazamento de
!calquera coordenada �, obviamente, o tama�o da celda: "a". Os bucles corren no n�mero de celda co-
!mezando en 0 para facilitar a asignaci�n de coordenadas. En consecuencia a �ltima celda debe ser a
!k-1. As coordenadas das celdas deben ir de 0 ata L-a (�ltima celda) porque a celda situada en x=L,
!y=L ou z=L colocar�a part�culas f�ra da caixa.

!Defino o �ndice s que recorrer� t�dolos elementos dos vectores posici�n. P��oo a cero para
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

!Introduzo unha condici�n de erro por se non se colocaron t�dalas part�culas.

if (s/=N) then
  write (*,*) 'Non se colocaron t�dalas part�culas'
  write (*,*) '   s=',s,'e   N=',N
endif

!Agora pretendo mover t�dalas part�culas unha distancia "infinitesimal" da s�a posici�n inicial de
!rede fcc. Para iso creo uns vectores de n�meros pseudo-aleatorios entre 0 e 1.

call random_number(difx)
call random_number(dify)
call random_number(difz)

!Agora reescalo estes vectores para que vaian entre -0.001 e 0.001 e desprazo as coordenadas iniciais
!nesa cantidade. Canto m�is grande sexa a caixa menos significativo ser�.

rx=rx+0.001d00*(2*difx-1)
ry=ry+0.001d00*(2*dify-1)
rz=rz+0.001d00*(2*difz-1)

!A continuaci�n,calculo a enerx�a potencial. Para iso chamo � subrutina potencialLJ que ten como
!variables de entrada N,ldp,ldpi,L33,Rc2,Rc6,rx,ry,rz,cor e devolve U,phiv,phivv,fx,fy,fz. A ex-
!plicaci�n incl�ese na subrutina.

call potencialLJ(N,ldp,ldpi,L33,Rc2,Rc6,rx,ry,rz,cor,U,phiv,phivv,fx,fy,fz)

!A continuaci�n realizo a asignaci�n de velocidades aleatorias a cada part�cula. Para iso xenero
!tres vectores de n�meros pseudo-aleatorios entre 0 e 1 que son as compo�entes das velocidades para
!cada part�cula.

call random_number(vx)
call random_number(vy)
call random_number(vz)

!Coas seguintes ordes reescalo as velocidades para que vaian de -1 a 1.

vx=2*vx-1
vy=2*vy-1
vz=2*vz-1

!Calculo o momento total que como son variables reducidas non � m�is que as sumas de compo�entes da
!velocidade.

Px=sum(vx)
Py=sum(vy)
Pz=sum(vz)

!Como non ser� cero (e debe selo), desprazamos as compo�entes das velocidades nunha cantidade igual
!� s�a media (coincide por ser variables reducidas) para que o momento total sexa nulo. Debe selo
!ao ser a colectividade microcan�nica e, polo tanto, un sistema illado.

vx=vx-Px/N
vy=vy-Py/N
vz=vz-Pz/N

!Calculo a enerx�a cin�tica "aleatoria" resultante

T0=0.5*sum(vx*vx+vy*vy+vz*vz)

!Amoso en pantalla os valores obtidos para a enerx�a cin�tica e potencial e pido que se introduza
!a enerx�a total.

write (*,*) '   U=',U
write (*,*) '   T=',T0
write (*,*) 'Insertar enerx�a total'
read (*,*) E

!Calculo a nova enerx�a cin�tica e se a enerx�a introducida non � adecuada dar� un erro.

T=E-U

if(T<0) then
  write (*,*) 'Erro: enerx�a cin�tica negativa'
  stop
endif

!Reescalo as velocidades mediante as seguintes li�as. Obviamente o momento total segue sendo nulo.

factor=sqrt(T/T0)

vx=factor*vx
vy=factor*vy
vz=factor*vz

!Calculo o momento (de novo) e a aceleraci�n total:

Px=sum(vx)
Py=sum(vy)
Pz=sum(vz)

Ax=sum(fx)
Ay=sum(fy)
Az=sum(fz)

!Amoso en pantalla t�dalas variables interesantes da simulaci�n:

write (*,*) ' N=',N           !N�mero de part�culas
write (*,*) ' a=',a           !Lado de cada celda
write (*,*) ' k=',k           !N�mero de celdas por eixo
write (*,*) ' rho=',rho       !Densidade de part�culas
write (*,*) ' L=',L           !Lado da caixa de simulaci�n
write (*,*) ' V=',L3          !Volume da caixa de simulaci�n
write (*,*) ' U=',U           !Enerx�a potencial total
write (*,*) ' dU/dV=',phiv    !Primeira derivada da enerx�a potencial
write (*,*) ' d2U/dV2=',phivv !Segunda derivada da enerx�a potencial
write (*,*) ' T=',T           !Enerx�a cin�tica
write (*,*) ' E=',E           !Enerx�a total
write (*,*)
write (*,*) ' P=',Px,Py,Pz    !Momento total
write (*,*) ' F=',Ax,Ay,Az    !Aceleraci�n total

!Defino os formatos a usar no gravado de arquivos.

9000	format(1pe13.6,8(1x,e13.6))
8000	format(1pe19.12,8(1x,e19.12))
8001	format(1pe19.12,4(1x,e19.12))
8002	format(4(i8))

!Gravo en ASCII os vectores de posici�ns, velocidades e aceleraci�ns
!As tres primeiras columnas son x,y,z, as tres seguintes vx,vy,vz e as tres �ltimas ax,ay,az.

open(10,file='Matrix.dat')
	do i=1,N
	  write(10,9000) rx(i),ry(i),rz(i),vx(i),vy(i),vz(i),fx(i),fy(i),fz(i)
    enddo
close(10)

!Gravo en binario os vectores de posici�ns, velocidades e aceleraci�ns.

open(20,file='Matrixbin.dat',form='unformatted')
	write(20) rx,ry,rz,vx,vy,vz,fx,fy,fz
close(20)

!Gravo en ASCII as variables escalares de inter�s:

open(30,file='N�meros.dat')
	write (30,8002) N,k,L,Rc                           !Enteiros
    write (30,8000) ndp,ldp,ldpi,a,rho,L3,Rc2,Rc3,Rc6  !Reais 1
	write (30,8001) U,phiv,phivv,T,E                   !Reais 2
close(30)

!Fin.

stop

end program rede