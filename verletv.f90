subroutine verletv(N,ldp,ldpi,L33,Rc2,Rc6,cor,rx,ry,rz,vx,vy,vz,fx,fy,fz,dt,U,phiv,phivv)
!
use prec
implicit none
!
!Defino as variables de entrada e sa�da:
!-rx,ry,rz: compo�entes do vector posici�n de cada part�cula
!-vx,vy,vz: compo�entes do vector velocidade de cada part�cula
!-fx,fy,fz: compo�entes do vector aceleraci�n de cada part�cula
!
real (kind=dp),dimension(N),intent(inout) :: rx,ry,rz,vx,vy,vz,fx,fy,fz
!
!Agora as variables s� de entrada:
!-N: 	n�mero de part�culas da simulaci�n
!-ldp: 	lado do cubo de simulaci�n en dobre precisi�n
!-ldpi: inverso de ldp
!-L33: 	inverso do triple do volume
!-Rc2: 	raio de corte ao cadrado
!-Rc6: 	inverso do raio de corte � sexta potencia
!-cor: 	factor para correcci�ns � enerx�a potencial e derivadas.
!-dt:	paso de tempo da simulaci�n
!
integer (kind=int),intent(in) :: N
real (kind=dp),intent(in) :: ldp,ldpi,L33,Rc2,Rc6,cor,dt
!
!E as variables s� de sa�da:
!-U: 	 enerx�a potencial
!-phiv:  primeira derivada da enerx�a potencial
!-phivv: segunda derivada da enerx�a potencial
!
real (kind=dp),intent(out) :: U,phiv,phivv
!
!Agora defino variables auxiliares:
!-t22: 0.5*dt*dt
!-t12: 0.5*dt
!
real (kind=dp) :: t22,t12
!
t12=0.5d00*dt
t22=t12*dt
!
!Calculo as posici�ns.
!
rx=rx+vx*dt+fx*t22
ry=ry+vy*dt+fy*t22
rz=rz+vz*dt+fz*t22
!
!Calculo a primeira parte da velocidade.
!
vx=vx+fx*t12
vy=vy+fy*t12
vz=vz+fz*t12
!
!Chamo a subrutina da enerx�a potencial que me da as novas aceleraci�ns, a enerx�a potencial e as
!s�as derivadas.
!
call potencialLJ(N,ldp,ldpi,L33,Rc2,Rc6,rx,ry,rz,cor,U,phiv,phivv,fx,fy,fz)
!
!Uso estas novas aceleraci�ns para obter o resto das velocidades.
!
vx=vx+fx*t12
vy=vy+fy*t12
vz=vz+fz*t12
!
!Fin
!
return
end subroutine verletv