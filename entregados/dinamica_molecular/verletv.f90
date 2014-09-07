subroutine verletv(N,ldp,ldpi,L33,Rc2,Rc6,cor,rx,ry,rz,vx,vy,vz,fx,fy,fz,dt,U,phiv,phivv)
!
use prec
implicit none
!
!Defino as variables de entrada e saída:
!-rx,ry,rz: compoñentes do vector posición de cada partícula
!-vx,vy,vz: compoñentes do vector velocidade de cada partícula
!-fx,fy,fz: compoñentes do vector aceleración de cada partícula
!
real (kind=dp),dimension(N),intent(inout) :: rx,ry,rz,vx,vy,vz,fx,fy,fz
!
!Agora as variables só de entrada:
!-N: 	número de partículas da simulación
!-ldp: 	lado do cubo de simulación en dobre precisión
!-ldpi: inverso de ldp
!-L33: 	inverso do triple do volume
!-Rc2: 	raio de corte ao cadrado
!-Rc6: 	inverso do raio de corte á sexta potencia
!-cor: 	factor para correccións á enerxía potencial e derivadas.
!-dt:	paso de tempo da simulación
!
integer (kind=int),intent(in) :: N
real (kind=dp),intent(in) :: ldp,ldpi,L33,Rc2,Rc6,cor,dt
!
!E as variables só de saída:
!-U: 	 enerxía potencial
!-phiv:  primeira derivada da enerxía potencial
!-phivv: segunda derivada da enerxía potencial
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
!Calculo as posicións.
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
!Chamo a subrutina da enerxía potencial que me da as novas aceleracións, a enerxía potencial e as
!súas derivadas.
!
call potencialLJ(N,ldp,ldpi,L33,Rc2,Rc6,rx,ry,rz,cor,U,phiv,phivv,fx,fy,fz)
!
!Uso estas novas aceleracións para obter o resto das velocidades.
!
vx=vx+fx*t12
vy=vy+fy*t12
vz=vz+fz*t12
!
!Fin
!
return
end subroutine verletv