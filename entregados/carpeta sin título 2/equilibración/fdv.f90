program fdv
!
!O programa fdv calcula a distribución das velocidades dadas polo programa equilibración.
!
!Defino enteiro e real e afirmo que non hai variables implícitas.
!
use prec
implicit none
!
!Declaro variables:
!
integer (kind=int),parameter :: N=500,pasos=500000,gravar=100,dimm=pasos/(2*gravar)
integer (kind=int) :: i,j,mx,my,mz,m2
real (kind=dp), dimension(dimm,N) :: velx,vely,velz,vel2
real (kind=dp) :: vxmax,vymax,vzmax,v2max,p,pinv,dvx,dvy,dvz,dv2
integer (kind=int), allocatable, dimension(:) :: distrx,distry,distrz,distr2
real (kind=dp), allocatable, dimension(:) :: fdvx,fdvy,fdvz,fdv2
!
!Defino a precisión para as velocidades e a súa inversa.
!
p=0.01d00;
pinv=1.d00/p;
!
!Recollo tódalas velocidades e asigno tódalas velocidades dunha mesma compoñente á unha soa variable.
!
open (10,file='uves.dat',form='unformatted')
	do i=1,dimm
	    read (10) velx(i,:),vely(i,:),velz(i,:)
    enddo
close (10)
!
!Calculo tamén o módulo das velocidades.
!
vel2=(velx*velx+vely*vely+velz*velz)**(0.5d00)
!
!Chamo a subrutina big que me calcula o maior elemento dun vector.
!
call big(velx,dimm,N,vxmax)
call big(vely,dimm,N,vymax)
call big(velz,dimm,N,vzmax)
call big(vel2,dimm,N,v2max)
!
!Defino o intervalo no cal calquera velocidade intermedia adopta o mesmo valor.
!
dvx=p*ceiling(pinv*vxmax)
dvy=p*ceiling(pinv*vymax)
dvz=p*ceiling(pinv*vzmax)
dv2=p*ceiling(pinv*v2max)
!
!Dimensiono as variables da distribución
!
allocate (distrx(floor(-pinv*dvx):ceiling(pinv*dvx)),&
		 &distry(floor(-pinv*dvy):ceiling(pinv*dvy)),&
         &distrz(floor(-pinv*dvz):ceiling(pinv*dvz)),&
         &distr2(0:ceiling(pinv*dv2)),&
         &fdvx(floor(-pinv*dvx):ceiling(pinv*dvx)),&
		 &fdvy(floor(-pinv*dvy):ceiling(pinv*dvy)),&
         &fdvz(floor(-pinv*dvz):ceiling(pinv*dvz)),&
         &fdv2(0:ceiling(pinv*dv2)))
!
!Anulo os acumuladores da distribución.
!
distrx=0
distry=0
distrz=0
distr2=0
!
!Comeza o reconto, as variables mx,my,mz,m2 recollen a que intervalo pertence cada valor. Despois
!aumento a conta do intervalo correspondente en unha unidade.
!
do i=1,dimm
  do j=1,N
      mx=nint(pinv*velx(i,j))
      my=nint(pinv*vely(i,j))
      mz=nint(pinv*velz(i,j))
      m2=nint(pinv*vel2(i,j))
      distrx(mx)=distrx(mx)+1
      distry(my)=distry(my)+1
      distrz(mz)=distrz(mz)+1
      distr2(m2)=distr2(m2)+1
  enddo
enddo
!
!Renormalizo dividindo polo número de contas totais.
!
fdvx=dble(distrx)/dble(dimm*N)
fdvy=dble(distry)/dble(dimm*N)
fdvz=dble(distrz)/dble(dimm*N)
fdv2=dble(distr2)/dble(dimm*N)
!
!Gravado de datos:
!
open (20,file='fdvx.dat')
	do i=floor(-pinv*dvx),ceiling(pinv*dvx)
    	write (20,7000) i,fdvx(i)
    enddo
close (20)
!
open (30,file='fdvy.dat')
    do i=floor(-pinv*dvy),ceiling(pinv*dvy)
    	write (30,7000) i,fdvy(i)
    enddo
close (30)
!
open (40,file='fdvz.dat')
    do i=floor(-pinv*dvz),ceiling(pinv*dvz)
    	write (40,7000) i,fdvz(i)
    enddo
close (40)
!
open (50,file='fdv2.dat')
    do i=0,ceiling(pinv*dv2)
    	write (50,7000) i,fdv2(i)
    enddo
close (50)
!
!Defino o formato.
!
7000	format(i8,1x,1pe13.6)
!
!Pecho as variables abertas.
!
deallocate (distrx,distry,distrz,distr2,fdvx,fdvy,fdvz,fdv2)
!
!Fin.
!
stop
end program fdv