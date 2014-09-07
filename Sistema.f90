program sistema
!
!Este programa o �nico que fai � reunir t�dolos datos relevantes dun sistema e po�elos
!nun mesmo arquivo, para que quede m�is bonito.
!
!Defino as precisi�ns e afirmo que non hai variables impl�citas.
!
integer,parameter :: int = SELECTED_INT_KIND(9)       !Define un enteiro
integer,parameter ::  dp = SELECTED_REAL_KIND(15,307) !Define un real
!
implicit none
!
!Declaro variables:
!
!
!
integer(kind=int) :: N,k,L,Rc
real(kind=dp) E,Eeq,U,Ueq,phiv,phivv,T,Teq,rho,a,Tempmedia,sTemp,rTemp,Presionmedia,sPresion,rPresion,&
&alphaEmedia,salphaE,ralphaE,Cvmedia,sCv,rCv,gammamedia,sgamma,rgamma,kSmedia,skS,rkS,alphaS,salphaS,&
&ralphaS,kT,skT,rkT,Cp,sCp,rCp,alphaP,salphaP,ralphaP
!
!Defino formatos a usar:
7000	format(4(i8))
7001	format(1pe19.12)
7002	format(1pe19.12,4(1x,e19.12))
8000	format(3(1pe19.12))
!
!Collo datos do primeiro programa, obvio alg�ns, qu�dome so con N, k, L, Rc, a, rho, U, phiv, phivv, T e E
!
open(10,file='N�meros.dat')
	read (10,7000) N,k,L,Rc                           !Enteiros
    read (10,7001) a								  !Paso temporal, que obvio
    read (10,7002) a,a,a,a,rho						  !Reais 1
	read (10,7002) U,phiv,phivv,T,E                   !Reais 2
close(10)
!
!Collo datos do terceito programa, collo so os primeiros datos porque supostamente xa est� en equilibrio.
!
open(20,file='Resultados.dat')
	read (20,7000) Eeq,Eeq,Teq,Ueq !Aqu� estaba gravado o tempo pero obvieino.
close(20)
!
!Collo tam�n as propiedades termodin�micas do sistema:
!
open(30,file='coeficientes.dat')
	read(30,8000) Tempmedia,sTemp,rTemp
    read(30,8000) Presionmedia,sPresion,rPresion
    read(30,8000) alphaEmedia,salphaE,ralphaE
    read(30,8000) Cvmedia,sCv,rCv
    read(30,8000) gammamedia,sgamma,rgamma
    read(30,8000) kSmedia,skS,rkS
    read(30,8000) alphaS,salphaS,ralphaS
    read(30,8000) kT,skT,rkT
    read(30,8000) Cp,sCp,rCp
    read(30,8000) alphaP,salphaP,ralphaP
close(30)
!
!O �ltimo programa da o coeficiente de difusi�n
open(40,file='Coef_difusion.dat')
	write(40,7001) D
close(40)
!
!Agora gardo todo nun mesmo arquivo:
!
9000	format('N=      ',i8 ,5x,'N�mero de part�culas')
9001	format('V=      ',f.0,5x,'Volume do sistema')
9002	format('E=      ',f.0,5x,'Enerx�a total do sistema',5x,'Eeq=      ',f.1,5x,'Enerx�a total no equilibrio')
9003	format('k=      ',i8 ,5x,'N�mero de celdas por lado')
9004	format('L=      ',i8 ,5x,'Lado da caixa de simulaci�n')
9005	format('Rc=     ',i8 ,5x,'Raio de corte do potencial')
9006	format('a=      ',f.2,5x,'Lado de cada celda')
9007	format('rho=    ',f.2,5x,'Densidade')
9008	format('U=      ',f.1,5x,'Enerx�a potencial',5x,'Ueq=      ',f.1,5x,'Enerx�a potencial no equilibrio')
9009	format('T=      ',f.2,5x,'Enerx�a cin�tica ',5x,'Teq=      ',f.1,5x,'Enerx�a potencial no equilibrio')
9010	format('phi_v=  ',f.2,5x,'Primeira derivada da enerx�a potencial respecto ao volume')
9011	format('phi_vv= ',f.2,5x,'Segunda  derivada da enerx�a potencial respecto ao volume')
9012	format('T=      ',f.4,
!
open(50,file='Sistema.dat')
	write(50,) N
!
stop
end program sistema