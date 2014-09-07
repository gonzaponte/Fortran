program medias
!
!O programa medias calcula os valores medios das constantes do sistema calculadas polo programa
!dinamica. Tamén calcula outras constantes así como o erro de tódalas constantes.
!
!O de sempre, defino precisión e afirmo que non hai variables implícitas.
!
use prec
implicit none
!
!Declaro variables:
!-N:  		  número de partículas do sistema
!-L:		  lado da caixa de simulación
!-i:		  índice para percorrer bucles
!-x:  		  número de datos gravados en dobre precisión
!-xi: 		  inversa de x
!-y:		  factor para incertezas
!-L3inv:	  inversa do volume da caixa de simulación
!-Temp:		  temperatura
!-Presion:	  presión
!-alphaE:	  coeficiente de dilatación a enerxía constante
!-Cv:		  capacidade calorífica a volume constante
!-gamma:	  gamma de Gruneisen
!-kS:		  coeficiente de compresibilidade adiabático
!-###media:	  media da variable ###
!-###2:		  cadrado da media da variable ###
!-s###:		  incerteza da variable ###
!-s###2:	  cadrado da incerteza da variable ###
!-alphaS:	  coecifiente de dilatación isentrópico
!-kT:		  coeficiente de compresibilidade isotérmico
!-Cp:		  capacidade calorífica a presión constante
!-alphaP:	  coeficiente de dilatación isobárico
!-r##:		  erro relativo da variable ###
!
integer(kind=int) :: N,L,i
real(kind=dp),parameter :: x=20.d00,xi=1.d00/x,y=1.d00/(x*(x-1))
real(kind=dp) :: L3inv
real(kind=dp),dimension(nint(x)) :: Temp,Presion,alphaE,Cv,gamma,kS
real(kind=dp) :: Tempmedia,Presionmedia,alphaEmedia,Cvmedia,gammamedia,kSmedia
real(kind=dp) :: Temp2,Cv2,gamma2,kS2
real(kind=dp) :: sTemp,sPresion,salphaE,sCv,sgamma,skS
real(kind=dp) :: sTemp2,sCv2,sgamma2,skS2
real(kind=dp) :: alphaS,kT,Cp,alphaP
real(kind=dp) :: salphaS,skT,sCp,salphaP
real(kind=dp) :: rTemp,rPresion,ralphaE,rCv,rgamma,rkS,ralphaS,rkT,rCp,ralphaP
!
!Defino formatos:
!
8003	format(i8,i8,i8)
8004	format(1pe19.12,5(1x,e19.12))
9000	format('Temperatura:',1x,1pe19.12,1x,'±',e19.12,5x,e9.2,1x,'%')
9001	format('Presión:',5x,1pe19.12,1x,'±',e19.12,5x,e9.2,1x,'%')
9002	format('alphaE:',6x,1pe19.12,1x,'±',e19.12,5x,e9.2,1x,'%')
9003	format('Cv:',10x,1pe19.12,1x,'±',e19.12,5x,e9.2,1x,'%')
9004	format('gamma:',7x,1pe19.12,1x,'±',e19.12,5x,e9.2,1x,'%')
9005	format('kS:',10x,1pe19.12,1x,'±',e19.12,5x,e9.2,1x,'%')
9006	format('alphaS:',6x,1pe19.12,1x,'±',e19.12,5x,e9.2,1x,'%')
9007	format('kT:',10x,1pe19.12,1x,'±',e19.12,5x,e9.2,1x,'%')
9008	format('Cp:',10x,1pe19.12,1x,'±',e19.12,5x,e9.2,1x,'%')
9009	format('alphaP:',6x,1pe19.12,1x,'±',e19.12,5x,e9.2,1x,'%')
!
!Leo o número de partículas e o lado da caixa do arquivo de datos da simulación:
!
open(10,file='Números.dat')
	read(10,8003) N,L,L
close(10)
!
!Calculo o valor de L3inv:
!
L3inv=1.d00/dble(L*L*L)
!
!Leo o ficheiro de constantes do sistema e asigno os valores ás variables do programa
!
open(20,file='Constantes.dat')
	read (20,8004) Temp(1),Presion(1),alphaE(1),Cv(1),gamma(1),kS(1) !Para saltarme o primeiro dato, que da mal resultado
	do i=1,nint(x)
    	read (20,8004) Temp(i),Presion(i),alphaE(i),Cv(i),gamma(i),kS(i)
    enddo
close(20)
!
!Calculo medias dos datos que acabo de ler sumando e dividindo polo número de termos.
!
Tempmedia    = sum(Temp)   * xi
Presionmedia = sum(Presion)* xi
alphaEmedia  = sum(alphaE) * xi
Cvmedia      = sum(Cv)	   * xi
gammamedia   = sum(gamma)  * xi
kSmedia      = sum(kS)     * xi
!
!Calculo os cadrados das medias que serán útiles posteriormente.
!
Temp2    = Tempmedia  * Tempmedia
Cv2      = Cvmedia	  * Cvmedia
gamma2   = gammamedia * gammamedia
kS2      = kSmedia    * kSmedia
!
!Calculo as incertezas de tódalas medias como a raíz cadrada da suma dos cadrados das desviacións á media
!sobre a diferenza entre o cadrado do número de termos e o número de termos.
!
sTemp    = dsqrt(sum((Temp-Tempmedia)*(Temp-Tempmedia))*y)
sPresion = dsqrt(sum((Presion-Presionmedia)*(Presion-Presionmedia))*y)
salphaE  = dsqrt(sum((alphaE-alphaEmedia)*(alphaE-alphaEmedia))*y)
sCv      = dsqrt(sum((Cv-Cvmedia)*(Cv-Cvmedia))*y)
sgamma   = dsqrt(sum((gamma-gammamedia)*(gamma-gammamedia))*y)
skS      = dsqrt(sum((kS-kSmedia)*(kS-kSmedia))*y)
!
!Calculo os cadrados das incertezas das medias que son útiles posteriormente.
!
sTemp2    = sTemp  * sTemp
sCv2      = sCv    * sCv
sgamma2   = sgamma * sgamma
skS2      = skS    * skS
!
!Calculo o resto de constantes do sistema a partir das anteriores.
!
alphaS = -gammamedia*Tempmedia ; alphaS=1.d00/alphaS
kT     = 1.d00/kSmedia-Tempmedia*Cvmedia*gamma2*L3inv ; kT=1.d00/kT
Cp     = Cvmedia*kT/kSmedia
alphaP = Cvmedia*gammamedia*kT*L3inv
!
!Calculo as incertezas destas últimas cantidades por propagación de erros.
!
salphaS = -alphaS*dsqrt(sgamma2/gamma2+sTemp2/Temp2)
skT     = kT*kT*dsqrt(skS2/(kS2*kS2)+L3inv*L3inv*(sTemp2*Cv2*gamma2*gamma2+sCv2*Temp2*gamma2*gamma2+4*sgamma2*gamma2*Temp2*Cv2))
sCp     = Cp*dsqrt(sCv2/Cv2+skT*skT/(kT*kT)+skS2/kS2)
salphaP = alphaP*dsqrt(sCv2/Cv2+sgamma2/gamma2+skT*skT/(kT*kT))
!
!Calculo os erros relativos en tanto por cen como dúas veces o erro sobre o valor por cen.
!
rTemp	 = 200 * sTemp    / Tempmedia
rPresion = 200 * sPresion / Presionmedia
ralphaE	 = 200 * salphaE  / alphaEmedia
rCv		 = 200 * sCv      / Cvmedia
rgamma	 = 200 * sgamma   / gammamedia
rkS		 = 200 * skS      / kSmedia
ralphaS	 = 200 * salphaS  / alphaS
rkT  	 = 200 * skT      / kT
rCp		 = 200 * sCp      / Cp
ralphaP  = 200 * salphaP  / alphaP
!
!Gardo os valores  nun arquivo:
!
open(30,file='Coef_termo.dat')
	write(30,9000) Tempmedia,2*sTemp,rTemp
    write(30,9001) Presionmedia,2*sPresion,rPresion
    write(30,9002) alphaEmedia,2*salphaE,(-ralphaE)
    write(30,9003) Cvmedia,2*sCv,rCv
    write(30,9004) gammamedia,2*sgamma,rgamma
    write(30,9005) kSmedia,2*skS,rkS
    write(30,9006) alphaS,2*salphaS,(-ralphaS)
    write(30,9007) kT,2*skT,rkT
    write(30,9008) Cp,2*sCp,rCp
    write(30,9009) alphaP,2*salphaP,ralphaP
close(30)
!
!E noutro que non leva os nomes das constantes:
!
9010	format(3(1pe19.12))
!
open(40,file='coeficientes.dat')
	write(40,9010) Tempmedia,2*sTemp,rTemp
    write(40,9010) Presionmedia,2*sPresion,rPresion
    write(40,9010) alphaEmedia,2*salphaE,(-ralphaE)
    write(40,9010) Cvmedia,2*sCv,rCv
    write(40,9010) gammamedia,2*sgamma,rgamma
    write(40,9010) kSmedia,2*skS,rkS
    write(40,9010) alphaS,2*salphaS,(-ralphaS)
    write(40,9010) kT,2*skT,rkT
    write(40,9010) Cp,2*sCp,rCp
    write(40,9010) alphaP,2*salphaP,ralphaP
close(40)
!
!Fin.
!
stop
end program medias