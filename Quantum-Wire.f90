program TMDCCircular
parameter (nr_steps = 10000)
implicit real*8 (a,b,d-h,o-z)
implicit complex*16 (c)
implicit integer (i-n)
real*8 muj, mujp1
DOUBLE PRECISION r(nr_steps),D(nr_steps),E(nr_steps),vhj(nr_steps),potencial(nr_steps),rstep(3)
DOUBLE PRECISION,ALLOCATABLE:: ZLPK(:,:)
INTEGER IFAIL(NR_STEPS),IWORK(5*NR_STEPS)
DOUBLE PRECISION    W(NR_STEPS),WORK(5*NR_STEPS),WF(NR_STEPS),PSI(nr_steps,10)
 CHARACTER JOBZ,RANGE

common / limits1 / rstep,r_max,a0
	
ALLOCATE (ZLPK(NR_STEPS,NR_STEPS)) !WINDOWS

open(unit = 11, file = 'potential.dat')
open(unit = 12, file = 'waves.dat')
open(unit = 13, file = 'Energies.dat')

pi = 4.D0*atan(1.D0)
r_min=1d-8
rstep(1)=1.0d-4 !short range
rstep(2)=5.0d-2 !middle range
rstep(3)=9.0d-2 !long range

am1 = 0.067D0 !Effective mass

a0=0.52917721D0 / am1         
ry=13.6056923D+3 * am1 

l = 0 !angular momentum
l2 = l*l

RR = 40.D0/a0 !Quantum wire radius
V0 = 245.D0/Ry !External potential, or band offset

rdim = r_min
do i = 1,nr_steps
 if(i.le.1000)then
  rdim = rdim + rstep(1) 
 elseif(i.gt.1000.and.i.le.5000)then
  rdim = rdim + rstep(2) 
 elseif(i.gt.5000)then
  rdim = rdim + rstep(3) 
 endif
 r(i) = rdim/a0

 if(r(i).le.RR)then
  potencial(i) = 0.D0
 else
  potencial(i) = V0
 endif

 write(11,*)r(i)*a0,potencial(i)*ry !writting potential into a file
enddo

r_max = r(nr_steps)

!-------- Defining a tri-diagonal matrix -------------------------------
do i=1,nr_steps
	   
 if(i.ne.1)then
  rojm1=r(i-1)
 else
  rojm1=r_min
 endif
 roj=r(i)
 if(i.ne.nr_steps)then
  rojp1=r(i+1)
 else
  rojp1=r(i) 
 endif
        
 rojp12=(rojp1+roj)/2.
 if(i.eq.NR_STEPS)rojp12=roj
 rojm12=(roj+rojm1)/2.
 if(i.eq.1)rojp12=roj
	   
 hj=dsqrt((rojp12*rojp12-rojm12*rojm12)/2.)
	   
 vhj(i)=hj
 
 if(i.le.nr_steps-2)then 
  rojp2=r(i+2)
 else
  rojp2=r(i)+2*rstep(3)
 endif
	 
 rojp12p1=(rojp2+rojp1)/2.
 rojm12p1=(rojp1+roj)/2.

 hjp1=dsqrt((rojp12p1*rojp12p1-rojm12p1*rojm12p1)/2.)
	  
	   
 muj=(roj+rojm1)/(2.*(roj-rojm1))
 if(i.eq.1)muj=0
 if(i.eq.NR_STEPS)muj=0
   
 if(i.ne.NR_STEPS)mujp1=(rojp1+roj)/(2.*(rojp1-roj))
 if(i.eq.NR_STEPS-1)mujp1=0
 if(i.eq.NR_STEPS)mujp1=0

 D(i)= l2/(r(i)**2) + potencial(i) + (mujp1+muj)/(hj*hj)
 E(i)=-mujp1/(hj*hjp1)
 	
enddo
!----------------------------------------------------------------------


!       ------- LAPACK ----------
ABSTOL=0.0d0  ! TolerÃ¢ncia do erro -> deixar 0.0d0
IL=1       ! IL e IU sÃ£o as posiÃ§Ãµes dos autoestados
IU=10       !maximum calculated level
VL=0.d0    ! VL e VU sÃ£o os limites do intervalo numÃ©rico
VU=0.d0
RANGE='I'  ! =I calcula do IL-Ã©simo ao IU-Ã©simo estado
	           ! =V calcula autovalores dentro do intervalo (VL,VU]
!       -------------------------

JOBZ='V' !calcula funÃ§Ã£o de onda e energia

!       Input:  Vetores D e E
 call DSTEVX(JOBZ,RANGE,NR_STEPS,D,E,VL,VU,IL,IU,ABSTOL,IM,W,ZLPK,NR_STEPS,WORK,IWORK,IFAIL,INFO)	   
!       Output: Vetor W dÃ¡ os autovalores, ZLPK dÃ¡ os autovetores

do nivel = 1,IU
 do i=1,nr_steps
  WF(i)=ZLPK(i,nivel)/vhj(i)
 enddo

 call Norma(nr_steps,r,WF,S)

 do i=1,nr_steps
  PSI(i,nivel)=WF(i)/S  !funcao normalizada
 enddo
 
 if(nivel.eq.1)then
  do i=1,nr_steps
   write(12,*)r(i)*a0,PSI(i,nivel)
  enddo
 endif

enddo

20         format(3e16.6)
 
do nivel=1,IU
 write(13,*)nivel, W(nivel)*Ry !deixa a energia em meV
 write(*,*)nivel, W(nivel)*Ry !deixa a energia em meV    
 do i=1,nr_steps
  WF(i)=PSI(i,nivel)
 enddo
enddo

!write(*,*)'Results should be', Ry*(2.4048D0/RR)**2, Ry*(5.5201D0/RR)**2, Ry*(8.6537D0/RR)**2
!write(*,*)'Results (l = 1) should be', Ry*(3.8317D0/RR)**2, Ry*(7.0156D0/RR)**2, Ry*(10.1735D0/RR)**2

end program

!=======================================================================
subroutine Norma(nx,r,W,S) 
!                                            Wave function normalization
!=======================================================================
real*8 a,b,h,S,S0,S1,pi
integer nx,i
real*8 rstep(3),r_max,r(nx)
real*8 W(nx)    
	
common / limits1 / rstep,r_max,a0


pi=4.*atan(1.)
S=0.0
S1=0.0 
S0=0.0
do i=1,nx
   if(i.le.1000)then
      S0=rstep(1)*W(i)*W(i)*r(i)
   elseif(i.gt.1000.and.i.le.5000)then
      S0=rstep(2)*W(i)*W(i)*r(i)
   elseif(i.gt.5000)then
      S0=rstep(3)*W(i)*W(i)*r(i)
   endif
   
   S1=S1+S0
   
enddo
	
S=dsqrt(2.*pi*S1)
      
return
end



