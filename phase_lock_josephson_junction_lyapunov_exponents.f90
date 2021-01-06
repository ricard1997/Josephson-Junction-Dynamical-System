program josephson_junction_simulations
!implicit none
integer,parameter::n=12,nstep=310000,ndim=3
real*8::htry,t,pi
real*8,dimension(n)::dvardt,var,varscal,varout
real*8,dimension(ndim,ndim)::Mat
real*8,dimension(ndim)::modulos,lya
real*8,dimension(nstep,2)::func
real*8::dx,i0,a,b,s,Per,si,ra,ran
integer::i,ncy,nprom,reg,clock
		
reg=400!promedios calculados
pi=dacos(-1.d0)  
dx=0.4d0/reg!grilla de i_0
frec=2d0
Per=2.d0*pi/frec
ncy=320
nprom=256    !durante cuentos periodos se hace la integracion
htry= Per/float(ncy)

open(file='exponentesdelpromediosubida.dat',unit=14,status='unknown')
open(file='exponentesdelpromediobajada.dat',unit=15,status='unknown')
open(file='curvaiv.dat',unit=16,status='unknown')
open(file='curvaivbajada.dat',unit=17,status='unknown')
open(file='variacionesdecondiciones1.dat',unit=18,status='unknown')


!Generating random numbers
call system_clock(clock)
call random_seed(clock)
call random_number(ra)
call random_number(ran)
i0=0.6d0-dx
var(1)=pi*ra
var(2)=pi*ran
var(3)=pi*ra
t=0.d0	
si=1.d0




do i=1,2*reg
	t=0.0d0 

	if (i.gt.reg) si=-1.d0
	!if (i0.lt.0.85d0) si=1.d0
	i0=i0+si*dx
	print*, i0,i
	
	call exponentesprom(var,dvardt,n,t,htry,i0,nstep,func,lya)
	a=func(1,2)
	b=func(ncy*nprom,2)
	call integ(func(:,1),a,b,s,ncy*nprom)
	if(i.lt.reg) write(16,*) i0, (s/(b-a))/frec
	if(i.ge.reg) write(17,*) i0, (s/(b-a))/frec
	if(i.lt.reg) write(14,*) i0, lya
	if(i.ge.reg) write(15,*) i0, lya
	if(int(100*i0).gt.85.and.int(100*i0).lt.100) then
		!print*,i0,var(1),var(2),var(3), 'valor'
		!print*,lya, (s/(b-a))/frec
		write(18,*) i0,var(1),var(2),var(3)
	endif
	
enddo
  
close(14)
close(16)

contains
!________________subrutina promedio_____________________

subroutine integ(func,a,b,s,nit)
integer::n
real*8::a,b,s,al
integer::nit,j
real*8:: del,suma,tnm,x
real*8,dimension(nit)::func
	tnm=nit
	del=(b-a)/nit
	x=a+0.5*del
	suma=0.
	do j=1,nit
		suma=suma+func(j)
		x=x+del
	enddo 
	s=((b-a)*suma/tnm)
end subroutine integ









!________________________-subrutina para la integracion-___________


subroutine exponentesprom(var,dvardt,n,t,htry,i0,nstep,func,lya)
implicit none
integer,parameter::ndim=3
integer::nstep,n
real*8::htry,t,pi,dx
real*8,dimension(n)::dvardt,var,varscal,varout
real*8,dimension(ndim,ndim)::Mat
real*8,dimension(ndim)::modulos,lya
real*8::i0,may
real*8,dimension(nstep,2)::func
integer::j,k,co,con,nort,ncy,nper,tau,nprom
pi=dacos(-1.d0)
ncy=320   !cuantas divisiones hay en unperiodo
nper=256  !rechaza nper cantidad de pediodos
nprom=256
tau=1

!Cond iniciales de la matriz M
!La matriz M tiene vectores columna como vectores iniciales
Mat=0.d0
do k=1,ndim
	Mat(k,k)=1.d0
enddo
!Se evolucionan como variables cada espacio de la matriz
do k=1,ndim
	do co=1,ndim
		var(k*ndim+co)=Mat(k,co)
	enddo
enddo
nort=1
lya=0.d0
co=1
con=0

!Algorithm for computing Lyapunov exponents
do while (co.eq.1)
	con=con+1
	call derivs(t,var,dvardt,i0)
	call rkck(var,dvardt,n,t,htry,varout,i0)
	t=t+htry
	var=varout
	if(var(1).gt.pi) var(1)=var(1)-2*pi
	if(var(1).lt.-pi) var(1)=var(1)+2*pi
	if(var(3).gt.2*pi) var(3)=var(3)-2*pi
	if(con.gt.nper*ncy) then
		if(con-nper*ncy.gt.nprom*ncy) co=0
		func(con-nper*ncy,1)=var(2)
		func(con-nper*ncy,2)=t
	endif

	!Ortonormalizacion de Gramschimith
	if(int(t).eq.nort*tau) then
		!Lleno la matriz Mat
		do k=1,ndim	
			do j=1,ndim
				Mat(k,j)=var(ndim*k+j)
			enddo
		enddo
		call ortogo(Mat,ndim,modulos)
		do j=1,ndim
			do k=1,ndim
				var(j*ndim+k)=Mat(j,k)
			enddo
		enddo
		do j=1,ndim
			lya(j)=lya(j)+dlog(modulos(j))
		enddo
		if(con.gt.nstep) co=0
		nort=nort+1
			
	endif
enddo
do j=1,ndim
	lya(j)=(1.d0/(nort*tau))*lya(j)
enddo
end subroutine exponentesprom





!_________________________________________rkck_________________________________


subroutine rkck(var,dvardt,n,t,dt,varout,i0)

!Runge-Kutta Cash-Karp (fifth order)
implicit none
integer::n
real*8:: dt,t
real*8,dimension(n)::dvardt,var,varout
integer,parameter::NMAX=50
integer::i
real*8,dimension(NMAX)::ak2,ak3,ak4,ak5,ak6,vartemp
real*8,parameter::CA2=.2d0,CA3=.3d0,CA4=.6d0,CA5=1.d0,CA6=.875d0,AB21=.2d0,AB31=3.d0/40.d0
real*8,parameter::AB32=9.d0/40.d0,AB41=.3,AB42=-.9d0,AB43=1.2d0,AB51=-11.d0/54.d0,AB52=2.5d0
real*8,parameter::AB53=-70.d0/27.d0,AB54=35.d0/27.d0,AB61=1631.d0/55296.d0,AB62=175.d0/512.d0
real*8,parameter::AB63=575.d0/13824.d0,AB64=44275.d0/110592.d0,AB65=253.d0/4096.d0
real*8,parameter::BC1=37.d0/378.d0,BC3=250.d0/621.d0,BC4=125.d0/594.d0,BC6=512.d0/1771.d0
real*8,parameter::DBC1=BC1-2825.d0/27648.d0,DBC3=BC3-18575.d0/48384.d0
real*8,parameter::DBC4=BC4-13525.d0/55296.d0,DBC5=-277.d0/14336.d0,DBC6=BC6-.25d0
real*8::i0
do  i=1,n        !primer paso
	vartemp(i)=var(i)+AB21*dt*dvardt(i)
enddo 
call derivs(t+CA2*dt,vartemp,ak2,i0) !aqui se obtiene k2
do i=1,n        !segundo paso
	vartemp(i)=var(i)+dt*(AB31*dvardt(i)+AB32*ak2(i))
enddo 
call derivs(t+CA3*dt,vartemp,ak3,i0)  !aqui se obtiene k3
do i=1,n        !tercer paso
	vartemp(i)=var(i)+dt*(AB41*dvardt(i)+AB42*ak2(i)+AB43*ak3(i))
enddo
call derivs(t+CA4*dt,vartemp,ak4,i0)   !aqui se obtiene k4
do i=1,n        !cuarto paso
	vartemp(i)=var(i)+dt*(AB51*dvardt(i)+AB52*ak2(i)+AB53*ak3(i)+AB54*ak4(i))
enddo
call derivs(t+CA5*dt,vartemp,ak5,i0)   !Aqui se obtiene k5
do i=1,n        !quinto paso
	vartemp(i)=var(i)+dt*(AB61*dvardt(i)+AB62*ak2(i)+AB63*ak3(i)+AB64*ak4(i)+AB65*ak5(i))
enddo
call derivs(t+CA6*dt,vartemp,ak6,i0)  !aqui se obtiene k6
do i=1,n  !aqui se obtiene el vector de salida
	varout(i)=var(i)+dt*(BC1*dvardt(i)+BC3*ak3(i)+BC4*ak4(i)+BC6*ak6(i))
enddo 
end subroutine rkck

!_______________________________________derivs__________________________________________


subroutine derivs(t,var,dvardt,i0)
real*8::t,var(*),dvardt(*)
real*8::i0,frec,i1,alpha
integer::i,con,conja
integer,parameter::ndim=3
real,dimension(n,n)::Jac

!Josephson junction parameters
frec=2.d0
i1=80
alpha=0.8d0


!Dynamical system - Josephson Junction
dvardt(1)=var(2)
dvardt(2)=i0+i1*dsin(var(3))-alpha*var(2)-dsin(var(1))
dvardt(3)=frec  

!Jacobian
Jac(1,1)=0.d0
Jac(1,2)=1.d0
Jac(1,3)=0.d0
Jac(2,1)=-dcos(var(1))
Jac(2,2)=-alpha
Jac(2,3)=-i1*dcos(var(3))
Jac(3,1)=0.d0
Jac(3,2)=0.d0
Jac(3,3)=0.d0
do j=1,ndim*ndim
	dvardt(ndim+j)=0.d0
enddo  
con=0
conja=1

!Completing equations to compute Lyapunov exponents
do i=1,ndim*ndim
	con=con+1      
	do j=1,ndim
		if(con.gt.ndim) conja=conja+1      
		if(con.gt.ndim) con=1
		dvardt(ndim+i)=dvardt(ndim+i)+Jac(conja,j)*var((ndim+con)+(j-1)*ndim)  		
  	enddo
enddo
end subroutine derivs




!__________________________________Subrutinas necesarias para la ortonormalizacion_____________________

subroutine mult(mat1,mat2,respu,n)
	integer::i,j,k,n
	real*8,dimension(n,n)::mat1,mat2,respu
	real*8::suma
	
	do i=1,n
		do j=1,n
		        suma=0.d0
			do k=1,n
			        
				suma=suma+mat1(i,k)*mat2(k,j)
				!if (i.eq.1.and.j.eq.1) print*,suma
			enddo
			respu(i,j)=suma
		enddo
	enddo
	
end subroutine

subroutine ortogo(mat,n,modulos)
	integer::n,i,j,k
	real*8,dimension(n,n)::mat
	real*8,dimension(n)::modulos,vec
	real*8::modu
	call norm(mat(:,1),n,modu)
	!print*,'matriz inical-------------------------------'
	!do i=1,n
		!print*,mat(:,i)
        !enddo
	!print*,'matriz inical-------------------------------'
	
	
	modulos(1)=modu
	!print*,mat(:,1),'vec 1'
	do i=2,n
		vec=0.d0
		do j=1,i-1	 
			vec=vec+prod(mat(:,i),mat(:,j),n)*mat(:,j)
			!print*,prod(mat(:,i),mat(:,j),n),'::',i,j
			!print*,vec,'::',i,j
		enddo
	!	print*,vec
		mat(:,i)=mat(:,i)-vec
		!print*,vec,'::',i,j,'matriz'
		call norm(mat(:,i),n,modu)
		!print*,mat(:,i),'vec',i
		modulos(i)=modu
	enddo

end subroutine


double precision function prod(v1,v2,n)
	integer::n,i
	real*8,dimension(n)::v1,v2
	real*8::suma
	suma=0.d0
	do i=1,n
		suma=suma+v1(i)*v2(i)
	enddo
	prod=suma

end function

subroutine norm(v,n,modu)
	integer::n,i
	real*8,dimension(n)::v
	real*8::modu
	modu=prod(v,v,n)
	modu=dsqrt(modu)
	v=(1.d0/modu)*v
end subroutine


end program
