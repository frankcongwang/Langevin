module init
  implicit none
  real(8), parameter :: kT = 1d0                       !needs to be modified
  real(8), parameter :: h = 0.002d0                      !needs to be modified
  real(8), parameter :: w = 1d0
  real(8), parameter :: m = 1d0
  real(8), parameter :: Mex = 50d0
  real(8), parameter :: pressure_ex = 1d0
  integer, parameter :: eqstep=1d4/h
  integer, parameter :: tsstep=1d5/h
  integer, parameter :: sample=2
!  integer, parameter :: Mtb=4   !or 6 the length of the NHC
!  real(8), parameter :: mQ(Mtb) = 1d0
  real(8), parameter :: pi=3.14159265358979d0
  real(8) :: gamma = 1000.0d0
  real(8) :: gammaV = 5000.0d0
end module init

subroutine calForcex(fn, x, V)
  use init
  implicit none
  real(8) :: fn, x, V
  fn = -m*w**2*V/2/pi*sin(2*pi*x/V)              !needs to be modified
end subroutine calForcex

subroutine calForceV(fv, x, V)
  use init
  implicit none
  real(8) :: fv, x, V
  fv = -m*w*w/2.0d0/pi*((1.0d0-cos(2*pi*x/V))*V/pi-x*sin(2*pi*x/V))     !needs to be modified
end subroutine calForceV

subroutine calPressure(Pins,Vol,enek,Fn,rn,Fv)
  use init
  implicit none
  real(8) :: Pins, Vol, enek, Fn, rn, Fv
!  real(8) :: 
  Pins=1d0/Vol*2.0d0*enek+(rn-(floor(rn/Vol+0.5d0)*Vol))*Fn/Vol+Fv
end subroutine calPressure

subroutine kEnergy(Ek,pn)
        use init
        implicit none
        real(8) :: Ek, pn
        Ek=pn*pn/2.0d0/m
end subroutine

subroutine restrictCoord(qn,qv)
        use init
        implicit none
        real(8) :: qn,qv
        integer :: i,j
                if (qn>=qv) qn=qn-qv
                if (qn<0.0d0) qn=qn+qv
end subroutine

subroutine MolecularDynamics(qn,pn,qv,pv,fn,fv,Pressure,enek)
        use init
        use random
        implicit none
        real(8) :: qn,pn,qv,pv,fn,fv,Pressure,enek
        real(8) :: qv0,fp0,fn0,a1,b1,a2,b2,c1,c2,d1,d2,beta1,beta2
        real(8) :: rand
        c1=gammaV*h/2/Mex
        a1=(1-c1)/(1+c1)
        b1=1/(1+c1)
        c2=gamma*h/2/m
        a2=(1-c2)/(1+c2)
        b2=1/(1+c2)
        d1=exp(-h*gammaV/Mex)
        d2=exp(-h*gamma/m)
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call kEnergy(enek,pn)
          call calPressure(Pressure,qv,enek,fn,qn,fv)
          fp0=Pressure-Pressure_ex
          call random_normal(rand)
          beta1=sqrt(2*gammaV*kT*h)*rand
          qv0=qv
          qv=qv+h*b1*pv/Mex+0.5d0*h*h*b1/Mex*(Pressure-Pressure_ex)+0.5d0*h*b1/Mex*beta1
          call calForcex(fn,qn,qv)
          fn0=fn
          call random_normal(rand)
          beta2=sqrt(2*gamma*kT*h)*rand
          qn=qn*qv/qv0+h*2*qv/(qv+qv0)*b2*(pn/m+0.5d0*h/m*fn+0.5d0/m*beta2)
          call restrictCoord(qn,qv)
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call kEnergy(enek,pn)
          call calPressure(Pressure,qv,enek,fn,qn,fv)
          pv=a1*pv+0.5d0*h*(a1*fp0+Pressure-Pressure_ex)+b1*beta1
          pn=a2*pn+0.5d0*h*(a2*fn0+fn)+b2*beta2
!write(*,*) qn,pn,qv,pv
!read(*,*)
          call kEnergy(enek,pn)
          call calPressure(Pressure,qv,enek,fn,qn,fv)
  end subroutine

!  subroutine tbStat(qn,pn,qv,pv,fn,fv,Pressure,enek)
!        use random
!        use init
!        implicit none
!        real(8) :: qn,pn,qv,pv,fn,fv,Pressure,enek
!        real(8) :: rand,a,b
!        integer :: i
!        call kEnergy(enek,pn)
!        call calForcex(fn,qn,qv)
!        call calForceV(fv,qn,qv)
!        call calPressure(Pressure,qv,enek,fn,qn,fv)
!          pv = pv + 0.5d0*h*(Pressure-Pressure_ex)
!!write(*,*) "**",qn,pn,qv,pv,fn
!!read(*,*)
!          qn = qn + 0.5d0*h*pn/m
!          call restrictCoord(qn,qv)
!          qn = qn*exp(0.5d0*h*pv/Mex/qv)
!          call restrictCoord(qn,qv)
!          qv = qv + 0.5d0*h*pv/Mex
!        a=exp(-h*gammaV/Mex)
!       call random_normal(rand)
!        pv=pv*a+sqrt(m*kT*(1-a*a))*rand
!        b=exp(-h*gamma/m)
!       call random_normal(rand)
!        pn=pn*b+sqrt(Mex*kT*(1-b*b))*rand
!          qv = qv + 0.5d0*h*pv/Mex
!          qn = qn*exp(0.5d0*h*pv/Mex/qv)
!          call restrictCoord(qn,qv)
!          qn = qn + 0.5d0*h*pn/m
!          call restrictCoord(qn,qv)
!        call kEnergy(enek,pn)
!        call calForcex(fn,qn,qv)
!        call calForceV(fv,qn,qv)
!        call calPressure(Pressure,qv,enek,fn,qn,fv)
!          pv = pv + 0.5d0*h*(Pressure-Pressure_ex)
!end subroutine




program main
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qv, pv, fn, fv, Pressure
!  real(8) :: qt(Mtb), pt(Mtb)
  integer :: i, j
  real*8 :: eptmp, ektmp, ep(sample), ek(sample),pres(sample),ep_ave,ek_ave,pres_ave,ep_std,ek_std,pres_std
  real*8 :: cor_ep(tsstep/2),cor_ek(tsstep/2)
  real*8 :: t, cortimep, cortimek,cortimep_std,cortimek_std
  integer :: n, vcount(100)
  character(30) :: c
  vcount=0
  open(22,file='resultJensen.maindat')
  ep(:)=0
  ek(:)=0
    write(*,*) 'gamma=',gamma, 'dt=', h
pres=0.0d0
    do j=1, sample
       write(c,'(I2)') j 
       write(*,*) 'Sample=', j
       open(33,file=trim('traj_'//adjustl(c)))
       open(999,file=trim('note_'//adjustl(c)))
       call random_normal(rand)
       pn = sqrt(m*kT/2)*rand
       call random_normal(rand)
       pv = 0.2d0*rand
       call random_number(rand)
       qn = rand
       call random_normal(rand)
       qv = 3.0d0 !1.5d0+0.05d0*rand
!       do i=1,Mtb
!         call random_normal(rand)
!         qt(i)=0.2d0*rand
!         call random_normal(rand)
!         pt(i)=0.2d0*rand
!       end do
!write(*,*) qn,pn,qv,pv
!read(*,*)
    !   write(*,*) 'sample=', j
    !      pause
          call calForcex(fn,qn,qv)
       do i=1, eqstep
         call MolecularDynamics(qn,pn,qv,pv,fn,fv,Pressure,ektmp)
         eptmp = m*w**2*qv**2/4/pi**2*(1-cos(2*pi*qn/qv))       !needs to be modified
!         ep(j)  = ep(j) + eptmp/tsstep
!         ek(j)  = ek(j) + ektmp/tsstep
!         pres(j)=pres(j)+ Pressure/tsstep
         write(999,*) i,pn,qn,pv,qv
!       p  pause
!                    write(33,'(I16,F16.8,F16.8,F16.8,F16.8)') i,eptmp,ektmp,Pressure,eptmp+ektmp+Pressure*qv
           
       end do

       do i=1, tsstep
         call MolecularDynamics(qn,pn,qv,pv,fn,fv,Pressure,ektmp)

         eptmp = m*w**2*qv**2/4/pi**2*(1-cos(2*pi*qn/qv))       !needs to be modified
         ep(j)  = ep(j) + eptmp/tsstep
         ek(j)  = ek(j) + ektmp/tsstep
         pres(j)=pres(j)+ Pressure/tsstep
!   write(*,*) pres(j)
!   read(*,*)
         write(999,*) i,pn,qn,pv,qv
         n=floor(10.0d0*qv)+1
         if(n<=100) vcount(n)=vcount(n)+1
!         pause
                     write(33,'(I16,F16.8,F16.8,F16.8,F16.8)') i,eptmp,ektmp,Pressure,eptmp+ektmp+Pressure*qv
           
        
         if (mod(i, tsstep/10+1) .eq. 0) then
             write(*,*) real(i)/real(tsstep)*100, '%'
             write(*,*) qn, pn,Pressure
         end if
        enddo
      close(33)
      close(999)
   enddo
   ep_ave = sum(ep)/sample
   ep_std = sqrt(sum((ep-ep_ave)**2)/(sample-1)/sample)
   ek_ave = sum(ek)/sample
   ek_std = sqrt(sum((ek-ek_ave)**2)/(sample-1)/sample)
   pres_ave=sum(pres)/sample
   write(*,*) pres_ave
   read(*,*)
   pres_std=sqrt(sum((pres-pres_ave)**2)/(sample-1)/sample)
   write(22,'(F8.2,F8.2,F16.8,F16.8,F16.8,F16.8)') h,gamma, ep_ave,ep_std,ek_ave,ek_std,pres_ave,pres_std
open(998,file="countv.txt")
   do i=1,100
     write(998,*) i,vcount(i)
     end do
close(998)
end program

