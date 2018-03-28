module init
  implicit none
  real(8), parameter :: kT = 1d0                       !needs to be modified
  real(8), parameter :: h = 0.02d0                      !needs to be modified
       real(8), parameter :: w = 1d0
       real(8), parameter :: m = 1d0
       real(8), parameter :: Mex = 18d0
       real(8), parameter :: pressure_ex = 1d0
       integer, parameter :: eqstep=1d5/h
       integer, parameter :: tsstep=1d5/h
       integer, parameter :: sample=20
       real(8), parameter :: pi=3.14159265358979d0
       real(8), parameter :: gamma=1.0d0
       real(8), parameter :: gammaV=1.0d0
        real(8), parameter :: a = exp(-gamma*h)
        real(8), parameter :: b = exp(-gammaV*h)
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
        real(8) :: coe,rand
!          call tbStat(qn,pn,qv,pv,fn,fv,Pressure,enek,0.5d0*h)
        
        coe=exp(0.25d0*pv/Mex*h)
          pn = pn/(coe)**4 + 0.5d0*h*fn/(coe)**2
          call calForceV(fv,qn,qv)
          call kEnergy(enek,pn)
          call calPressure(Pressure,qv,enek,fn,qn,fv)
        pv=pv+0.5d0*h*(qv*(Pressure-Pressure_ex)+2.0d0*enek) !unchanged
        coe=exp(0.25d0*pv/Mex*h)
        qn=qn*coe*coe+0.5d0*h*(pn/Mex)*coe
          call restrictCoord(qn,qv)
        qv=qv*coe*coe
          call random_normal(rand)
          pv = b*pv + sqrt((1-b*b)*kT)*sqrt(Mex)*rand
          call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 
        !  qn = qn + 0.5d0*h*(pn/m)
        !  call restrictCoord(qn,qv)
        !  call tbStat(qn,pn,qv,pv,fn,fv,Pressure,enek,h)
        !  qn = qn + 0.5d0*h*(pn/m)
        !  call restrictCoord(qn,qv)
        !  call calForcex(fn,qn,qv)
        !  call calForceV(fv,qn,qv)
        coe=exp(0.25d0*pv/Mex*h)
        qv=qv*coe*coe
        qn=qn*coe*coe+0.5d0*h*(pn/Mex)*coe
          call restrictCoord(qn,qv)
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call kEnergy(enek,pn)
          call calPressure(Pressure,qv,enek,fn,qn,fv)
        pv=pv+0.5d0*h*(qv*(Pressure-Pressure_ex)+2.0d0*enek) !unchanged
        coe=exp(0.25d0*pv/Mex*h)
          pn = pn/(coe)**4 + 0.5d0*h*fn/(coe)**2
          call kEnergy(enek,pn)
          call calPressure(Pressure,qv,enek,fn,qn,fv)
  end subroutine

!  subroutine tbStat(qn,pn,qv,pv,fn,fv,Pressure,enek,dtstat)
!         use init
!         use NHCtherm
!         implicit none
!         real(8) :: qn,pn,qv,pv,fn,fv,Pressure,enek,dtstat
!         integer :: i
!         call kEnergy(enek,pn)
!         call calForcex(fn,qn,qv)
!         call calForceV(fv,qn,qv)
!         call calPressure(Pressure,qv,enek,fn,qn,fv)
!         pv=pv+0.5d0*dtstat*(qv*(Pressure-Pressure_ex)+2.0d0*enek) !unchanged
!         qv=
!         pv=pv*exp(-0.5d0*dtstat*pt(1)/Mq(1))
!         pn=pn*exp(-dtstat*(2.0d0*pv/Mex+pt(1)/Mq(1)))
!         
!         qn=qn*exp(dtstat*pv/Mex)
!           call restrictCoord(qn,qv)
!         qv=qv*exp(dtstat*pv/Mex)
!         do i=1,Mtb
!           qt(i)=qt(i)+dtstat*pt(i)
!         end do
!
!         pv=pv*exp(-0.5d0*dtstat*pt(1)/Mq(1))
!         call kEnergy(enek,pn)
!         call calForcex(fn,qn,qv)
!         call calForceV(fv,qn,qv)
!         call calPressure(Pressure,qv,enek,fn,qn,fv)
!         pv=pv+0.5d0*dtstat*(qv*(Pressure-Pressure_ex)+2.0d0*enek)
!!         pv=pv*exp(-0.25d0*dtstat*pt(1)/Mq(1))
!
!         gt(1)=2.0d0*enek+pv*pv/Mex-2.0d0*kT
!         do i=1,Mtb-1
!           pt(i)=pt(i)*exp(-0.25d0*dtstat*pt(i+1)/Mq(i+1))
!
!           pt(i)=pt(i)+0.5d0*dtstat*gt(i)
!!           if (i==1) gt(2)= pt(2)*pt(2)/Mq(2)-kT
!!           if (i>1) 
!gt (i+1)=pt(i)*pt(i)/Mq(i)-kT
!
!           pt(i)=pt(i)*exp(-0.25d0*dtstat*pt(i+1)/Mq(i+1))
!!        if (i==1) gt(2)= pt(2)*pt(2)/Mq(2)-kT
!!        if (i>1)  gt(i+1)=pt(i+1)*pt(i+1)/Mq(i+1)-kT
!
!         end do
!         pt(Mtb)=pt(Mtb)+0.5d0*dtstat*gt(Mtb)
!en d subroutine




program main
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qv, pv, fn, fv, Pressure
  integer :: i, j
  real*8 :: eptmp, ektmp, ep(sample), ek(sample),pres(sample),ep_ave,ek_ave,pres_ave,ep_std,ek_std,pres_std
  real*8 :: cor_ep(tsstep/2),cor_ek(tsstep/2)
  real*8 :: t, cortimep, cortimek,cortimep_std,cortimek_std
  integer :: n, vcount(100)
  character(30) :: c
  vcount=0
  open(22,file='result.maindat')
  ep(:)=0
  ek(:)=0
    write(*,*) 'Middle', 'dt=', h
pres=0.0d0

    do j=1, sample
       write(c,'(I2)') j 
       write(*,*) 'Sample=', j
       open(33,file=trim('traj_'//adjustl(c)))
       open(999,file=trim('note_'//adjustl(c)))

       call random_normal(rand)
       pn = rand
       call random_normal(rand)
       pv = 0.05d0*rand
       call random_normal(rand)
       qn = rand
       call random_normal(rand)
       qv = 10.0d0 !1.5d0+0.05d0*rand

   !    write(*,*) 'sample=', j
    !      pause
          call calForcex(fn,qn,qv)
          fv=0.0d0
          Pressure=0.0d0
          ektmp=0.0d0
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
         write(999,*) i,pn,qn,pv,qv
         n=floor(10.0d0*qv)+1
         if(n<=100) vcount(n)=vcount(n)+1
!         pause
                     write(33,'(I16,F16.8,F16.8,F16.8,F16.8)') i,eptmp,ektmp,Pressure,eptmp+ektmp+Pressure*qv
           
        
         if (mod(i, tsstep/10+1) .eq. 0) then
             write(*,*) real(i)/real(tsstep)*100, '%'
             write(*,*) qn, pn
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
   pres_std=sqrt(sum((pres-pres_ave)**2)/(sample-1)/sample)
   write(22,'(F8.2,F16.8,F16.8,F16.8,F16.8)') h, ep_ave,ep_std,ek_ave,ek_std,pres_ave,pres_std
open(998,file="countv.txt")
   do i=1,100
     write(998,*) i,vcount(i)
     end do
close(998)
end program

