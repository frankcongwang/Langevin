module init
  implicit none
  real(8), parameter :: kT = 1d0                       !needs to be modified
  real(8), parameter :: h = 0.01d0                      !needs to be modified
  real(8), parameter :: w = 1d0
  real(8), parameter :: m = 1d0
  real(8), parameter :: Mex = 18d0
  real(8), parameter :: mQ = 1d0
  real(8), parameter :: pressure_ex = 2d0
  integer, parameter :: eqstep=1d2/h
  integer, parameter :: tsstep=1d3/h
  integer, parameter :: sample=20
  real(8), parameter :: pi=3.14159265358979d0
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
  fv = -m*w**2*V/2/pi/pi*(1-cos(2*pi*x/V)-V*x*pi*log(V)*sin(2*pi*x/V))                        !needs to be modified
end subroutine calForceV

subroutine calPressure(Pins,Vol,pl,Fl,rl,Fv)
  use init
  implicit none
  real(8) :: Pins, Vol, pl, Fl, rl, Fv
!  real(8) :: 
  Pins=1d0/Vol*(pl**2/m)+Fv !+rl*Fl
end subroutine calPressure



program molphys
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qv, pv, qt1, pt1, qt2, pt2, fn, fv, a, b, coe1, coe2, coe3, Pressure
  real(8) :: gamma = 0.8d0
  real(8) :: gammav = 1.0d0
  integer :: i, j
  real*8 :: eptmp, ektmp, ep(sample), ek(sample),pres(sample),ep_ave,ek_ave,pres_ave,ep_std,ek_std,pres_std
  real*8 :: cor_ep(tsstep/2),cor_ek(tsstep/2)
  real*8 :: t, cortimep, cortimek,cortimep_std,cortimek_std
  integer :: n!, ndt
  character(30) :: c
  open(22,file='result.maindat')
  ep(:)=0
  ek(:)=0
    write(*,*) 'gamma=',gamma, 'dt=', h

    do j=1, sample
       write(c,'(I2)') j
       write(*,*) 'Sample=', j
       open(33,file=trim('traj_'//adjustl(c)))

       call random_normal(rand)
       pn = rand-0.5d0
       call random_normal(rand)
       pv = rand-0.5d0
       call random_normal(rand)
       pt1= rand-0.5d0
       call random_normal(rand)
       pt2= rand-0.5d0
       call random_number(rand)
       qn = 1d0*(rand-0.5d0)
       call random_normal(rand)
       qv = 1d0+0.2d0*rand
       call random_normal(rand)
       qt1= 1d0*(rand-0.5d0)
       call random_normal(rand)
       qt2= 1d0*(rand-0.5d0)

    !   write(*,*) 'sample=', j
    !      pause
          coe3=exp(0.5d0*h*pv/Mex)
          call calForcex(fn,qn,qv)
       do i=1, eqstep
          pn = pn + 0.5d0*h*fn
          qn = qn + 0.5d0*h*(pn/m)
          qv=qv*coe3
         
          pt2= pt2+0.5d0*h*(pt1*pt1/mQ-kT)
          coe2=exp(-0.25d0*h*pt2/mQ)
          pt1= pt1*coe2
          pt1= pt1+0.5d0*h*(pn*pn/m+pv*pv/Mex-2*kT)
          pt1= pt1*coe2

          coe1=exp( -0.25d0*h*pt1/mQ)
          pv = pv* coe1
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
  !        write(*,*) Pressure
          pv = pv+0.5d0*h*((pressure_ex-Pressure)*qv+pn**2/m)
          pv = pv*coe1


          qt1= qt1+ h*pt1/mQ
          qt2= qt2+ h*pt2/mQ

          coe3=exp(0.5d0*h*pv/Mex)
          pn = pn*coe1**4/(coe3**4)
          qn = qn*coe3*coe3

          coe1=exp( -0.25d0*h*pt1/mQ)
          pv = pv* coe1
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
  !        write(*,*) Pressure
          pv = pv+0.5d0*h*((pressure_ex-Pressure)*qv+pn**2/m)
          pv = pv*coe1
          pt1= pt1*coe2
          pt1= pt1+0.5d0*h*(pn*pn/m+pv*pv/Mex-2*kT)
          pt1= pt1*coe2
          pt2= pt2+0.25d0*h*(pt1*pt1/mQ-kT)


          coe3=exp(0.5d0*h*pv/Mex)
          qv=qv*coe3
          qn = qn + 0.5d0*h*(pn/m)
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
          pn = pn + 0.5d0*h*fn

         eptmp = m*w**2*qv**2/4/pi**2*(1-cos(2*pi*qn/qv))       !needs to be modified
         ektmp = 0.5*pn**2/m
         ep(j)  = ep(j) + eptmp/tsstep
         ek(j)  = ek(j) + ektmp/tsstep
         pres(j)=pres(j)+ Pressure/tsstep
!         write(*,*) pn,qn,pv,qv,Pressure
!         pause
                     write(33,'(I16,F16.8,F16.8,F16.8,F16.8)') i,eptmp,ektmp,Pressure,qv
           
       end do

       do i=1, tsstep
          pn = pn + 0.5d0*h*fn
          qn = qn + 0.5d0*h*(pn/m)
          qv=qv*coe3
         
          pt2= pt2+0.5d0*h*(pt1*pt1/mQ-kT)
          coe2=exp(-0.25d0*h*pt2/mQ)
          pt1= pt1*coe2
          pt1= pt1+0.5d0*h*(pn*pn/m+pv*pv/Mex-2*kT)
          pt1= pt1*coe2

          coe1=exp( -0.25d0*h*pt1/mQ)
          pv = pv* coe1
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
  !        write(*,*) Pressure
          pv = pv+0.5d0*h*((pressure_ex-Pressure)*qv+pn**2/m)
          pv = pv*coe1


          qt1= qt1+ h*pt1/mQ
          qt2= qt2+ h*pt2/mQ

          coe3=exp(0.5d0*h*pv/Mex)
          pn = pn*coe1**4/(coe3**4)
          qn = qn*coe3*coe3

          coe1=exp( -0.25d0*h*pt1/mQ)
          pv = pv* coe1
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
  !        write(*,*) Pressure
          pv = pv+0.5d0*h*((pressure_ex-Pressure)*qv+pn**2/m)
          pv = pv*coe1
          pt1= pt1*coe2
          pt1= pt1+0.5d0*h*(pn*pn/m+pv*pv/Mex-2*kT)
          pt1= pt1*coe2
          pt2= pt2+0.25d0*h*(pt1*pt1/mQ-kT)


          coe3=exp(0.5d0*h*pv/Mex)
          qv=qv*coe3
          qn = qn + 0.5d0*h*(pn/m)
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
          pn = pn + 0.5d0*h*fn

         eptmp = m*w**2*qv**2/4/pi**2*(1-cos(2*pi*qn/qv))       !needs to be modified
         ektmp = 0.5*pn**2/m
         ep(j)  = ep(j) + eptmp/tsstep
         ek(j)  = ek(j) + ektmp/tsstep
         pres(j)=pres(j)+ Pressure/tsstep
!         write(*,*) pn,qn,pv,qv,Pressure
!         pause
                     write(33,'(I16,F16.8,F16.8,F16.8,F16.8)') i,eptmp,ektmp,Pressure,qv

        
         if (mod(i, tsstep/10+1) .eq. 0) then
             write(*,*) real(i)/real(tsstep)*100, '%'
             write(*,*) qn, pn
         end if
         eptmp = m*w**2*qv**2/4/pi**2*(1-cos(2*pi*qn/qv))       !needs to be modified
         ektmp = 0.5*pn**2/m
         ep(j)  = ep(j) + eptmp/tsstep
         ek(j)  = ek(j) + ektmp/tsstep
         pres(j)=pres(j)+ Pressure/tsstep
                     write(33,'(I16,F16.8,F16.8,F16.8,F16.8)') i,eptmp,ektmp,Pressure,qv
        enddo
      close(33)
   enddo
   ep_ave = sum(ep)/sample
   ep_std = sqrt(sum((ep-ep_ave)**2)/(sample-1)/sample)
   ek_ave = sum(ek)/sample
   ek_std = sqrt(sum((ek-ek_ave)**2)/(sample-1)/sample)
   pres_ave=sum(pres)/sample
   pres_std=sqrt(sum((pres-pres_ave)**2)/(sample-1)/sample)
   write(22,'(F8.2,F8.2,F16.8,F16.8,F16.8,F16.8)') h,gamma, ep_ave,ep_std,ek_ave,ek_std,pres_ave,pres_std
end program molphys

