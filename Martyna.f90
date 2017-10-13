module init
  implicit none
  real(8), parameter :: kT = 1d0                       !needs to be modified
  real(8), parameter :: h = 0.01d0                      !needs to be modified
  real(8), parameter :: w = 1d0
  real(8), parameter :: m = 1d0
  real(8), parameter :: Mex = 18d0
  real(8), parameter :: pressure_ex = 1d0
  integer, parameter :: eqstep=1d5/h
  integer, parameter :: tsstep=1d6/h
  integer, parameter :: sample=20
  real(8), parameter :: pi=3.14159265358979d0
end module init

program main
  call molphys
end program main
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

subroutine molphys
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qv, pv, qnp1, pnp1, fn, fv, a, b, Pressure
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
    a = exp(-gamma*h)
    b = exp(-gammaV*h)
    write(*,*) 'gamma=',gamma, 'dt=', h

    do j=1, sample
       write(c,'(I2)') j
       write(*,*) 'Sample=', j
       open(33,file=trim('traj_'//adjustl(c)))

       call random_normal(rand)
       pn = rand-0.5d0
       call random_normal(rand)
       pv = rand-0.5d0
       call random_number(rand)
       qn = 4d0*(rand-0.5d0)
       call random_normal(rand)
       qv = 1d0+0.5d0*rand

    !   write(*,*) 'sample=', j
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
   !       write(*,*) fn,fv
          call calPressure(Pressure,qv,pn,fn,qn,fv)
   !       write(*,*) Pressure
    !      pause
       do i=1, eqstep
          pv = pv + 0.5*h*((pressure_ex-Pressure)*qv+pn**2/m)
          pn = pn*exp(-0.5*h*pv/Mex) + 0.5*h*fn
          qv = qv*exp( 0.5*h*pv/Mex)
          qn = qn*exp( 0.5*h*pv/Mex) + 0.5*h*(pn/m)
          
          call random_normal(rand)
          pv = b*pv + sqrt((1-b*b)*kT)*sqrt(Mex)*rand
          call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 
                                
          qn = qn*exp( 0.5*h*pv/Mex) + 0.5*h*(pn/m)
          qv = qv*exp( 0.5*h*pv/Mex)
          call calForcex(fn,qn,qv)
          call calForceV(fn,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
          pn = pn*exp(-0.5*h*pv/Mex) + 0.5*h*fn
          pv = pv + 0.5*h*((pressure_ex-Pressure)*qv+pn**2/m)
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)

         eptmp = m*w**2*qv**2/4/pi**2*(1-cos(2*pi*qn/qv))       !needs to be modified
         ektmp = 0.5*pn**2/m
         ep(j)  = ep(j) + eptmp/tsstep
         ek(j)  = ek(j) + ektmp/tsstep
         pres(j)=pres(j)+ Pressure/tsstep
  !       write(*,*) pn,qn,pv,qv,Pressure
  !       pause
                     write(33,'(I16,F16.8,F16.8,F16.8,F16.8)') i,eptmp,ektmp,Pressure,qv
           
       end do

       do i=1, tsstep
          pv = pv + 0.5*h*((pressure_ex-Pressure)*qv+pn**2/m)
          pn = pn*exp(-0.5*h*pv/Mex) + 0.5*h*fn
          qv = qv*exp( 0.5*h*pv/Mex)
          qn = qn*exp( 0.5*h*pv/Mex) + 0.5*h*(pn/m)
          
          call random_normal(rand)
          pv = b*pv + sqrt((1-b*b)*kT)*sqrt(Mex)*rand
          call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 
                                
          qn = qn*exp( 0.5*h*pv/Mex) + 0.5*h*(pn/m)
          qv = qv*exp( 0.5*h*pv/Mex)
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
          pn = pn*exp(-0.5*h*pv/Mex) + 0.5*h*fn
          pv = pv + 0.5*h*((pressure_ex-Pressure)*qv+pn**2/m)
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)

        
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
end subroutine molphys

