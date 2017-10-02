module init
  implicit none
  real(8), parameter :: kT = 1d0                       !needs to be modified
  real(8), parameter :: h = 0.1d0                      !needs to be modified
  integer, parameter :: num = 3000
  real(8), parameter :: w = 1d0
  real(8), parameter :: m = 1d0
  integer, parameter :: eqstep=1d7/h
  integer, parameter :: tsstep=1d7/h
  integer, parameter :: sample=20
  real(8), parameter :: pi=3.14159265358979d0
end module init

program main
  call molphys
end program main
subroutine calForcex(fn, x, V)
  use init
  implicit none
  real(8) :: fn, x
  fn = -m*w**2*V/2/pi*sin(2*pi*x/V)              !needs to be modified
end subroutine calForcex

subroutine calForceV(fn, x, V)
  use init
  implicit none
  real(8) :: fn, x
  fn = -m*w**2*V/2/pi/pi*(1-cos(2*pi*x/V)-V*x*pi*ln(V)*sin(2*pi*x/V))                        !needs to be modified
end subroutine calForceV

subroutine calpressure(Pins,pi,Fi,ri)
  use init
  implicit none
  real(8) :: Pins, pi, Fi, ri
  Pins=

subroutine molphys
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qc, pc, qnp1, pnp1, fn, fnp1, a, b
  real(8) :: gamma = 0.8d0
  integer :: i, j
  real*8 :: eptmp, ektmp, ep(sample), ek(sample),ep_ave,ek_ave,ep_std,ek_std
  real*8 :: cor_ep(tsstep/2),cor_ek(tsstep/2)
  real*8 :: t, cortimep, cortimek,cortimep_std,cortimek_std
  integer :: n!, ndt
  character(30) :: c
  open(22,file='result.maindat')
  ep(:)=0
  ek(:)=0
    a = exp(-gamma*h)
    write(*,*) 'gamma=',gamma, 'dt=', h

    do j=1, sample
       write(c,'(I2)') j
       write(*,*) 'Sample=', j
       open(33,file=trim('traj_'//adjustl(c)))

       call random_normal(rand)
       pn = rand
       call random_normal(rand)
       pc = rand
       call random_number(rand)
       qn = 4d0*(rand-0.5d0)
       call random_normal(rand)
       pc = 1d0*(rand-0.5d0)
       call calForce(fn, qn)

    !   write(*,*) 'sample=', j
       do i=1, eqstep
          pn = pn + 0.5*h*fn

          qn = qn + 0.5*h*pn/m
          
          call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 
                                
          qn = qn + 0.5*h*pn/m                                 
          call calForce(fn, qn)

          pn = pn + 0.5*h*fn

           
       end do

       do i=1, tsstep
         pn = pn + 0.5*h*fn

          qn = qn + 0.5*h*pn/m
          
          call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 
                                
          qn = qn + 0.5*h*pn/m                                 
          call calForce(fn, qn)

          pn = pn + 0.5*h*fn
         if (mod(i, tsstep/10+1) .eq. 0) then
             write(*,*) real(i)/real(tsstep)*100, '%'
             write(*,*) qn, pn
         end if
         eptmp = 0.5*k*qn**2                                                               !needs to be modified
         ektmp = 0.5*pn**2/m
         ep(j) = ep(j)+eptmp/tsstep
         ek(j) = ek(j) + ektmp/tsstep
         if (i>tsstep-2000000) then
                     write(33,'(I16,F16.8,F16.8)') i,eptmp,ektmp
           endif
        enddo
      close(33)
   enddo
   ep_ave = sum(ep)/sample
   ep_std = sqrt(sum((ep-ep_ave)**2)/(sample-1)/sample)
   ek_ave = sum(ek)/sample
   ek_std = sqrt(sum((ek-ek_ave)**2)/(sample-1)/sample)
   write(22,'(F8.2,F8.2,F16.8,F16.8,F16.8,F16.8)') h,gamma, ep_ave,ep_std,ek_ave,ek_std
end subroutine molphys

