module init
  implicit none
  real(8), parameter :: unitl=1.0d0
  real(8), parameter :: kT = 1d0                                                  !needs to be modified
  real(8), parameter :: h = 0.1d0                                                       !needs to be modified
  integer, parameter :: num = 30
  real(8), parameter :: w = 1d0
  real(8), parameter :: m = 1d0
  integer, parameter :: eqstep=1d7/h
  integer, parameter :: tsstep=1d7/h
  integer, parameter :: sample=20
  integer, parameter :: limitofcut=1000
  real(8), parameter :: cutofforce=0.0000001d0
  real(8), parameter :: pi=3.14159265358979d0
  real(8) :: pn(num), qn(num)
end module init

program main
  call molphys
end program main

subroutine calV2(Epot,x,V)
  use init
  implicit none
  real(8) :: Epot, x, V
  Epot=k/x+1/2*log(r)
end subroutine calV2

subroutine calPotential(pot, V, n)
  use init
  implicit none
  real(8) :: pot, V
  integer :: n
  real(8), optinial :: pot
  real(8) :: energy1,force2,distance,energy1,energy2
  integer :: a,b
  Fn=0
  do i=1, limitofcut
    a=int(i/num)
    b=mod(a,num)
    distance=a*unitl+q(n+b)-q(n)
    call calImpact(force1,distance,V)
    distance=-a*unitl+q(n-b)-q(n)
    call calImpact(force2,distance,V)
    Fn=Fn+force1+force2
    if ((force1 < cutofforce) .and. (force2<cutofforce)) then
      exit
    end if
  end do
end subroutine calPotential

subroutine calImpact(fn, x, V)
  use init
  implicit none
  real(8) :: fn, x, V
  if x>0 then
    fn = k/x/x-1/2/x!-m*w**2*V/2/pi*sin(2*pi*x/V)                                                                              !needs to be modified
  else
    fn = -k/x/x-1/2/x
end subroutine calImpact

subroutine calForce(Fn, V, n)
  use init
  implicit none
  real(8) :: Fn, V
  integer :: n
  real(8) :: force1,force2,distance
  integer :: a,b
  Fn=0
  do i=1, limitofcut
    a=int(i/num)
    b=mod(a,num)
    distance=a*unitl+q(n+b)-q(n)
    call calImpact(force1,distance,V)
    distance=-a*unitl+q(n-b)-q(n)
    call calImpact(force2,distance,V)
    Fn=Fn+force1+force2
    if ((force1 < cutofforce) .and. (force2<cutofforce)) then
      exit
    end if
  end do
end subroutine calForce  


subroutine molphys
  use init
  use random
  implicit none
  real(8) :: rand, fn, fnp1, a, b
!  real(8) :: pn(num), qn(num)
  real(8) :: gamma = 0.8d0
  integer :: i, j, l
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

!       call random_normal(rand)
!       pn = rand
!       call random_number(rand)
!       qn = 4d0*(rand-0.5d0)
!       call calForce(fn, qn)

    !   write(*,*) 'sample=', j

!initial
       do i=1,num 
         qn(i)=unitl*i
         call random_normal(rand)
         pn(i)=rand
       end do

       do i=1, eqstep
          do l=1, num
          call calForce(fn, V, l)
          pn(l) = pn(l) + 0.5*h*fn

          qn(l) = qn(l) + 0.5*h*pn(l)/m
          if (qn(l)>unitl*num) then
            qn(l)=qn(l)-unitl*num
          else if(qn(l)<0) then
            qn(l)=qn(l)+unitl*num
          end if
          call random_normal(rand)
          pn(l) = a*pn(l) + sqrt((1-a*a)*kT)*sqrt(m)*rand 
                                
          qn(l) = qn(l) + 0.5*h*pn(l)/m                                 
          call calForce(fn, qn)

          pn(l) = pn(l) + 0.5*h*fn
          end do
           
       end do

       do i=1, tsstep
          do l=1, num
          call calForce(fn, V, l)
          pn(l) = pn(l) + 0.5*h*fn

          qn(l) = qn(l) + 0.5*h*pn(l)/m
          if (qn(l)>unitl*num) then
            qn(l)=qn(l)-unitl*num
          else if(qn(l)<0) then
            qn(l)=qn(l)+unitl*num
          end if
          call random_normal(rand)
          pn(l) = a*pn(l) + sqrt((1-a*a)*kT)*sqrt(m)*rand 
                                
          qn(l) = qn(l) + 0.5*h*pn(l)/m                                 
          call calForce(fn, qn)

          pn(l) = pn(l) + 0.5*h*fn
          end do
         if (mod(i, tsstep/10+1) .eq. 0) then
             write(*,*) real(i)/real(tsstep)*100, '%'
             write(*,*) qn, pn
         end if
         eptmp = !!!!!!!!!!!!!!!!!!1                                                               !needs to be modified
         ektmp = 0.5*sum(pn**2/m)
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

