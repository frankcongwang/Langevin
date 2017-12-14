module init
  implicit none
  real(8), parameter :: kT = 227.0d0/3.1577464d5                                                  !needs to be modified
  real(8), parameter :: h = 5.0d0                      !needs to be modified
  real(8), parameter :: sigma = 3.405d0/0.5291772d0
  real(8), parameter :: m = 39.948d0
  real(8), parameter :: Mex = 180d0*40d0
  real(8), parameter :: pressure_ex = 1d0
  real(8), parameter :: initiallength=34.9826d0
  integer, parameter :: Noa=864
  integer, parameter :: eqstep=1d2/h
  integer, parameter :: tsstep=1d3/h
  integer, parameter :: sample=2
  integer, parameter :: Mtb=4   !or 6 the length of the NHC
  real(8), parameter :: mQ(Mtb) = 1d0
  real(8), parameter :: pi=3.14159265358979d0
  real(8), parameter :: rcut=2.5*sigma
end module init

subroutine Initial(qx,qy,qz)
        use init
        implicit none
        real(8) :: qx(Noa), qy(Noa), qz(Noa)
        integer :: period=6
        integer :: i,j,k,m
        real(8) :: l=initialength/period
        do i=1,Noa
          j=floor(Noa/4/period/period)
          qz(Noa)=l*j
          k=Noa-j*4*period*period
          j=floor(k/4/period)
          qy(Noa)=l*j
          k=k-j*4*period
          j=floor(k/4)
          qx(Noa)=l*j
          k=k-j*4
          select case(k)
          case(1)
                  qx(Noa)=qx(Noa)+l/2.0d0
                  qy(Noa)=qy(Noa)+l/2.0d0
          case(2)
                  qy(Noa)=qy(Noa)+l/2.0d0
                  qz(Noa)=qz(Noa)+l/2.0d0
          case(3)
                  qx(Noa)=qx(Noa)+l/2.0d0
                  qz(Noa)=qz(Noa)+l/2.0d0
          end select
  end subroutine



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
  fv = -m*w**2*V/2/pi/pi*(1-cos(2*pi*x/V)+V*x*pi*log(V)*sin(2*pi*x/V))                        !needs to be modified
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

subroutine MolecularDynamics(qn,pn,qv,pv,fn,fv,Pressure,enek)
        use init
        implicit none
        real(8) :: qn,pn,qv,pv,fn,fv,Pressure,enek
          pn = pn + 0.5d0*h*fn
          qn = qn + 0.5d0*h*(pn/m)
!          call tbStat(qn,pn,qv,pv,Pressure,enek)
          qn = qn + 0.5d0*h*(pn/m)
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          pn = pn + 0.5d0*h*fn
          call kEnergy(enek,pn)
          call calPressure(Pressure,qv,enek,fn,qn,fv)
  end subroutine

  subroutine tbStat(qn,pn,qv,pv,Pressure,enek)
        use init
        implicit none
        real(8) :: qn,pn,qv,pv,fn,fv,Pressure,enek
        real(8) :: qt(Mtb),pt(Mtb)
        real(8) :: gt(Mtb)
        integer :: i
        gt(Mtb)= pt(Mtb-1)*pt(Mtb-1)/Mq(Mtb-1)-kT

        pt(Mtb)=pt(Mtb)+0.5d0*h*gt(Mtb)

        do i=1,Mtb-1
          pt(Mtb-i)=pt(Mtb-i)*exp(-0.25d0*h*pt(Mtb-i+1)/Mq(Mtb-i+1))

          if (Mtb-i==1) gt(Mtb-i)=2.0d0*enek+pv*pv/Mex-kT !unchanged
          if (Mtb-i>1)  gt(Mtb-i)=pt(Mtb-i)*pt(Mtb-i)/Mq(Mtb-i)-kT

          pt(Mtb-i)=pt(Mtb-i)*exp(-0.25d0*h*pt(Mtb-i+1)/Mq(Mtb-i+1))

        end do
        pv=pv*exp(-0.25d0*pt(1)/Mq(1))
        call kEnergy(enek,pn)
        call calForcex(fn,qn,qv)
        call calForceV(fv,qn,qv)
        call calPressure(Pressure,qv,enek,fn,qn,fv)
        pv=pv+0.5d0*h*(qv*(Pressure-Pressure_ex)+2.0d0*enek) !unchanged
        pv=pv*exp(-0.25d0*pt(1)/Mq(1))
        pn=pn*exp(-h*(2.0d0*pv/Mex+pt(1)/Mq(1)))
        
        qn=qn*exp(h*pv/Mex)
        qv=qv*exp(h*pv/Mex)
        do i=1,Mtb
          qt(i)=qt(i)+h*pt(i)
        end do

        pv=pv*exp(-0.25d0*pt(1)/Mq(1))
        call kEnergy(enek,pn)
        call calForcex(fn,qn,qv)
        call calForceV(fv,qn,qv)
        call calPressure(Pressure,qv,enek,fn,qn,fv)
        pv=pv+0.5d0*h*(qv*(Pressure-Pressure_ex)+2.0d0*enek)
        pv=pv*exp(-0.25d0*pt(1)/Mq(1))

        gt(1)=2.0d0*enek+pv*pv/Mex-kT
        do i=1,Mtb-1
          pt(i)=pt(i)*exp(-0.25d0*h*pt(i+1)/Mq(i+1))

          if (i==1) gt(Mtb)= pt(Mtb-1)*pt(Mtb-1)/Mq(Mtb-1)-kT
          if (Mtb-i>1)  gt(Mtb-i)=pt(Mtb-i)*pt(Mtb-i)/Mq(Mtb-i)-kT

          pt(i)=pt(i)*exp(-0.25d0*h*pt(i+1)/Mq(i+1))

        end do

end subroutine




program main
  use init
  use random
  implicit none
  real(8) :: rand, fn, fv, a, b, coe1, coe2, coe3, Pressure
  real(8) :: qx(Noa), qy(Noa), qz(Noa), px(Noa), py(Noa), pz(Noa)
  real(8) :: qv(3), pv(3)
  real(8) :: qt(Mtb), pt(Mtb)
  real(8) :: gamma = 0.8d0
!  real(8) :: gammav = 1.0d0
  integer :: i, j
  real*8 :: eptmp, ektmp, ep(sample), ek(sample),pres(sample),ep_ave,ek_ave,pres_ave,ep_std,ek_std,pres_std
  real*8 :: cor_ep(tsstep/2),cor_ek(tsstep/2)
  real*8 :: t, cortimep, cortimek,cortimep_std,cortimek_std
  integer :: n!, ndt
  character(30) :: c
  qv=34.9826d0
  open(22,file='result.maindat')
  ep(:)=0
  ek(:)=0
    write(*,*) 'gamma=',gamma, 'dt=', h

    do j=1, sample
       write(c,'(I2)') j 
       write(*,*) 'Sample=', j
       open(33,file=trim('traj_'//adjustl(c)))
       open(999,file=trim('note_'//adjustl(c)))

       call initial(qx,qy,qz)
       call random_normal(rand)
       pv = 0.0d0 !0.05d0*(rand-0.5d0)
       call random_normal(rand)
       pn = 1.0d0 !1.5d0+0.05d0*rand
       qv=initiallength
       do i=1,Mtb
         call random_normal(rand)
         qt(i)=0.2d0*(rand-0.5d0)
         call random_normal(rand)
         pt(i)=0.2d0*(rand-0.5d0)
       end do

    !   write(*,*) 'sample=', j
    !      pause
          call calForcex(fn,qn,qv)
       do i=1, eqstep
         call MolecularDynamics(qn,pn,qv,pv,fn,fv,Pressure,ektmp)
         eptmp = m*w**2*qv**2/4/pi**2*(1-cos(2*pi*qn/qv))       !needs to be modified
         ep(j)  = ep(j) + eptmp/tsstep
         ek(j)  = ek(j) + ektmp/tsstep
         pres(j)=pres(j)+ Pressure/tsstep
         write(999,*) i,pn,qn,pv,qv
!       p  pause
                     write(33,'(I16,F16.8,F16.8,F16.8,F16.8)') i,eptmp,ektmp,Pressure,qv
           
       end do

       do i=1, tsstep
         call MolecularDynamics(qn,pn,qv,pv,fn,fv,Pressure,ektmp)

         eptmp = m*w**2*qv**2/4/pi**2*(1-cos(2*pi*qn/qv))       !needs to be modified
         ep(j)  = ep(j) + eptmp/tsstep
         ek(j)  = ek(j) + ektmp/tsstep
         pres(j)=pres(j)+ Pressure/tsstep
         write(999,*) i,pn,qn,pv,qv
!         pause
                     write(33,'(I16,F16.8,F16.8,F16.8,F16.8)') i,eptmp,ektmp,Pressure,qv
           
        
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
   write(22,'(F8.2,F8.2,F16.8,F16.8,F16.8,F16.8)') h,gamma, ep_ave,ep_std,ek_ave,ek_std,pres_ave,pres_std
end program

