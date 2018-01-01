module init
  implicit none
  real(8), parameter :: kT = 227.0d3/3.1577464d5                                                  !needs to be modified
  real(8), parameter :: h = 0.50d0                      !needs to be modified
  real(8), parameter :: sigma = 3.405d0/0.5291772d0
  real(8), parameter :: epsl = 119.8d0/3.1577464d5
  real(8), parameter :: m = 39.948d0
  real(8), parameter :: Mex = 18d4
  real(8), parameter :: pressure_ex = 1d5
  real(8), parameter :: initiallength=34.9826d0
  integer, parameter :: Noa=12
  integer, parameter :: eqstep=1d2/h
  integer, parameter :: tsstep=1d5/h
  integer, parameter :: sample=2
  integer, parameter :: Mtb=4   !or 6 the length of the NHC
  real(8), parameter :: mQ(Mtb) = 1d0
  real(8), parameter :: pi=3.14159265358979d0
  real(8), parameter :: rcut2=(2.5*sigma)**2
end module init

subroutine Initial(qx,qy,qz)
        use init
        implicit none
        real(8) :: qx(Noa), qy(Noa), qz(Noa)
        integer,parameter :: period=2
        integer :: i,j,k
        real(8) :: l=initiallength/period
        do i=1,Noa
          j=i/4/period/period
          qz(i)=l*j
          k=i-j*4*period*period
          j=k/4/period
          qy(i)=l*j
          k=k-j*4*period
          j=k/4
          qx(i)=l*j
          k=k-j*4
          select case(k)
          case(1)
                  qx(i)=qx(i)+l/2.0d0
                  qy(i)=qy(i)+l/2.0d0
          case(2)
                  qy(i)=qy(i)+l/2.0d0
                  qz(i)=qz(i)+l/2.0d0
          case(3)
                  qx(i)=qx(i)+l/2.0d0
                  qz(i)=qz(i)+l/2.0d0
          end select
        end do
  end subroutine

subroutine restrictCoord(qn,qv)
        use init
        implicit none
        real(8) :: qn(Noa,3),qv
        integer :: i,j
        do j=1,3
        do i=1,Noa
                if (qn(i,j)>=qv) qn(i,j)=qn(i,j)-qv
                if (qn(i,j)<0.0d0) qn(i,j)=qn(i,j)+qv
        end do
        end do
end subroutine


subroutine calPotential(Ep,qn,qv)
        use init
        implicit none
        real(8) :: Ep, qn(Noa,3),qv
        real(8) :: dist2
        integer :: i,j
        Ep=0.0d0
        do i=1,Noa
        do j=1,Noa
          if (j<i) then
                  dist2=sum((qn(i,:)-qn(j,:))**2)
                  if (dist2>=3*qv**2) write(*,*) "mmp"
                  if (dist2<rcut2)   Ep=Ep+4.0d0*epsl*(sigma**2/dist2)**6-(sigma**2/dist2)**3!;write(22,*) 1,j,Ep;read(*,*)
          end if
        end do
        end do
end subroutine

subroutine calKinetic(Ek,pn)
        use init
        implicit none
        real(8) :: Ek, pn(Noa,3)
        Ek=0.0d0
        Ek=sum(pn(:,:)**2)/2.0d0/m
end subroutine

subroutine calForcex(fn, qn)
        use init
        implicit none
        real(8) :: fn(Noa,3), qn(Noa,3)
        real(8) :: dist2, force(3)
        integer :: i,j
        fn=0.0d0
        do i=1,Noa
        do j=1,Noa
          if (j>=i) exit
          dist2=sum((qn(i,:)-qn(j,:))**2)
          if (dist2<rcut2) then
                  force=24.0d0*epsl*((2*sigma**12/dist2**7)-(sigma**6/dist2**4))*(qn(i,:)-qn(j,:))
                  fn(i,:)=fn(i,:)+force
                  fn(j,:)=fn(j,:)-force
        !          write(*,*) fn(:,:)
        !          read(*,*)
          end if
        end do
        end do
end subroutine calForcex

subroutine calPressure(Pins,Vol,enek,Fn,qn)
  use init
  implicit none
  real(8) :: Pins, Vol, enek, Fn(Noa,3), qn(Noa,3)
  real(8) :: qc(3),virial
  integer :: i
  do i=1,3
    qc(i)=sum(qn(:,i))/Noa
  end do
  virial=0.0d0
  do i=1,Noa
    virial=virial+sum((qn(i,:)-qc(:))*Fn(i,:))
  end do
  Pins=3.0d0/Vol*(2.0d0*enek+virial)
end subroutine calPressure

subroutine MolecularDynamics(qn,pn,qv,pv,fn,Pressure,enek)
        use init
        implicit none
        real(8) :: qn(Noa,3),pn(Noa,3),qv,pv,fn(Noa,3),Pressure,enek
          pn = pn + 0.5d0*h*fn
          qn = qn + 0.5d0*h*(pn/m)
          call restrictCoord(qn,qv)
!          call tbStat(qn,pn,qv,pv,fn,Pressure,enek)
          qn = qn + 0.5d0*h*(pn/m)
          call calForcex(fn,qn)
!          call calForceV(fv,qn,qv)
          call restrictCoord(qn,qv)
          pn = pn + 0.5d0*h*fn
          call calKinetic(enek,pn)
          call calPressure(Pressure,qv,enek,fn,qn)
  end subroutine

  subroutine tbStat(qn,pn,qv,pv,fn,Pressure,enek)
        use init
        implicit none
        real(8) :: qn(Noa,3),pn(Noa,3),qv,pv,fn(Noa,3),Pressure,enek
        real(8) :: qt(Mtb),pt(Mtb)
        real(8) :: gt(Mtb)
        integer :: i
        gt(Mtb)= pt(Mtb-1)*pt(Mtb-1)/Mq(Mtb-1)-kT

        pt(Mtb)=pt(Mtb)+0.5d0*h*gt(Mtb)

        do i=1,Mtb-1
          pt(Mtb-i)=pt(Mtb-i)*exp(-0.25d0*h*pt(Mtb-i+1)/Mq(Mtb-i+1))

          if (Mtb-i==1) gt(Mtb-i)=2.0d0*enek+pv*pv/Mex-kT !unchanged
          if (Mtb-i>1)  gt(Mtb-i)=pt(Mtb-i)*pt(Mtb-i)/Mq(Mtb-i)-kT
          pt(Mtb-i)=pt(Mtb-i)+0.5d0*h*gt(Mtb-i)

          pt(Mtb-i)=pt(Mtb-i)*exp(-0.25d0*h*pt(Mtb-i+1)/Mq(Mtb-i+1))

        end do
        pv=pv*exp(-0.25d0*pt(1)/Mq(1))
        call calKinetic(enek,pn)
        call calForcex(fn,qn)
!        call calForceV(fv,qn,qv)
        call calPressure(Pressure,qv,enek,fn,qn)
        pv=pv+0.5d0*h*(qv*(Pressure-Pressure_ex)+2.0d0*enek) !unchanged
        pv=pv*exp(-0.25d0*pt(1)/Mq(1))
        pn=pn*exp(-h*(2.0d0*pv/Mex+pt(1)/Mq(1)))
        
        qn=qn*exp(h*pv/Mex)
          call restrictCoord(qn,qv)
        qv=qv*exp(h*pv/Mex)
        do i=1,Mtb
          qt(i)=qt(i)+h*pt(i)
        end do

        pv=pv*exp(-0.25d0*pt(1)/Mq(1))
        call calKinetic(enek,pn)
        call calForcex(fn,qn)
!        call calForceV(fv,qn,qv)
        call calPressure(Pressure,qv,enek,fn,qn)
        pv=pv+0.5d0*h*(qv*(Pressure-Pressure_ex)+2.0d0*enek)
        pv=pv*exp(-0.25d0*pt(1)/Mq(1))

        gt(1)=2.0d0*enek+pv*pv/Mex-kT
        do i=1,Mtb-1
          pt(i)=pt(i)*exp(-0.25d0*h*pt(i+1)/Mq(i+1))

          pt(i)=pt(i)+0.5d0*h*gt(i)
          if (i==1) gt(Mtb)= pt(Mtb-1)*pt(Mtb-1)/Mq(Mtb-1)-kT
          if (Mtb-i>1)  gt(Mtb-i)=pt(Mtb-i)*pt(Mtb-i)/Mq(Mtb-i)-kT

          pt(i)=pt(i)*exp(-0.25d0*h*pt(i+1)/Mq(i+1))

        end do

end subroutine




program main
  use init
  use random
  implicit none
  real(8) :: rand, fn(Noa,3), Pressure
  real(8) :: qn(Noa,3),pn(Noa,3)
  real(8) :: qv, pv
  real(8) :: qt(Mtb), pt(Mtb)
  integer :: i, j, l
  real*8 :: eptmp, ektmp, ep(sample), ek(sample),pres(sample),ep_ave,ek_ave,pres_ave,ep_std,ek_std,pres_std
  real*8 :: cor_ep(tsstep/2),cor_ek(tsstep/2)
  real*8 :: t, cortimep, cortimek,cortimep_std,cortimek_std
  integer :: n!, ndt
  character(30) :: c
  open(22,file='result.maindat')
  ep(:)=0
  ek(:)=0
    write(*,*) 'dt=', h

    do j=1, sample
       write(c,'(I2)') j 
       write(*,*) 'Sample=', j
       open(33,file=trim('traj_'//adjustl(c)))
       open(999,file=trim('note_'//adjustl(c)))

       call initial(qn(1,1),qn(1,2),qn(1,3))
       call random_normal(rand)
       pv = 0.0d0 !sqrt(Mex*kT)*(rand-0.5d0)
       do i=1,Noa
         do l=1,3
           call random_normal(rand)
           pn = sqrt(m*kT)*(rand-0.5d0)
         end do
       end do
       qv=initiallength
!       do i=1,Noa
!       write(22,*) pn(i,1),pn(i,2),pn(i,3)
!       end do
!      read(*,*) 
!       do i=1,Mtb
!         call random_normal(rand)
!         qt(i)=sqrt(mQ(i))*(rand-0.5d0)
!         call random_normal(rand)
!         pt(i)=sqrt(mQ(i))*(rand-0.5d0)
!       end do
!
    !   write(*,*) 'sample=', j
    !      pause
          call calForcex(fn,qn)
!       do i=1,Noa
!       write(22,*) fn(i,1),fn(i,2),fn(i,3)
!       end do
!      read(*,*) 
       do i=1, eqstep
         call MolecularDynamics(qn,pn,qv,pv,fn,Pressure,ektmp)
         call calPotential(eptmp,qn,qv)
         ep(j)  = ep(j) + eptmp/tsstep
         ek(j)  = ek(j) + ektmp/tsstep
         pres(j)=pres(j)+ Pressure/tsstep
         write(999,*) i,pn,qn,pv,qv
!       p  pause
                     write(33,'(I16,E16.8,E16.8,E16.8,E16.8)') i,eptmp,ektmp,Pressure,ektmp+eptmp
           
       end do

       do i=1, tsstep
         call MolecularDynamics(qn,pn,qv,pv,fn,Pressure,ektmp)

         call calPotential(eptmp,qn,qv)
         ep(j)  = ep(j) + eptmp/tsstep
         ek(j)  = ek(j) + ektmp/tsstep
         pres(j)=pres(j)+ Pressure/tsstep
         write(999,*) i,pn,qn,pv,qv
!         pause
                     write(33,'(I16,E16.8,E16.8,E16.8,E16.8)') i,eptmp,ektmp,Pressure,ektmp+eptmp
           
        
         if (mod(i, tsstep/10+1) .eq. 0) then
             write(*,*) real(i)/real(tsstep)*100, '%'
   !          write(*,*) qn, pn
         end if
        end do
      close(33)
      close(999)
   end do
   ep_ave = sum(ep)/sample
   ep_std = sqrt(sum((ep-ep_ave)**2)/(sample-1)/sample)
   ek_ave = sum(ek)/sample
   ek_std = sqrt(sum((ek-ek_ave)**2)/(sample-1)/sample)
   pres_ave=sum(pres)/sample
   pres_std=sqrt(sum((pres-pres_ave)**2)/(sample-1)/sample)
!   write(22,'(F8.2,F8.2,F16.8,F16.8,F16.8,F16.8)') h, ep_ave,ep_std,ek_ave,ek_std,pres_ave,pres_std
       do i=1,Noa
       write(22,*) qn(i,1),qn(i,2),qn(i,3)
       end do

end program

