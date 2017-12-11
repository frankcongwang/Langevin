
          pn = pn + 0.5d0*h*fn
          qn = qn + 0.5d0*h*(pn/m)
          qv=qv*coe3
         
          pt2= pt2+0.25d0*h*(pt1*pt1/mQ-kT)
          coe2=exp(-0.125d0*h*pt2/mQ)
          pt1= pt1*coe2
          pt1= pt1+0.25d0*h*(pn*pn/m+pv*pv/Mex-2*kT)
          pt1= pt1*coe2
          qt1= qt1+ 0.5d0+h*pt1/mQ
          qt2= qt2+ 0.5d0+h*pt2/mQ
          coe1=exp( -0.25d0*h*pt1/mQ)
          pv = pv* coe1
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
  !        write(*,*) Pressure
          pv = pv+0.5d0*h*((pressure_ex-Pressure)*qv+pn**2/m)
          pv = pv*coe1
          pt1= pt1*coe2
          pt1= pt1+0.25d0*h*(pn*pn/m+pv*pv/Mex-2*kT)
          pt1= pt1*coe2
          pt2= pt2+0.25d0*h*(pt1*pt1/mQ-kT)

          coe1=exp( -0.25d0*h*pt1/mQ)
          coe3=exp(0.5d0*h*pv/Mex)
          pn = pn*coe1**4/(coe3**4)
          qn = qn*coe3*coe3

          pt2= pt2+0.25d0*h*(pt1*pt1/mQ-kT)
          coe2=exp(-0.125d0*h*pt2/mQ)
          pt1= pt1*coe2
          pt1= pt1+0.25d0*h*(pn*pn/m+pv*pv/Mex-2*kT)
          pt1= pt1*coe2
          qt1= qt1+ 0.5d0+h*pt1/mQ
          qt2= qt2+ 0.5d0+h*pt2/mQ
          coe1=exp( -0.25d0*h*pt1/mQ)
          pv = pv* coe1
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
  !        write(*,*) Pressure
          pv = pv+0.5d0*h*((pressure_ex-Pressure)*qv+pn**2/m)
          pv = pv*coe1
          pt1= pt1*coe2
          pt1= pt1+0.25d0*h*(pn*pn/m+pv*pv/Mex-2*kT)
          pt1= pt1*coe2
          pt2= pt2+0.25d0*h*(pt1*pt1/mQ-kT)


          coe3=exp(0.5d0*h*pv/Mex)
          qv=qv*coe3
          qn = qn + 0.5d0*h*(pn/m)
          call calForcex(fn,qn,qv)
          call calForceV(fv,qn,qv)
          call calPressure(Pressure,qv,pn,fn,qn,fv)
          pn = pn + 0.5d0*h*fn
