      program test_clout
      include 'implno.dek'

! declare
      double precision x,theta,ans(9),exact(9),fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta


! in physical terms x is chemical potantial dividedd by kerg*temp
! theta is the test mass energy divided by kerg*temp
      x     = 0.4d0
      theta = 0.1d0


! evaluate the fermi integrals using cloutman's method
      call clout(-0.5d0,x,theta,ans(1))
      call clout(0.5d0,x,theta,ans(2))
      call clout(1.0d0,x,theta,ans(3))
      call clout(1.5d0,x,theta,ans(4))
      call clout(2.0d0,x,theta,ans(5))
      call clout(2.5d0,x,theta,ans(6))
      call clout(3.0d0,x,theta,ans(7))
      call clout(4.0d0,x,theta,ans(8))
      call clout(5.0d0,x,theta,ans(9))



! evaluate the exact fermi integrals by using the more general
      call dfermi(-0.5d0,x,theta,exact(1),fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)
      call dfermi(0.5d0,x,theta,exact(2),fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)
      call dfermi(1.0d0,x,theta,exact(3),fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)
      call dfermi(1.5d0,x,theta,exact(4),fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)
      call dfermi(2.0d0,x,theta,exact(5),fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)
      call dfermi(2.5d0,x,theta,exact(6),fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)
      call dfermi(3.0d0,x,theta,exact(7),fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)
      call dfermi(4.0d0,x,theta,exact(8),fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)
      call dfermi(5.0d0,x,theta,exact(9),fdeta,fdtheta,fdeta2,fdtheta2,fdetadtheta)



! say what we got
      write(6,110) 'fermi integrals at x =',x,' theta =',theta
 110  format(1x,a,1pe12.4,a,1pe12.4)

      write(6,111) 'simpson extrapolation','exact quadrature'
 111  format(1x,t16,a,t42,a)

      write(6,112) ans(1),exact(1),ans(2),exact(2),ans(3),exact(3), &
                   ans(4),exact(4),ans(5),exact(5),ans(6),exact(6), &
                   ans(7),exact(7),ans(8),exact(8),ans(9),exact(9)

 112  format(1x,'order -1/2 :',1p2e24.16,/, &
             1x,'order  1/2 :',1p2e24.16,/, &
             1x,'order    1 :',1p2e24.16,/, &
             1x,'order  3/2 :',1p2e24.16,/, &
             1x,'order    2 :',1p2e24.16,/, &
             1x,'order  5/2 :',1p2e24.16,/, &
             1x,'order    3 :',1p2e24.16,/, &
             1x,'order    4 :',1p2e24.16,/, &
             1x,'order    5 :',1p2e24.16)

      stop 'normal termination'
      end




! routines for cloutman's method
      include 'cloutman.f90'



! routines for the brutally efficient quadrature integrations
      include 'fermi_dirac_quadrature.f90'
