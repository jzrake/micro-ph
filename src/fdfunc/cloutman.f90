      subroutine clout(dk,eta,theta,fer)
      include 'implno.dek'

!..evaluates the fermi-dirac integrals via direct integration.
!..input is the order dk, the degeneracy parameter eta,
!..and the relativity parameter theta.
!..output is the value of the fermi function.

!..implements cloutman's method (apjs 71, 677, 1989) for the fermi-dirac functions.
!..in essense, this method does a simpson's integration for three transforms
!..on three grids. the final answer is derived by aitken extrapolation of
!..the grid's answers in tandem with the "b and g" transforms.


!..declare the pass
      double precision dk,eta,theta,fer

!..locals
!..mptmin sets the number of grid points on the coarsest of the three grids
!..increasing mptmin increases the accuracy but slows the code down

      integer          mpts,mptmin,msize,itrap
      parameter        (mptmin = 101, msize=4*mptmin, itrap=0)

      integer          nfail,i,ngb,nm,m2,m4,j,ii
      double precision b,bad,bgt,bkt, &
                       denom,extrap,f(msize),fmax,fsum,gkt, &
                       gxfrm,h,rgb,rob,rog,clfunc,xalue(msize),x(msize), &
                       xmax,relerr,digits,pextr,cextr,xdum


!..a is the lower integration limit,
!..binc is the increment in the upper integration limit
!..bratio is the maximum value of f(b)/max(f)
!..gk and bk are the increments in the b and g transforms

      double precision third,a,binc,bratio,gk,bk
      parameter        (third  = 1.0d0/3.0d0, &
                        a      = 0.0d0, &
                        binc   = 0.5d0, &
                        bratio = 1.0d-12, &
                        gk     = 1.0d0, &
                        bk     = 1.1d0)


!..initialize
      bad  = 1.0d0
      mpts = mptmin
      nfail = 0


!..this loop increments the b & g transforms
       do ngb=1,3
        b = bad
        if (ngb .eq. 2) b = bad + gk
        if (ngb .eq. 3) b = bad + gk


!..start of mesh size loop
!..increments the number of mesh points and upper integration limit
        do nm=1,3
60       continue
         m2 = (mpts - 3)/2
         m4 = m2 + 1
         h  = (b-a)/float(mpts - 1)
         fmax = -1.0d100
         do j=1,mpts
          x(j) = a + h * dfloat(j-1)
          f(j) = clfunc(dk,x(j),eta,theta)
          if (f(j) .ge. fmax) then
           xmax = x(j)
           fmax = f(j)
          end if
         enddo
         if (ngb .ne. 1  .or. nm.ne.1 .or. &
             f(mpts).lt.bratio*fmax .or. b .gt. 100.0d0) go to 80
         b = b + binc
         bad = b
         go to 60
80       continue


!..use trapezoidal rule
         if (itrap .eq. 1) then
          fsum = 0.0d0
          do i=1,mpts
           fsum = fsum + f(i)
          enddo
          xalue(nm) = h*(fsum - 0.5d0*(f(1) + f(mpts)))

!..or use simpsons rule
         else
          fsum = 0.0d0
          do ii = 1,m2
           i = ii
           i = m2 - ii + 1
           j = m4 - i
           fsum = fsum + f(2*j+1) + f(2*j) * 2.0d0
          enddo
          xalue(nm) = third*h*(2.0d0*fsum+f(1)+f(mpts)+4.0d0*f(mpts-1))
         end if

!..end of mesh size loop
         if (nm .lt. 3) mpts = 2*(mpts-1) + 1
        enddo


!..now do the aitken extrapolation for the three mesh sizes
        denom = xalue(3) + xalue(1) - 2.0d0 * xalue(2)
        if (denom .ne. 0.0d0) then
         extrap = xalue(3) - (xalue(3) - xalue(2))**2 / denom
        else
         extrap = xalue(3)
         nfail = nfail + 1
        end if

!..store the intermediate b & g xform results and back for another set
        if (ngb .eq. 1) then
         bgt = extrap
         rgb = f(mpts)
        end if
        if (ngb .eq. 2) then
         gkt = extrap
         rog = f(mpts)
        end if
        if (ngb .eq. 3) then
         bkt = extrap
         rob = f(mpts)
        end if
        mpts = mptmin

!..end of b & g transform loop
       enddo



!..calculate the b and g transforms
       xdum = 1.0d0/rgb
       rog  = rog * xdum
       rob  = bk * rob * xdum


!..here is the value of the extrapolated integral
      fer = (gkt - rog*bgt) / (1.0d0 - rog)
      return
      end





      double precision function clfunc(dk,x,eta,theta)
      include 'implno.dek'
      double precision dk,x,eta,theta

!..this is the fermi-dirac integral after the substitution z**2 = x
!..see cloutman equation 7, modified for a non-zero theta

      clfunc = x**(2.0d0*dk + 1.0d0) * sqrt(1.0d0 + 0.5d0*theta*x*x) &
               * 2.0d0/(dexp(x*x - eta) + 1.0d0)

      return
      end
