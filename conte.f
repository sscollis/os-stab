c***********************************************************************
      program conte_chan
c***********************************************************************
c 
c     Purpose:  Solve the Orr-Sommerfeld equation using 4th order Runge-
c               Kutta explicit integration with pseudoorthogonalization
c
c               This version assumes a parabolic channel profile
c
c     Citation: Follows the method of S.D. Conte in 
c               SIAM Review, v8, n3, July 1966
c               "The Numerical Solution of Linear Boundary Value
c                Problems"
c
c     Output:   Complex valued eigenfunction is phi 
c
c               fort.10     Mean velocity profile
c               fort.11     phi(y)
c               fort.12     phi'(y)
c               fort.13     phi''(y)
c               fort.14     phi'''(y)
c
c     Author:   Scott Collis
c
c     Date:     6-30-1992
c
c     Revised:  12-29-2019
c
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (imax=10000)
      integer     n
      complex     alpha, c, rgamma
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax)

      common      /setup/ n, alpha, u, d2u, Re, ymin, ymax, h, 
     &                    rgamma, uspl, d2uspl, ydat
     
      common      /eig/   c
c***********************************************************************
      parameter (neq=4)
      complex   bc1(neq), bc2(neq)
      character input*1
c
c     Defaults
c      
      n = 2000
      Re = 10000.0
      alphar = 1.0
      alphai = 0.0
      cr =  0.23753
      ci =  0.00374
      ymin = -1.0
      ymax =  1.0
c
c     Input parameters
c
      write (*,10)
  10  format (/,10x,'Solve Orr-Sommerfeld for Channel using shooting',/)
      write (*,15)
  15  format(/,1x,'Read from keyboard, file, default (k,d) ==> ',$)
      read (*,'(a)') input
      if (input .eq. 'k' .or. input .eq. 'K') then
        write (*,20)
  20    format (/,1x,'Enter number of steps ==> ',$)
        read (*,*) n
        write (*,30)
  30    format (/,1x,'Enter Reynolds number ==> ',$)
        read (*,*) Re
        write (*,40) 
  40    format (/,1x,'Enter alpha_r ==> ',$)
        read (*,*) alphar
        write (*,45) 
  45    format (/,1x,'Enter alpha_i ==> ',$)
        read (*,*) alphai
        write (*,50)
  50    format (/,1x,'Enter guess for (c_r, c_i) ==> ',$)
        read (*,*) cr, ci
      end if
c
c     Echo input
c
      write (*,100) n
 100  format (/,1x,'n = ',i5)
      write (*,110) Re
 110  format (1x,'Re = ',e20.10)
      write (*,120) alphar, alphai
 120  format (1x,'alpha = (',e20.10,', ',e20.10,')')
      write (*,130) cr, ci
 130  format (1x,'c = (',e20.10,', ',e20.10,')',/)
c
c     Set constants
c
      alpha = cmplx(alphar,alphai)
      c = cmplx(cr,ci)
      h = (ymax-ymin)/float(n)
c
c     Check the problem size
c
      if (n .gt. imax) then
        write (*,300)
 300    format (/,/,1x,'N > Idim...Stopping',/,/)
        goto 210
      endif
c
c     This really just makes a parabolic profile
c
      call SOLVE_BL
c
c     set the boundary conditions
c
      do i = 1, neq
        bc1(i) = 0.0
        bc2(i) = 0.0
      end do
      
      call CONTE(n, neq, 2, bc1, bc2, ymin, ymax)

 210  continue      
      stop
      end

C***********************************************************************
      subroutine CONTE(nstep, n, r, yo, yf, to, tf)
C***********************************************************************
C
C     First order linear boundary value problem solver using Conte's
C     method.  Fourth order Runge-Kutta is used for time advancement
C
C     nstep = number of integration steps between [t0,tf]
C     n     = number of ODEs to integrate
C     r     = number of particular solutions 
C     yo    = starting values
C     yf    = ending values
C     to    = Interval start value
C     tf    = Interval end value
C
c***********************************************************************
      complex     c
      common      /eig/ c
C***********************************************************************
      integer     i, m, q, r, s, mi, mj, qend, IPVT(n-r), icount
      real        t, tq(0:nstep), to, tf, h
      complex     yo(n), yf(n), B(n-r,0:nstep), err, cold, ctemp
      complex     U(n,n-r,0:nstep), P(n-r,n-r,0:nstep), z(n,n-r)
      complex     y(n,0:nstep), v(n,0:nstep), w(n), omega(n,0:nstep)
      complex     eta(n),A(n,n-r),x(n),FAC(n-r,n-r), det1, det
      complex     olderr, INPROD, ut(n,n-r), fd, vmax

      complex     cm1, cm2, errm1, errm2, qt, At, Bt, Ct, Utemp(n)

      real        test, testalpha, pi, det2, AI(n,n-r), angle
      real        aa, bb, cc, fdxr, fdxi, fdyr, fdyi
      logical     norm, eigfun 
      external    INPROD, FHOMO, FPART

      real        errtol, maxcount, eigtol
c
c     initialize variables
c
      errtol = 1.0e-15
      maxcount = 20
      eigtol = 1.0e-15
      eigfun = .true.
c
c     set the normalization constraint
c
      pi = ACOS(-1.0)
c
c     Orthonormalize when angle is greater than testalpha
c     angle = [0,90] degrees (10 deg seems good)
c
      angle = 10.0
      testalpha = angle*pi/180.0
c
c     compute the step size
c      
      h = (tf-to)/nstep
c
c     Begin the eigenvalue iteration loop
c
      icount = 1
      err = 1.0
      cm1 = c + 1.0
      do while ((abs(err) .ge. errtol .and. icount .le. maxcount) .and.
     &          (abs(c-cm1) .ge. eigtol) )
c
c       q = number of normalizations
c       qend = total number of normalizations
c       tq = time at normalization q
c
        q = 0
        tq(0) = to
c
c       Set the Homogeneous initial conditions
c
        do m = 1, n-r
          do i = 1, n
            if (i .eq. r+m) then
              u(i,m,0) = 1.0
            else
              u(i,m,0) = 0.0
            end if
          end do
        end do
c
c       Set the particular initial conditions
c
        do i = 1, r
          v(i,0) = yo(i)
        end do
        do i = r+1, n
          v(i,0) = 0.0
        end do
c
c       Integrate the homo. and particular equations
c      
        do k = 1, nstep
          t = to + h*k
c
c         Loop thru all homogeneous solutions
c
          do m = 1, n-r
            call RK4(n, U(1,m,k-1), U(1,m,k), t-h, h, FHOMO)
          end do
c
c         Now the particular solution
c
          call RK4(n, v(1,k-1), v(1,k), t-h, h, FPART)
c
c         Check to see it orthonomalization is required
c       
          norm = .false.
          do mi = 1, n-r
            do mj = 1, n-r
              if (mi .ne. mj) then
                aa = ABS(inprod(n, U(1,mi,k), U(1,mi,k)))
                bb = ABS(inprod(n, U(1,mj,k), U(1,mj,k)))
                cc = ABS(inprod(n, U(1,mi,k), U(1,mj,k)))
#ifdef USE_ANALYTIC_INPROD
                test = ACOS(MIN(1.0,ABS(cc)/SQRT(aa*bb)))
#else
                test = ACOS(ABS(cc)/SQRT(aa*bb))
#endif
                if (test .le. testalpha) norm = .true.
#ifdef VERBOSE
                write (*,*) k, mi, mj, 180.0*test/pi, 
     &                      180.0*testalpha/pi, norm
#endif
              end if
            end do
          end do
c
c         Perform normalization
c       
c         if (norm .or. (k .eq. nstep) .or. (mod(k,10) .eq. 0) ) then
          if (norm .or. (k .eq. nstep)) then
            q = q + 1
            tq(q) = t
            if (k .eq. nstep) qend = q
c
c           Gram-Schmidt
c
            w(1) = SQRT(inprod(n, U(1,1,k), U(1,1,k)))
            do i = 1, n
              z(i,1) = U(i,1,k)/w(1)
            end do
            do mi = 2, (n-r)
              do i = 1, n
                eta(i) = U(i,mi,k)
              end do
              do mj = mi-1, 1, -1
                do i = 1, n
                  eta(i) = eta(i) - inprod(n, U(1,mi,k), z(1,mj))*z(i,mj)
                end do
              end do
              w(mi) = SQRT(inprod(n, eta, eta))
              do i = 1, n
                z(i,mi) = eta(i)/w(mi)
              end do
            end do
c
c           Now I have the orthonormal basis in z and 
c           the norms in w so I can compute the P orthonormalization 
c           matrix
c
            do j = 1, n-r
              do i = 1, j
                if (j .eq. i) then
                  P(j,i,q) = 1.0/w(j)
                else
                  P(j,i,q) = 0.0
                  do s = i, j-1
                    P(j,i,q) = P(j,i,q)-inprod(n,U(1,j,k),z(1,s))/w(j)*
     &                         P(s,i,q)
                  end do
                end if
              end do
            end do
c
c           Check the P matrix
c
            if (.false.) then
              do i = 1, n
                do m = 1, n-r
                ut(i,m) = 0.0
                  do j = 1, n-r
                    ut(i,m) = ut(i,m) + U(i,j,k)*P(m,j,q)
                  end do
                end do
              end do
              do i = 1,n
                write (*,*) i,(ut(i,m), m = 1, n-r)
              end do
              write (*,*)
              do i = 1,n
                write (*,*) i,( z(i,m), m = 1, n-r)
              end do
              write (*,*)
              write (*,*)
            end if
c
c           Now update the U matrix with the orthonormal values
c       
            do i = 1, n
              do m = 1, n-r
                U(i,m,k) = z(i,m)
              end do
            end do
c
c           Form the orthogonal complement of the particular solution
c
            do j = 1, n-r
              omega(j,q) = inprod(n, v(1,k), U(1,j,k))
            end do
            do i = 1, n
              do j = 1, n-r
                v(i,k) = v(i,k) - U(i,j,k)*omega(j,q)
              end do
            end do
  
          end if
        end do
c
c       write out the current solutions
c
#ifdef VERBOSE_SOLUTIONS
        do k = 1, nstep
          t = to + h*k
          write (20,10) t,REAL(U(1,1,k)),AIMAG(U(1,1,k)),REAL(U(1,2,k)),
     &                    AIMAG(U(1,2,k)),REAL(v(1,k)),AIMAG(v(1,k))
 10       format (1x,7(e12.4,1x))
        end do
#endif
c
c       check boundary conditions
c
        if (.true.) then
c
c         strictly enforce the zero BC
c
          B(1,qend) = -U(1,2,nstep)/U(1,1,nstep)
          B(2,qend) = 1.0
        else
c
c         strictly enforce the zero slope BC
c
          B(1,qend) =  1.0
          B(2,qend) = -U(2,1,nstep)/U(2,2,nstep)
        end if
#ifdef CHECK_DETERMINATE
c
c       Check the determinate (requires IMSL)
c
        do i = 1, n 
          x(i) = yf(i) - v(i,nstep)
          do m = 1, n-r
            A(i,m) = U(i,m,nstep)
          end do
        end do
        CALL WRCRN ('A', n, n-r, A, n, 0)
c
c       compute the determinate
c      
        call LFTCG( n-r, A, n, fac, n-r, ipvt)
        call LFDCG( n-r, fac, n-r, ipvt, det1, det2)   
        det = det1*10**det2
        write (*,*) 'det = ', det
#endif
        ctemp = c
        if (icount .eq. 1) then
          cm2 = c
          err = U(2,1,nstep)*B(1,qend) + U(2,2,nstep)*B(2,qend)
          errm2 = err
          c = CMPLX( REAL(c)*.9999, AIMAG(c) )
        else if (icount .eq. 2) then
          cm1 = c
          err = U(2,1,nstep)*B(1,qend) + U(2,2,nstep)*B(2,qend)
          errm1 = err
c         c = CMPLX( REAL(c), AIMAG(c)*.9999 )
          fdxr =  REAL(err-errm2) / REAL(c-cm2)
          fdxi = AIMAG(err-errm2) / REAL(c-cm2)
          fd = CMPLX( fdxr, fdxi )
          c = c - err/fd
        else
          err = U(2,1,nstep)*B(1,qend) + U(2,2,nstep)*B(2,qend)
          qt = (c-cm1)/(cm1-cm2)
          At = qt*err-qt*(1.+qt)*errm1+qt**2*errm2
          Bt = (2.*qt+1.)*err-(1.+qt)**2*errm1+qt**2*errm2
          Ct = (1.+qt)*err
          if ( ABS(Bt+SQRT(Bt**2-4.*At*Ct)) .gt.
     &         ABS(Bt-SQRT(Bt**2-4.*At*Ct)) )  then
             c = ctemp-(ctemp-cm1)*2.*Ct/(Bt+SQRT(Bt**2-4.*At*Ct))
          else
             c = ctemp-(ctemp-cm1)*2.*Ct/(Bt-SQRT(Bt**2-4.*At*Ct))
          end if
          cm2 = cm1
          cm1 = ctemp
          errm2 = errm1
          errm1 = err
        end if

        write (*,30) icount, real(c), aimag(c), real(err), aimag(err)
  30    format (1x,i4,2(e17.8,e17.8,3x))

        icount = icount + 1
      end do
c
c.... Output Eigenvalue
c      
      write (*,40) real(ctemp), aimag(ctemp), qend
  40  format (/,'Eigenvalue = ',e17.8,1x,e17.8,2x,i5/)
c
c     Second Pass to compute the eigenfunctions
c
      if (eigfun) then
        vmax = 0.0
        k = nstep
        do i = 1, n
          y(i,k) = v(i,k)
          do m = 1, n-r
            y(i,k) = y(i,k) + U(i,m,k)*B(m,q)
          end do
        end do
        do m = 1, n-r
          B(m,q-1) = 0.0
          do j = 1, n-r
            B(m,q-1) = B(m,q-1) + P(j,m,q)*(B(j,q)-omega(j,q))
          end do
        end do
        do k = nstep-1, 0, -1
          t = to + h*k
          if ( t .lt. tq(q-1) ) then
            q = q - 1
            do m = 1, n-r
              B(m,q-1) = 0.0
              do j = 1, n-r
                B(m,q-1) = B(m,q-1) + P(j,m,q)*(B(j,q)-omega(j,q))
              end do
            end do
          end if
          do i = 1, n
            y(i,k) = v(i,k)
            do m = 1, n-r
              y(i,k) = y(i,k) + U(i,m,k)*B(m,q-1)
            end do
          end do
          if ( ABS(y(i,k)) .gt. ABS(vmax) ) then
            vmax = y(i,k)
          end if
        end do
c
c       Output eigenfunction and derivatives: u, u', u'', u'''
c
        do k = 0, nstep
          t = to + h*k
          write (11,20) t, REAL(y(1,k)/vmax),AIMAG(y(1,k)/vmax)
          write (12,20) t, REAL(y(2,k)/vmax),AIMAG(y(2,k)/vmax)
          write (13,20) t, REAL(y(3,k)/vmax),AIMAG(y(3,k)/vmax)
          write (14,20) t, REAL(y(4,k)/vmax),AIMAG(y(4,k)/vmax)
  20      format (1x,3(ES16.8E3,1x))
        end do
      end if

      return
      end

C***********************************************************************
      function INPROD(n, v1, v2)
C***********************************************************************
C
C     Perform and inner product on two comples vectors, v1 and v2
C
c***********************************************************************
      complex     v1(n), v2(n)
      complex     INPROD
      integer     n, i
#if USE_ANALYTIC_INPROD
c
c     Lets define a different complex inner product so that the
c     iteration is analytic
c
      INPROD = 0.0
      do i = 1, n
        INPROD = INPROD + v1(i)*v2(i)
      end do
#else
      INPROD = 0.0
      do i = 1, n
        INPROD = INPROD + v1(i)*conjg(v2(i))
      end do
#endif 
      return
      end
C***********************************************************************
      subroutine RK4(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step using fourth order (complex) Runge-Kutta
C
c***********************************************************************
      external    FUNC
      integer     neq
      real        to, h
      complex     yo(neq), yf(neq)
      complex     f(neq), k1(neq), k2(neq), k3(neq), k4(neq), q(neq)
      
      call FUNC(neq, yo, to, f)
      do j = 1 , 4
        k1(j) = h*f(j)
        q(j) = yo(j) + 0.5*k1(j)
      end do
      call FUNC(neq, q, to+0.5*h, f)
      do j = 1 , 4
        k2(j) = h*f(j)
        q(j) = yo(j) + 0.5*k2(j)
      end do
      call FUNC(neq, q, to+0.5*h, f)
      do j = 1 , 4
        k3(j) = h*f(j)
        q(j) = yo(j) + k3(j)
      end do
      call FUNC(neq, q, to+h, f)
      do j = 1 , 4
        k4(j) = h*f(j)
        yf(j) = yo(j)+k1(j)/6.+(k2(j)+k3(j))/3.+k4(j)/6.
      end do

      return
      end

C***********************************************************************
      subroutine FHOMO(neq,yo,t,yf)
C***********************************************************************
C
C     Function evaluation for the Orr-Sommerfeld equation
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (imax=10000)
      integer     n
      complex     alpha, c, rgamma
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax)

      common      /setup/ n, alpha, u, d2u, Re, ymin, ymax, h, 
     .                    rgamma, uspl, d2uspl, ydat
     
      common      /eig/   c
c***********************************************************************
      complex     yo(4), yf(4)
      real        t
      integer     neq
      real        UU, dUU, d2UU
c
c     Get the velocity field
c
#ifdef USE_SPLINE
      call SPEVAL(n+1,ydat,U,Uspl,t,UU)
      call SPEVAL(n+1,ydat,d2U,d2Uspl,t,d2UU)
#else
      call profile(t,UU,dUU,d2UU)
#endif
      do j = 1 , 3
        yf(j) = yo(j+1)
      end do
      yf(4) = 2.*alpha**2*yo(3)-alpha**4*yo(1) + 
     .        (0,1.)*alpha*Re*((UU-c)*(yo(3)-alpha**2*yo(1))-
     .        d2UU*yo(1)) 

      return
      end

C***********************************************************************
      subroutine FPART(neq, yo,t,yf)
C***********************************************************************
C
C     Function evaluation for the Orr-Sommerfeld equation
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (imax=10000)
      integer     n
      complex     alpha, c, rgamma
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax)

      common      /setup/ n, alpha, u, d2u, Re, ymin, ymax, h, 
     .                    rgamma, uspl, d2uspl, ydat
     
      common      /eig/   c
c***********************************************************************
      complex     yo(4), yf(4)
      real        t
      integer     neq
      real        UU, dUU, d2UU
c
c     Get the velocity field
c
#ifdef USE_SPLINE
      call SPEVAL(n+1,ydat,U,Uspl,t,UU)
      call SPEVAL(n+1,ydat,d2U,d2Uspl,t,d2UU)
#else
      call PROFILE(t,UU,dUU,d2UU)
#endif

      do j = 1 , 3
        yf(j) = yo(j+1)
      end do
      yf(4) = 2.*alpha**2*yo(3)-alpha**4*yo(1) + 
     &        (0,1.)*alpha*Re*((UU-c)*(yo(3)-alpha**2*yo(1))-
     &        d2UU*yo(1)) 

      return
      end

C***********************************************************************
      subroutine SOLVE_BL
C***********************************************************************
C
C     Setup/solve for the velocity profile. In this case, the profile
C     for an incompressible channel flow is analytic
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (imax=10000)
      integer     n
      complex     alpha, c, rgamma
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax)

      common      /setup/ n, alpha, u, d2u, Re, ymin, ymax, h, 
     &                    rgamma, uspl, d2uspl, ydat
     
      common      /eig/   c
c***********************************************************************
      integer     i, j, k, p
      real        xi(3,imax), f(3), eta(3), y, pi, us, dus, d2us
      real        k1(3), k2(3), k3(3), k4(3), err, x2old, f1old

      pi = acos(-1.0)
#ifdef USE_SPLINE
c     Now assemble the velocity field
c
      do i = 0, n
        y = ymin + i*h
        u(i) = 1-y**2
      end do
c
c     Compute 2nd order finite difference approximation to 
c     second derivative
c
      do i = 1, n-1
        d2u(i) = (u(i+1)-2*u(i)+u(i-1))/(h)**2
      end do
      d2u(0) = (-u(4)+4*u(3)-5*u(2)+2*u(1))/(h)**2
      d2u(n) = (2*u(n)-5*u(n-1)+4*u(n-2)-u(n-3))/(h)**2
c
c     Need to interpolate the velocity profile to evaluate it at
c     arbitrary y
c
      do i = 0, n
        ydat(i) = ymin + i*h
      end do
      call SPLINE(n+1,ydat,u,uspl)
      call SPLINE(n+1,ydat,d2u,d2uspl)
c
c     Write profile on 2x resolve mesh
c
c     open (10, FILE='chan.dat', ACTION='WRITE', FORM='FORMATTED')
      do i = 0, 2*n
        y = ymin + i*h/2.0
        call SPEVAL(n+1,ydat,u,uspl,y,us)
        call SPEVAL(n+1,ydat,d2u,d2uspl,y,d2us)
        write (10,10) y, us, d2us
  10    format (1x,5(ES16.8E3,1x))
      end do
c     close(10)
#else
      do i = 0, n
        y = ymin + i*h
        call profile(y,u(i),dudy,d2u(i))
      end do
      do i = 0, 2*n
        y = ymin + i*h/2.0
        call profile(y,us,dus,d2us)
        write (10,10) y, us, d2us
  10    format (1x,5(ES16.8E3,1x))
      end do
#endif
      write (*,20)
  20  format (1x,'Channel velocity profile completed...',/)

      return
      end

C***********************************************************************
      subroutine profile(y,u,dudy,d2udy2)
C***********************************************************************
      real y, u, dudy, d2udy2
C***********************************************************************
      u      = 1.0 - y**2
      dudy   = -2.0*y
      d2udy2 = -2.0
      return
      end

C***********************************************************************
      subroutine set_profile(n,y,u,dudy,d2udy2)
C***********************************************************************
      integer     n
      real        y(0:n), u(0:n), dudy(0:n), d2udy2(0:n)
C***********************************************************************
      do i = 0, n
        call profile(y(i),u(i),dudy(i),d2udy2(i))
      end do
      return
      end

C***********************************************************************
      SUBROUTINE SPLINE (N,X,Y,FDP)      
C***********************************************************************
C-----THIS SUBROUTINE COMPUTES THE SECOND DERIVATIVES NEEDED 
C-----IN CUBIC SPLINE INTERPOLATION.  THE INPUT DATA ARE:    
C-----N = NUMBER OF DATA POINTS          
C-----X = ARRAY CONTAINING THE VALUES OF THE INDEPENDENT VARIABLE      
C-----    (ASSUMED TO BE IN ASCENDING ORDER)       
C-----Y = ARRAY CONTAINING THE VALUES OF THE FUNCTION AT THE 
C-----    DATA POINTS GIVEN IN THE X ARRAY         
C-----THE OUTPUT IS THE ARRAY FDP WHICH CONTAINS THE SECOND  
C-----DERIVATIVES OF THE INTERPOLATING CUBIC SPLINE.         
      DIMENSION X(N),Y(N),A(N),B(N),C(N),R(N),FDP(N)  
C-----COMPUTE THE COEFFICIENTS AND THE RHS OF THE EQUATIONS. 
C-----THIS ROUTINE USES THE CANTILEVER CONDITION.  THE PARAMETER       
C-----ALAMDA (LAMBDA) IS SET TO 1. BUT THIS CAN BE USER-MODIFIED.      
C-----A,B,C ARE THE THREE DIAGONALS OF THE TRIDIAGONAL SYSTEM;         
C-----R IS THE RIGHT HAND SIDE.  THESE ARE NOW ASSEMBLED.    
      ALAMDA = 1.    
      NM2 = N - 2    
      NM1 = N - 1    
      C(1) = X(2) - X(1)       
      DO 1 I=2,NM1   
      C(I) = X(I+1) - X(I)     
      A(I) = C(I-1)  
      B(I) = 2.*(A(I) + C(I))  
      R(I) = 6.*((Y(I+1) - Y(I))/C(I) - (Y(I) - Y(I-1))/C(I-1))        
    1 CONTINUE       
      B(2) = B(2) + ALAMDA * C(1)        
      B(NM1) = B(NM1) + ALAMDA * C(NM1)  
C-----AT THIS POINT WE COULD CALL A TRIDIAGONAL SOLVER SUBROUTINE      
C-----BUT THE NOTATION IS CLUMSY SO WE WILL SOLVE DIRECTLY.  THE       
C-----NEXT SECTION SOLVES THE SYSTEM WE HAVE JUST SET UP.    
      DO 2 I=3,NM1   
      T = A(I)/B(I-1)          
      B(I) = B(I) - T * C(I-1) 
      R(I) = R(I) - T * R(I-1) 
    2 CONTINUE       
      FDP(NM1) = R(NM1)/B(NM1) 
      DO 3 I=2,NM2   
      NMI = N - I    
      FDP(NMI) = (R(NMI) - C(NMI)*FDP(NMI+1))/B(NMI)         
    3 CONTINUE       
      FDP(1) = ALAMDA * FDP(2) 
      FDP(N) = ALAMDA * FDP(NM1)         
C-----WE NOW HAVE THE DESIRED DERIVATIVES SO WE RETURN TO THE          
C-----MAIN PROGRAM.  
      RETURN         
      END  
C***********************************************************************
      SUBROUTINE SPEVAL (N,X,Y,FDP,XX,F) 
C***********************************************************************
C-----THIS SUBROUTINE EVALUATES THE CUBIC SPLINE GIVEN       
C-----THE DERIVATIVE COMPUTED BY SUBROUTINE SPLINE.          
C-----THE INPUT PARAMETERS N,X,Y,FDP HAVE THE SAME 
C-----MEANING AS IN SPLINE.    
C-----XX = VALUE OF INDEPENDENT VARIABLE FOR WHICH 
C-----     AN INTERPOLATED VALUE IS REQUESTED      
C-----F =  THE INTERPOLATED RESULT       
      DIMENSION X(20),Y(20),FDP(20)      
C-----THE FIRST JOB IS TO FIND THE PROPER INTERVAL.          
      NM1 = N - 1    
      DO 1 I=1,NM1   
      IF (XX.LE.X(I+1)) GO TO 10         
    1 CONTINUE       
C-----NOW EVALUATE THE CUBIC   
   10 DXM = XX - X(I)          
      DXP = X(I+1) - XX        
      DEL = X(I+1) - X(I)      
      F = FDP(I)*DXP*(DXP*DXP/DEL - DEL)/6.        
     1   +FDP(I+1)*DXM*(DXM*DXM/DEL - DEL)/6.     
     2   +Y(I)*DXP/DEL + Y(I+1)*DXM/DEL 
      RETURN        
      END 
