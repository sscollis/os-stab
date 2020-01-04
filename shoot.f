C***********************************************************************
      subroutine CONTE(nstep, n, r, yo, yf, to, tf, c, eigfun,
     &                 FHOMO, FPART, INPROD)
C***********************************************************************
C
C     First order linear boundary value problem solver using Conte's
C     method.  Fourth order Runge-Kutta is used for time advancement
C
C     nstep  = number of integration steps between [t0,tf]
C     n      = number of ODEs to integrate
C     r      = number of particular solutions 
C     yo     = starting values
C     yf     = ending values
C     to     = Interval start value
C     tf     = Interval end value
C     c      = Current eigenvalue
C     eigfun = Whether to output eigenfunction at the end
C     FHOMO  = Subroutine for homogeneous equation
C     FPART  = Subroutin for particular equation
C     INPROD = Function to compute inner product
C
c***********************************************************************
      complex     c
      integer     i, m, q, r, s, mi, mj, qend, IPVT(n-r), icount
      real        t, tq(0:nstep), to, tf, h
      complex     yo(n), yf(n), B(n-r,0:nstep), err, cold, ctemp
      complex     U(n,n-r,0:nstep), P(n-r,n-r,0:nstep), z(n,n-r)
      complex     y(n,0:nstep), v(n,0:nstep), w(n), omega(n,0:nstep)
      complex     eta(n),A(n,n-r),x(n),FAC(n-r,n-r), det1, det
      complex     olderr, ut(n,n-r), fd, vmax

      complex     cm1, cm2, errm1, errm2, qt, At, Bt, Ct, Utemp(n)

      real        test, testalpha, pi, det2, AI(n,n-r), angle
      real        aa, bb, cc, fdxr, fdxi, fdyr, fdyi
      logical     norm, eigfun 
      complex     INPROD
      external    INPROD, FHOMO, FPART

      real        errtol, maxcount, eigtol
c
c     initialize variables
c
      errtol = 1.0e-15
      maxcount = 20
      eigtol = 1.0e-15
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
#ifdef USE_RKCK45
            call CRKCK45(n, U(1,m,k-1), U(1,m,k), t-h, h, FHOMO)
#else
            call CRK4(n, U(1,m,k-1), U(1,m,k), t-h, h, FHOMO)
#endif
          end do
c
c         Now the particular solution
c
#ifdef USE_RKCK45
          call CRKCK45(n, v(1,k-1), v(1,k), t-h, h, FPART)
#else
          call CRK4(n, v(1,k-1), v(1,k), t-h, h, FPART)
#endif
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
      write (*,45) icount, abs(err), abs(c-cm1)
  45  format ('  Eigenvalue iterations = ', i5/,'  |error| = ',
     &        es17.8/,'  |c-cm1| = ', es17.8/)
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
      subroutine CONTE_IC(nstep, testalpha, n, r, yo, yf, to, tf, c, 
     &                    eigfun, FIC, FHOMO, FPART, INPROD)
C***********************************************************************
C
C     First order linear boundary value problem solver using Conte's
C     method with a non Dirac function initial condition.  
C     Runge-Kutta is used for time advancement but there are
C     other options available via compile defines
C
C     nstep  = number of integration steps between [t0,tf]
C     n      = number of ODEs to integrate
C     r      = number of particular solutions 
C     yo     = starting values
C     yf     = ending values
C     to     = Interval start value
C     tf     = Interval end value
C     c      = Current eigenvalue
C     eigfun = Whether to output eigenfunction at the end
C     FIC    = Subroutine to compute the initial condition
C     FHOMO  = Subroutine for homogeneous equation
C     FPART  = Subroutin for particular equation
C     INPROD = Function to compute inner product
C
C***********************************************************************
      integer     i, m, q, r, s, mi, mj, qend, IPVT(n-r), icount, north
      real        t, tq(0:nstep), to, tf, h 
      complex     yo(n), yf(n), B(n-r,0:nstep), err, c, cold, ctemp
      complex     U(n,n-r,0:nstep), P(n-r,n-r,0:nstep), z(n,n-r)
      complex     y(n,0:nstep), v(n,0:nstep), w(n-r), omega(n,0:nstep)
      complex     eta(n),A(n,n-r),x(n),FAC(n-r,n-r), det1, det
      complex     olderr, ut(n,n-r), fd, vmax
      complex     cm1, cm2, errm1, errm2, qt, At, Bt, Ct, Utemp(n)
      real        aa, bb, cc
      real        test, testalpha, pi, det2, AI(n,n-r)
      real        fdxr, fdxi, fdyr, fdyi
      logical     norm, eigfun

      complex     INPROD
      external    INPROD, FIC, FHOMO, FPART

      real        errtol, maxcount, eigtol
#ifdef USE_NR_ODEINT
      external    RKQC
#endif
c
c     Setup for ZVODE
c
#ifdef USE_VODE
      real RWORK(20+(n*(n-r)))
      complex ZWORK(15*(n*(n-r)))
      integer IWORK(30+(n*(n-r)))
      external FHOMO_VODE
      real RTOL, ATOL, RPAR
      integer ITOL, ITASK, ISTATE, IOPT, LZW, LRW, LIW, MF, IPAR(2)
      ITOL = 1
      RTOL = 1.0D-9
      ATOL = 1.0D-8
      ITASK = 1
      ISTATE = 1
      IOPT = 0
      LZW = 15*(n*(n-r))
      LRW = 20+(n*(n-r))
      LIW = 30+(n*(n-r))
      MF = 10
      RPAR = 0.0
      IPAR(1) = n 
      IPAR(2) = n-r
#endif
c
c     initialize variables
c
      errtol = 1.0e-14
      maxcount = 20
      eigtol = 1.0e-14
c
c     set the normalization constraint
c
      pi = ACOS(-1.0)
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
      do while ( ((abs(err) .ge. errtol) .or. (abs(c-cm1) .ge. eigtol)) 
     &           .and. (icount .le. maxcount) ) 
        q = 0
        tq(0) = tf
c
c       Set the initial conditions
c
        k = 0
        call FIC(n,r,U(1,1,k))
c
c       Gram-Schmidt
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
c       Now update the U matrix with the orthonormal values
c       
        do i = 1, n
          do m = 1, n-r
            U(i,m,k) = z(i,m)
          end do
        end do
#ifdef CHECK_ORTHOGONALITY
        aa = ABS(inprod(n, U(1,1,k), U(1,1,k)))
        bb = ABS(inprod(n, U(1,2,k), U(1,2,k)))
        cc = ABS(inprod(n, U(1,1,k), U(1,2,k)))
        test = ACOS(cc/SQRT(aa*bb))*180./pi
        write (*,*) k, mi, mj, test
#endif
c
c       Integrate the homo. and particular equations
c      
        do k = 1, nstep
          t = tf - h*k
c
c         Loop thru all homogeneous solutions
c
          do m = 1, n-r
#ifdef USE_NR_ODEINT
            do i = 1, n
              Utemp(i) = U(i,m,k-1)
            end do
            call ODEINT(Utemp,n,t+h,t,1.E-6,-h,1.e-20,nok,nbad,
     &                  FHOMO,RKQC)
            do i = 1, n
              U(i,m,k) = Utemp(i)
            end do
c#elif USE_VODE
c            ISTATE=1
c            rt0 = t+h
c            rt1 = t
c            do i = 1, n
c              Utemp(i) = U(i,m,k-1)
c            end do
cc           write(*,*) "Calling ZVODE...", rt0, rt1
c            call ZVODE(FHOMO_VODE,n,Utemp,rt0,rt1,ITOL,RTOL,ATOL,ITASK,
c     &                 ISTATE,IOPT,ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,MF,
c     &                 RPAR,IPAR)
c           write(*,*) "Finixhed ZVODE...", rt0, rt1
c            if (ISTATE .lt. 0) then
c              write(*,*) "Error in ZVODE:  ISTATE = ", ISTATE
c              call exit(1)
c            end if
c            do i = 1, n
c              U(i,m,k) = Utemp(i)
c            end do
#elif USE_LSRK14
            call CLSRK14(n, U(1,m,k-1), U(1,m,k), t+h, -h, FHOMO)
#elif USE_RKCK45
            call CRKCK45( n, U(1,m,k-1), U(1,m,k), t+h, -h, FHOMO )
#else
            call CRK4(n, U(1,m,k-1), U(1,m,k), t+h, -h, FHOMO)
#endif
          end do
#if 0
          if (norm .or. k.eq.1) ISTATE=1
          ISTATE=1
          rt0 = t+h
          rt1 = t
          do m = 1, n-r
            do i = 1, n
              U(i,m,k) = U(i,m,k-1)
            end do
          end do
c         write(*,*) "Calling ZVODE...", rt0, rt1, n*(n-r)
          call ZVODE(FHOMO_VODE,n*(n-r),U(1,1,k),rt0,rt1,ITOL,RTOL,ATOL,
     &               ITASK,ISTATE,IOPT,ZWORK,LZW,RWORK,LRW,IWORK,LIW,
     &               JEX,MF,RPAR,IPAR)
          if (ISTATE .lt. 0) then
            write(*,*) "Error in ZVODE:  ISTATE = ", ISTATE
            call exit(1)
          end if
#else
c         write(*,*) "Hack...fix it"
c         stop
#endif
c
c         Test to see if normalization is required
c
          norm = .false.
          ta = pi*testalpha/180.0
          do mi = 1, n-r
            do mj = 1, n-r
              if (mi .ne. mj) then
                aa = ABS(inprod(n, U(1,mi,k), U(1,mi,k)))
                bb = ABS(inprod(n, U(1,mj,k), U(1,mj,k)))
                cc = ABS(inprod(n, U(1,mi,k), U(1,mj,k)))
                test = acos(cc/SQRT(aa*bb))
                if (test .lt. ta) norm = .true.
c               write (*,*) k, mi, mj, test, testalpha, norm
              end if
            end do
          end do
c
c         Perform normalization
c       
          if ( norm .or. (k .eq. nstep) ) then
            q = q + 1
c           write (*,*) 'Normalization = ', q
            tq(q) = t
            if (k .eq. nstep) then
              qend = q
c             write (*,*) 'qend = ',qend
            end if
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
                  eta(i) = eta(i) - inprod(n, U(1,mi,k), 
     &                     z(1,mj))*z(i,mj)
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
  
          end if
        end do
c
c       check boundary conditions
c
        if (.true.) then
c
c         strictly enforce the zero BC
c
          B(1,qend) = -U(1,2,nstep)/U(1,1,nstep)
          B(2,qend) = 1.0
c         olderr = err
c         err = U(2,1,nstep)*B(1,qend) + U(2,2,nstep)*B(2,qend)
        else
c
c         strictly enforce the zero slope BC
c
          B(1,qend) =  1.0
          B(2,qend) = -U(2,1,nstep)/U(2,2,nstep)
c         olderr = err
c         err = U(1,1,nstep)*B(1,qend) + U(1,2,nstep)*B(2,qend)
        end if
        
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
        write (*,30) icount,real(ctemp),aimag(ctemp),real(err),
     &               aimag(err)
  30    format (1x,i4,2(e17.8,e17.8,3x))
        icount = icount + 1
      end do
c
c.... Output Eigenvalue
c      
      icount = icount - 1
      write (*,40) real(ctemp), aimag(ctemp), qend
  40  format (/,'Eigenvalue = ',e17.8,1x,e17.8,2x,i5/)
      write (*,45) icount, abs(err), abs(c-cm1)
  45  format ('  Eigenvalue iterations = ', i5/,'  |error| = ',
     &        es17.8/,'  |c-cm1| = ', es17.8/)
c
c     Second Pass to compute the eigenfunctions
c
      if (eigfun) then
        q  = qend
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
            B(m,q-1) = B(m,q-1) + P(j,m,q)*B(j,q) 
          end do
        end do
        do k = nstep-1, 0, -1
          t = tf - h*k
          if ( t .gt. tq(q-1) ) then
            q = q - 1
            do m = 1, n-r
              B(m,q-1) = 0.0
              do j = 1, n-r
                B(m,q-1) = B(m,q-1) + P(j,m,q)*B(j,q)
              end do
            end do
          end if
          do i = 1, n
            y(i,k) = 0.0
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
          t = tf - h*k
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
      subroutine SRK4(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step using fourth order (real) Runge-Kutta
C
c***********************************************************************
      external    FUNC
      integer     neq
      real        to, h
      real        yo(neq), yf(neq)
      real        f(neq), k1(neq), k2(neq), k3(neq), k4(neq), q(neq)
      
      call FUNC(neq, yo, to, f)
      do j = 1 , neq
        k1(j) = h*f(j)
        q(j) = yo(j) + 0.5*k1(j)
      end do
      call FUNC(neq, q, to+0.5*h, f)
      do j = 1 , neq
        k2(j) = h*f(j)
        q(j) = yo(j) + 0.5*k2(j)
      end do
      call FUNC(neq, q, to+0.5*h, f)
      do j = 1 , neq
        k3(j) = h*f(j)
        q(j) = yo(j) + k3(j)
      end do
      call FUNC(neq, q, to+h, f)
      do j = 1 , neq
        k4(j) = h*f(j)
        yf(j) = yo(j)+k1(j)/6.+(k2(j)+k3(j))/3.+k4(j)/6.
      end do

      return
      end

C***********************************************************************
      subroutine CRK4(neq, yo, yf, to, h, FUNC)
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
      do j = 1 , neq
        k1(j) = h*f(j)
        q(j) = yo(j) + 0.5*k1(j)
      end do
      call FUNC(neq, q, to+0.5*h, f)
      do j = 1 , neq
        k2(j) = h*f(j)
        q(j) = yo(j) + 0.5*k2(j)
      end do
      call FUNC(neq, q, to+0.5*h, f)
      do j = 1 , neq
        k3(j) = h*f(j)
        q(j) = yo(j) + k3(j)
      end do
      call FUNC(neq, q, to+h, f)
      do j = 1 , neq
        k4(j) = h*f(j)
        yf(j) = yo(j)+k1(j)/6.+(k2(j)+k3(j))/3.+k4(j)/6.
      end do

      return
      end

C***********************************************************************
      SUBROUTINE SPLINE (N,X,Y,FDP)      
C***********************************************************************
C
C.... Note: this routine is in the public domain and available
C     at https://web.stanford.edu/class/me200c/
C
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
C
C.... Note: this routine is in the public domain and available
C     at https://web.stanford.edu/class/me200c/
C
C-----THIS SUBROUTINE EVALUATES THE CUBIC SPLINE GIVEN       
C-----THE 2ND DERIVATIVE COMPUTED BY SUBROUTINE SPLINE.          
C-----THE INPUT PARAMETERS N,X,Y,FDP HAVE THE SAME 
C-----MEANING AS IN SPLINE.    
C-----XX = VALUE OF INDEPENDENT VARIABLE FOR WHICH 
C-----     AN INTERPOLATED VALUE IS REQUESTED      
C-----F =  THE INTERPOLATED RESULT       
      DIMENSION X(N),Y(N),FDP(N)      
C-----THE FIRST JOB IS TO FIND THE PROPER INTERVAL.          
#if USE_NR_HUNT
c
c     Search using bisection with a good guess
c
      I = IOLD
      IF (XX.EQ.X(1)) THEN
        I = 1
      ELSE IF (XX.EQ.X(N)) THEN
        I = N
      ELSE
        call HUNT (X,N,XX,I)
      END IF
      IOLD = I
#elif 1
      I = IOLD
      IF (XX.EQ.X(1)) THEN
        I = 1
      ELSE IF (XX.EQ.X(N)) THEN
        I = N
      ELSE
        call BISECT (X,N,XX,I)
      ENDiF
      IOLD = I
#else
c
c     This is really a slow way of searching
c
      NM1 = N - 1
      DO 1 I=1,NM1
      IF (XX.LE.X(I+1)) GO TO 10
    1 CONTINUE   
#endif
C-----NOW EVALUATE THE CUBIC   
   10 DXM = XX - X(I)          
      DXP = X(I+1) - XX        
      DEL = X(I+1) - X(I)      
      F = FDP(I)*DXP*(DXP*DXP/DEL - DEL)/6.        
     1   +FDP(I+1)*DXM*(DXM*DXM/DEL - DEL)/6.     
     2   +Y(I)*DXP/DEL + Y(I+1)*DXM/DEL 
      RETURN        
      END 

C***********************************************************************
      subroutine BISECT(X,N,XX,I)
C***********************************************************************
      dimension X(N)
C***********************************************************************
      il = I-1
      ir = N-1
      do while (ir-il .gt. 1) 
        im = ISHFT(ir+il,-1) 
        if ( X(im+1) > xx ) then
          ir = im
        else
          il = im
        end if
      end do
      I = il+1
      return
      end

C***********************************************************************
      SUBROUTINE SPDER(N,X,Y,FDP,XX,F,FP,FPP)
C***********************************************************************
C
C.... Note: this routine is in the public domain and available
C     at https://web.stanford.edu/class/me200c/
C
C-----THIS SUBROUTINE EVALUATES THE CUBIC SPLINE GIVEN       
C-----THE 2ND DERIVATIVE COMPUTED BY SUBROUTINE SPLINE.          
C-----THE INPUT PARAMETERS N,X,Y,FDP HAVE THE SAME 
C-----MEANING AS IN SPLINE.    
C-----XX = VALUE OF INDEPENDENT VARIABLE FOR WHICH 
C-----     AN INTERPOLATED VALUE IS REQUESTED      
C-----F =  THE INTERPOLATED RESULT       
C-----FP = THE INTERPOLATED DERIVATIVE RESULT       
      DIMENSION X(N),Y(N),FDP(N)
C-----THE FIRST JOB IS TO FIND THE PROPER INTERVAL.          
#if USE_NR_HUNT
c
c     Search using bisection with a good guess
c
      I = IOLD
      IF (XX.EQ.X(1)) THEN
        I = 1
      ELSE IF (XX.EQ.X(N)) THEN
        I = N
      ELSE
        call HUNT (X,N,XX,I)
      END IF
      IOLD = I
c
c     Try out a bubble sort
c
#elif 1
      I = IOLD
      IF (XX.EQ.X(1)) THEN
        I = 1
      ELSE IF (XX.EQ.X(N)) THEN
        I = N
      ELSE
        call BISECT (X,N,XX,I)
      ENDiF
      IOLD = I
#else
c
c     This is really a slow way of searching
c
      NM1 = N - 1
      DO 1 I=1,NM1
      IF (XX.LE.X(I+1)) GO TO 10
    1 CONTINUE   
#endif
C-----NOW EVALUATE THE CUBIC   
   10 continue
C     write(*,*) I, X(I), XX, X(I+1)
      DXM = XX - X(I)
      DXP = X(I+1) - XX
      DEL = X(I+1) - X(I)
      F = FDP(I)*DXP*(DXP*DXP/DEL - DEL)/6.0
     1   +FDP(I+1)*DXM*(DXM*DXM/DEL - DEL)/6.0
     2   +Y(I)*DXP/DEL + Y(I+1)*DXM/DEL
      FP= FDP(I)*(-3.0*DXP*DXP/DEL + DEL)/6.0
     1   +FDP(I+1)*(3.0*DXM*DXM/DEL - DEL)/6.0
     2   -Y(I)/DEL + Y(I+1)/DEL
      FPP=FDP(I)*DXP/DEL+FDP(I+1)*DXM/DEL
      RETURN     
      END

C***********************************************************************
      subroutine SLSRK14(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step with low-storage Runge Kutta 14
C
C***********************************************************************
      external FUNC
      integer  neq
      real     to, h
      real     yo(neq), yf(neq)
      real     f(neq), yt(neq)
C***********************************************************************
      parameter (nsubstep=14)
      real A(0:14), B(0:13), c(0:13), w(0:14)
      data A / 0.0, 
     &         -0.7188012108672410, 
     &         -0.7785331173421570,
     &         -0.0053282796654044, 
     &         -0.8552979934029281, 
     &         -3.9564138245774565, 
     &         -1.5780575380587385,
     &         -2.0837094552574054, 
     &         -0.7483334182761610,
     &         -0.7032861106563359,  
     &          0.0013917096117681,
     &         -0.0932075369637460, 
     &         -0.9514200470875948,
     &         -7.1151571693922548, 
     &         0.0/
      data B / 0.0367762454319673,
     &         0.3136296607553959,
     &         0.1531848691869027,
     &         0.0030097086818182,
     &         0.3326293790646110,
     &         0.2440251405350864,
     &         0.3718879239592277,
     &         0.6204126221582444,
     &         0.1524043173028741,
     &         0.0760894927419266,
     &         0.0077604214040978,
     &         0.0024647284755382,
     &         0.0780348340049386,
     &         5.5059777270269628 /
      data c / 0.0,
     &         0.0367762454319673,
     &         0.1249685262725025,
     &         0.2446177702277698,
     &         0.2476149531070420,
     &         0.2969311120382472,
     &         0.3978149645802642,
     &         0.5270854589440328,
     &         0.6981269994175695,
     &         0.8190890835352128,
     &         0.8527059887098624,
     &         0.8604711817462826,
     &         0.8627060376969976,
     &         0.8734213127600976 /
      data w / -0.116683473041717417,
     &         0.213493962104674251,
     &         0.128620987881127052,
     &         4.610096100109887907,
     &         -5.386527768056724064,
     &         1.445540684241274576,
     &         -0.761388932107154526,
     &         0.543874700576422732,
     &         0.102277834602298279,
     &         0.07127466608688701188,
     &         -3.459648919807762457,
     &         37.20095449534884580,
     &         -39.09786206496502814,
     &         5.505977727026962754,
     &         0.0 /
      do j = 1, neq
        yf(j) = yo(j)
      end do
      do i = 0, nsubstep-1
        t = to + c(i)*h
        call FUNC(neq, yf, t, f)
        do j = 1, neq
          yt(j) = A(i)*yt(j) + h*f(j)
        end do
        do j = 1, neq
          yf(j) = yf(j) + B(i)*yt(j)
        end do
        if (i+1 .lt. nsubstep) then
          t = to + c(i+1)*h
        else
          t = to + h
        end if
      end do
      return
      end

C***********************************************************************
      subroutine CLSRK14(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step with low-storage Runge Kutta 14
C
C***********************************************************************
      external FUNC
      integer  neq
      real     to, h
      complex  yo(neq), yf(neq)
      complex  f(neq), yt(neq)
C***********************************************************************
      parameter (nsubstep=14)
      real A(0:14), B(0:13), c(0:13), w(0:14)
      data A / 0.0, 
     &         -0.7188012108672410, 
     &         -0.7785331173421570,
     &         -0.0053282796654044, 
     &         -0.8552979934029281, 
     &         -3.9564138245774565, 
     &         -1.5780575380587385,
     &         -2.0837094552574054, 
     &         -0.7483334182761610,
     &         -0.7032861106563359,  
     &          0.0013917096117681,
     &         -0.0932075369637460, 
     &         -0.9514200470875948,
     &         -7.1151571693922548, 
     &         0.0/
      data B / 0.0367762454319673,
     &         0.3136296607553959,
     &         0.1531848691869027,
     &         0.0030097086818182,
     &         0.3326293790646110,
     &         0.2440251405350864,
     &         0.3718879239592277,
     &         0.6204126221582444,
     &         0.1524043173028741,
     &         0.0760894927419266,
     &         0.0077604214040978,
     &         0.0024647284755382,
     &         0.0780348340049386,
     &         5.5059777270269628 /
      data c / 0.0,
     &         0.0367762454319673,
     &         0.1249685262725025,
     &         0.2446177702277698,
     &         0.2476149531070420,
     &         0.2969311120382472,
     &         0.3978149645802642,
     &         0.5270854589440328,
     &         0.6981269994175695,
     &         0.8190890835352128,
     &         0.8527059887098624,
     &         0.8604711817462826,
     &         0.8627060376969976,
     &         0.8734213127600976 /
      data w / -0.116683473041717417,
     &         0.213493962104674251,
     &         0.128620987881127052,
     &         4.610096100109887907,
     &         -5.386527768056724064,
     &         1.445540684241274576,
     &         -0.761388932107154526,
     &         0.543874700576422732,
     &         0.102277834602298279,
     &         0.07127466608688701188,
     &         -3.459648919807762457,
     &         37.20095449534884580,
     &         -39.09786206496502814,
     &         5.505977727026962754,
     &         0.0 /
      do j = 1, neq
        yf(j) = yo(j)
      end do
      do i = 0, nsubstep-1
        t = to + c(i)*h
        call FUNC(neq, yf, t, f)
        do j = 1, neq
          yt(j) = A(i)*yt(j) + h*f(j)
        end do
        do j = 1, neq
          yf(j) = yf(j) + B(i)*yt(j)
        end do
        if (i+1 .lt. nsubstep) then
          t = to + c(i+1)*h
        else
          t = to + h
        end if
      end do
      return
      end

C***********************************************************************
      subroutine SRKCK45(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step Runge-Kutta Cash-Karp method
C
C***********************************************************************
      external FUNC
      integer  neq
      real     to, h, t
      real     yo(neq), yf(neq), yt(neq)
      real     yk(neq,6), ye(neq)
C***********************************************************************
      real b(6,5)
      real a(6), c(6), d(6)
      data a / 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 /
      data b / 0.0, 0.2, 0.075, 0.3, -0.2037037037037037,
     &         0.029495804398148147,
     &         0.0, 0.0, 0.225, -0.9, 2.5, 0.341796875,
     &         0.0, 0.0, 0.0, 1.2, -2.5925925925925926,
     &         0.041594328703703706,
     &         0.0, 0.0, 0.0, 0.0, 1.2962962962962963,
     &         0.40034541377314814,
     &         0.0, 0.0, 0.0, 0.0, 0.0, 0.061767578125 /
      data c / 0.09788359788359788, 0.0, 0.4025764895330113,
     &         0.21043771043771045, 0.0, 0.2891022021456804 /
      data d / -0.004293774801587311, 0.0, 0.018668586093857853,
     &         -0.034155026830808066, -0.019321986607142856,
     &         0.03910220214568039 /
c
c     Test data
c
#ifdef OS_DEBUG
      do i = 1, 6
        do j = 1, 5
          write(*,*) i, j, b(i,j)
        end do
      end do
      stop
#endif
c
c     Stage 1 - 6
c
      do m = 1, 6
        t = to + a(m)*h
        do n = 1, neq
          yt(n) = yo(n)
        end do
        do k = 1, m-1
          do n = 1, neq
            yt(n) = yt(n) + b(m,k)*yk(n,k)
          end do
        end do
        call FUNC(neq, yt, t, yk(1,m))
        do n = 1, neq
          yk(n,m) = h * yk(n,m)
        end do
      end do
c
c     Final solution and error
c
      do n = 1, neq
        yf(n) = yo(n)
        ye(n) = 0.0
      end do
      do k = 1, 6
        do n = 1, neq
          yf(n) = yf(n) + c(k)*yk(n,k)
          ye(n) = ye(n) + d(k)*yk(n,k)
        end do
      end do

      return
      end

C***********************************************************************
      subroutine CRKCK45(neq, yo, yf, to, h, FUNC)
C***********************************************************************
C
C     Advance one time step Runge-Kutta Cash-Karp method
C
C***********************************************************************
      external FUNC
      integer  neq
      real     to, h, t
      complex  yo(neq), yf(neq), yt(neq)
      complex  yk(neq,6), ye(neq)
C***********************************************************************
      real b(6,5)
      real a(6), c(6), d(6)
      data a / 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 /
      data b / 0.0, 0.2, 0.075, 0.3, -0.2037037037037037, 
     &         0.029495804398148147,
     &         0.0, 0.0, 0.225, -0.9, 2.5, 0.341796875,
     &         0.0, 0.0, 0.0, 1.2, -2.5925925925925926, 
     &         0.041594328703703706,
     &         0.0, 0.0, 0.0, 0.0, 1.2962962962962963,
     &         0.40034541377314814,
     &         0.0, 0.0, 0.0, 0.0, 0.0, 0.061767578125 /
      data c / 0.09788359788359788, 0.0, 0.4025764895330113,
     &         0.21043771043771045, 0.0, 0.2891022021456804 /
      data d / -0.004293774801587311, 0.0, 0.018668586093857853,
     &         -0.034155026830808066, -0.019321986607142856,
     &         0.03910220214568039 / 
c
c     Test data
c
#ifdef OS_DEBUG
      do i = 1, 6
        do j = 1, 5
          write(*,*) i, j, b(i,j)
        end do
      end do
      stop
#endif
c
c     Stage 1 - 6
c
      do m = 1, 6 
        t = to + a(m)*h
        do n = 1, neq
          yt(n) = yo(n)
        end do
        do k = 1, m-1
          do n = 1, neq
            yt(n) = yt(n) + b(m,k)*yk(n,k)
          end do
        end do 
        call FUNC(neq, yt, t, yk(1,m))
        do n = 1, neq
          yk(n,m) = h * yk(n,m)
        end do
      end do
c
c     Final solution and error
c
      do n = 1, neq
        yf(n) = yo(n)
        ye(n) = 0.0
      end do
      do k = 1, 6
        do n = 1, neq
          yf(n) = yf(n) + c(k)*yk(n,k)
          ye(n) = ye(n) + d(k)*yk(n,k)
        end do
      end do

      return
      end 
