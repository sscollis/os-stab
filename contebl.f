c***********************************************************************
      program conte_bl
c***********************************************************************
c 
c     Purpose:  Solve the Orr-Sommerfeld equation using 4th order Runge-
c               Kutta explicit integration with orthogonalization using
c               the method of Conte.
c
c               This routine finds a single eigensolution to the Orr-
c               Sommerfeld equation for the boundary layer.
c
c               This routine works rather well but to get good
c               eigenvalues you must really resolve the boundary-layer
c               profile well.  The OS shooting solution appears to be 
c               quite sensitive to the profile.
c
c     Output:   Complex valued eignenfunction is phi
c
c               Set eigfun=1 if you want eigenfunction
c
c               fort.10     Mean velocity profile
c               fort.11     phi(y)
c               fort.12     phi'(y)
c               fort.13     phi''(y)
c               fort.14     phi'''(y)
c
c     Author:   Scott Collis
c
c     Date:     7-17-1992
c
c     Revised:  12-29-2019 
c
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (imax=50000)
      integer     nbl
      complex     alpha, c, rgamma
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   rgamma, beta, uspl, d2uspl, ydat
     
      common      /eig/  c, alpha, Re
c***********************************************************************
      parameter (neq=4, inunit=20)
      complex   bc1(neq), bc2(neq)
      real      testalpha
      integer   north, nstep
      logical   eigfun
      character input*1
      character infile*20
      character string*80
c
c     Defaults
c
      nbl = 1000
      nstep = 1000
      north = 50
c
c     testalpha = [0,90] degrees (10 deg seems good)
c
      testalpha = 10.0 
      Re = 580.
      beta = 0.0
      alphar = 0.179
      alphai = 0.0
      cr =  0.36413124E+00
      ci =  0.79527387E-02
      ymin = 0.0
      ymax = 17.0
c
c     ymin = 20.0 leads to an underflow and denormal, 17 is okay
c
      eigfun = .true.
      f2p = 0.5
c
c     Input parameters
c
      write (*,10)
  10  format (/,10x,'Solve Orr-Sommerfeld for Boundary Layer using ',
     &        'shooting',/)
      write (*,15)
  15  format(/,1x,'Read from keyboard, file, default (k,f,d) ==> ',$)
      read (*,'(a)') input
      if (input .eq. 'k' .or. input .eq. 'K') then
        write (*,20)
  20    format (/,1x,'Enter number of steps for BL solve ==> ',$)
        read (*,*) nbl
        write (*,25)
  25    format (/,1x,'Enter number of steps for shooting ==> ',$)
        read (*,*) nstep
        write (*,30)
  30    format (/,1x,'Enter Reynolds number ==> ',$)
        read (*,*) Re
        write (*,35)
  35    format (/,1x,'Enter beta ==> ',$)
        read (*,*) beta
        write (*,40) 
  40    format (/,1x,'Enter alpha_r ==> ',$)
        read (*,*) alphar
        write (*,45) 
  45    format (/,1x,'Enter alpha_i ==> ',$)
        read (*,*) alphai
        write (*,50)
  50    format (/,1x,'Enter guess for (c_r, c_i) ==> ',$)
        read (*,*) cr, ci
        write (*,55)
  55    format (/,1x,'Enter Ymin and Ymax ==> ',$)
        read (*,*) Ymin, Ymax
      else if (input .eq. 'f' .or. input .eq. 'F') then
        write (*,60)
  60    format (/,1x,'Enter filename ==> ',$)
        read (*,'(a)') infile
        open (unit=inunit, file=infile, form='formatted')
        string(1:1) = '#'
        do while (string(1:1) .eq. '#')
          read(inunit,'(a)',err=99) string
        end do
        call readi (string,nbl)
        read(inunit,'(a)') string
        call readi (string,nstep)
        read(inunit,'(a)') string
        call readr (string,testalpha)
        read(inunit,'(a)') string
        call readr (string,Re)
        read(inunit,'(a)') string
        call readr (string,beta)
        read(inunit,'(a)') string
        call readr (string,f2p)
        read(inunit,'(a)') string
        call readr (string,alphar)
        read(inunit,'(a)') string
        call readr (string,alphai)
        read(inunit,'(a)') string
        call readr (string,cr)
        read(inunit,'(a)') string
        call readr (string,ci)
        read(inunit,'(a)') string
        call readr (string,Ymin)
        read(inunit,'(a)') string
        call readr (string,Ymax)
      end if
c
c     Echo input
c
      write (*,100) nbl
 100  format (/,1x,'nbl = ',i5)
      write (*,102) nstep
 102  format (1x,'nstep = ',i5)
      write (*,105) testalpha
 105  format (1x,'Test Alpha = ',e20.10)
      write (*,110) Re
 110  format (1x,'Re = ',e20.10)
      write (*,115) beta
 115  format (1x,'Beta = ',e20.10)
      write (*,117) f2p
 117  format (1x,'f2p = ',e20.10)
      write (*,140) ymin, ymax
 140  format (1x,'Ymin = ',e20.10,'  Ymax = ',e20.10)
      write (*,120) alphar, alphai
 120  format (1x,'alpha = (',e20.10,', ',e20.10,')')
      write (*,130) cr, ci
 130  format (1x,'c = (',e20.10,', ',e20.10,')',/)
c
c     Set constants
c
      alpha = cmplx(alphar,alphai)
      c = cmplx(cr,ci)
      h = (ymax-ymin)/float(nbl)
c
c     Fix to make my nondimensionalization match Mack's
c
      Re = Re*SQRT(2.)
      alpha = alpha*SQRT(2.)
      ymin = ymin
      ymax = ymax
c
c     Check problem size
c      
      if (nbl .gt. imax) then
        write (*,300) 
 300    format (/,/,1x,'N > Idim...Stopping',/,/)
        goto 210
      end if
c
c     Solve the Blasius boundary layer equations
c 
      call SOLVE_BL
c
c     set the boundary conditions
c
      do i = 1, neq
        bc1(i) = 0.0
        bc2(i) = 0.0
      end do
      
      call CONTE(nstep,testalpha,neq,2,bc1,bc2,ymin,ymax,eigfun)
c
c     Read error
c      
      goto 210
  99  write (*,200) 
 200  format(/,/1x,'Error in input file...',/,/)
 210  continue

      stop
      end

C***********************************************************************
      subroutine CONTE(nstep, testalpha, n, r, yo, yf, to, tf, eigfun)
C***********************************************************************
C
C     First order linear boundary value problem solver using Conte's
C     method.  Fourth order Runge-Kutta is used for time advancement
C
c***********************************************************************
      complex     c, alpha
      real        Re
      common      /eig/ c, alpha, Re
C***********************************************************************
      integer     i, m, q, r, s, mi, mj, qend, IPVT(n-r), icount, north
      real        t, tq(0:nstep), to, tf, h 
      complex     yo(n), yf(n), B(n-r,0:nstep), err, cold, ctemp
      complex     U(n,n-r,0:nstep), P(n-r,n-r,0:nstep), z(n,n-r)
      complex     y(n,0:nstep), v(n,0:nstep), w(n-r), omega(n,0:nstep)
      complex     eta(n),A(n,n-r),x(n),FAC(n-r,n-r), det1, det
      complex     olderr, INPROD, ut(n,n-r), fd, vmax, rgamma
      complex     cm1, cm2, errm1, errm2, qt, At, Bt, Ct, Utemp(n)
      real        aa, bb, cc
      real        test, testalpha, pi, det2, AI(n,n-r)
      real        fdxr, fdxi, fdyr, fdyi
      logical     norm, eigfun
      external    INPROD, FHOMO, FPART, RKQC

      real        errtol, maxcount, eigtol
c
c     initialize variables
c
      errtol = 1.0e-14
      maxcount = 20
      eigtol = 1.0e-14
      eigfun = .true.
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
        U(1,1,0) = CEXP(-alpha*tf)
        U(2,1,0) = (-alpha)*CEXP(-alpha*tf)
        U(3,1,0) = (alpha**2)*CEXP(-alpha*tf)
        U(4,1,0) = (-alpha**3)*CEXP(-alpha*tf)
        
        rgamma = SQRT(alpha**2+(0.,1.)*alpha*Re*(1.0-c))
  
        U(1,2,0) = CEXP(-rgamma*tf)
        U(2,2,0) = (-rgamma)*CEXP(-rgamma*tf)
        U(3,2,0) = (rgamma**2)*CEXP(-rgamma*tf)
        U(4,2,0) = (-rgamma**3)*CEXP(-rgamma*tf)
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
#elif USE_LSRK14
            call CLSRK14(n, U(1,m,k-1), U(1,m,k), t+h, -h, FHOMO)
#else
            call CRK4(n, U(1,m,k-1), U(1,m,k), t+h, -h, FHOMO)
#endif
          end do
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
      subroutine FHOMO(neq, yo,t,yf)
C***********************************************************************
C
C     Function evaluation for the Orr-Sommerfeld equation
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (imax=50000)
      integer     n
      complex     alpha, c, rgamma
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   rgamma, beta, uspl, d2uspl, ydat
     
      common      /eig/  c, alpha, Re
c***********************************************************************
      integer     neq
      complex     yo(neq), yf(neq)
      real        t
      real        UU, d2UU
      real        Pi
      pi = acos(-1.0)
c
c     Get the velocity field
c
#ifdef USE_SPINE_DERIVATIVE
      call SPEVAL(nbl+1,ydat,U,Uspl,t,UU)
      call SPEVAL(nbl+1,ydat,d2U,d2Uspl,t,d2UU)
#else
      call SPDER(nbl+1,ydat,U,Uspl,t,UU,dUU,d2UU)
#endif
      do j = 1 , neq-1
        yf(j) = yo(j+1)
      end do
      yf(neq) = 2.*alpha**2*yo(3)-alpha**4*yo(1) + 
     &          (0.,1.)*alpha*Re*((UU-c)*(yo(3)-alpha**2*yo(1))-
     &          d2UU*yo(1)) 

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
      parameter   (imax=50000)
      integer     n
      complex     alpha, c, rgamma
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   rgamma, beta, uspl, d2uspl, ydat
     
      common      /eig/  c, alpha, Re
c***********************************************************************
      integer     neq
      complex     yo(neq), yf(neq)
      real        t
      real        UU, d2UU
      real        Pi
      pi = acos(-1.0)
c
c     Get the velocity field
c
#ifdef USE_SPINE_DERIVATIVE
      call SPEVAL(nbl+1,ydat,U,Uspl,t,UU)
      call SPEVAL(nbl+1,ydat,d2U,d2Uspl,t,d2UU)
#else
      call SPDER(nbl+1,ydat,U,Uspl,t,UU,dUU,d2UU)
#endif
      do j = 1 , neq-1
        yf(j) = yo(j+1)
      end do
      yf(neq) = 2.*alpha**2*yo(3)-alpha**4*yo(1) + 
     &          (0.,1.)*alpha*Re*((UU-c)*(yo(3)-alpha**2*yo(1))-
     &          d2UU*yo(1)) 

      return
      end

C***********************************************************************
      subroutine SOLVE_BL
C***********************************************************************
C
C     Integrate the boundary layer similarity equation to get velocity
C     profile.
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (imax=50000)
      integer     n
      complex     alpha, c, rgamma
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   rgamma, beta, uspl, d2uspl, ydat
     
      common      /eig/  c, alpha, Re
c***********************************************************************
      integer     i, j, k, p
      real        xi(3,0:imax), f(3), eta(3), y
      real        k1(3), k2(3), k3(3), k4(3), err, x2old, f1old
      external    BLASIUS, RKQCR, BL
#ifdef USE_RKF45
      real t, tout, work(imax), abserr, relerr
      integer iflag, iwork(5)
#endif
      do j = 1,3
        do i = 0,nbl
          xi(j,i) = 0.0
        end do
        f(j) = 0.0
        eta(j) = 0.0
      end do
c
c     Set the boundary conditions including guess for f"(0)
c      
      xi(1,0) = 0
      xi(2,0) = 0
      xi(3,0) = f2p
      err = 1.0
      p = 1
      do while ( abs(err) .gt. 1.e-10)

        do i = 1, nbl
          y = ymin + float(i)*h
#ifdef USE_NR_ODEINT
          do j = 1, 3
            eta(j) = xi(j,i-1)
          end do
          call ODEINTR(eta,3,y-h,y,1.E-7,h/2.,1.e-20,nok,nbad,
     &                 BLASIUS,RKQCR)
          do j = 1, 3
            xi(j,i) = eta(j)
          end do
#elif USE_RKF45
          t = y-h
          tout = y
          do j = 1, 3
            eta(j) = xi(j,i-1)
          end do
          relerr = 1.0e-9
          abserr = 0.0
          iflag = 1
          call RKF45(BL,3,eta,t,tout,relerr,abserr,iflag,work,iwork)
          if (iflag .ne. 2) call exit(1)
          do j = 1, 3
            xi(j,i) = eta(j)
          end do
#elif USE_LSRK14
          call SLSRK14(3, xi(1,i-1), xi(1,i), y-h, h, BLASIUS)
#else
          call SRK4(3, xi(1,i-1), xi(1,i), y-h, h, BLASIUS)
#endif
        end do
c
c       Check f'(ymax)
c      
        if (p .eq. 1) then
          xi3old = xi(3,0)
          xi(3,0) = xi(3,0)*.99
        else
          xi3temp = xi(3,0)
          xi(3,0) = xi(3,0)+((xi3old-xi(3,0))/(xi2old-xi(2,nbl)))*
     &              (1.0 - xi(2,nbl))
          xi3old = xi3temp
        end if
        p = p + 1
        xi2old = xi(2,nbl)
        err = 1.0 - xi2old
c       write (*,*) p, err
      end do
c
c     Now assemble the velocity field
c
      do i = 0, nbl
        u(i) = xi(2,i)
      end do
c
c     Compute 2nd order finite difference approximation to 
c     second derivative
c
      do i = 1, nbl-1
        d2u(i) = (u(i+1)-2*u(i)+u(i-1))/(h)**2
      end do
      d2u(0) = (-u(4)+4*u(3)-5*u(2)+2*u(1))/(h)**2
      d2u(nbl) = (2*u(nbl)-5*u(nbl-1)+4*u(nbl-2)-u(nbl-3))/(h)**2
c
c     Need to interpolate the velocity profile to evaluate it at
c     arbitrary y
c
      do i = 0, nbl
        ydat(i) = ymin + i*h
      end do
      call SPLINE(nbl+1,ydat,u,uspl)
c
c     Use the spline result for the second derivative instead of FD
c
      do i = 0, nbl
        d2u(i) = uspl(i)
      end do
      call SPLINE(nbl+1,ydat,d2u,d2uspl)

c     open (10, FILE='bl.dat', ACTION='WRITE', FORM='FORMATTED')
      do i = 0, nbl
        y = ymin + i*h
        call SPDER(nbl+1,ydat,u,uspl,y,us,dus,ddus)
        call SPEVAL(nbl+1,ydat,d2u,d2uspl,y,d2us)
        write (10,10) y, us, d2us, u(i), d2u(i)
c       write (10,10) y, us, dus, ddus, d2us, u(i), d2u(i)
  10    format (1x,7(ES16.8E3,1x))
      end do
c     close(10)
      
      write (*,20) 
  20  format (1x,'Blasius velocity profile completed...',/)

      return
      end

C***********************************************************************
      subroutine BL(y,xi,f)
C***********************************************************************
      real y, xi(3), f(3)
c     write(*,*) "In BL..."
      call BLASIUS(3,xi,y,f)
      return
      end

C***********************************************************************
      subroutine BLASIUS(neq,xi,y,f)
C***********************************************************************
C
C     Function evaluation for the boundary-layer similarity equation
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (imax=50000)
      integer     nbl
      complex     alpha, c, rgamma
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   rgamma, beta, uspl, d2uspl, ydat
     
      common      /eig/  c, alpha, Re
c***********************************************************************
      integer     neq
      real        xi(neq), f(neq), y

      do j = 1 , 2
        f(j) = xi(j+1)
      end do
      f(3) = -xi(1)*xi(3)-beta*(1.-xi(2)**2) 

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
      subroutine READI(string, I)
C***********************************************************************
C
C     Read an integer from an input file
C
c***********************************************************************

      integer   I, iloc, floc
      character string*80

      iloc = index (string,'=')
      if (iloc .ne. 0) then
        iloc = iloc + 1
        floc = index (string,'#')
        if (floc .eq. 0) then 
          floc = 80
        else
          floc = floc - 1
        end if
        read (string(iloc:floc),'(I10)',err=99) I
      else
        write (*,10)
  10    format (/,/,1x,'ERROR in input file...',/,/)
        stop
      end if
      goto 100
c
c     Error 
c
  99  continue
      write (*,10)
      stop

 100  return
      end

C***********************************************************************
      subroutine READR(string, R)
C***********************************************************************
C
C     Read an integer from an input file
C
c***********************************************************************
      integer    iloc, floc
      real       R
      character  string*80

      iloc = index (string,'=')
      if (iloc .ne. 0) then
        iloc = iloc + 1
        floc = index (string,'#')
        if (floc .eq. 0) then
          floc = 80
        else
          floc = floc - 1
        end if
        read (string(iloc:floc),'(g20.10)',err=99) R
      else
        write (*,10)
 10     format (/,/,1x,'ERROR in input file...',/,/)
        stop
      end if
      goto 100
c
c     Error 
c
  99  continue
      write (*,10)
      stop

 100  return
      end

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
