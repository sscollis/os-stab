c***********************************************************************
c> \file orrsom.f
c> \brief Solves the Orr-Sommerfeld equation for incompressible boundary
c>        layers solving first the Blasius equation 
c> \author S. Scott Collis
c***********************************************************************
      program OrrSom
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
c     Author:   Scott Collis
c
c     Date:     7-17-92
c
c     Revised:  10-29-92     changed tabs to spaces <-- big deal
c     Revised:  1-11-2020    code cleanup
c
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=20000)
      integer     n
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p

      common      /data/  n, u, d2u, ymin, ymax, h, f2p,
     &                    gamma, beta, uspl, d2uspl, ydat
     
      common      /eig/   c, alpha, Re
c***********************************************************************
      parameter    (neq=4)
      integer       inunit
      complex       bc1(neq), bc2(neq), omega
      real          testalpha
      integer       north, nstep
      logical       eigfun
      character*1   input
      character*20  infile
      character*80  string

      inunit = 20
      
      write (*,10)
  10  format (/,10x,'Solve Orr-Sommerfeld (Shooting)',/)
#ifdef USE_FILE_INPUT
      write (*,60)
  60  format (/,1x,'Enter filename ==> ',$)
      read (*,'(a)') infile
      open (unit=inunit, file=infile, form='formatted')
#else
      inunit = 5
#endif
      string(1:1) = '#'
      do while (string(1:1) .eq. '#')
        read(inunit,'(a)',err=99) string
      end do
      call readi (string,n)
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
c
c     Echo input
c
      write (*,100) n
 100  format (/,1x,'n = ',i5)
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
      h = (ymax-ymin)/float(n)
c
c     Fix to make my nondimensionalization match Mack's
c
#ifdef USE_MACK_NORMALIZATION
c
c     SSC: 1/7/2020 Why using the 1.7207877 ?
c
      Re = Re*SQRT(2.)/1.7207877
      alpha = alpha*SQRT(2.)/1.7207877
#else
      Re = Re*SQRT(2.)
      alpha = alpha*SQRT(2.)
#endif
      ymin = ymin
      ymax = ymax
      
      if (n .gt. idim) then
        write (*,300) 
 300    format (/,/,1x,'N > Idim...Stopping',/,/)
        goto 210
      end if
      
      call SOLVE_BL
c
c     set the boundary conditions
c
      do i = 1, neq
        bc1(i) = 0.0
        bc2(i) = 0.0
      end do
      
      eigfun = .true.
      call CONTE(nstep,testalpha,neq,2,bc1,bc2,ymin,ymax,eigfun)

      goto 210
c
c     Read error
c      
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

      common      /eig/   c, alpha, Re
C***********************************************************************
      integer     i, m, q, r, s, mi, mj, qend, IPVT(n-r), icount, north
      real        t, tq(0:nstep), to, tf, h 
      complex     yo(n), yf(n), B(n-r,0:nstep), err, cold, ctemp
      complex     U(n,n-r,0:nstep), P(n-r,n-r,0:nstep), z(n,n-r)
      complex     y(n,0:nstep), v(n,0:nstep), w(n-r), omega(n,0:nstep)
      complex     eta(n),A(n,n-r),x(n),FAC(n-r,n-r), det1, det
      complex     olderr, INPROD, ut(n,n-r), fd, max, gamma
      complex     cm1, cm2, errm1, errm2, qt, At, Bt, Ct, Utemp(n)
      real        aa, bb, cc
      real        test, testalpha, pi, det2, AI(n,n-r)
      real        fdxr, fdxi, fdyr, fdyi
      logical     norm, eigfun
      external    INPROD, FHOMO, FPART, RKQC
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
      err = 1.
      do while ((abs(err) .ge. 1.0e-8) .and. (icount .le. 20) .and.
     &          (abs(c-cm1) .ge. 1.0e-12) )
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
        
        gamma = SQRT(alpha**2+(0.,1.)*alpha*Re*(1.0-c))

        U(1,2,0) = CEXP(-gamma*tf)
        U(2,2,0) = (-gamma)*CEXP(-gamma*tf)
        U(3,2,0) = (gamma**2)*CEXP(-gamma*tf)
        U(4,2,0) = (-gamma**3)*CEXP(-gamma*tf)
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
c       aa = ABS(inprod(n, U(1,1,k), U(1,1,k)))
c       bb = ABS(inprod(n, U(1,2,k), U(1,2,k)))
c       cc = ABS(inprod(n, U(1,1,k), U(1,2,k)))
c       test = ACOS(cc/SQRT(aa*bb))*180./pi
c       write (*,*) k, mi, mj, test
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
            call ODEINT(Utemp,n,t+h,t,1.E-5,-h/2.,1.e-20,nok,nbad,
     &                  FHOMO,RKQC)
c           write (*,*) k, nok, nbad
            do i = 1, n
              U(i,m,k) = Utemp(i)
            end do
#else
            call CRK4(n, U(1,m,k-1), U(1,m,k), t+h, -h, FHOMO)
#endif
          end do
c
c         Test to see if normalization is required
c
          norm = .false.
          do mi = 1, n-r
            do mj = 1, n-r
              if (mi .ne. mj) then
                aa = ABS(inprod(n, U(1,mi,k), U(1,mi,k)))
                bb = ABS(inprod(n, U(1,mj,k), U(1,mj,k)))
                cc = ABS(inprod(n, U(1,mi,k), U(1,mj,k)))
#ifdef USE_ANALYTIC_INPROD
                test = 180.0*ACOS(MIN(1.0,ABS(cc)/SQRT(aa*bb)))/pi
#else
                test = 180.0*ACOS(ABS(cc)/SQRT(aa*bb))/pi
#endif
c               write (*,*) k, mi, mj, test
c		                
c               test = ACOS( ABS ( inprod(n, U(1,mi,k), U(1,mj,k)) / 
c      &               SQRT( inprod(n, U(1,mi,k), U(1,mi,k))*
c      &                     inprod(n, U(1,mj,k), U(1,mj,k)) ) ) )
c     
c               test = cc/SQRT(aa*bb)
c               if (test .gt. testalpha) norm = .true.
                if (test .le. testalpha) norm = .true.
              end if
            end do
          end do
c
c         Perform normalization
c	
          if ( norm .or. (k .eq. nstep) ) then
            q = q + 1
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
                  eta(i) = eta(i)-inprod(n,U(1,mi,k),z(1,mj))*z(i,mj)
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
     &              P(s,i,q)
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
      
      write (*,40) real(ctemp),aimag(ctemp),qend
  40  format (/,'Eigenvalue = ',e17.8,1x,e17.8,2x,i5/)

c
c     Second Pass
c

      if (eigfun) then
        write(*,50) 
  50    format('Second pass: computing and writing eigenfunction')
        max = 0.0
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
          if ( ABS(y(i,k)) .gt. ABS(max) ) then
            max = y(i,k)
          end if
        end do
        do k = 0, nstep
          t = tf - h*k
          write (11,20) t, REAL(y(1,k)/max),AIMAG(y(1,k)/max)
          write (12,20) t, REAL(y(2,k)/max),AIMAG(y(2,k)/max)
          write (13,20) t, REAL(y(3,k)/max),AIMAG(y(3,k)/max)
          write (14,20) t, REAL(y(4,k)/max),AIMAG(y(4,k)/max)
          write (15,20) t, REAL(y(2,k)/max),REAL((0.,-1.)*ALPHA*y(1,k)/max)
  20      format (1x,3(ES17.8E3,1x) )
        end do
      end if

      return
      end

C***********************************************************************
      function INPROD(n, v1, v2)
C***********************************************************************
C
C     Perform and inner product on two complex vectors, v1 and v2
C
c***********************************************************************
      complex     v1(n), v2(n)
      complex     INPROD
      integer     n, i
 
      INPROD = 0.0
      do i = 1, n
        INPROD = INPROD + v1(i)*conjg(v2(i))
c       INPROD = INPROD + v1(i)*v2(i)
      end do
      
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
      parameter   (idim=20000)
      integer     n
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p

      common      /data/  n, u, d2u, ymin, ymax, h, f2p,
     &                    gamma, beta, uspl, d2uspl, ydat
     
      common      /eig/   c, alpha, Re
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
      call SPEVAL(n+1,ydat,U,Uspl,t,UU)
      call SPEVAL(n+1,ydat,d2U,d2Uspl,t,d2UU)

      do j = 1 , neq-1
        yf(j) = yo(j+1)
      end do
      yf(neq) = (1./alpha/Re*(2.*alpha**2*yo(3)-alpha**4*yo(1)) + 
     &          (0.,1.)*((UU-c)*(yo(3)-alpha**2*yo(1))-
     &          d2UU*yo(1)))*alpha*Re 

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
      parameter   (idim=20000)
      integer     n
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p

      common      /data/  n, u, d2u, ymin, ymax, h, f2p,
     &                    gamma, beta, uspl, d2uspl, ydat
     
      common      /eig/   c, alpha, Re
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
      call SPEVAL(n+1,ydat,U,Uspl,t,UU)
      call SPEVAL(n+1,ydat,d2U,d2Uspl,t,d2UU)

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
      parameter   (idim=20000)
      integer     n
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p

      common      /data/  n, u, d2u, ymin, ymax, h, f2p,
     &                    gamma, beta, uspl, d2uspl, ydat
     
      common      /eig/   c, alpha, Re
c***********************************************************************
      integer     i, j, k, p
      real        xi(3,0:idim), f(3), eta(3), y
      real        k1(3), k2(3), k3(3), k4(3), err, x2old, f1old
      external    BLASIUS, RKQCR
      
      do j = 1,3
        do i = 0,n
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

        do i = 1, n
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
          xi(3,0) = xi(3,0)+((xi3old-xi(3,0))/(xi2old-xi(2,n)))*
     &              (1.0 - xi(2,n))
          xi3old = xi3temp
        end if
        p = p + 1
        xi2old = xi(2,n)
        err = 1.0 - xi2old
c       write (*,*) p, err
      end do
c
c     Now assemble the velocity field
c
      do i = 0, n
        u(i) = xi(2,i)
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
c
c     Use the spline result for the second derivative
c
      do i = 0, n
        d2u(i) = uspl(i)
      end do
      call SPLINE(n+1,ydat,d2u,d2uspl)

      do i = 0, n
        y = ymin + i*h
          call SPEVAL(n+1,ydat,u,uspl,y,us)
          call SPEVAL(n+1,ydat,d2u,d2uspl,y,d2us)
        write (10,10) y, us, d2us, u(i), d2u(i)
  10    format (1x,5(1PE20.13E2,1x))
      end do
      
      write (*,20) 
  20  format (1x,'Velocity Profile completed...',/)

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
      parameter   (idim=20000)
      integer     n
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p

      common      /data/  n, u, d2u, ymin, ymax, h, f2p,
     &                    gamma, beta, uspl, d2uspl, ydat
     
      common      /eig/   c, alpha, Re
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
      DIMENSION X(N),Y(N),FDP(N)      
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
C***********************************************************************
      subroutine READI(string, I)
C***********************************************************************
C
C     Read an integer from an input file
C
c***********************************************************************
      integer       I, iloc, floc
      character*80  string

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
      integer       iloc, floc
      real          R
      character*80  string

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
