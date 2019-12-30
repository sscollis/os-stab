      program Solve_Similar_BL
c******************************************************************************
c
c     Purpose:  This program solves the Boundary Layer similarity equation
c               using a Chebyshev collocation method.  Written originally
c               for ME 308 and ME 351B final projects.
c
c               This code works reasonably well but the convergence is rather
c               slow.  I think it can be made faster but I am not sure that
c               it is worth it.  I now have an RK4 boundary layer equation
c               integrator that works quite well.  See CONTEBL.F for example.
c
c     Author:   S. Scott Collis
c
c     Date:     2-22-92
c
c     Revised:  9-18-92
c
c******************************************************************************
      integer     idim
      parameter   (idim=128)
      real        u(0:idim), v(0:idim), xi(0:idim), th(0:idim)
      real        maxerr
      integer     i, j, k, n, nbig, maxcount
      
      real        gamma, etaout, phi1(0:idim), phi2(0:idim), beta, alpha
      common      /metrics/  gamma, etaout, phi1, phi2, beta, alpha
      
      write (*,10)
 10   format (/,/,10x,'Solve Boundary Layer Similarity Equation')
      write (*,20)
 20   format (/,1x,'Enter the number of modes ==> ',$)
      read (*,*) n
      write (*,60)
 60   format (/,1x,'Enter Beta ==> ',$)
      read (*,*) beta
      write (*,30)
 30   format (/,1x,'Enter alpha (alpha < 0) ==> ',$)
      read (*,*) alpha
      write (*,40)
 40   format (/,1x,'Enter the max allowable error ==> ',$)
      read (*,*) maxerr
      write (*,70)
 70   format (/,1x,'Enter the max number of iterations ==> ',$)
      read (*,*) maxcount
      write (*,50)
 50   format (/,1x,'Interpolate to mesh size ==> ',$)
      read (*,*) nbig
      
      gamma = 1.2
      etaout = 15.0

      call MAKE_MESH (n, th, xi)
      call MAKE_METRICS (n, xi)
      call INIT_PROFILE (n, u, v, th, xi)
      call SOLVE_BL (n, u, v, th, xi, maxerr, maxcount)
      call SAVE_VELOCITY(n,u,v,xi)
      call CHEBYINT(n,u,v,nbig)
      
      stop
      end
C******************************************************************************
      subroutine MAKE_MESH(n, th, xi)
C******************************************************************************
C
C     Setup the collocation points in the mapped coordinate, xi and in 
C     chebyshev space, th.
C
C******************************************************************************
      integer    n
      real       th(0:n), xi(0:n)
      real       pi, dth
      integer    i
      
      pi = Acos (-1.0)
      dth = pi/float(n)
      
      do i = 0, n
        th(i) = i*dth
        xi(i) = cos(th(i))
      end do
      
      return
      end
C******************************************************************************
      subroutine MAKE_METRICS(n, xi)
C******************************************************************************
C
C     Compute the first and second metrics for the coordinate transformation
C     I am using the truncated mapping of Street et. al.
C
C******************************************************************************
      integer    n
      real       xi(0:n), A, B
      integer    i
      
      parameter   (idim=128)
      real        gamma, etaout, phi1(0:idim), phi2(0:idim), beta, alpha
      common      /metrics/  gamma, etaout, phi1, phi2, beta, alpha

      write (*,*)
      write (*,*) 'Mesh and Metrics'
      write (*,*)
      do i = 0, n
        A = etaout/2.*(1-TANH(gamma))
        B = gamma/2.
        eta = etaout*(1-tanh(gamma))*.5*(xi(i)+1)/
     .               (1-tanh(gamma/2.*(xi(i)+1)))
        phi1(i) = 1./(A/2.*(1+EXP(2.*B*(1+xi(i))))+A*B*
     .            EXP(2.*B*(1+xi(i)))*(1+xi(i)))
        phi2(i) = -2.*A*B*(1+B*(xi(i)+1))*EXP(2.*B*(xi(i)+1))*phi1(i)**3
        write (*,10) eta, xi(i), phi1(i), phi2(i)
  10    format (4(es12.5,1x))
      end do
        
      return
      end
C******************************************************************************
      subroutine INIT_PROFILE(n, u, v, th, xi)
C******************************************************************************
C
C     Setup the initial boundary layer profiles
C
C******************************************************************************
      integer     n
      real        u(0:n), v(0:n), th(0:n), xi(0:n)
      real        eta
      integer     i, ans

      parameter   (idim=128)
      real        gamma, etaout, phi1(0:idim), phi2(0:idim), beta, alpha
      common      /metrics/  gamma, etaout, phi1, phi2, beta, alpha

      write (*,10)
  10  format (/,1x,'Read in a profile (1,0) ==> ',$)
      read (*,*) ans
      
      if (ans.eq.1) then
        call READ_VELOCITY(n,u,v)
      else
        do i = 0, n
          eta = etaout*(1-tanh(gamma))*.5*(xi(i)+1)/
     .                 (1-tanh(gamma/2.*(xi(i)+1)))
          u(i) = .5*(xi(i)+1) 
        end do
c
c       Now given the initial u, v must be set to satisfy continuity
c
        call COMPUTE_V (n, u, v, xi)
      end if
      
      write (*,*) 'Initial profiles'
      write (*,*)
      do i = 0, n
        eta = etaout*(1-tanh(gamma))*.5*(xi(i)+1)/
     .               (1-tanh(gamma/2.*(xi(i)+1)))
        write (*,12) eta, u(i), v(i)
  12    format (1x,3(e12.5,2x))
      end do
      write (*,*)

      return
      end
C******************************************************************************
      subroutine SPECTRAL(n, u, v, th, xi, L)
C******************************************************************************
C
C     Compute the spectral operator given a velocity field using Fast
C     Chebyshev transform
C
C******************************************************************************
      integer     n
      real        u(0:n), v(0:n), th(0:n), xi(0:n), L(1:n-1)
      integer     i, j
      real        eta

      parameter   (idim=128)
      real        gamma, etaout, phi1(0:idim), phi2(0:idim), beta, alpha
      common      /metrics/  gamma, etaout, phi1, phi2, beta, alpha

      real        u1(0:idim), u2(0:idim)
      
      if (n.gt.idim) then
        write (*,*) 'Error.  N > idim in SPECTRAL'
        stop
      end if
c
c     Compute the first derivative
c
      call DCHEBYSHEV (u, u1, u2, n)
      
c      write (*,*)
c      do i = 0, n
c        eta = etaout*(1-tanh(gamma))*.5*(xi(i)+1)/
c     .               (1-tanh(gamma/2.*(xi(i)+1)))
c        write (*,10) eta,u(i),phi1(i)*u1(i),
c     .             (phi1(i)**2*u2(i)+phi2(i)*u1(i)),phi1(i),phi2(i)
c   10   format (1x,6(e11.5,1x))
c      end do
c
C     Now form the spectral operator, L
c      
      do i = 1, n-1
        L(i) = phi1(i)**2*u2(i)+(phi2(i)-v(i)*phi1(i))*u1(i)
      end do
      
      return
      end

C******************************************************************************
      subroutine COMPUTE_V (n, u, v, xi)
C******************************************************************************
C
C     Compute v(i) from continuity given u(i).
C
C******************************************************************************
      integer     n, i
      real        u(0:n), v(0:n), xi(0:n)
      
      parameter   (idim=128)
      real        gamma, etaout, phi1(0:idim), phi2(0:idim), beta, alpha
      common      /metrics/  gamma, etaout, phi1, phi2, beta, alpha

      real        utemp(0:idim)
      
      do i = 0, n
        utemp(i) = -1./phi1(i)*u(i)
      end do
      call ICHEBYSHEV (utemp,v,n)
      
      return
      end

C******************************************************************************
      subroutine FINITE_DIFF(n, u, v, xi, Ha, Hb, Hc, b, L)
C******************************************************************************
C
C     Generate the second order finite difference preconditioning matrix.
C     Note that central difference is used which results in a tridiagonal
C     matrix.  Note that phi1 and phi2 are the first and 
C     second metrics of the coordinate transformation.
C
C******************************************************************************
      integer     n
      real        u(0:n), v(0:n), xi(0:n), Ha(n-1), Hb(n-1), Hc(n-1)
      real        b(1:n-1), L(1:n-1)
      real        x10, x11, x01
      integer     i, j

      parameter   (idim=128)
      real        gamma, etaout, phi1(0:idim), phi2(0:idim), beta, alpha
      common      /metrics/  gamma, etaout, phi1, phi2, beta, alpha
c
c     Finite difference operator must also have beta(1-u**2) term. maybe
c      
      do i = 1, n-1
        x10 = xi(i+1)-xi(i)
        x11 = xi(i+1)-xi(i-1)
        x01 = xi(i)-xi(i-1)
c
c       Richardson Iteration
c
        Ha(i) = (phi1(i)**2/x01-(phi2(i)-v(i)*phi1(i)))/x11
        Hb(i) = -1.0*phi1(i)**2*(1./x10 + 1./x01)/x11
        Hc(i) = (phi1(i)**2/x10+(phi2(i)-v(i)*phi1(i)))/x11
c
c       Simple pseudo time-step
c
c       Ha(i) = -(phi1(i)**2/x01-(phi2(i)-v(i)*phi1(i)))/x11
c       Hb(i) =  alpha+phi1(i)**2*(1./x10 + 1./x01)/x11
c       Hc(i) = -(phi1(i)**2/x10+(phi2(i)-v(i)*phi1(i)))/x11
      end do
      Ha(1)   = 0.0
      Hc(n-1) = 0.0
c
c     Form the R.H.S. and include the boundary conditions.  Note that 
c     u(0) = 1.0 and u(n) = 0.0
c
      do i = 1, n-1
        b(i) = alpha*(L(i)+beta*(1.-u(i)**2))
      end do
c      i = 1
c      x10 = xi(i+1)-xi(i)
c      x11 = xi(i+1)-xi(i-1)
c      x01 = xi(i)-xi(i-1)
c      b(i) = b(i)-(phi1(i)**2/x01-(phi2(i)-v(i)*phi1(i)))/x11*u(0)
c      i = n
c      x10 = xi(i+1)-xi(i)
c      x11 = xi(i+1)-xi(i-1)
c      x01 = xi(i)-xi(i-1)
c      b(i) = b(i)-(phi1(i)**2/x10+(phi2(i)-v(i)*phi1(i)))/x11*u(n)

      return
      end

C******************************************************************************
      subroutine SOLVE_BL(n, u, v, th, xi, maxerr, maxcount)
C******************************************************************************
C
C     Solve the similar boundary layer equations using a Chebyshev collocation
C     method.  On input u, and v contain the initial 'guess' and the non-linear
C     problem is solved iteratively starting from these inital profiles.  On
C     output u and v contain the solution
C
C******************************************************************************
      integer     n
      real        u(0:n), v(0:n), th(0:n), xi(0:n), maxerr
      parameter   (idim=128)
      real        du(idim), b(idim), L(idim)
      real        Ha(idim), Hb(idim), Hc(idim) 
      real        residual, eta, oldresidual
      integer     i, j, icount, maxcount

      real        gamma, etaout, phi1(0:idim), phi2(0:idim), beta, alpha
      common      /metrics/  gamma, etaout, phi1, phi2, beta, alpha

      real        H(idim,idim), work(idim), cond, z(idim), oa(idim,idim)
      integer     ipvt(idim)
      real        Ldiff
      real        u1(0:idim), u2(0:idim)

      if (n.gt.idim) then
        write (*,*) 'Error:  N > idim in SOLVE_BL'
        stop
      end if
     
      write(*,*) "Beginning Iterative Solve"
      write(*,*) 

      residual = 100.
      oldresidual = 1000.
      icount = 0
      do while ((residual .gt. maxerr).and.(icount.lt.maxcount))
        icount = icount + 1
c        do i = 1, idim
c          du(i) = 0.0
c          Ha(i) = 0.0
c          Hb(i) = 0.0
c          Hc(i) = 0.0
c          L(i) = 0.0
c          b(i) = 0.0
c        end do
c
c       Setup the spectral R.H.S.
c
        call SPECTRAL (n, u, v, th, xi, L)
c
c       Setup the finite difference opperator
c
        call FINITE_DIFF (n, u, v, xi, Ha, Hb, Hc, b, L)
c
c       Build the matrix
c
c        do i = 2, n-2
c          H(i,i-1) = Ha(i)
c          H(i,i)   = Hb(i)
c          H(i,i+1) = Hc(i)
c        end do
c        H(1,1) = Hb(1)
c        H(1,2) = Hc(1)
c        H(n-1,n-2) = Ha(n-1)
c        H(n-1,n-1) = Hb(n-1)
c        ndim = idim
c
c       H(u) should about equal L
c
c        write (*,*)
c        write (*,*) 'Spectral Finite Difference Comparison'
c        write (*,*)
c        do i = 1, n-1
c          Ldiff = 0.0
c          do j = 1, n-1
c            Ldiff = Ldiff + H(i,j)*u(j)
c          end do
c          if (i.eq.1) then 
c            Ldiff = Ldiff - b(1) + L(i)+beta*(1.-u(i)**2)
c          end if
c          write (*,*) xi(i), Ldiff, L(i)
c        end do
c        write (*,*)
                  
c        call DECOMP(ndim,n-1,H,COND,IPVT,WORK,OA,Z)                         
c        write (*,*) 'Condition = ',cond                                 
c        call SOLVE(ndim,n-1,H,B,IPVT,du)    
c        
        call TRIDIAG (n-1,du,Ha,Hb,Hc,b)
c
c        try without preconditioning
c
c         do i = 1, n-1
c           du(i) = b(i)
c         end do         
c
c       Update u vector, of course the b.c. don't change
c
        do i = 1, n-1
          u(i) = u(i) + du(i)
        end do
c
c       Now Chebyshev integrate to get the new v
c
        call COMPUTE_V (n, u, v, xi)
c
c       Check the residual to seen if it is small enough
c
        oldresidual = residual
        residual = 0.0
        do i = 1, n-1
          residual = residual + (du(i))**2
        end do
        residual = SQRT(residual/(n-1))
        write(*,*) icount, residual
        if ((mod(icount,50).eq.0)) then
          open(UNIT=10)
          call DCHEBYSHEV (u, u1, u2, n)
          do i = 0, n
            eta = etaout*(1-tanh(gamma))*.5*(xi(i)+1)/
     .                   (1-tanh(gamma/2.*(xi(i)+1)))
            write (10,10) eta,u(i),phi1(i)*u1(i),
     .                   (phi1(i)**2*u2(i)+phi2(i)*u1(i))
  10        format (1x,4(e11.5,2x))
          end do
          close(UNIT=10)
        end if
      end do
      write (*,*)
      write (*,*) "Completed solve with interations, residual"
      write (*,*) icount, residual

      return
      end
C******************************************************************************
C******************************************************************************
      SUBROUTINE DCHEBYSHEV(Y, Y1, Y2, N)
C******************************************************************************
C
C     Calculate the Chebyshev transform of a set Y of N+1 real-valued data 
C     points and compute the first and second derivatives in Chebyshev space.  
C     Then inverse transform and return the derivatives in Y1 and Y2 in real
C     space.  Note that Y is returned unscathed. 
C
C******************************************************************************
      REAL Y(0:N), Y1(0:N), Y2(0:N)
      PARAMETER (idim=128)
      REAL WORK(0:IDIM)

      IF (N.GT.IDIM) THEN
        WRITE (*,*) 'ERROR.  N > IDIM in DCHEBYSHEV'
        STOP
      END IF
C
C     SAVE THE INPUT VECTOR
C
      DO I = 0, N
        WORK(I) = Y(I)
      END DO
C
C     COMPUTE THE CHEBYSHEV TRANSFORM
C
      CALL chebyshev(Y,N,1)
C
C     NOW USE THE RECURSIVE RELATION TO TAKE THE FIRST DERIVATIVE
C
      Y1(N) = 0.0
      Y1(N-1) = 2.*FLOAT(N)*Y(N)
      DO K = N-2, 0, -1
        Y1(K) = Y1(K+2) + 2.*(K+1.)*Y(K+1)
      END DO
      Y1(0) = Y1(0)/2.
C
C     NOW REPEAT THE RECURSIVE RELATION TO TAKE THE SECOND DERIVATIVE
C
      Y2(N) = 0.0
      Y2(N-1) = 2.*FLOAT(N)*Y1(N)
      DO K = N-2, 0, -1
        Y2(K) = Y2(K+2) + 2.*(K+1.)*Y1(K+1)
      END DO
      Y2(0) = Y2(0)/2.
C
C     INVERSE TRANSFORM TO GET BACK TO REAL SPACE
C      
      CALL chebyshev(Y1,N,-1)
      CALL chebyshev(Y2,N,-1)
      DO I = 0, N
        Y(I) = WORK(I)
      END DO
      
      RETURN
      END

C***********************************************************************
      subroutine CHEBYSHEV (Y, N, ISIGN)
C***********************************************************************
C
C     Take the Chebyshev transform and normalize
C
C***********************************************************************
      integer  n, isign
      real     y(0:n), t(0:256)

      if (isign .ne. 1) then
        if (isign .ne. -1) then
          write (*,*) 'ERROR:  Invalid Isign in CHEBYSHEV'
          stop
        end if
      end if

      do i = 1, n
        t(i) = y(i)
      end do
            
      call COSFT3 (y,n,isign)

c     call COSFT (y,n,isign)
    
      if (isign .eq. 1) then
        do i = 0,n
          y(i) = y(i)/float(n)*2.
        end do
      end if
      
      return
      end

C******************************************************************************
      SUBROUTINE COSFT3 (Y,N,ISIGN)
C******************************************************************************
C
C     Do a brute force transform.  ISIGN = 1 does a forward transform.
C     ISIGN = -1 a backward transform.  This is a hybrid double/single
C     precision routine.
C
C******************************************************************************
      integer     N, ISIGN
      real        Y(0:N), T(0:N)
      real*16     YY(0:N), TT(0:N), PI

      PI = DACOS(-1.0D0)

      DO I = 0, N
        TT(I) = DBLE(Y(I))
      END DO

      if (isign .eq. 1) then
        DO I = 0, N
          YY(I) = TT(0)/2.D0 + TT(N)*COS(DFLOAT(I)*PI)/2.D0
          DO M = 1, N-1
            YY(I) = YY(I) + TT(M)*COS(DFLOAT(M)*DFLOAT(I)*PI/DFLOAT(N))
          END DO
          YY(I) = YY(I)*2./DFLOAT(N)
        end do
        YY(0) = YY(0)/2.D0
        YY(N) = YY(N)/2.D0
      else 
        DO I = 0, N
          YY(I) = 0.0
          DO M = 0, N
            YY(I) = YY(I) + TT(M)*COS(DFLOAT(M)*DFLOAT(I)*PI/DFLOAT(N))
          END DO
        END DO
      end if
C
C     GO BACK TO SINGLE PRECISION
C
      DO I = 0, N
        if (Y(i).le.1e-14) then
          Y(I) = 0.0
        end if
        Y(I) = SNGL(YY(I))
      END DO
      
      RETURN
      END

C*****************************************************************************
      subroutine TRIDIAG(num,u,a,b,c,f)
C*****************************************************************************
c
c     Solves the tridiagonal system B[a,b,c]u = f
c
c*****************************************************************************
      real    u(num), a(num), b(num), c(num), f(num)

      a(num) = a(num)/b(num)
      do i = num-1,1,-1
        b(i) = b(i)-c(i)*a(i+1)
        a(i) = a(i)/b(i)
      end do
      
      u(num) = f(num)/b(num)
      do i = num-1,1,-1
        u(i) = (f(i)-c(i)*u(i+1))/b(i)
      end do
      
      do i = 2, num
        u(i)=u(i)-a(i)*u(i-1)
      end do

      return
      end

C******************************************************************************
      SUBROUTINE ICHEBYSHEV(Y,YI,N)
C******************************************************************************
C
C     Calculate the Chebyshev transform of a set Y of N+1 real-valued data 
C     points and integrate in Chebyshev space.  The integral is returned
C     in real space in YI and Y is unscathed.
C
C******************************************************************************
      REAL Y(0:N), YI(0:N+1)
      PARAMETER (idim=128)
      REAL WORK(0:IDIM)
      
      DO I = 0, N
        WORK(I) = Y(I)
      END DO
C
C     COMPUTE THE CHEBYSHEV TRANSFORM
C
      CALL Chebyshev(Y,N,1)
C
C     NOW INTEGRATE
C     NOTE I AM ASSUMING THAT Y(N) IS SMALL SUCH THAT YI(N+1) = 0.0     
C
      YI(N+1) = Y(N)/(2.*(N+1.))
      YI(N) = Y(N-1)/(2.*N)
      YI(1) = 1./2.*(2.*Y(0)-Y(2))
      YI(0) = YI(1)
      DO K = 2, N-1
        YI(K) = 1./(2.*K)*(Y(K-1)-Y(K+1))
        YI(0) = YI(0) + (-1.)**(K-1)*YI(K) 
      END DO
      
      CALL Chebyshev(YI,N,-1)
      DO I = 0, N
        Y(I) = WORK(I)
      END DO

      RETURN
      END
C******************************************************************************
      SUBROUTINE CHEBYINT(N,U,V,NBIG)
C******************************************************************************
C
C
C******************************************************************************
      INTEGER   N, NBIG
      REAL      U(0:N), V(0:N)
      REAL      IU, IV, X, PI
      
      parameter   (idim=128)
      real        gamma, etaout, phi1(0:idim), phi2(0:idim), beta, alpha
      common      /metrics/  gamma, etaout, phi1, phi2, beta, alpha

      PI = ACOS(-1.0)

      call Chebyshev (u,n,1)
      call Chebyshev (v,n,1)

      open(10)
      do i = 0, NBIG
        X = I*PI/NBIG
        IU = 0.0
        IV = 0.0
        DO M = 0, N
          IU = IU + U(M)*COS(FLOAT(M)*X)
          IV = IV + V(M)*COS(FLOAT(M)*X)
        END DO
        IU = FLOAT(N)*IU*0.5
        IV = FLOAT(N)*IV*0.5
        eta = etaout*(1-tanh(gamma))*.5*(cos(X)+1)/
     .               (1-tanh(gamma/2.*(cos(X)+1)))
        write (10,10) eta, IU, IV
      end do
      close(10)
  10  format (1x,3(e16.8,4x))

      RETURN
      END
      
C******************************************************************************
      SUBROUTINE SAVE_VELOCITY(N,U,V, xi)
C******************************************************************************
C
C
C******************************************************************************
      INTEGER   N, i, j
      REAL      U(0:N), V(0:N), xi(0:N), eta
      character*15 filename
      
      parameter   (idim=128)
      real        gamma, etaout, phi1(0:idim), phi2(0:idim), beta, alpha
      common      /metrics/  gamma, etaout, phi1, phi2, beta, alpha

      write (*,5)
   5  format (/,1x,'Write file')
      write (*,10)
  10  format (1x,'Enter filename ==> ',$)
      read (*,'(a)') filename
      open (unit=11,file=filename,status='unknown')
      write (11,*) N
      do i = 0, n
        eta = etaout*(1.-tanh(gamma))*.5*(xi(i)+1.)/
     .               (1.-tanh(gamma/2.*(xi(i)+1.)))
        write (11,*) eta,u(i),v(i)
      end do
      close (11)
      write (*,*)
      write (*,*) 'Wrote file = ', filename
      
      RETURN
      END
C******************************************************************************
      SUBROUTINE READ_VELOCITY(N,U,V)
C******************************************************************************
C
C
C******************************************************************************
      INTEGER   N, i, j, nmode
      REAL      U(0:N), V(0:N),eta
      character*15 filename
      
      write (*,5)
   5  format (/,1x,'Read file',/)
      write (*,10)
  10  format (1x,'Enter filename ==> ',$)
      read (*,'(a)') filename
      open (unit=11,file=filename,status='unknown')
      read (11,*) nmode
      do i = 0, n
        read (11,*) eta,u(i),v(i)
      end do
      close (11)
      
      RETURN
      END
