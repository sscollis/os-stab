      program Orr_col
c***********************************************************************
c
c     Purpose:  This program solves the Orr-Sommerfeld equation using a 
c               Chebyshev-collocation method for channel flow.
c
c     Author:   S. Scott Collis
c
c     Date:     2-29-92
c
c     Revision: 9-18-92
c
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=1024)
      integer     n
      real        eta(0:idim), th(0:idim), lmap
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim)
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     &                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      real        alphar, alphai
c
c     Setup IMSL workspace
c
c      REAL             RWKSP(200000)
c      COMMON /WORKSP/  RWKSP
c      call IWKIN(200000)
c
c     User input
c      
      write (*,10)
  10  format (/,/,10x,'Solve Orr-Sommerfeld (Collocation)')
      write (*,20)
  20  format (/,1x,'Enter the number of modes ==> ',$)
      read (*,*) n
      write (*,30)
  30  format (/,1x,'Enter Reynolds number ==> ',$)
      read (*,*) Re
      write (*,40) 
  40  format (/,1x,'Enter alpha (R,I) ==> ',$)
      read (*,*) alphar, alphai
      
      alpha = cmplx (alphar, alphai)
      
      call MAKE_DERIVATIVES
      call MAKE_CHANNEL_METRICS

      call INIT_CHANNEL_PROFILE

      call MAKE_MATRIX
      call SOLVE_ORR_SOM
      
      stop
      end

C***********************************************************************
      subroutine MAKE_DERIVATIVES
C***********************************************************************
C
C     Make the required matrices that take derivatives in Chebyshev 
C     space
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=1024)
      integer     n
      real        eta(0:idim), th(0:idim), lmap
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim)
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     &                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      integer  m, p, LDD, i, j, k
      real     Identity(0:idim,0:idim), D1hat(0:idim,0:idim)
      
      LDD = idim
      call CHEBYD (D1hat, LDD, n)
      call CHEBYD (D1   , LDD, n)
c
c     We enforce the Neumann boundary conditions indirectly by setting
c     v' = 0 at 0 and n.  NOTE THAT D1-D4 SHOULD ONLY BE APPLIED TO THE
c     v FIELD ONLY.
c
      do i = 0, n
        D1(0,i) = 0.0
        D1(N,i) = 0.0
      end do
C
C     To get higher derivatives just do matrix multiplication
C      
      do i = 0, n
        do j = 0, n
          D2(i,j) = 0.0
          do k = 0, n
            D2(i,j) = D2(i,j) +  D1hat(i,k)*D1(k,j)
          end do
        end do
      end do
      do i = 0, n
        do j = 0, n
          D3(i,j) = 0.0
          do k = 0, n
            D3(i,j) = D3(i,j) +  D1hat(i,k)*D2(k,j)
          end do
        end do
      end do
      do i = 0, n
        do j = 0, n
          D4(i,j) = 0.0
          do k = 0, n
            D4(i,j) = D4(i,j) +  D1hat(i,k)*D3(k,j)
          end do
        end do
      end do

c      CALL WRRRN ('D1', N+1, N+1, D1, LDd+1, 0)
c      CALL WRRRN ('D2', N+1, N+1, D2, LDd+1, 0)
c      CALL WRRRN ('D3', N+1, N+1, D3, LDd+1, 0)
c      CALL WRRRN ('D4', N+1, N+1, D4, LDd+1, 0)
      
      return
      end

C***********************************************************************
      subroutine MAKE_CHANNEL_METRICS
C***********************************************************************
C
C     Setup the collocation points in the mapped coordinate, eta and in 
C     chebyshev space, th.  Also compute the transformation metrics.
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=1024)
      integer     n
      real        eta(0:idim), th(0:idim), lmap
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim)
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     &                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      real       pi, dth
      integer    i
      
      pi  = ACOS(-1.0)
      dth = pi/FLOAT(n)
c
c     Make mesh in transformed, eta, and Chebyshev space, th
c
      do i = 0, n
        th(i) = FLOAT(i)*dth
        eta(i) = COS(th(i))
      end do

      return
      end

C***********************************************************************
      subroutine INIT_CHANNEL_PROFILE
C***********************************************************************
C
C     Setup the initial channel profile
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=1024)
      integer     n
      real        eta(0:idim), th(0:idim), lmap
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim)
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     &                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      integer     i, LDD
      real        D1hat(0:idim,0:idim), D2hat(0:idim,0:idim)  
                
      LDD = idim
     
      open(10) 
      do i = 0, n
        u(i) = (1.-eta(i)**2)
        d1u(i) = -2.*eta(i)
        d2u(i) = -2.0
        write(10,20) eta(i), u(i), d1u(i), d2u(i)
 20     format (1x,4(es16.8e3,1x))
      end do
      close(10)
#if 0
c
c     Compute the collocation derivatives
c
      call CHEBYD (D1hat, LDD, n)
c     CALL WRRRN ('D1hat', N+1, N+1, D1hat, LDd+1, 0)
      do i = 0, n
        do j = 0, n
          D2hat(i,j) = 0.0
          do k = 0, n
            D2hat(i,j) = D2hat(i,j) +  D1hat(i,k)*D1hat(k,j)
          end do
        end do
      end do
      do i = 0, n
        d1u(i) = 0.0
        d2u(i) = 0.0
        do k = 0, n
            d1u(i) = d1u(i) + D1hat(i,k)*u(k)
            d2u(i) = d2u(i) + D2hat(i,k)*u(k)
        end do
c       write (*,10) eta(i), u(i), d1u(i), d2u(i)
c  10   format (1x,4(e12.4,1x))
      end do
#endif
      return
      end

C***********************************************************************
      subroutine MAKE_MATRIX
C***********************************************************************
C
C     This routine generates the matrices which are combined to make
C     the generalized eigenvalue problem.
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=1024)
      integer     n
      real        eta(0:idim), th(0:idim), lmap
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim)
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     &                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      real      rwork(0:idim), cwork(0:idim), identity(0:idim,0:idim)
      complex   work
      
      do i = 0, n
        do j = 0, n
          A4(i,j) = 0.0
          A3(i,j) = 0.0
          A2(i,j) = 0.0
          A1(i,j) = 0.0
          B2(i,j) = 0.0
          B1(i,j) = 0.0
          B0(i,j) = 0.0
          identity(i,j) = 0.0
        end do
        identity(i,i) = 1.0
      end do
      
      do i = 0, n
        do j = 0, n
          A4(i,j) = CMPLX( D4(i,j), 0.0 )
        end do
      end do

      do i = 0, n
        work = -2.0 * alpha**2 - cmplx(0.,1.)*alpha*Re*U(i)
        do j = 0, n
          A2(i,j) = work*D2(i,j)
        end do
      end do

      do i = 0, n
        work = alpha**4 + cmplx(0.,1.)*alpha**3*Re*U(i) + 
     &                    cmplx(0.,1.)*alpha*Re*d2u(i)
        do j = 0, n
          A0(i,j) = work*identity(i,j)
        end do
      end do

      do i = 0, n
        work = -1.0*cmplx(0.,1.)*Re
        do j = 0, n
          B2(i,j) = work*D2(i,j)
        end do
      end do

      do i = 0, n
        work = alpha**2*cmplx(0.,1.)*Re
        do j = 0, n
          B0(i,j) = work*identity(i,j)
        end do
      end do
                         
      return
      end

C***********************************************************************
      subroutine SOLVE_ORR_SOM
C***********************************************************************
C
C     This routine generates the discrete eigenvalue problem in 
C     Chebyshev space for the Orr-Sommerfeld equation.
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=1024)
      integer     n
      real        eta(0:idim), th(0:idim), lmap
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim)
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     &                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      complex     A(0:idim,0:idim), B(0:idim,0:idim), alp(0:idim)
      complex     T1(0:idim,0:idim), T2(0:idim,0:idim)
      complex     T3(0:idim,0:idim), T4(0:idim,0:idim)
      complex     beta(0:idim), evec(0:idim,0:idim), temp3(0:idim)
      complex     avec(0:idim,0:idim)
      complex     evec2(0:idim,0:idim), eval(0:idim)
      complex     P1(0:idim,0:idim), P2(0:idim,0:idim)
      real        temp1(0:idim), temp2(0:idim), CHECKEIG, nw
      integer     lda, ldb, ldevec, p, i, j, which, k, l, m
      complex     eigenvalue, eigenvector(0:idim), tvec(0:idim)
      integer     index(0:idim)
      logical     first, print
      external    CHECKEIG

      integer     info, lwork, ipvt(0:idim)
      complex     work(16*(idim+1))
      real        rwork(8*(idim+1))

      complex     zdotc, escale
      external    zdotc
C***********************************************************************
      lwork = 16*(idim+1)

      lda = idim+1
      ldb = idim+1
      ldevec = idim+1
c
c.... Build the matrix system
c      
      do i = 0, n
        do j = 0, n
          A(i,j) = A4(i,j)+A3(i,j)+A2(i,j)+A1(i,j)+A0(i,j)
          B(i,j) = B2(i,j)+B1(i,j)+B0(i,j)
        end do
      end do
c
c     Now enforce the v = 0 @ +_ 1 boundary condition instead of the
c     Orr-Sommerfeld equation at n = 0 and n = N.  The Neumann 
c     boundary condition is implicitly enforced.
c
      do j = 0, n
        B(0,j) = CMPLX(0.0, 0.0)
        B(j,0) = CMPLX(0.0, 0.0)
        B(n,j) = CMPLX(0.0, 0.0)
        B(j,n) = CMPLX(0.0, 0.0)
        A(0,j) = CMPLX(0.0, 0.0)
        A(j,0) = CMPLX(0.0, 0.0)
        A(n,j) = CMPLX(0.0, 0.0)
        A(j,n) = CMPLX(0.0, 0.0)
      end do
      A(0,0) = CMPLX(1.0, 0.0) 
      A(n,n) = CMPLX(1.0, 0.0)
c
c.... Solve the eigenvalue problem
c
#ifdef USE_ISML
c
c.... The top and bottom rows are trivial so that they could be removed
c
      CALL LINCG (N-1, B, LDA, T1, LDA)
      CALL MCRCR (N-1, N-1, T1, LDA, N-1, N-1, A, LDA,
     &            N-1, N-1, T4, LDA)

      CALL EVCCG (N-1, T4, LDA, eval, evec, ldevec)
#else
      if (.false.) then
c
c.... Cannot use this as B is singular 
c
         call ZGETRF(N+1, N+1, B, lda, ipvt, info)
         call ZGETRS('N', N+1, N+1, B, lda, ipvt, A, lda, info)
         call ZGEEV('V', 'V', N+1, A, lda, eval, avec,
     &              lda, evec, lda, work, lwork, rwork, info)
      else
c
c.... Use generalized engenvalue solver
c
         call ZGEGV( 'V', 'V', N+1, A, LDA, B, LDA, alp, beta, avec, 
     &               LDA, evec, LDA, work, lwork, rwork, info)
c
c.... compute the eigenvalues as ratios of alp and beta
c
         do i = 0, N
           if (beta(i).ne.0) then
             eval(i) = alp(i)/beta(i)
           else
             eval(i) = 0.0
           end if
         end do
      end if
#endif
c
c     Sort the eigenvalues.
c      
      do i = 0, N
        temp2(i) = IMAG(eval(i))
        index(i) = i
      end do

c     Need to issolate the most unstable eigenvalue and eigenvector
c     Must watch out, this routine isn't very robust.

      call PIKSR2(n+1,temp2,index)
      do i = 0, N
         temp1(i) = real(eval(index(i)))
         do j = 0, N
            A(j,i) = evec(j,index(i))
            B(j,i) = avec(j,index(i))
         end do
      end do
      
      do i = 0, N
         eval(i) = cmplx(temp1(i),temp2(i))
         do j = 0, N
            evec(j,i) = A(j,i)
            avec(j,i) = B(j,i)
         end do
      end do

      write (*,50) N+1
  50  format (/,1x,'Eigenvalues for N+1 = ',i4,/)
c
c     Print eigenvalues
c        
      do i = 0, N
        omega(i) = eval(i)
        write (*,35) i, REAL(omega(i)), IMAG(omega(i)), 
     &               REAL(omega(i)/alpha), IMAG(omega(i)/alpha)
  35    format (1x,i5,2x,4(1pe17.10,1x))
      end do
     
      escale = 0.0 
      write (*,40)
      read (*,*) which
      do while (which .ne. -1)

        eigenvalue = eval(which)
        do i = 0, N
          eigenvector(i) = evec(i,which)
        end do
c       write (*,*) CHECKEIG(N-1,T4,lda,eigenvalue,eigenvector)
c
c.... normalize the eigenvector
c
        do i = 0, N
          if (abs(evec(i,which)) .gt. abs(escale)) then
            escale = evec(i,which)
          endif
        end do
        escale = 1.0/(escale)
        call zscal( N+1, escale, evec(0,which), 1)
c
c.... check the orthogonality of the adjoint and regular eigenvectors
c
        escale = zdotc( N+1, avec(0,which), 1, evec(0,which), 1 )
        write(*,*) escale
c
c.... scale the adjoint by the inner product
c
        escale = 1.0/conjg(escale)
        call zscal( N+1, escale, avec(0,which), 1)
        escale = zdotc( N+1, avec(0,which), 1, evec(0,which), 1 )
        write(*,*) escale

c       call CHEBYSHEV(temp1,n,1)
c       call CHEBYSHEV(temp2,n,1)
        open(10)
        do i = 0, N
          write (11,999) eta(i), 
     &                   real(evec(i,which)), aimag(evec(i,which)), 
     &                   real(avec(i,which)), aimag(avec(i,which))
 999      format (1x,6(ES16.8E3,1x))
        end do
        close(11)
c
c.... Interpolate to a finer mesh
c       
        if (.false.) then
         do i = 1, n
           temp1(i) = REAL(avec(i,which))
           temp2(i) = IMAG(avec(i,which))
         end do
         call CHEBYINT(n, temp1, temp2, 256)
         end if
 
        write (*,40)
        read (*,*) which
      end do
 40   format (/,1x,'Eigenvector for which eigenvalue (-1 quits) ==> ',$)
      
      return
      end

C***********************************************************************
      function CHECKEIG(N,A,LDA,EVAL,EVEC)
C***********************************************************************
C
C     Check an eigenvalue and eigenvector
C
C***********************************************************************
      PARAMETER (idim=1024)
      
      integer    N
      complex    A(LDA,N), EVAL, EVEC(N)
      complex    X(idim), Y(idim)
      real       CHECKEIG
#ifdef USE_IMSL
      CALL MUCRV (N, N, A, LDA, N, EVEC, 1, N, X)
#else
      CALL ZGEMV ('N', N, N, 1.0, A, LDA, EVEC, 1, 0.0, X, 1)
#endif
      CHECKEIG = 0.0
      DO I = 1, N
        CHECKEIG = CHECKEIG + ABS(X(I)-EVAL*EVEC(I))
      END DO
      CHECKEIG = CHECKEIG/FLOAT(N)
      
      RETURN
      END

C***********************************************************************
      SUBROUTINE DCHEBYSHEV(N, Y, DY)
C***********************************************************************
C
C     Calculate the Chebyshev transform of a set Y of N+1 real-valued
C     data points and compute the derivative in Chebyshev space. 
C     Then inverse transform and return the derivative in DY in real 
C     space.  Note that Y is returned unscathed. 
C
C***********************************************************************
      INTEGER N
      REAL Y(0:N), DY(0:N)

      PARAMETER (idim=1024)
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
      CALL CHEBYSHEV(Y,N,1)
C
C     NOW USE THE RECURSIVE RELATION TO TAKE THE DERIVATIVE
C
      DY(N) = 0.0
      DY(N-1) = 2.*FLOAT(N)*Y(N)
      DO K = N-2, 0, -1
        DY(K) = DY(K+2) + 2.*(K+1.)*Y(K+1)
      END DO
      DY(0) = DY(0)/2.
C
C     INVERSE TRANSFORM TO GET BACK TO REAL SPACE
C      
      CALL CHEBYSHEV(DY,N,-1)

      DO I = 0, N
        Y(I) = WORK(I)
      END DO
      
      RETURN
      END

C***********************************************************************
      SUBROUTINE ICHEBYSHEV(Y,YI,N)
C***********************************************************************
C
C     Calculate the Chebyshev transform of a set Y of N+1 real-valued
C     data points and integrate in Chebyshev space.  The integral is 
C     returned in real space in YI and Y is unscathed.
C
C***********************************************************************
      REAL Y(0:N), YI(0:N+1)
      PARAMETER (idim=1024)
      REAL WORK(0:IDIM)
      
      DO I = 0, N
        WORK(I) = Y(I)
      END DO
C
C     COMPUTE THE CHEBYSHEV TRANSFORM
C
      CALL CHEBYSHEV(Y,N,1)
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
      
      CALL CHEBYSHEV(YI,N,-1)
      DO I = 0, N
        Y(I) = WORK(I)
      END DO

      RETURN
      END

C***********************************************************************
      SUBROUTINE CHEBYINT(nmode,U,V,NBIG)
C***********************************************************************
C
C     Interpolate two functions to a finer mesh
C
c***********************************************************************
      INTEGER   NBIG, nmode
      REAL      U(0:Nmode), V(0:Nmode)
      REAL      IU, IV, X, PI
      
      PI = ACOS(-1.0)

      call CHEBYSHEV (U,nmode,1)
      call CHEBYSHEV (V,nmode,1)

      open(30)
      do i = 0, NBIG
        X = I*PI/NBIG
        IU = 0.0
        IV = 0.0
        DO M = 0, Nmode
          IU = IU + U(M)*COS(FLOAT(M)*X)
          IV = IV + V(M)*COS(FLOAT(M)*X)
        END DO
        write (30,10) cos(x), Iu, Iv
      end do
      close(30)
  10  format (1x,3(es16.8e3,1x))

      RETURN
      END

C************************************************************************
      SUBROUTINE COSFT3 (Y,N,ISIGN)
C************************************************************************
C
C     Do a brute force transform.  ISIGN = 1 does a forward transform.
C     ISIGN = -1 a backward transform.  This is a hybrid double/single
C     precision routine.
C
C************************************************************************
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

C     call COSFT (y,n,isign)
      
C     if (isign .eq. 1) then
C       do i = 0,n
C         y(i) = y(i)/float(n)*2.
c          if (abs(y(i)) .lt. 1.e-6) then
c              y(i) = 0.0
c          end if
c           write (*,*) i, y(i)
C       end do
C    end if
      
      return
      end

C***********************************************************************
      subroutine CONVOLVE (NDIM, A, B, LDB, ES, LDE, C)
C***********************************************************************
C
C     Form the Chebyshev convolution of real vector A with real matrix B 
C     and put the result in real matrix E.
c
c     A  double
c     B  double
c     ES single
C
C***********************************************************************
      INTEGER NDIM, LDB, LDE
      REAL    ES(0:LDE,0:NDIM), C(0:NDIM)
      REAL*16 A(0:NDIM),  B(0:LDB,0:NDIM),  E(0:LDE,0:NDIM)

      call COSFT3D (A,NDIM,1)

      DO N = 0, NDIM
        DO J = 0, NDIM
          E(N,J) = 0.0D0
        END DO
        DO M = -NDIM, NDIM
          IF (ABS(N-M).LE.NDIM) THEN
            DO J = 0, NDIM
              E(N,J) = E(N,J)+ DBLE(C(ABS(N-M))) * A(ABS(N-M)) *
     &                         DBLE(C(ABS(M)))   * B(ABS(M),J)
            END DO
          END IF
        END DO
        DO J = 0, NDIM
          E(N,J) = E(N,J)/2.D0
        END DO
      END DO
      DO J = 0, NDIM
        E(0,J) = E(0,J)/2.D0
      END DO

      DO N = 0, NDIM
        DO J = 0, NDIM
          ES(N,J) = SNGL(E(N,J))
        END DO
      END DO

      RETURN
      END
            
C**************************************************************************
      SUBROUTINE COSFT3D (YY,N,ISIGN)
C**************************************************************************
C
C     Do a brute force transform. ISIGN = 1 does a forward transform.
C     ISIGN = -1 a backward transform. This one is fully double precision.
C
C**************************************************************************
      integer     N, ISIGN
      real*16     YY(0:N), TT(0:N), PI

      PI = DACOS(-1.0D0)

      DO I = 0, N
        TT(I) = YY(I)
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
      
      RETURN
      END
      
C**************************************************************************
      SUBROUTINE CHEBYD (D, LDD1, N)
C**************************************************************************
C
C     Calculation the Chebyshev collocation derivative matrix
C
C**************************************************************************
      PARAMETER (IDIM=1024)
      REAL      D(0:LDD1,0:N), PI
      REAL      C(0:IDIM)

      IF (N .GT. IDIM) THEN
        WRITE (*,*) 'ERROR:  N > IDIM in CHEBYD'
        STOP
      END IF
       
      PI = ACOS(-1.0)
      
      DO I = 1, N-1
        C(I) = 1.0
      END DO
      C(0) = 2.0
      C(N) = 2.0
            
      DO J = 0, N
        DO K = 0, N
          IF ( J.EQ.0 .AND. K.EQ.0) THEN
            D(J,K) = (2.0*FLOAT(N)**2+1.0)/6.0
          ELSE IF ( J.EQ.N .AND. K.EQ.N) THEN
            D(J,K) = -1.0*(2.*FLOAT(N)**2+1.0)/6.0
          ELSE IF (J.EQ.K) THEN
            D(J,K) =  -1.0*COS(FLOAT(J)*PI/FLOAT(N))/
     &                (2.*(1.-(COS(FLOAT(J)*PI/FLOAT(N)))**2))
          ELSE
            D(J,K) = C(J)*(-1.0)**(J+K)/(C(K)*(COS(FLOAT(J)*PI/
     &               FLOAT(N))-COS(FLOAT(K)*PI/FLOAT(N))))
          END IF
        END DO
      END DO

      RETURN
      END
C*************************************************************************
C
C     These routines are for READ_BL only!
C
C*************************************************************************
      SUBROUTINE DCHEBYSHEVF(Y,YP,N)
C*************************************************************************
C
C     Calculate the Chebyshev transform of a set Y of N+1 real-valued data 
C     points and compute the derivative in Chebyshev space.  Then inverse
C     transform and return the derivative in YP and the chebyshev 
C     coefficients of the derivative in Y.
C
C************************************************************************
      DIMENSION Y(0:N), YP(0:N)
C
C     NOW USE THE RECURSIVE RELATION TO TAKE THE DERIVATIVE
C
      YP(N) = 0.0
      YP(N-1) = 2.*FLOAT(N)*Y(N)
      DO K = N-2, 0, -1
        YP(K) = YP(K+2) + 2.*(K+1.)*Y(K+1)
      END DO
      YP(0) = YP(0)/2.

c     DO K = 0, N
c       Y(K) = YP(K)
c     END DO
      
      RETURN
      END

C***********************************************************************
      SUBROUTINE COSFT (Y,N,ISIGN)
C***********************************************************************
C
C     Do a brute force transform.  ISIGN = 1 does a forward transform.
C     ISIGN = -1 a backward transform.  This is a hybrid double/single
C     precision routine.
C
C***********************************************************************
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

C*************************************************************************
      SUBROUTINE CHEBYINTF(N,INTERP,X,Y)
C*************************************************************************
C
C     Do a Chebyshev interpolation based on the Chebyshev coefficients Y.
C     The interpolated result is given in INTERP.
C
C*************************************************************************
      INTEGER   N
      REAL      Y(0:N)
      REAL      INTERP, X

      INTERP = 0.0
      DO M = 0, N
        INTERP = INTERP + Y(M)*COS(FLOAT(M)*X)
      END DO
      
      RETURN
      END
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
