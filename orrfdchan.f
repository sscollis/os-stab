      program orr_fd_chan
c***********************************************************************
c
c     Purpose:  This program solves the Orr-Sommerfeld equation using a 
c               finite-difference method for channel profiles.
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
      parameter   (idim=512)
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
     .                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      integer     nar, nai, ians
      real        alphar(100), alphai(100)
      real        minar, maxar, incar, minai, maxai, incai, dar, dai
      complex     eigenvalue(100,100), eigenvector(0:idim), ctemp
      character*9 filename
      logical     print
c***********************************************************************
c
c     Setup IMSL workspace
c
#ifdef USE_IMSL
      REAL             RWKSP(150000)
      COMMON /WORKSP/  RWKSP
      call IWKIN(150000)
#endif
c
c     User input
c      
      write (*,10)
  10  format (/,/,10x,'Solve Orr-Sommerfeld (Finite-difference)')
      write (*,20)
  20  format (/,1x,'Enter the number of mesh points ==> ',$)
      read (*,*) n
      write (*,30)
  30  format (/,1x,'Enter Reynolds number ==> ',$)
      read (*,*) Re
c     write (*,65)
c 65  format (/,1x,'Enter mapping metric ==> ',$)
c     read (*,*) Lmap
      lmap = 1.0
      write (*,40) 
  40  format (/,1x,'Enter alpha_r (min,max,inc) ==> ',$)
      read (*,*) minar,maxar,incar
      write (*,45) 
  45  format (/,1x,'Enter alpha_i (min,max,inc) ==> ',$)
      read (*,*) minai,maxai,incai
      write (*,50) 
  50  format (/,1x,'Enter root filename ==> ',$)
      read (*,'(a)') filename
      write (*,52) 
  52  format (/,1x,'Print eigenvectors (1,0) ==> ',$)
      read (*,*) ians
      if (ians.eq.1) then
         print = .true.
      else
         print = .false.
      end if

      nar = aint((maxar-minar)/incar)+1
      nai = aint((maxai-minai)/incai)+1
      write (*,55) nar, nai
 55   format (/,1x,' nar = ',i5,'  nai = ',i5)
      do i = 1, nar
        alphar(i) = minar + (i-1.)*incar
      end do
      do i = 1, nai
        alphai(i) = minai + (i-1.)*incai
      end do

      nx = nar
      ny = nai
      iloc = index(filename,' ')-1
      open (unit=10,file=filename(1:iloc)//'.g',form='formatted',
     .      status='unknown')
      write (10,100) nx, ny
 100  format(5(i5,1x))
      write (10,110) ((alphar(i), i=1,nx), j=1,ny),
     .               ((alphai(j), i=1,nx), j=1,ny)
 110  format (5(ES16.8E3,1x))
      close (10)
      open (unit=12,file=filename(1:iloc)//'.dat',form='formatted',
     .      status='unknown')
      call MAKE_CHANNEL_METRICS
      call MAKE_DERIVATIVES
      call INIT_CHANNEL_PROFILE
      do j = 1, ny
        do i = 1, nx
          alpha = cmplx (alphar(i), alphai(j))
          call MAKE_MATRIX
          call SOLVE_ORR_SOM(ctemp, eigenvector, print)
          eigenvalue(i,j) = ctemp
          write (12,60) real(alpha),aimag(alpha),real(eigenvalue(i,j)),
     .                  aimag(eigenvalue(i,j))
          write (*,60)  real(alpha),aimag(alpha),real(eigenvalue(i,j)),
     .                  aimag(eigenvalue(i,j))
  60      format (1x,4(ES16.8E3,1x))
        end do
      end do

      open (unit=11,file=filename(1:iloc)//'.q',form='formatted',
     .      status='unknown')
      write (11,100) nx, ny
      write (11,110) 0.0,0.0,0.0,0.0
      write (11,110) ((1.0                   , i=1,nx), j=1,ny),
     .               ((real(eigenvalue(i,j)) , i=1,nx), j=1,ny),
     .               ((aimag(eigenvalue(i,j)), i=1,nx), j=1,ny),
     .               ((1.0                   , i=1,nx), j=1,ny)
      close (11)
      
      stop
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
      parameter   (idim=512)
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
     .                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      real       pi, deta
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
c
c     Make transformation metrics
c
      do i = 0, n
        m1(i) = 1.0
        m2(i) = 0.
        m3(i) = 0.
        m4(i) = 0.
        if (i.eq.0) then
          c(i) = 2.
        else
          c(i) = 1.
        end if
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
      parameter   (idim=512)
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
     .                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      integer     i, LDD
      real        D1hat(0:idim,0:idim), D2hat(0:idim,0:idim)  
                
      LDD = idim
      
      do i = 0, n
        u(i) = (1.-eta(i)**2)
      end do
c
c     Compute the finite difference derivatives
c
      call FiniteD1 (D1hat, LDD, n, eta)
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
      end do

      open(10)
      do i = 0, n
        write (10,10) eta(i), u(i), d1u(i), d2u(i)
      end do 
 10   format (1x,4(ES16.8E3,4x))
      close(10)

      return
      end

C***********************************************************************
      subroutine MAKE_DERIVATIVES
C***********************************************************************
C
C     Make the required matrices that take derivatives in Finite diff
C     space
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=512)
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
     .                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      integer  m, p, LDD, i, j, k
      real     Identity(0:idim,0:idim), D1hat(0:idim,0:idim)
      
      LDD = idim
      call FiniteD1 (D1hat, LDD, n, eta)
      call FiniteD1 (D1   , LDD, n, eta)
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

c     CALL WRRRN ('D1', N+1, N+1, D1, LDd+1, 0)
c     CALL WRRRN ('D2', N+1, N+1, D2, LDd+1, 0)
c     CALL WRRRN ('D3', N+1, N+1, D3, LDd+1, 0)
c     CALL WRRRN ('D4', N+1, N+1, D4, LDd+1, 0)
      
      return
      end

C***********************************************************************
      subroutine MAKE_BL_METRICS
C***********************************************************************
C
C     Setup the collocation points in the mapped coordinate, eta and in 
C     chebyshev space, th.  Also compute the transformation metrics.
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=512)
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
     .                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      real       pi, deta
      integer    i
      
      pi  = ACOS(-1.0)
      deta = 2./float(n)
c
c     Make mesh in transformed, eta
c
      do i = 0, n
        eta(i) = 1. - float(i)*deta
      end do
c
c     Make transformation metrics
c
      do i = 0, n
        m1(i) = (eta(i)-1.)**2/(2.*Lmap)
        m2(i) = (eta(i)-1.)**3/(2.*Lmap**2)
        m3(i) = 3.*(eta(i)-1.)**4/(4.*Lmap**3)
        m4(i) = 3.*(eta(i)-1.)**5/(2.*Lmap**4)
        if (i.eq.0) then
          c(i) = 2.
        else
          c(i) = 1.
        end if
      end do

      return
      end

C***********************************************************************
      SUBROUTINE FUNCD(X,F,DF)
C***********************************************************************
      REAL          X, ETA, GAMMA, ETAOUT
      COMMON        /map/  ETA

      GAMMA = 1.2
      ETAOUT  = 15.0

      F = ETAOUT*(1.-TANH(GAMMA))/2.*(X+1)/
     .           (1.-TANH(GAMMA/2.*(X+1)))-ETA
      A = ETAOUT/2.*(1-TANH(GAMMA))
      B = GAMMA/2.
      DF = A/2.*(1.+EXP(2.*B*(1.+X)) + 2.*B*EXP(2.*B*(1.+X)) +
     .           2.*B*EXP(2.*B*(1.+X))*X)
      RETURN
      END

C***********************************************************************
      subroutine INIT_BL_PROFILE
C***********************************************************************
C
C     Setup the initial BL profile
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=512)
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
     .                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      parameter     (nu = 64)
      real          utemp(0:nu),utemp1(0:nu),utemp2(0:nu),junk, pi
      real          gamma, yout, y, x, xi, IU, IV, IW
      integer       nmode, i, merge
      character*15  filename
      integer       LDD
      real          D1hat(0:idim,0:idim), D2hat(0:idim,0:idim)
      
      real          RTNEWT
      external      RTNEWT  

      common        /map/  x
             
      pi = ACOS(-1.0)   
      LDD = idim     
      gamma = 1.2
      xout  = 15.0

      write (*,5)
   5  format (/,1x,'Read Mean Profile',/)
      write (*,9)
   9  format (1x,'Enter filename ==> ',$)
      read (*,'(a)') filename

      open (unit=11,file=filename,status='unknown')

      read (11,*) nmode
      do i = 0, nmode
        read (11,*) y,utemp(i),junk
      end do
      close (11)
c
c     I need to spectrally interpolate this profile onto the new grid.
c 
      call CHEBYSHEV (utemp,nmode,1)
c
c     y is nondimensionalized by dr = sqrt(2*x*nu/u_inf)
c
      i = n
      y = Lmap*(1.+eta(i))/(1.-eta(i))
      do while (y/sqrt(2.) .le. xout)
        u(i) = 0.0
        X = Y/SQRT(2.)
        xi = RTNEWT(-2.,2.2,1e-12,FUNCD)
        do m = 0, nmode
          u(i) = u(i)+utemp(m)*COS(float(m)*ACOS(xi))
        end do
        i = i - 1
        y = Lmap*(1.+eta(i))/(1.-eta(i))
      end do
      merge = i
      write (*,11) merge
  11  format(/,1x,'Merging at i = ',i4)
      do j = merge, 0, -1
        u(j) = 1.0
      end do
c
c     Compute the collocation derivatives
c
      call FiniteD1 (D1hat, LDD, n, eta)
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
      end do
c
c     Ensure that the derivatives are zero when they are supposed to be.
c
      do j = merge-2, 0, -1
        d1u(j) = 0.0
        d2u(j) = 0.0
      end do

      do i = 0, n
        write (*,10) eta(i), u(i), d1u(i), d2u(i)
      end do 
 10   format (1x,4(e16.8,4x))
 
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
      parameter   (idim=512)
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
     .                     alpha,omega,type
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
C
C     Include the independent variable transformation metrics
C      
      do i = 0, n
        work = m1(i)**4
        do j = 0, n
          A4(i,j) = work*D4(i,j)
        end do
      end do

      do i = 0, n
        work = 6.*m1(i)**2*m2(i)
        do j = 0, n
          A3(i,j) = work*D3(i,j)
        end do
      end do

      do i = 0, n
        work = 3.*m2(i)**2+4.*m1(i)*m3(i)-2.*alpha**2*m1(i)**2-
     .         cmplx(0.,1.)*alpha*Re*U(i)*m1(i)**2
        do j = 0, n
          A2(i,j) = work*D2(i,j)
        end do
      end do

      do i = 0, n
        work = m4(i)-2.*alpha**2*m2(i)-cmplx(0.,1.)*alpha*Re*U(i)*m2(i)
        do j = 0, n
          A1(i,j) = work*D1(i,j)
        end do
      end do

      do i = 0, n
        work = alpha**4+cmplx(0.,1.)*alpha**3*Re*U(i)+
     .         cmplx(0.,1.)*alpha*Re*(d2u(i)*m1(i)**2+d1u(i)*m2(i))
        do j = 0, n
          A0(i,j) = work*identity(i,j)
        end do
      end do

      do i = 0, n
        work = -1.0*cmplx(0.,1.)*Re*m1(i)**2
        do j = 0, n
          B2(i,j) = work*D2(i,j)
        end do
      end do

      do i = 0, n
        work = -1.0*cmplx(0.,1.)*Re*m2(i)
        do j = 0, n
          B1(i,j) = work*D1(i,j)
        end do
      end do

      do i = 0, n
        work = cmplx(0.,1.)*alpha**2*Re
        do j = 0, n
          B0(i,j) = work*identity(i,j)
        end do
      end do
                         
      return
      end

C***********************************************************************
      function CHECKEIG(N,A,LDA,EVAL,EVEC)
C***********************************************************************
C
C     Check an eigenvalue and eigenvector
C
C***********************************************************************
      integer    N
      complex    A(LDA,N), EVAL, EVEC(N)
      complex    X(N), Y(N)
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
      subroutine SOLVE_ORR_SOM(eigenvalue, eigenvector, print)
C***********************************************************************
C
C     This routine generates the discrete eigenvalue problem in 
C     Chebyshev space for the Orr-Sommerfeld equation.
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=512)
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
     .                     alpha,omega,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0
c***********************************************************************
      complex     A(0:idim,0:idim), B(0:idim,0:idim), alp(0:idim)
      complex     T1(0:idim,0:idim), T2(0:idim,0:idim)
      complex     T3(0:idim,0:idim), T4(0:idim,0:idim)
      complex     beta(0:idim), evec(0:idim,0:idim), temp3(0:idim)
      complex     evec2(0:idim,0:idim), eval(0:idim)
      complex     avec(0:idim,0:idim)
      complex     P1(0:idim,0:idim), P2(0:idim,0:idim)
      complex     eigenvalue, eigenvector(0:idim)
      real        temp1(0:idim), temp2(0:idim), residual
      integer     lda, ldb, ldevec, ldavec, p, i, j, which, k, l, m
      real        rind(0:idim)
      integer     ind(0:idim)
      logical     first, print
      external    CHECKEIG
      complex     escale

      integer     INFO, LWORK, IPIV(0:idim)
      complex     WORK(16*(idim+1))
      real        RWORK(8*(idim+1))
C***********************************************************************
      LDA = idim+1
      LDB = idim+1
      LDEVEC = idim+1
      LDAVEC = idim+1
      
      do i = 1, n-1
        do j = 1, n-1
          A(i-1,j-1) = A4(i,j)+A3(i,j)+A2(i,j)+A1(i,j)+A0(i,j)
          B(i-1,j-1) = B2(i,j)+B1(i,j)+B0(i,j)
        end do
      end do
c
c     Now enforce the v = 0 @ +_ 1 boundary condition instead of the
c     Orr-Sommerfeld equation at n = 0 and n = N.  The Neumann 
c     boundary condition is implicitly enforced.
c
c      do j = 0, n
c        B(0,j) = CMPLX(0.0, 0.0)
c        B(n,j) = CMPLX(0.0, 0.0)
c        A(0,j) = CMPLX(0.0, 0.0)
c        A(n,j) = CMPLX(0.0, 0.0)
c      end do
c      A(0,0) = CMPLX(1.0, 0.0) 
c      A(n,n) = CMPLX(1.0, 0.0)
c
c     But the top and bottom rows are trivial so that they can be
c     removed
c
#ifdef USE_IMSL
c     write(*,*) "Solving Eigenvalue Problem with ISML interface"
      CALL LINCG (N-1, B, LDA, T1, LDA)
      CALL MCRCR (N-1, N-1, T1, LDA, N-1, N-1, A, LDA,
     .            N-1, N-1, T4, LDA)
      CALL EVCCG (N-1, T4, LDA, eval, evec, ldevec)
#else
#ifdef OS_USE_GENERALIZED_EIGENSOLVER
c
c.... use generalized eigenvalue solver
c
      LWORK = 16*(idim+1)
      INFO = 0
      call ZGEGV( 'V', 'V', N-1, A, LDA, B, LDA, alp, beta, avec,
     .             LDA, evec, LDA, work, lwork, rwork, info)
c
c.... compute the eigenvalues (given generalized eigenvalues)
c
      do i = 0, N-2
        if (beta(i).ne.0) then
          eval(i) = alp(i)/beta(i)
        else
          eval(i) = 0.0
        end if
      end do
#else
      LWORK = 2*(N-1)
      INFO = 0
      CALL ZLACPY('A', N-1, N-1, B, LDA, T1, LDA)
      CALL ZGETRF(N-1, N-1, T1, LDA, IPIV, INFO)
      CALL ZGETRI(N-1, T1, LDA, IPIV, WORK, LWORK, INFO)
      CALL ZGEMM('N','N',N-1,N-1,N-1,1.0,T1,LDA,A,LDA,0.0,T4,LDA)
      CALL ZGEEV('V','V',N-1, T4, LDA, eval, avec, LDAVEC, 
     &           evec, LDEVEC, WORK, LWORK, RWORK, INFO)
#endif
#endif
      do i = 0, N-2 
        temp1(i) = REAL(eval(i))
        temp2(i) = AIMAG(eval(i))
      end do

c     Need to issolate the most unstable eigenvalue and eigenvector
c     Must watch out, this routine isn't very robust.

      do i = 0,N-2
        rind(i) = REAL(i)
      end do
      call SORT(n-1,temp2,rind)
      do i = 0,N-2
        ind(i) = INT(rind(i))
      end do
      first = .true.
      do i = n-2,0,-1
        if ((abs(temp2(i)).lt. 10.0).and.(first).and.
     &      (int(temp1(ind(i))).ne.999)) then
          first = .false.
          eigenvalue = cmplx(temp1(ind(i)),temp2(i))
          iloc = i
        end if
      end do
#ifdef ORR_CHECK_EIG 
      do j = 0, n-2
        eigenvector(j) = evec(j,ind(iloc))
      end do
      residual = CHECKEIG(N-1,T4,lda,eigenvalue,eigenvector)
      if (residual .gt. .01) then
        write (*,*) 'WARNING eigenvalue not converged!'
      end if
#endif
      if (print) then
        write (*,100) residual
 100    format (/,1x,'Residual = ',f17.10)
        write (*,31) Re, real(alpha), aimag(alpha)
        write (*,32)
        do i = 1, N
          write (*,35) REAL(eval(i)), AIMAG(eval(i)), i
          write (20,36) REAL(eval(i))/alpha, AIMAG(eval(i))/alpha
        end do
        write (*,40)
        read (*,*) which
        do while (which .ne. 0)
#if 1
c
c.... normalize the eigenvector for output
c
          escale = 0.0
          do i = 1, N-1
            if (abs(evec(i-1,ind(iloc))) .gt. abs(escale)) then
              escale = evec(i-1,ind(iloc))
            endif
          end do
          write(*,*) "escale = ", escale
          escale = 1.0/(escale)
          call zscal( N-1, escale, evec(0,ind(iloc)), 1)
#endif
          temp1(0) = 0.0
          temp2(0) = 0.0
          temp1(n) = 0.0
          temp2(n) = 0.0
          do j = 1, n-1
            temp1(j) = REAL(evec(j-1,ind(iloc)))
            temp2(j) = AIMAG(evec(j-1,ind(iloc)))
          end do
          open(11)
          do j = 0, N
            write (11,60) eta(j), temp1(j), temp2(j), 
     &           REAL(avec(j,ind(iloc))), AIMAG(avec(j,ind(iloc)))
          end do
          close(11)
          write (*,40)
          read (*,*) which
        end do
      end if

 31   format(/,1x,3(e17.10,4x))
 32   format(/,1x,'       w_r                w_i              count',/)
 35   format(1x,2(e17.10,4x),i5)
 36   format(1x,2(e17.10,4x))
 40   format(/,1x,'Eigenvector for which eigenvalue (0 quits) ==> ',$)
 55   format(/,1x,'Performance Index = ',g12.5,/)
 60   format(1x,6(ES16.8E3,1x))
     
      return
      end

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

      PARAMETER (idim=512)
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
      PARAMETER (idim=512)
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

      do i = 0, NBIG
        X = I*PI/NBIG
        IU = 0.0
        IV = 0.0
        DO M = 0, Nmode
          IU = IU + U(M)*COS(FLOAT(M)*X)
          IV = IV + V(M)*COS(FLOAT(M)*X)
        END DO
        write (33,10) cos(x),Iu,Iv
      end do
  10  format (1x,3(e16.8,4x))

      RETURN
      END

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
          YY(I) = TT(0)/2.D0 + TT(N)*DCOS(DFLOAT(I)*PI)/2.D0
          DO M = 1, N-1
            YY(I) = YY(I) + TT(M)*DCOS(DFLOAT(M)*DFLOAT(I)*PI/DFLOAT(N))
          END DO
          YY(I) = YY(I)*2./DFLOAT(N)
        end do
        YY(0) = YY(0)/2.D0
        YY(N) = YY(N)/2.D0
      else 
        DO I = 0, N
          YY(I) = 0.0
          DO M = 0, N
            YY(I) = YY(I) + TT(M)*DCOS(DFLOAT(M)*DFLOAT(I)*PI/DFLOAT(N))
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

c     call COSFT (y,n,isign)
      
c     if (isign .eq. 1) then
c       do i = 0,n
c         y(i) = y(i)/float(n)*2.
c       end do
c     end if
      
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
     .                         DBLE(C(ABS(M)))   * B(ABS(M),J)
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
            
C******************************************************************************
      SUBROUTINE COSFT3D (YY,N,ISIGN)
C******************************************************************************
C
C     Do a brute force transform.  ISIGN = 1 does a forward transform.
C     ISIGN = -1 a backward transform.  This one is fully double precision.
C
C******************************************************************************
      integer     N, ISIGN
      real*16     YY(0:N), TT(0:N), PI

      PI = DACOS(-1.0D0)

      DO I = 0, N
        TT(I) = YY(I)
      END DO

      if (isign .eq. 1) then
        DO I = 0, N
          YY(I) = TT(0)/2.D0 + TT(N)*DCOS(DFLOAT(I)*PI)/2.D0
          DO M = 1, N-1
            YY(I) = YY(I) + TT(M)*DCOS(DFLOAT(M)*DFLOAT(I)*PI/DFLOAT(N))
          END DO
          YY(I) = YY(I)*2./DFLOAT(N)
        end do
        YY(0) = YY(0)/2.D0
        YY(N) = YY(N)/2.D0
      else 
        DO I = 0, N
          YY(I) = 0.0
          DO M = 0, N
            YY(I) = YY(I) + TT(M)*DCOS(DFLOAT(M)*DFLOAT(I)*PI/DFLOAT(N))
          END DO
        END DO
      end if
      
      RETURN
      END
      
C******************************************************************************
      SUBROUTINE FiniteD1 (D, LDD1, N, eta)
C******************************************************************************
C
C     Calculation the finite difference derivative matrix
C
C******************************************************************************
      PARAMETER (IDIM=512)
      REAL      D(0:LDD1,0:N), PI, eta(0:n)
      REAL      C(0:IDIM), DX

      IF (N .GT. IDIM) THEN
        WRITE (*,*) 'ERROR:  N > IDIM in FiniteD'
        STOP
      END IF
      
      DX = 2./FLOAT(N)
      
      DO I = 0,N
        DO J = 0,N
          D(I,J) = 0.0
        END DO
      END DO
      D(0,0) =  1./(eta(0)-eta(1))
      D(0,1) =  -1./(eta(0)-eta(1))
      D(N,N-1) =  1./(eta(N-1)-eta(N))
      D(N,N) =  -1./(eta(N-1)-eta(N))
      DO I = 1, N-1
        D(I,I-1) =  1./(eta(I-1)-eta(i+1))
        D(I,I+1) = -1./(eta(I-1)-eta(i+1))
      END DO

      RETURN
      END

C******************************************************************************
      SUBROUTINE CHEBYD (D, LDD1, N)
C******************************************************************************
C
C     Calculation the Chebyshev collocation derivative matrix
C
C******************************************************************************
      PARAMETER (IDIM=512)
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
     .                (2.*(1.-(COS(FLOAT(J)*PI/FLOAT(N)))**2))
          ELSE
            D(J,K) = C(J)*(-1.0)**(J+K)/(C(K)*(COS(FLOAT(J)*PI/
     .               FLOAT(N))-COS(FLOAT(K)*PI/FLOAT(N))))
          END IF
        END DO
      END DO

      RETURN
      END
