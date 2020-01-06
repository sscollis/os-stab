C> \file orrucbl.f
C> \brief Solves Orr-Sommerfeld using method of Deissler 
C> \author S. Scott Collis
c***********************************************************************
      program orr_uc_bl
c***********************************************************************
c
c     Purpose:  This program solves the Orr-Sommerfeld equation using a 
c               Chebyshev-collocation method.
c
c     Version:  This version solves the boudary layer problem with
c               algebraic mapping and allows for moving reference frame
c               so as to determine absolute/convective instability like
c               Deissler.
c
c     Author:   S. Scott Collis
c
c     Date:     3-14-92
c
c     Revision: 9-18-92
c
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=256)
      integer     n
      real        eta(0:idim), th(0:idim), lmap, uorg(0:idim)
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim), uc
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     dA2(0:idim,0:idim), dA1(0:idim,0:idim)
      complex     dA0(0:idim,0:idim), dB0(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,uc,uorg,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0, dA2, dA1,
     .                     dA0, dB0
c***********************************************************************
      integer     nar, nai, ians
      real        alphar(100), alphai(100), delta
      real        minar, maxar, incar, minai, maxai, incai, dar, dai
      complex     eigenvalue(100), eigenvector(0:idim), ctemp, A, B
      complex     alp(3), as, oldeig(3,100), fd
      real        oldeigr1, oldeigi1
      real        oldeigr2, oldeigi2
      character*9 filename
      logical     print
c***********************************************************************
c
c     Setup IMSL workspace
c
#ifdef USE_IMSL
      REAL             RWKSP(100000)
      COMMON /WORKSP/  RWKSP
      call IWKIN(100000)
#endif
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
      write (*,65)
  65  format (/,1x,'Enter reference frame velocity u_c ==> ',$)
      read (*,*) uc
      Lmap = 2.0
      write (*,40) 
  40  format (/,1x,'Enter alpha_1 (r,i) ==> ',$)
      read (*,*) alphar(1),alphai(1)
      write (*,45)
  45  format (/,1x,'Enter delta ==> ',$)
      read (*,*) delta
      write (*,46)
  46  format (/,1x,'Enter deltau ==> ',$)
      read (*,*) deltau
      write (*,50)
  50  format (/,1x,'New run (1/0) ==> ',$)
      read (*,*) irun
      if (irun.eq.0) then
        write (*,70) 
  70    format (/,1x,'Enter oldeig1 (r,i) ==> ',$)
        read (*,*) oldeigr1,oldeigi1
        write (*,85) 
  85    format (/,1x,'Enter oldeig2 (r,i) ==> ',$)
        read (*,*) oldeigr2,oldeigi2
        do i = 1, 3
          oldeig(i,1) = CMPLX(oldeigr1,oldeigi1)
          oldeig(i,2) = CMPLX(oldeigr2,oldeigi2)
        end do
        istart = 3
      else
        istart = 1
      end if
      
      alphar(2) = alphar(1) + delta
      alphar(3) = alphar(1) - delta
      alphai(2) = alphai(1) - delta
      alphai(3) = alphai(1) - delta
      ians = 0
      if (ians.eq.1) then
         print = .true.
      else
         print = .false.
      end if

      call MAKE_DERIVATIVES
      call MAKE_BL_METRICS
      call INIT_BL_PROFILE
      as = cmplx(10,10)
      
      do k = istart, 40+istart-1

        resid = 10.
        if (k .le. 2)  then
          do i = 1, 3
            oldeig(i,k) = 999.
          end do
        else
          do i = 1, 3
            fd = (oldeig(i,k-1)-oldeig(i,k-2))/deltau
            oldeig(i,k) = oldeig(i,k-1) + fd*deltau
          end do
        end if
        call SHIFT_BL
        icount = 0
        do while (abs(resid) .gt. 1e-7 )
          icount = icount + 1
          do i = 1, 3
            aso = as
            alpha = cmplx (alphar(i), alphai(i))
            alp(i) = alpha
            call MAKE_MATRIX
            ctemp = oldeig(i,k)
            call SOLVE(ctemp, eigenvector, print)
            eigenvalue(i) = ctemp
            oldeig(i,k) = ctemp
c           write (*,60)  real(alpha),aimag(alpha),real(eigenvalue(i)),
c    .                    aimag(eigenvalue(i))
  60        format (1x,4(e17.10,4x))
          end do
          A  = cmplx(0.,-1.)*(eigenvalue(2)-eigenvalue(1))
          B  = cmplx(0.,-1.)*(eigenvalue(3)-eigenvalue(1))
          as = (A*(alp(3)**2-alp(1)**2)-B*(alp(2)**2-alp(1)**2))/
     .         (2.*B*(alp(1)-alp(2))-2.*A*(alp(1)-alp(3)))
          resid = as-aso
          if (k .eq. 1 .or. icount .gt. 5) then
            write (*,80) real(as),aimag(as), abs(resid)
  80        format (1x,'==> ',3(e17.10,4x))
          end if
          alphar(1) = real(as)
          alphai(1) = aimag(as)
          alphar(2) = alphar(1) + delta
          alphar(3) = alphar(1) - delta
          alphai(2) = alphai(1) - delta
          alphai(3) = alphai(1) - delta
        end do
        write (*,90)  uc, real(eigenvalue(1)), 
     .                aimag(eigenvalue(1)),
     .                real(as),aimag(as)
 90     format (1x,f6.4,2x,4(e17.10,4x))
        uc = uc - deltau
      end do

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
      parameter   (idim=256)
      integer     n
      real        eta(0:idim), th(0:idim), lmap, uorg(0:idim)
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim), uc
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     dA2(0:idim,0:idim), dA1(0:idim,0:idim)
      complex     dA0(0:idim,0:idim), dB0(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,uc,uorg,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0, dA2, dA1,
     .                     dA0, dB0
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
      subroutine MAKE_BL_METRICS
C***********************************************************************
C
C     Setup the collocation points in the mapped coordinate, eta and in 
C     chebyshev space, th.  Also compute the transformation metrics.
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=256)
      integer     n
      real        eta(0:idim), th(0:idim), lmap, uorg(0:idim)
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim), uc
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     dA2(0:idim,0:idim), dA1(0:idim,0:idim)
      complex     dA0(0:idim,0:idim), dB0(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,uc,uorg,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0, dA2, dA1,
     .                     dA0, dB0
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
      parameter   (idim=256)
      integer     n
      real        eta(0:idim), th(0:idim), lmap, uorg(0:idim)
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim), uc
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     dA2(0:idim,0:idim), dA1(0:idim,0:idim)
      complex     dA0(0:idim,0:idim), dB0(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,uc,uorg,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0, dA2, dA1,
     .                     dA0, dB0
c***********************************************************************
      parameter     (nu = 64)
      real          utemp(0:nu),utemp1(0:nu),utemp2(0:nu),junk, pi
      real          gamma, yout, y, x, xi, IU, IV, IW
      integer       nmode, i, merge
      character*15  filename
      integer       LDD
      real          D1hat(0:idim,0:idim), D2hat(0:idim,0:idim)  

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
c      write (*,11) merge
  11  format(/,1x,'Merging at i = ',i4)
      do j = merge, 0, -1
        u(j) = 1.0
      end do
c
c     Compute the collocation derivatives
c
      call CHEBYD (D1hat, LDD, n)
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
        Uorg(i) = u(i)
      end do
c
c     Ensure that the derivatives are zero when they are supposed to be.
c
      do j = merge-2, 0, -1
        d1u(j) = 0.0
        d2u(j) = 0.0
      end do
c
c      call CHEBYSHEV (U,n,1)
c      call CHEBYSHEV (d1u,n,1)
c      call CHEBYSHEV (d2u,n,1)
c
c      do i = 0, 256
c        X = float(I)*PI/256
c        IU = 0.0
c        IV = 0.0
c        IW = 0.0
c        DO M = 0, N
c          IU = IU +   U(M)*COS(FLOAT(M)*X)
c          IV = IV + d1u(M)*COS(FLOAT(M)*X)
c          IW = IW + d2u(M)*COS(FLOAT(M)*X)
c        END DO
c        write (*,10) cos(x),Iu,Iv,Iw
c      end do
c
c      call CHEBYSHEV (U,n,-1)
c      call CHEBYSHEV (d1u,n,-1)
c      call CHEBYSHEV (d2u,n,-1)

      return
      end
C***********************************************************************
      subroutine SHIFT_BL
C***********************************************************************
C
C     Shift the boundary layer profile to a new frame of reference.
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=256)
      integer     n
      real        eta(0:idim), th(0:idim), lmap, uorg(0:idim)
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim), uc
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     dA2(0:idim,0:idim), dA1(0:idim,0:idim)
      complex     dA0(0:idim,0:idim), dB0(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,uc,uorg,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0, dA2, dA1,
     .                     dA0, dB0
c***********************************************************************
c
c     Need to interpolate the velocity profile to evaluate it at
c     arbitrary y
c
      do i = 0, n
        U(i) = Uorg(i) - Uc
      end do

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
      parameter   (idim=256)
      integer     n
      real        eta(0:idim), th(0:idim), lmap, uorg(0:idim)
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim), uc
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     dA2(0:idim,0:idim), dA1(0:idim,0:idim)
      complex     dA0(0:idim,0:idim), dB0(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,uc,uorg,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0, dA2, dA1,
     .                     dA0, dB0
c***********************************************************************
      real      identity(0:idim,0:idim)
      complex   work, ai
      
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
      subroutine HTRAN(N,A,B,LD)
C***********************************************************************
C
C     Take the complex conjugate transpose of A and put it in B
C
C***********************************************************************
      integer    N, LD
      complex    A(LD,N), B(LD,N)
      
      do i = 1, N
        do j = 1, N
          B(j,i) = CONJG(A(i,j))
        end do
      end do

      RETURN
      END
C***********************************************************************
      subroutine CXDOTY (N,X,Y,C)
C***********************************************************************
C
C     Take the complex conjugate dot product of vector x and y
C
C***********************************************************************
      integer    N
      complex    X(N), Y(N), C
      
      C = CMPLX (0.0,0.0)
      do I = 1, N
        C = C + CONJG(X(I))*Y(I)
      end do

      RETURN
      END
C***********************************************************************
C               S O L V E   O R R   S O M M E R F E L D 
C***********************************************************************
      subroutine SOLVE(eigenvalue,eigenvector,print)
C***********************************************************************
C
C     This routine generates the discrete eigenvalue problem in 
C     Chebyshev space for the Orr-Sommerfeld equation.
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=256)
      integer     n
      real        eta(0:idim), th(0:idim), lmap, uorg(0:idim)
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim), uc
      real        c(0:idim), u(0:idim), d1u(0:idim), d2u(0:idim), Re
      real        D1(0:idim,0:idim), D2(0:idim,0:idim)
      real        D3(0:idim,0:idim), D4(0:idim,0:idim)
      complex     A1(0:idim,0:idim), A2(0:idim,0:idim)
      complex     A3(0:idim,0:idim), A4(0:idim,0:idim)
      complex     A0(0:idim,0:idim), B0(0:idim,0:idim)
      complex     B1(0:idim,0:idim), B2(0:idim,0:idim)
      complex     dA2(0:idim,0:idim), dA1(0:idim,0:idim)
      complex     dA0(0:idim,0:idim), dB0(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,uc,uorg,type
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0, dA2, dA1,
     .                     dA0, dB0
c***********************************************************************
      complex     A(0:idim,0:idim), B(0:idim,0:idim)
      complex     dA(0:idim,0:idim), dB(0:idim,0:idim)
      complex     T1(0:idim,0:idim), T2(0:idim,0:idim)
      complex     T3(0:idim,0:idim), T4(0:idim,0:idim)

      complex     eval(0:idim), evec(0:idim,0:idim), tvec(0:idim)
      complex     aeval(0:idim), aevec(0:idim,0:idim), tavec(0:idim) 

      complex     eigenvalue, eigenvector(0:idim)
      complex     dwda, dw, dalp, y(0:idim), prod

      real        temp1(0:idim), temp2(0:idim), residual, CHECKEIG, min
      real        temp3(0:idim), temp4(0:idim), dalpha, diff(0:idim)

      integer     lda, ldb, ldevec, p, i, j, which, k, l, m
      integer     index(0:idim), index2(0:idim)
      logical     first, print
      external    CHECKEIG
C***********************************************************************
      lda = idim+1
      ldb = idim+1
      ldevec = idim+1
      
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
      write(*,*) "Solving Eigenvalue Problem with IMSL interface"
      CALL LINCG (N-1, B, LDA, T1, LDA)
      CALL MCRCR (N-1, N-1, T1, LDA, N-1, N-1, A, LDA,
     .            N-1, N-1, T4, LDA)

      CALL EVCCG (N-1, T4, LDA, eval, evec, ldevec)
#else
      write(*,*) "Solving Eigenvalue Problem with new LAPACK interface"
      LWORK = 2*(N-1)
      INFO = 0
      CALL ZLACPY('A', N-1, N-1, B, LDA, T1, LDA)
      CALL ZGETRF(N-1, N-1, T1, LDA, IPIV, INFO)
      CALL ZGETRI(N-1, T1, LDA, IPIV, WORK, LWORK, INFO)
      CALL ZGEMM ('N','N',N-1,N-1,N-1,1.0,T1,LDA,A,LDA,0.0,T4,LDA)
      CALL ZGEEV ('N','V',N-1, T4, LDA, eval, VL, LDVL, evec, ldevec,
     .            WORK, LWORK, RWORK, INFO)
#endif
c
      do i = 0, N-2 
        temp1(i) = REAL(eval(i))
        temp2(i) = AIMAG(eval(i))
      end do

c     Need to issolate the most unstable eigenvalue and eigenvector
c     Must watch out, this routine isn't very robust.

      do i = 0,N-2
        index(i) = i
      end do
      call PIKSR2(n-1,temp2,index)

      min = 10000.
      do i = 0, n-2
        temp1(i) = REAL(eval(index(i)))
        if ( ABS(ABS(eigenvalue)-ABS(cmplx(temp1(i),temp2(i))))
     .      .le. min ) then
          iloc = i
          min = ABS(ABS(eigenvalue)-ABS(cmplx(temp1(i),temp2(i))))
        end if
      end do
      
      if (ABS(eigenvalue).gt.99.) then
        eigenvalue = cmplx(temp1(n-2),temp2(n-2))
        iloc = n-2
      else
        eigenvalue = cmplx(temp1(iloc),temp2(iloc))
      end if
      
c     write (*,*) 'Eigenvalue at ',iloc
      eigenvector(0) = cmplx(0.0,0.0)
      eigenvector(n) = cmplx(0.0,0.0)
      
      do j = 1, n-1
        eigenvector(j) = evec(j-1,index(iloc))
        tvec(j-1)  = evec(j-1,index(iloc))
      end do
      
      residual = CHECKEIG (N-1,T4,lda,eigenvalue,tvec)
      if (residual .gt. .01) then
        write (*,*) 'WARNING eigenvalue not converged!'
      end if
      
      if (print) then
        write (*,100) residual
 100    format (/,1x,'Residual = ',e17.10)
        write (*,31) Re, real(alpha), aimag(alpha)
        write (*,32)
        do i = n-2,0,-1
          write (*,37) REAL(eval(index(i))/alpha),
     .                 AIMAG(eval(index(i))/alpha),index(i)
c          write (20,37) REAL(eval(index(i))),AIMAG(eval(index(i))),
c     .                  REAL(aeval(index(i))),AIMAG(aeval(index(i))),
c     .                  index(i)
        end do
        write (*,40)
        read (*,*) which
          do while (which .ne. 0)
          temp1(0) = 0.0
          temp2(0) = 0.0
          temp1(n) = 0.0
          temp2(n) = 0.0
          temp3(0) = 0.0
          temp4(0) = 0.0
          temp3(n) = 0.0
          temp4(n) = 0.0
          do j = 1, n-1
            temp1(j) = REAL(evec(j-1,index(iloc)))
            temp2(j) = AIMAG(evec(j-1,index(iloc)))
            temp3(j) = REAL(aevec(j-1,index(iloc)))
            temp4(j) = AIMAG(aevec(j-1,index(iloc)))
          end do
c          do j = 0, n
c            write (*,38) eta(j),temp1(j),temp2(j),temp3(j),temp4(j)
c          end do
          write (*,*)
          call CHEBYINT (n, temp1, temp2, 128)
c          write (*,*)
c          call CHEBYINT (n, temp3, temp4, 128)
          write (*,40)
          read (*,*) which
        end do
      end if

 31   format(/,1x,3(e17.10,4x))
 32   format(/,1x,'       w_r                w_i              count',/)
 35   format(1x,2(e17.10,4x),i5)
 36   format(1x,2(e17.10,4x))
 37   format(1x,2(e17.10,2x),i5)
 38   format(1x,5(e15.8,1x))
 40   format(/,1x,'Eigenvector for which eigenvalue (0 quits) ==> ',$)
 55   format(/,1x,'Performance Index = ',g12.5,/)
      
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

      PARAMETER (idim=256)
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
      PARAMETER (idim=256)
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
        write (*,10) cos(x),Iu,Iv
      end do
  10  format (1x,3(e16.8,4x))

      RETURN
      END

C***********************************************************************
      SUBROUTINE COSFT3 (Y,N,ISIGN)
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
            
C***********************************************************************
      SUBROUTINE COSFT3D (YY,N,ISIGN)
C***********************************************************************
C
C     Do a brute force transform.  ISIGN = 1 does a forward transform.
C     ISIGN = -1 a backward transform.  This one is fully double
C     precision.
C
C***********************************************************************
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
      
C***********************************************************************
      SUBROUTINE CHEBYD (D, LDD1, N)
C***********************************************************************
C
C     Calculation the Chebyshev collocation derivative matrix
C
C***********************************************************************
      PARAMETER (IDIM=256)
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
