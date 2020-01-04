c***********************************************************************
      program conte_chan
c***********************************************************************
c 
c     Purpose:  Solve the Orr-Sommerfeld equation using Runge-
c               Kutta explicit integration with pseudo-orthogonalization
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

      common      /setup/ n, alpha, u, d2u, Re, ymin, ymax, h
     
      common      /eig/   c

      complex     INPROD
      external    OSHOMO, OSPART, INPROD
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
      call SET_MEAN_FLOW 
c
c     set the boundary conditions
c
      do i = 1, neq
        bc1(i) = 0.0
        bc2(i) = 0.0
      end do
      
      call CONTE(n, neq, 2, bc1, bc2, ymin, ymax, c, 
     &           OSHOMO, OSPART, INPROD)

 210  continue      
      stop
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
      subroutine OSHOMO(neq,yo,t,yf)
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

      common      /setup/ n, alpha, u, d2u, Re, ymin, ymax, h
     
      common      /eig/   c
c***********************************************************************
      complex     yo(4), yf(4)
      real        t
      integer     neq
      real        UU, dUU, d2UU
c
c     Get the velocity field
c
      call CHANNEL(t,UU,dUU,d2UU)

      do j = 1 , 3
        yf(j) = yo(j+1)
      end do
      yf(4) = 2.*alpha**2*yo(3)-alpha**4*yo(1) + 
     &        (0,1.)*alpha*Re*((UU-c)*(yo(3)-alpha**2*yo(1))-
     &        d2UU*yo(1)) 

      return
      end

C***********************************************************************
      subroutine OSPART(neq,yo,t,yf)
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

      common      /setup/ n, alpha, u, d2u, Re, ymin, ymax, h
     
      common      /eig/   c
c***********************************************************************
      complex     yo(4), yf(4)
      real        t
      integer     neq
      real        UU, dUU, d2UU
c
c     Get the velocity field
c
      call CHANNEL(t,UU,dUU,d2UU)

      do j = 1 , 3
        yf(j) = yo(j+1)
      end do
      yf(4) = 2.*alpha**2*yo(3)-alpha**4*yo(1) + 
     &        (0,1.)*alpha*Re*((UU-c)*(yo(3)-alpha**2*yo(1))-
     &        d2UU*yo(1)) 

      return
      end

C***********************************************************************
      subroutine SET_MEAN_FLOW 
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

      common      /setup/ n, alpha, u, d2u, Re, ymin, ymax, h
     
      common      /eig/   c
c***********************************************************************
      integer     i, j, k, p
      real        xi(3,imax), f(3), eta(3), y, pi, us, dus, d2us
      real        k1(3), k2(3), k3(3), k4(3), err, x2old, f1old

      pi = acos(-1.0)
c
c     Now assemble the velocity field
c
      do i = 0, n
        y = ymin + i*h
        call CHANNEL(y,u(i),dudy,d2u(i))
      end do
c
c     Write profile on 2x resolve mesh
c
c     open (10, FILE='chan.dat', ACTION='WRITE', FORM='FORMATTED')
      do i = 0, 2*n
        y = ymin + i*h/2.0
        call CHANNEL(y,us,dus,d2us)
        write (10,10) y, us, d2us
  10    format (1x,5(ES16.8E3,1x))
      end do
c     close(10)

      write (*,20)
  20  format (1x,'Channel velocity profile completed...',/)

      return
      end

C***********************************************************************
      subroutine CHANNEL(y,u,dudy,d2udy2)
C***********************************************************************
C
C     Incompressible channel flow profile
C
C***********************************************************************
      real y, u, dudy, d2udy2
C***********************************************************************
      u      = 1.0 - y**2
      dudy   = -2.0*y
      d2udy2 = -2.0
      return
      end
