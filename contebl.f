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
      complex     alpha, c
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   beta, uspl, d2uspl, ydat
     
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

      complex   INPROD
      external  OS_IC, OS_HOMO, OS_PART, INPROD
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
      
      call CONTE2(nstep,testalpha,neq,2,bc1,bc2,ymin,ymax,c,eigfun,
     &            OS_IC, OS_HOMO, OS_PART, INPROD)
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
      subroutine OS_IC(n,r,Uo)
C***********************************************************************
C
C     Set the initial condition at infinity for the Orr-Sommerfeld
C     equation for a boundary layer.   This works better than using
C     a simple Green's function
C
C***********************************************************************
      integer     n, r
      complex     Uo(n,n-r)
C***********************************************************************
      parameter   (imax=50000)
      integer     nbl
      complex     alpha, c
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   beta, uspl, d2uspl, ydat

      common      /eig/ c, alpha, Re
C***********************************************************************
      complex     rgamma
C***********************************************************************
      Uo(1,1) = CEXP(-alpha*ymax)
      Uo(2,1) = (-alpha)*CEXP(-alpha*ymax)
      Uo(3,1) = (alpha**2)*CEXP(-alpha*ymax)
      Uo(4,1) = (-alpha**3)*CEXP(-alpha*ymax)

      rgamma = SQRT(alpha**2+(0.,1.)*alpha*Re*(1.0-c))

      Uo(1,2) = CEXP(-rgamma*ymax)
      Uo(2,2) = (-rgamma)*CEXP(-rgamma*ymax)
      Uo(3,2) = (rgamma**2)*CEXP(-rgamma*ymax)
      Uo(4,2) = (-rgamma**3)*CEXP(-rgamma*ymax)

      return
      end

C***********************************************************************
      subroutine FHOMO_VODE(neq,t,yo,yf,rpar,ipar)
C***********************************************************************
      integer neq
      integer ipar(2)
      real rpar
      real t
      complex yo(IPAR(1),IPAR(2)), yf(IPAR(1),IPAR(2))
c     write(*,*) NEQ, IPAR(1), IPAR(2)
      do m = 1, IPAR(2) 
        call OS_HOMO(IPAR(1),yo(1,m),t,yf(1,m))
      end do
      return
      end

C***********************************************************************
      subroutine OS_HOMO(neq,yo,t,yf)
C***********************************************************************
C
C     Function evaluation for the Orr-Sommerfeld equation
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (imax=50000)
      integer     n
      complex     alpha, c
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   beta, uspl, d2uspl, ydat
     
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
      subroutine OS_PART(neq, yo,t,yf)
C***********************************************************************
C
C     Function evaluation for the Orr-Sommerfeld equation
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (imax=50000)
      integer     n
      complex     alpha, c
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   beta, uspl, d2uspl, ydat
     
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
      complex     alpha, c
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   beta, uspl, d2uspl, ydat
     
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
#elif USE_RKCK45
          call SRKCK45(3, xi(1,i-1), xi(1,i), y-h, h, BLASIUS)
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
#ifdef USE_SPLINE_DERIVATIVE
        call SPEVAL(nbl+1,ydat,u,uspl,y,us)
        call SPEVAL(nbl+1,ydat,d2u,d2uspl,y,d2us)
        write (10,10) y, us, d2us, u(i), d2u(i)
#else
        call SPDER(nbl+1,ydat,u,uspl,y,us,dus,d2us)
#ifdef USE_OUTPUT_DERIVATIVE
        write (10,10) y, us, dus, d2us, u(i), d2u(i)
#else
        write (10,10) y, us, d2us, u(i), d2u(i)
#endif
#endif
  10    format (1x,7(ES16.8E3,1x))
      end do
c     close(10)
      
      write (*,20) 
  20  format (1x,'Blasius velocity profile completed...',/)

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
      complex     alpha, c
      real        U(0:imax), d2U(0:imax), Re, ymin, ymax, h, beta
      real        Uspl(0:imax), d2Uspl(0:imax), ydat(0:imax), f2p

      common      /mean/ nbl, u, d2u, ymin, ymax, h, f2p,
     &                   beta, uspl, d2uspl, ydat
     
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
