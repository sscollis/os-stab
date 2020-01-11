c***********************************************************************
c> \file orrbfs.f
c> \brief Solves the Orr-Sommerfeld equation for incompressible
c>        backward-facing step with slip boundary on upper wall 
c> \author S. Scott Collis
c***********************************************************************
      program orr_bfs
c***********************************************************************
c
c     Purpose:  This program solves the Orr-Sommerfeld equation using a 
c               Chebyshev-collocation method for backward facing step
c               profiles provide from the DNS of Hung Le.  The
c               boundary conditions are basically that of a channel
c               with a slip wall at the upper boundary, similar to
c               the DNS boundary conditions.
c
c     Author:   S. Scott Collis
c
c     Date:     3-14-92
c
c     Revision: 1-5-2020 
c
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=256)
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
      complex     dA2(0:idim,0:idim), dA1(0:idim,0:idim)
      complex     dA0(0:idim,0:idim), dB0(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      character*20 filename
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,type,filename
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0, dA2, dA1,
     .                     dA0, dB0
c***********************************************************************
      integer     nar, nai, ians
      real        alphar(100), alphai(100)
      real        minar, maxar, incar, minai, maxai, incai, dar, dai
      complex     oldeig, dwda
      complex     eigenvalue(100,100), eigenvector(0:idim), ctemp
      logical     print
c***********************************************************************
c
c     Setup IMSL workspace
c
#ifdef USE_IMSL
      REAL             RWKSP(10000)
      COMMON /WORKSP/  RWKSP
      call IWKIN(10000)
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

      nar = AINT((maxar-minar)/incar)+1
      nai = AINT((maxai-minai)/incai)+1
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
      oldeig = cmplx(-999.,-999.)

      iloc = index(filename,' ')-1
      open (unit=10,file=filename(1:iloc)//'.g',form='unformatted',
     .      status='unknown')
      write (10) nx, ny
 100  format(i5,1x,$)
      write (10)((alphar(i), i=1,nx), j=1,ny),
     .          ((alphai(j), i=1,nx), j=1,ny)
 110  format (1x,e15.8,$)
      close (10)
      open (unit=12,file=filename(1:iloc)//'.dat',form='formatted',
     .      status='unknown')
      call MAKE_DERIVATIVES
      call INIT_BFS_PROFILE
      write (*,*)
      do j = 1, ny
        do i = 1, nx
          alpha = cmplx (alphar(i), alphai(j))
          call MAKE_MATRIX
          if (j.eq.2) then
            oldeig = eigenvalue(i,j-1)
          else if (j.gt.2) then
c            oldeig = eigenvalue(i,j-1)+(eigenvalue(i,j-1)-
c     .               eigenvalue(i,j-2))/(alphai(j-1)-alphai(j-2))*
c     .               (alphai(j)-alphai(j-1))
            oldeig = eigenvalue(i,j-1)
          else
            oldeig = cmplx(-9999.,-9999.)
          end if
          call SOLVE(ctemp, eigenvector, print, oldeig, dwda)
          eigenvalue(i,j) = ctemp
          write (12,60) real(alpha),aimag(alpha),real(eigenvalue(i,j)),
     .                  aimag(eigenvalue(i,j)),real(dwda),aimag(dwda)
          write (*,60)  real(alpha),aimag(alpha),real(eigenvalue(i,j)),
     .                  aimag(eigenvalue(i,j)),real(dwda),aimag(dwda)
  60      format (1x,6(e17.10,1x))
        end do
      end do

      open (unit=11,file=filename(1:iloc)//'.q',form='unformatted',
     .      status='unknown')
      write (11) nx, ny
      write (11) 0.0,0.0,0.0,0.0
      write (11)  ((1.0                   , i=1,nx), j=1,ny),
     .            ((real(eigenvalue(i,j)) , i=1,nx), j=1,ny),
     .            ((aimag(eigenvalue(i,j)), i=1,nx), j=1,ny),
     .            ((1.0                   , i=1,nx), j=1,ny)
      close (11)

      
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
      real        eta(0:idim), th(0:idim), lmap
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim)
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
      character*20 filename
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,type,filename
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
c     We enforce noslip indirectly by setting v' = 0 at -1.  
c     NOTE THAT D1-D4 SHOULD ONLY BE APPLIED TO THE v FIELD ONLY.
c
      do i = 0, n
c       D1(0,i) = 0.0
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
c
c     At the upper surface we allow the flow to slip by setting v'' = 0
c     at 1.  NOTE THAT D1-D4 SHOULD ONLY BE APPLIED TO THE v FIELD ONLY.
c
      do i = 0, n
        D2(0,i) = 0.0
c       D2(N,i) = 0.0
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
      real        eta(0:idim), th(0:idim), lmap
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim)
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
      character*20 filename
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,type,filename
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
        m1(i) = 2.0/Lmap
        m2(i) = 0.0
        m3(i) = 0.0
        m4(i) = 0.0
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
      subroutine INIT_BFS_PROFILE
C***********************************************************************
C
C     Setup the initial BL profile
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=256)
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
      complex     dA2(0:idim,0:idim), dA1(0:idim,0:idim)
      complex     dA0(0:idim,0:idim), dB0(0:idim,0:idim)
      complex     v(0:idim), alpha, omega(0:idim)
      character*1 type
      character*20 filename
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,type,filename
      common      /deriv/  D1,D2,D3,D4,lmap
      common      /matrix/ A4, A3, A2, A1, A0, B2, B1, B0, dA2, dA1,
     .                     dA0, dB0
c***********************************************************************
      parameter    (idim2=1000)
      integer      i, nbl, LDD
      character*20 profile
      real         junk, us(0:idim2), d2us(0:idim2), pi, xitemp
      real         xi(0:idim2), d1us(0:idim2), d1max, d2max
      real         ydat(0:idim2), ymin, ymax, h, ut(0:idim2)
      real         uspl(0:idim2), d2Uspl(0:idim2), metric, Lmap2
      real         D1hat(0:idim,0:idim), D2hat(0:idim,0:idim)  
      character*1  temp

      pi = ACOS(-1.0)
      LDD = idim

      if (.false.) then
        write (*,10)
  10    format (/,1x,'Read Mean Profile',/)
        write (*,20)
  20    format (1x,'Enter filename ==> ',$)
        read (*,'(a)') profile
        open (unit=11,file=profile,status='unknown')
  
        i = 0
        read (11,35) temp, nbl
  35    format (a1,i5)
        nbl = nbl - 1
        do i = 0, nbl
          read (11,*,end=100) ydat(i),Ut(i),junk
        end do
        close (11)
 100    continue
      else
        nbl = 100
        do i = 0, nbl
          ydat(i) = -1.0+i*1.0/nbl
          Ut(i) = (1.0-ydat(i)**2)
        end do
      end if
c
c     Make a mirror image
c
      do i = nbl-1, 0, -1
        ydat(2*nbl-i) = 2.*ydat(nbl)-ydat(i)
        Ut(2*nbl-i) = Ut(i)
      end do
      nbl2 = 2*nbl
c
c     Interpolate data onto a Chebyshev Grid
c
      ymin2 = ydat(0)
      ymax2 = ydat(nbl2)
      h2 = (ymax2-ymin2)/float(nbl)
      Lmap2 = ymax2-ymin2

      write (*,11) nbl2
  11  format (/,1x,'Nbl = ',i5)
      call SPLINE(nbl2+1,ydat,Ut,uspl)
      ilen = index(filename,' ')-1
      open(unit=8,file=filename(1:ilen)//'.org',form='formatted',
     .     status='unknown')
      do i = 0, nbl2
        write (8,40) ydat(i),ut(i),uspl(i)
      end do 
      close (8)     
      dth = pi/(nbl2)
      ilen = index(filename,' ')-1
      open(unit=9,file=filename(1:ilen)//'.spl',form='formatted',
     .     status='unknown')
      do i = 0, nbl2
        thtemp = i*dth
        xi(i) = COS(thtemp)*Lmap2/2. + Lmap2/2. + ymin2
        call SPEVAL(nbl2+1,ydat,Ut,uspl,xi(i),us(i))
        write (9,40) xi(i), us(i), thtemp
      end do
      close (9)
      
      write (*,66) nbl2
   66 format (/,1x,'nbl =  ',i5,'   Enter ncut ==> ',$)
      read (*,*) ifilter
      write (*,77) ifilter
   77 format (/,1x,'Low pass filter at n = ',i5)

      call COSFT(us,nbl2,1)
      do i = ifilter, nbl2
        us(i) = 0.0
      end do
      call DCHEBYSHEVF(us,d1us,nbl2)
      call DCHEBYSHEVF(d1us,d2us,nbl2)
c     call COSFT(us,nbl,-1)
c     call COSFT(d1us,nbl,-1)
c     call COSFT(d2us,nbl,-1)
c     do i = 0, nbl
c       write (10,50) xi(i),us(i),d1us(i)*metric,d2us(i)*metric**2
c  50   format (1x,4(e12.4,1x))
c     end do
c
c     Convert to theta grid
c
      ymin = ydat(0)
      ymax = ydat(nbl)
      h = (ymax-ymin)/float(nbl)
      Lmap = ymax-ymin

      call MAKE_BL_METRICS

      do i = 0, n
        xitemp = COS(th(i))*Lmap/2. + Lmap/2. + ymin
        xitemp2 = (xitemp-(ymin2 + Lmap2/2.))*2./Lmap2
        if (xitemp2 .ge. 1.) then
          xitemp2 = 1.0
        else if (xitemp2 .le. -1.) then
          xitemp2 = -1.0
        end if
        thtemp = ACOS( xitemp2 )
        call CHEBYINTF(nbl2,u(i),thtemp,us)
        call CHEBYINTF(nbl2,d1u(i),thtemp,d1us)
        call CHEBYINTF(nbl2,d2u(i),thtemp,d2us)
        d1u(i) = d1u(i)*Lmap/Lmap2
        d2u(i) = d2u(i)*(Lmap/Lmap2)**2
      end do
c
c     Compute the collocation derivatives
c
c      call CHEBYD (D1hat, LDD, n)
c      do i = 0, n
c        do j = 0, n
c          D2hat(i,j) = 0.0
c          do k = 0, n
c            D2hat(i,j) = D2hat(i,j) +  D1hat(i,k)*D1hat(k,j)
c          end do
c        end do
c      end do
c      do i = 0, n
c        d1u(i) = 0.0
c        d2u(i) = 0.0
c        do k = 0, n
c            d1u(i) = d1u(i) + D1hat(i,k)*u(k)
c            d2u(i) = d2u(i) + D2hat(i,k)*u(k)
c        end do
c      end do
      
      ilen = index(filename,' ')-1
      open(unit=10,file=filename(1:ilen)//'.vel',form='formatted',
     .     status='unknown')
      do i = 0, n
        xitemp = eta(i)*Lmap/2. + Lmap/2. + ymin
        write (10,12) xitemp,u(i),d1u(i)*2./Lmap,d2u(i)*(2./Lmap)**2
      end do 
  12  format (1x,4(e16.8,4x))
      close (10)

      write (*,30) 

  30  format (/,1x,'Velocity Profile completed...',/)
  40  format (1x,5(e12.5,2x))
  90  format (1x,i5,1x,e12.5)

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
      real        eta(0:idim), th(0:idim), lmap
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim)
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
      character*20 filename
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,type,filename
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
      
      ai = cmplx(0.,1.)

c     Make the derivative (wrt alpha) matrices

      do i = 0, n
        work = -4.*alpha*m1(i)**2 - ai*Re*U(i)*m1(i)**2 
        do j = 0, n
          dA2(i,j) = work*D2(i,j)
        end do
      end do

      do i = 0, n
        work = -4.*alpha*m2(i) - ai*Re*U(i)*m2(i)
        do j = 0, n
          dA1(i,j) = work*D1(i,j)
        end do
      end do

      do i = 0, n
        work =  4.*alpha**3 + 3.*ai*alpha**2*Re*U(i) + 
     .          ai*Re*( d2U(i)*m1(i)**2 + d1U(i)*m2(i) )
        do j = 0, n
          dA0(i,j) = work*identity(i,j)
        end do
      end do

      do i = 0, n
        work = ai*2.*alpha*Re
        do j = 0, n
          dB0(i,j) = work*identity(i,j)
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
#if USE_IMSL      
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
      
      C = CMPLX(0.0,0.0)
      do I = 1, N
        C = C + CONJG(X(I))*Y(I)
      end do

      RETURN
      END

C***********************************************************************
C               S O L V E   O R R   S O M M E R F E L D 
C***********************************************************************
      subroutine SOLVE(eigenvalue,eigenvector,print,oldeig,dwda)
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
      real        eta(0:idim), th(0:idim), lmap
      real        m1(0:idim), m2(0:idim), m3(0:idim), m4(0:idim)
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
      character*20 filename
      
      common      /data/   n,eta,th,m1,m2,m3,m4,c,u,d1u,d2u,Re,v,
     .                     alpha,omega,type,filename
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

      real        temp1(0:idim), temp2(0:idim), residual, CHECKEIG
      real        temp3(0:idim), temp4(0:idim), dalpha, diff(0:idim)
      complex     oldeig

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
          dA(i-1,j-1) = dA2(i,j)+dA1(i,j)+dA0(i,j)
          dB(i-1,j-1) = dB0(i,j)
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
      CALL LINCG (N-1, B, LDA, T1, LDA)
      CALL MCRCR (N-1, N-1, T1, LDA, N-1, N-1, A, LDA,
     .            N-1, N-1, T4, LDA)
      CALL EVCCG (N-1, T4, LDA, eval, evec, ldevec)
C
C     Solve adjoint problem too
C      
      CALL MCRCR (N-1, N-1, A, LDA, N-1, N-1, T1, LDA,
     .            N-1, N-1, T3, LDA)
      call HTRAN (N-1, T3, T2, LDA)
      call EVCCG (N-1, T2, LDA, aeval, aevec, ldevec)
#else
      write(*,*) "Need to implement LAPACK/BLAS"
#endif      
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
      do i = 0, n-2
        temp1(i) = REAL(eval(index(i)))
      end do
      
      eigenvalue = cmplx(temp1(n-2),temp2(n-2))
      
      iloc = n-2

c      first = .true.
c      if ( real(oldeig) .le. -999.) then
c        eigenvalue = cmplx(temp1(n-2),temp2(n-2))
c        iloc = n-2
c      else
c        do j = 0, n-2
c          diff(j) = abs(real(temp1(j)) - real(oldeig))
c          index2(j) = j
c        end do
c        call PIKSR2(n-1,diff,index2)
c        eigenvalue = cmplx(temp1(index2(1)),temp2(index2(1)))
c        iloc = index2(1)
c      end if

      eigenvector(0) = cmplx(0.0,0.0)
      eigenvector(n) = cmplx(0.0,0.0)
      
      do j = 1, n-1
        eigenvector(j) = evec(j-1,index(iloc))
        tvec(j-1)  = evec(j-1,index(iloc))
        tavec(j-1) = aevec(j-1,index(iloc))
      end do
      
c     write (*,120) real(eigenvalue),aimag(eigenvalue),
c    .              real(oldeig), aimag(oldeig)
c120  format (1x,'Eigs = ',4(e15.8,2x))
      
      residual = CHECKEIG (N-1,T4,lda,eigenvalue,tvec)
      if (residual .gt. .01) then
        write (*,*) 'WARNING eigenvalue not converged!'
      end if
      
c     Compute dw/da

       do i = 0, n-2
         do j = 0, n-2
           T1(i,j) = dB(i,j)*eigenvalue-dA(i,j)
         end do
       end do
#ifdef USE_IMSL
       call MUCRV  (N-1, N-1, T1, LDA, N-1, tvec, 1, N-1, temp1)
       call CXDOTY (N-1, tavec, temp1, dw)
       call MUCRV  (N-1, N-1, B, LDA, N-1, tvec, 1, N-1, temp1)
       call CXDOTY (N-1, tavec, temp1, dalp)
#else
       write(*,*) "Need to implement LAPACK/BLAS"
#endif
       dwda = -dw/dalp
       dalpha = -AIMAG(eigenvalue)/AIMAG(dwda)
       nalpha = alpha + dalpha
       nw = REAL(eigenvalue) + dalpha*REAL(dwda)
c      write (*,105) real(dwda), aimag(dwda), dalpha
  105  format (/,1x,'dw/da = ',e15.8,',  ',e15.8,' i   dalpha = ',
     .          e15.8)

      if (print) then
        write (*,100) residual
  100   format (/,1x,'Residual = ',e17.10)
c       write (*,105) real(dwda), aimag(dwda), dalpha
c 105   format (/,1x,'dw/da = ',e15.8,',  ',e15.8,' i   dalpha = ',
c    .          e15.8)
c       write (*,115) nalpha, nw
c 115   format (/,1x,'nalpha = ',e15.8,',   nw = ',e15.8)
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
        do while (which .ne. -1)
          temp1(0) = 0.0
          temp2(0) = 0.0
          temp1(n) = 0.0
          temp2(n) = 0.0
          temp3(0) = 0.0
          temp4(0) = 0.0
          temp3(n) = 0.0
          temp4(n) = 0.0
          do j = 1, n-1
            temp1(j) =  REAL( evec(j-1,which))
            temp2(j) = AIMAG( evec(j-1,which))
            temp3(j) =  REAL(aevec(j-1,which))
            temp4(j) = AIMAG(aevec(j-1,which))
          end do
c         do j = 0, n
c           write (*,38) eta(j),temp1(j),temp2(j),temp3(j),temp4(j)
c         end do
c
c         Write out normal solution
c
          write (*,*)
          call WRITEOUT (alpha,eval(which)/alpha,Re,n,temp1,temp2,
     .                   256,Lmap,ymin)
c
c         Write out Adjoint solution
c
c         write (*,*)
c         call CHEBYINT (n, temp3, temp4, 128)

          write (*,40)
          read (*,*) which
        end do
      end if

 31   format(/,1x,3(e17.10,4x))
 32   format(/,1x,'       c_r                c_i              count',/)
 35   format(1x,2(e17.10,4x),i5)
 36   format(1x,2(e17.10,4x))
 37   format(1x,2(e17.10,2x),i5)
 38   format(1x,5(e15.8,1x))
 40   format(/,1x,'Eigenvector for which eigenvalue (-1 quits) ==> ',$)
 55   format(/,1x,'Performance Index = ',g12.5,/)
      
      return
      end

C***********************************************************************
      SUBROUTINE WRITEOUT(ALPHA,C,RE,NMODE,EIGR,EIGI,NBIG,LMAP,YMIN)
C***********************************************************************
C
C     INTERPOLATE EIGENFUNCTIONS TO A FINER MESH AND OUTPUT
C
c***********************************************************************
      INTEGER   NBIG, NMODE
      REAL      EIGR(0:NMODE), EIGI(0:NMODE), LMAP, XI, YMIN
      REAL      TR, TI, X(0:NBIG), Y(0:NBIG), PI
      REAL      D1(0:NBIG,0:NBIG),TIME,U(0:NBIG,0:NBIG),V(0:NBIG,0:NBIG)
      REAL      PSIR(0:NBIG,0:NBIG), PSII(0:NBIG,0:NBIG)
      COMPLEX   T(0:NBIG), DT(0:NBIG), ALPHA, C
      
      PI = ACOS(-1.0)

      call CHEBYSHEV (EIGR,NMODE,1)
      call CHEBYSHEV (EIGI,NMODE,1)

      OPEN (UNIT=20,FILE='eigfun.dat',FORM='FORMATTED',STATUS='UNKNOWN')
      DO J = 0, NBIG
        Y(J) = J*PI/NBIG
        X(J) = J*2.*PI/NBIG
        TR = 0.0
        TI = 0.0
        DO M = 0, NMODE
          TR = TR + EIGR(M)*COS(FLOAT(M)*Y(J))
          TI = TI + EIGI(M)*COS(FLOAT(M)*Y(J))
        END DO
        T(J) = CMPLX(TR,TI)
        XI = COS(Y(J))*LMAP/2.+LMAP/2.+YMIN
        WRITE (20,10) XI, TR, TI
      END DO
      CLOSE (20)
C
C     NEED TO COMPUTE U AND V
C
      call CHEBYD (D1, NBIG, NBIG)
      DO I = 0, NBIG
        DT(I) = 0.0
        DO K = 0, NBIG
            DT(I) = DT(I) + D1(I,K)*T(K)
        END DO
        DT(I) = DT(I)*2./LMAP
      END DO

      OPEN (UNIT=21,FILE='eigvel.dat',FORM='FORMATTED',STATUS='UNKNOWN')
      TIME = 0.0
      DO I = 0, NBIG
        DO J = 0, NBIG
          U(I,J) = REAL(DT(J)*CEXP((0.,1.)*ALPHA*(X(I)-C*TIME)))
          V(I,J) = REAL((0.,-1.)*ALPHA*T(J)*
     .                  CEXP((0.,1.)*ALPHA*(X(I)-C*TIME)))
          PSIR(I,J) = REAL(T(J)*CEXP((0.,1.)*ALPHA*(X(I)-C*TIME)))
          PSII(I,J) = AIMAG(T(J)*CEXP((0.,1.)*ALPHA*(X(I)-C*TIME)))
          IF (I.EQ.0) THEN
            XI = COS(Y(J))*LMAP/2.+LMAP/2.+YMIN
            WRITE (21,20) XI, U(I,J), V(I,J), REAL(DT(J)),
     .                    AIMAG(DT(J)), REAL(T(J)), AIMAG(T(J))
          END IF
        END DO
      END DO
      CLOSE (21)
C
C     Write a Plot3d file
C
      open (unit=10,file='g.dat',form='unformatted',status='unknown')
      write (10) NBIG, NBIG
      write (10) (( X(I), i=1,NBIG), j=1,NBIG),
     .           (( COS(Y(J))*LMAP/2.+LMAP/2.+YMIN,i=1,NBIG),j=1,NBIG)

      open (unit=11,file='q.dat',form='unformatted',status='unknown')
      write (11) NBIG, NBIG
      write (11) 0.0,0.0,RE,TIME
      write (11)  ((PSIR(i,j), i=1,NBIG), j=1,NBIG),
     .            ((   U(i,j), i=1,NBIG), j=1,NBIG),
     .            ((   V(i,j), i=1,NBIG), j=1,NBIG),
     .            ((PSII(i,j), i=1,NBIG), j=1,NBIG)
          

  10  format (1x,3(e16.8,4x))
  20  format (1x,7(e16.8,4x))

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
            
C*************************************************************************
      SUBROUTINE COSFT3D (YY,N,ISIGN)
C*************************************************************************
C
C     Do a brute force transform.  ISIGN = 1 does a forward transform.
C     ISIGN = -1 a backward transform.  This one is fully double 
C     precision.
C
C*************************************************************************
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
      
C*************************************************************************
      SUBROUTINE CHEBYD (D, LDD1, N)
C*************************************************************************
C
C     Calculation the Chebyshev collocation derivative matrix
C
C*************************************************************************
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
