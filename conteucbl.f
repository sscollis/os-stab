c***********************************************************************
      program Conte_uc_bl
c***********************************************************************
c 
c     Purpose:  Solve the Orr-Sommerfeld equation using 4th order Runge-
c               Kutta explicit integration with orthogonalization using
c               the method of Conte.
c
c               This routine tracks the saddle point in the complex
c               wavenumber plane by shifting the reference frame
c               velocity, u_c.
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
c     Revised:  9-18-92         Tabs are 8 characters
c
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=20000)
      integer     nbl
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta, uc
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p
      real        Uorg(0:idim)

      common      /data/  nbl, u, d2u, ymin, ymax, h, f2p, uc,
     .                    gamma, beta, uspl, d2uspl, ydat, Uorg
     
      common      /eig/   c, alpha, Re
c***********************************************************************
      parameter     (neq=4, inunit=20)
      complex       bc1(neq), bc2(neq)
      real          testalpha 
      integer       north, nstep, num, omega, saddle, readbl
      logical       eigfun
      character*1   input
      character*20  infile
      character*80  string
c***********************************************************************
      complex       alphas,resid,eigenvalue(3),alphasold,A,B,alp(3)
      real          delta,alphar(3),alphai(3),deltauc
c***********************************************************************
c
c     Defaults
c
      nbl = 500
      nstep = 500
      north = 50
      testalpha = 0.0
      Re = 580.
      beta = 0.0
      alphar(1) = 0.179
      alphai(1) = 0.0
      cr =  0.36413124E+00
      ci =  0.79527387E-02
      ymin = 0.0
      ymax = 20.0
      f2p = 0.5
      delta = 0.0
      deltauc = 0.0
      num = 1
      Uc = 0.0
      omega = 0
      saddle = 0
      readbl = 1
c
c     Input from keyboard
c      
      write (*,10)
  10  format (/,10x,'Solve Orr-Sommerfeld (Shooting)',/)
      write (*,15)
  15  format(/,1x,'Read from keyboard, file, default (k,f,d) ==> ',$)
      read (*,'(a)') input
      if (input .eq. 'k' .or. input .eq. 'K') then
        write (*,20)
  20    format (/,1x,'Enter number of profile points ==> ',$)
        read (*,*) nbl
        write (*,25)
  25    format (/,1x,'Enter number of integration steps ==> ',$)
        read (*,*) nstep
        write (*,27)
  27    format (/,1x,'Enter Test Alpha ==> ',$)
        read (*,*) testalpha
        write (*,30)
  30    format (/,1x,'Enter Reynolds number ==> ',$)
        read (*,*) Re
        write (*,35)
  35    format (/,1x,'Enter beta ==> ',$)
        read (*,*) beta
        write (*,37)
  37    format (/,1x,'Enter guess for BL second deriv. at wall ==> ',$)
        read (*,*) f2p
        write (*,40) 
  40    format (/,1x,'Enter alpha_r ==> ',$)
        read (*,*) alphar(1)
        write (*,45) 
  45    format (/,1x,'Enter alpha_i ==> ',$)
        read (*,*) alphai(1)
        write (*,50)
  50    format (/,1x,'Enter guess for (c_r, c_i) ==> ',$)
        read (*,*) cr, ci
        write (*,55)
  55    format (/,1x,'Enter Ymin and Ymax ==> ',$)
        read (*,*) Ymin, Ymax
        write (*,65)
  65    format (/,1x,'Enter reference frame velocity u_c ==> ',$)
        read (*,*) uc
        write (*,70)
  70    format (/,1x,'Enter delta ==> ',$)
        read (*,*) delta
        write (*,71)
  71    format (/,1x,'Enter deltauc ==> ',$)
        read (*,*) deltauc
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
        call readr (string,alphar(1))
        read(inunit,'(a)') string
        call readr (string,alphai(1))
        read(inunit,'(a)') string
        call readr (string,cr)
        read(inunit,'(a)') string
        call readr (string,ci)
        read(inunit,'(a)') string
        call readr (string,Ymin)
        read(inunit,'(a)') string
        call readr (string,Ymax)
        read(inunit,'(a)') string
        call readr (string,Uc)
        read(inunit,'(a)') string
        call readr (string,delta)
        read(inunit,'(a)') string
        call readr (string,deltauc)
        read(inunit,'(a)') string
        call readi (string,num)
        read(inunit,'(a)') string
        call readi (string,omega)
        read(inunit,'(a)') string
        call readi (string,saddle)
      end if
c
c     Check fixed dimensions
c
      if (nbl .gt. idim) then
        write (*,300) 
 300    format (/,/,1x,'N > Idim...Stopping',/,/)
        goto 210
      end if
c
c     Set constants
c
      alpha = CMPLX (alphar(1), alphai(1))
      c = CMPLX (cr,ci)
      h = (ymax-ymin)/float(nbl)
      if (omega .eq. 1) c = c/alpha
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
      write (*,120) alphar(1), alphai(1)
 120  format (1x,'alpha = (',e20.10,', ',e20.10,')')
      write (*,130) REAL(c), AIMAG(c)
 130  format (1x,'c = (',e20.10,', ',e20.10,')')
      write (*,145) Uc
 145  format (1x,'Uc = ',e20.10)
      write (*,150) delta
 150  format (1x,'delta = ',e20.10)
      write (*,155) deltauc
 155  format (1x,'deltauc = ',e20.10)
      write (*,160) num
 160  format (1x,'num = ',i5,/)
c
c     set the boundary conditions
c
      do i = 1, neq
        bc1(i) = 0.0
        bc2(i) = 0.0
      end do
c
c     Fix to make my nondimensionalization match Mack's
c
      Re = Re*SQRT(2.)
      alphar(1) = alphar(1)
      alphai(1) = alphai(1)
      alphar(2) = alphar(1) + delta
      alphar(3) = alphar(1) - delta
      alphai(2) = alphai(1) - delta
      alphai(3) = alphai(1) - delta
      ymin = ymin
      ymax = ymax
c      
c     Loop for different convection speeds
c
      if (readbl.eq.1) then
        call READ_BL
      else
        call SOLVE_BL
      end if
      
      if (saddle .eq. 1) then
        eigfun = .false.
        alphas = cmplx(10.,10.)
        do k = 1, num
c         write (*,*) uc
          call SHIFT_BL
          resid = 10.
          icount = 0
          do while ( abs(resid) .gt. 1.E-5 )
            icount = icount + 1
            do i = 1, 3
              alphasold = alphas
              alpha = cmplx (alphar(i), alphai(i))*SQRT(2.)
              alp(i) = alpha/SQRT(2.)
              call CONTE(nstep,testalpha,neq,2,bc1,bc2,ymin,ymax,eigfun)
              eigenvalue(i) = c*alp(i)
c             write (*,61)  real(alp(i)),aimag(alp(i)),real(eigenvalue(i))
c    .                      ,aimag(eigenvalue(i))
  61          format (1x,4(e17.10,4x))
            end do
            A  = cmplx(0.,-1.)*(eigenvalue(2)-eigenvalue(1))
            B  = cmplx(0.,-1.)*(eigenvalue(3)-eigenvalue(1))
            alphas = (A*(alp(3)**2-alp(1)**2)-B*(alp(2)**2-alp(1)**2))/
     .               (2.*B*(alp(1)-alp(2))-2.*A*(alp(1)-alp(3)))
            resid = alphas-alphasold
            if (k .eq. 1 .or. icount .gt. 4) then
              write (*,80) real(alphas), aimag(alphas), abs(resid)
  80          format (1x,'==> ',3(e17.10,4x))
            end if
            alphar(1) = real(alphas)
            alphai(1) = aimag(alphas)
            alphar(2) = alphar(1) + delta
            alphar(3) = alphar(1) - delta
            alphai(2) = alphai(1) - delta
            alphai(3) = alphai(1) - delta
          end do
          write (*,90)  uc, real(eigenvalue(1)),aimag(eigenvalue(1)),
     .                  real(alphas),aimag(alphas)
 90       format (1x,f6.4,2x,4(e17.10,4x))
          uc = uc - deltauc
        end do
      else
        eigfun = .true.
        alpha = cmplx (alphar(1), alphai(1))*SQRT(2.)
        call CONTE(nstep,testalpha,neq,2,bc1,bc2,ymin,ymax,eigfun)
        write (*,61) real(alpha),aimag(alpha),real(c),
     .               aimag(c)
      end if
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
      cm1 = 1000.

c      do while ((abs(err) .ge. 1.0e-8) .and. (icount .le. 20) .and.
c     .          (abs(c-cm1) .ge. 1.0e-12) )

      do while ((abs(err) .ge. 1.0e-8) .and. (icount .le. 20) .and.
     .          (abs(c-cm1) .ge. 1.0e-12) )

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
             do i = 1, n
              Utemp(i) = U(i,m,k-1)
            end do
            call ODEINT(Utemp,n,t+h,t,1.E-6,-h,1.e-20,nok,nbad,
     .                  FHOMO,RKQC)
c           write (*,*) k, nok, nbad
            do i = 1, n
              U(i,m,k) = Utemp(i)
            end do
c           call RK4SCOTT(n, U(1,m,k-1), U(1,m,k), t+h, -h, FHOMO)
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
                test = cc/SQRT(aa*bb)
c               write (*,*) k, mi, mj, test
                if (test .ge. testalpha) norm = .true.
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
     .                         P(s,i,q)
                  end do
                end if
              end do
            end do
c
c           Check the P matrix
c
c            if (.false.) then
c             do i = 1, n
c               do m = 1, n-r
c               ut(i,m) = 0.0
c                 do j = 1, n-r
c                   ut(i,m) = ut(i,m) + U(i,j,k)*P(m,j,q)
c                 end do
c               end do
c             end do
c  
c             do i = 1,n
c               write (*,*) i,(ut(i,m), m = 1, n-r)
c             end do
c             write (*,*)
c             do i = 1,n
c               write (*,*) i,( z(i,m), m = 1, n-r)
c             end do
c             write (*,*)
c             write (*,*)
c            end if
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
c        if (.true.) then
c
c         strictly enforce the zero BC
c
          B(1,qend) = -U(1,2,nstep)/U(1,1,nstep)
          B(2,qend) = 1.0
c         olderr = err
c         err = U(2,1,nstep)*B(1,qend) + U(2,2,nstep)*B(2,qend)
c        else
c
c         strictly enforce the zero slope BC
c
c         B(1,qend) =  1.0
c         B(2,qend) = -U(2,1,nstep)/U(2,2,nstep)
c         olderr = err
c         err = U(1,1,nstep)*B(1,qend) + U(1,2,nstep)*B(2,qend)
c        end if
        
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
     .         ABS(Bt-SQRT(Bt**2-4.*At*Ct)) )  then
             c = ctemp-(ctemp-cm1)*2.*Ct/(Bt+SQRT(Bt**2-4.*At*Ct))
          else
             c = ctemp-(ctemp-cm1)*2.*Ct/(Bt-SQRT(Bt**2-4.*At*Ct))
          end if
          cm2 = cm1
          cm1 = ctemp
          errm2 = errm1
          errm1 = err
        end if
c        write (*,30) icount,real(ctemp),aimag(ctemp),real(err),
c     .               aimag(err)
  30    format (1x,i4,2(e17.8,e17.8,3x))

        icount = icount + 1
      end do
      
      if (eigfun) then
c
c       Second Pass
c
        write (*,40) real(c),aimag(c),qend
  40    format (/,'Eigenvalue = ',e17.8,1x,e17.8,2x,i5/)

        q = qend
        max = 0.0
        do m = 1, n-r
          B(m,q-1) = 0.0
          do j = 1, n-r
            B(m,q-1) = B(m,q-1) + P(j,m,q)*B(j,q) 
          end do
        end do
        do i = 1, n
          y(i,nstep) = 0.0
          do m = 1, n-r
            y(i,nstep) = y(i,nstep) + U(i,m,nstep)*B(m,q)
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
        end do

        arg = ATAN2( AIMAG(y(1,0)), REAL(y(1,0)) )

        do k = 0, nstep
          do i = 1, n
            y(i,k) = y(i,k)*CEXP( (0.,-1.)*arg )
          end do
          if ( ABS( REAL( y(1,k) ) ) .gt. ABS(REAL(max)) ) then
            max = y(1,k)
          end if
        end do

        do k = 0, nstep
          t = tf - h*k
          write (11,20) t, REAL(y(1,k)/max),AIMAG(y(1,k)/max)
          write (12,20) t, REAL(y(2,k)/max),AIMAG(y(2,k)/max)
          write (13,20) t, REAL(y(3,k)/max),AIMAG(y(3,k)/max)
          write (14,20) t, REAL(y(4,k)/max),AIMAG(y(4,k)/max)
  20      format ( 1x, 3(e17.8,2x) )
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
 
      INPROD = 0.0
      do i = 1, n
c       INPROD = INPROD + v1(i)*conjg(v2(i))
        INPROD = INPROD + v1(i)*v2(i)
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
      integer     nbl
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta, uc
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p
      real        Uorg(0:idim)

      common      /data/  nbl, u, d2u, ymin, ymax, h, f2p, uc,
     .                    gamma, beta, uspl, d2uspl, ydat, Uorg
     
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
      call SPEVAL(nbl+1,ydat,U,Uspl,t,UU)
      call SPEVAL(nbl+1,ydat,d2U,d2Uspl,t,d2UU)

      do j = 1 , neq-1
        yf(j) = yo(j+1)
      end do
      yf(neq) = (1./alpha/Re*(2.*alpha**2*yo(3)-alpha**4*yo(1)) + 
     .          (0.,1.)*((UU-c)*(yo(3)-alpha**2*yo(1))-
     .          d2UU*yo(1)))*alpha*Re 

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
      integer     nbl
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta, uc
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p
      real        Uorg(0:idim)

      common      /data/  nbl, u, d2u, ymin, ymax, h, f2p, uc,
     .                    gamma, beta, uspl, d2uspl, ydat, Uorg
     
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
      call SPEVAL(nbl+1,ydat,U,Uspl,t,UU)
      call SPEVAL(nbl+1,ydat,d2U,d2Uspl,t,d2UU)

      do j = 1 , neq-1
        yf(j) = yo(j+1)
      end do
      yf(neq) = 2.*alpha**2*yo(3)-alpha**4*yo(1) + 
     .          (0.,1.)*alpha*Re*((UU-c)*(yo(3)-alpha**2*yo(1))-
     .          d2UU*yo(1)) 

      return
      end

C***********************************************************************
      subroutine READ_BL
C***********************************************************************
C
C     Read an arbitrary profile from .
C
c***********************************************************************
c     Common variables
c***********************************************************************
      parameter   (idim=20000)
      integer     nbl
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta, uc
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p
      real        Uorg(0:idim)

      common      /data/  nbl, u, d2u, ymin, ymax, h, f2p, uc,
     .                    gamma, beta, uspl, d2uspl, ydat, Uorg
     
      common      /eig/   c, alpha, Re
c***********************************************************************
      integer      i
      character*20 filename
      
      write (*,10)
  10  format (/,1x,'Read Mean Profile',/)
      write (*,20)
  20  format (1x,'Enter filename ==> ',$)
      read (*,'(a)') filename
      open (unit=11,file=filename,status='unknown')

      i = -1
      do while (.true.)
        read (11,*,end=100) ydat(i),uorg(i),junk
        i = i + 1
      end do
      close (11)
 100  nbl = i
       
      call SHIFT_BL

      write (*,30) 
  30  format (1x,'Velocity Profile completed...',/)

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
      integer     nbl
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta, uc
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p
      real        Uorg(0:idim)

      common      /data/  nbl, u, d2u, ymin, ymax, h, f2p, uc,
     .                    gamma, beta, uspl, d2uspl, ydat, Uorg
     
      common      /eig/   c, alpha, Re
c***********************************************************************
      integer     i, j, k, p
      real        xi(3,0:idim), f(3), eta(3), y
      real        k1(3), k2(3), k3(3), k4(3), err, x2old, f1old
      external    BLASIUS, RKQCR
      
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
          do j = 1, 3
            eta(j) = xi(j,i-1)
          end do
          call ODEINTR(eta,3,y-h,y,1.E-7,h/2.,1.e-20,nok,nbad,BLASIUS,RKQCR)
          do j = 1, 3
            xi(j,i) = eta(j)
          end do
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
     .              (1.0 - xi(2,nbl))
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
        Uorg(i) = xi(2,i)
        ydat(i) = ymin + i*h
      end do
c
c     Compute 2nd order finite difference approximation to 
c     second derivative
c
c      do i = 1, nbl-1
c        d2u(i) = (u(i+1)-2*u(i)+u(i-1))/(h)**2
c      end do
c      d2u(0) = (-u(4)+4.*u(3)-5.*u(2)+2.*u(1))/(h)**2
c      d2u(n) = (2.*u(nbl)-5.*u(nbl-1)+4.*u(nbl-2)-u(nbl-3))/(h)**2
      
      call SHIFT_BL

      write (*,20) 
  20  format (1x,'Velocity Profile completed...',/)

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
      parameter   (idim=20000)
      integer     nbl
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta, uc
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p
      real        Uorg(0:idim)

      common      /data/  nbl, u, d2u, ymin, ymax, h, f2p, uc,
     .                    gamma, beta, uspl, d2uspl, ydat, Uorg
     
      common      /eig/   c, alpha, Re
c***********************************************************************
c
c     Need to interpolate the velocity profile to evaluate it at
c     arbitrary y
c
      do i = 0, nbl
        U(i) = Uorg(i) - Uc
      end do
      call SPLINE(nbl+1,ydat,u,uspl)
c
c     Use the spline result for the second derivative
c
      do i = 0, nbl
        d2u(i) = uspl(i)
      end do
      call SPLINE(nbl+1,ydat,d2u,d2uspl)

      do i = 0, nbl
        y = ymin + i*h
        call SPEVAL(nbl+1,ydat,u,uspl,y,us)
        call SPEVAL(nbl+1,ydat,d2u,d2uspl,y,d2us)
        write (10,10) y, us, d2us, u(i), d2u(i)
  10    format (1x,5(E16.5E4,1x))
      end do

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
      integer     nbl
      complex     alpha, c, gamma
      real        U(0:idim), d2U(0:idim), Re, ymin, ymax, h, beta, uc
      real        Uspl(0:idim), d2Uspl(0:idim), ydat(0:idim), f2p
      real        Uorg(0:idim)

      common      /data/  nbl, u, d2u, ymin, ymax, h, f2p, uc,
     .                    gamma, beta, uspl, d2uspl, ydat, Uorg
     
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
      DIMENSION    X(N),Y(N),FDP(N)
      INTEGER      IOLD
      COMMON /OLD/ IOLD     
c
c-----THE FIRST JOB IS TO FIND THE PROPER INTERVAL.          
c
c     This is really a stupid way of searching
c
c      NM1 = N - 1    
c      DO 1 I=1,NM1   
c      IF (XX.LE.X(I+1)) GO TO 10         
c    1 CONTINUE       
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
      SUBROUTINE HUNT(XX,N,X,JLO)
C***********************************************************************
C
C     TAKEN FROM NUMERICAL RECIPES (P. 91)
C
C***********************************************************************
      REAL XX(N)
      LOGICAL ASCND
      ASCND=XX(N).GT.XX(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
        JLO=0
        JHI=N+1
        GO TO 3
      END IF
      INC = 1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
 1      JHI=JLO+INC
        IF(JHI.GT.N)THEN
          JHI=N+1
        ELSE IF (X.GE.XX(JHI).EQV.ASCND) THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
        END IF
      ELSE
        JHI=JLO
 2      JLO=JHI-INC
        IF(JLO.LT.1)THEN
          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
        END IF
      END IF
      
 3    IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      END IF
      GO TO 3
      END 
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,
     .                  DERIVS,RKQC)
C***********************************************************************
      PARAMETER (MAXSTP=10000, NMAX=10, TWO=2.0, ZERO=0.0, TINY=1.E-30)
      
      COMPLEX YP(10,200)
      COMMON /PATH/ KMAX, KCOUNT, DXSAV, XP(200), YP
      
      COMPLEX YSTART(NVAR), YSCAL(NMAX), Y(NMAX), DYDX(NMAX)
      EXTERNAL DERIVS, RKQC

      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT=0
      DO I = 1, NVAR
        Y(I)=YSTART(I)
      END DO
      IF (KMAX.GT.0) XSAV=X-DXSAV*TWO
      DO NSTP=1,MAXSTP
        CALL DERIVS(NVAR, Y, X, DYDX)
        DO I = 1, NVAR
          YSCAL(I)=ABS(Y(I))+ABS(H*DYDX(I))+TINY
        END DO
        IF(KMAX.GT.0) THEN
          IF(ABS(X-XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX-1) THEN
              KOUNT=KOUNT+1
              XP(KOUNT)=X
              DO I = 1, NVAR
                YP(I,KOUNT)=Y(I)
              END DO
              XSAV = X
            END IF
          END IF
        END IF
        IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        END IF
        IF((X-X2)*(X2-X1).GE.ZERO)THEN
          DO I = 1, NVAR
            YSTART(I)=Y(I)
          END DO
          IF(KMAX.NE.0)THEN
            KOUNT=KOUNT+1
            XP(KOUNT)=X
            DO I = 1, NVAR
              YP(I,KOUNT)=Y(I)
            END DO
          END IF
          RETURN
        END IF
        IF(ABS(HNEXT).LT.HMIN) PAUSE 'Stepsize smaller than minimum'
        H=HNEXT
      END DO
      PAUSE 'Too many steps'
      RETURN
      END
C***********************************************************************
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
C***********************************************************************
      PARAMETER (NMAX=10,PGROW=-0.20,PSHRNK=-0.25,FCOR=1./15.,
     .           ONE=1.,SAFETY=0.9,ERRCON=6.E-4)
      EXTERNAL DERIVS
      COMPLEX  Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      
      XSAV=X
      DO I = 1,N
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
      END DO
      H=HTRY
  1   HH=0.5*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(N, YTEMP, X, DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
      IF(X.EQ.XSAV) PAUSE 'Stepsize not significant in RKQC.'
      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.0
      DO I = 1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX = MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
      END DO
      ERRMAX = ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON) THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.*H
        END IF
      END IF
      DO I = 1, N
        Y(I) = Y(I)+YTEMP(I)*FCOR
      END DO
      RETURN
      END
C***********************************************************************
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
C***********************************************************************
      PARAMETER (NMAX=10)
      COMPLEX Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      external derivs
      
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO I = 1, N
        YT(I)=Y(I)+HH*DYDX(I)
      END DO
      CALL DERIVS(N, YT, XH, DYT)
      DO I = 1, N
        YT(I)=Y(I)+HH*DYT(I)
      END DO
      CALL DERIVS(N, YT, XH, DYM)
      DO I = 1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
      END DO
      CALL DERIVS(N, YT, X+H, DYT)
      DO I = 1, N
        YOUT(I) = Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
      END DO
      RETURN
      END
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
      SUBROUTINE ODEINTR(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,
     .                   DERIVS,RKQCR)
C***********************************************************************
      PARAMETER (MAXSTP=10000, NMAX=10, TWO=2.0, ZERO=0.0, TINY=1.E-30)
      
      real    YP(10,200)
      COMMON /PATHR/ KMAX, KCOUNT, DXSAV, XP(200), YP
      
      real    YSTART(NVAR), YSCAL(NMAX), Y(NMAX), DYDX(NMAX)
      EXTERNAL DERIVS, RKQCR

      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT=0
      DO I = 1, NVAR
        Y(I)=YSTART(I)
      END DO
      IF (KMAX.GT.0) XSAV=X-DXSAV*TWO
      DO NSTP=1,MAXSTP
        CALL DERIVS(NVAR, Y, X, DYDX)
        DO I = 1, NVAR
          YSCAL(I)=ABS(Y(I))+ABS(H*DYDX(I))+TINY
        END DO
        IF(KMAX.GT.0) THEN
          IF(ABS(X-XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX-1) THEN
              KOUNT=KOUNT+1
              XP(KOUNT)=X
              DO I = 1, NVAR
                YP(I,KOUNT)=Y(I)
              END DO
              XSAV = X
            END IF
          END IF
        END IF
        IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
        CALL RKQCR(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        END IF
        IF((X-X2)*(X2-X1).GE.ZERO)THEN
          DO I = 1, NVAR
            YSTART(I)=Y(I)
          END DO
          IF(KMAX.NE.0)THEN
            KOUNT=KOUNT+1
            XP(KOUNT)=X
            DO I = 1, NVAR
              YP(I,KOUNT)=Y(I)
            END DO
          END IF
          RETURN
        END IF
        IF(ABS(HNEXT).LT.HMIN) PAUSE 'Stepsize smaller than minimum'
        H=HNEXT
      END DO
      PAUSE 'Too many steps'
      RETURN
      END
C***********************************************************************
      SUBROUTINE RKQCR(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
C***********************************************************************
      PARAMETER (NMAX=10,PGROW=-0.20,PSHRNK=-0.25,FCOR=1./15.,
     .           ONE=1.,SAFETY=0.9,ERRCON=6.E-4)
      EXTERNAL DERIVS
      real     Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      
      XSAV=X
      DO I = 1,N
        YSAV(I)=Y(I)
        DYSAV(I)=DYDX(I)
      END DO
      H=HTRY
  1   HH=0.5*H
      CALL RK4R(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(N, YTEMP, X, DYDX)
      CALL RK4R(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
      IF(X.EQ.XSAV) PAUSE 'Stepsize not significant in RKQC.'
      CALL RK4R(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.0
      DO I = 1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX = MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
      END DO
      ERRMAX = ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON) THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.*H
        END IF
      END IF
      DO I = 1, N
        Y(I) = Y(I)+YTEMP(I)*FCOR
      END DO
      RETURN
      END
C***********************************************************************
      SUBROUTINE RK4R(Y,DYDX,N,X,H,YOUT,DERIVS)
C***********************************************************************
      PARAMETER (NMAX=10)
      real       Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      external   derivs
      
      HH=H*0.5
      H6=H/6.
      XH=X+HH
      DO I = 1, N
        YT(I)=Y(I)+HH*DYDX(I)
      END DO
      CALL DERIVS(N, YT, XH, DYT)
      DO I = 1, N
        YT(I)=Y(I)+HH*DYT(I)
      END DO
      CALL DERIVS(N, YT, XH, DYM)
      DO I = 1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
      END DO
      CALL DERIVS(N, YT, X+H, DYT)
      DO I = 1, N
        YOUT(I) = Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
      END DO
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
C***********************************************************************
      subroutine RK4SCOTT(neq, yo, yf, to, h, FUNC)
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


      