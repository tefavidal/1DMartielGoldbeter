      subroutine ic(t,Nx,gamma,ro)
      
      implicit none
      double precision t
      integer i,Nx

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0


      double precision gamma(Nx),ro(Nx)
      double precision gamma0(10),ro0(10)
      double precision dke(Nx),dsigma(Nx)
      double precision gammanullc(274),ro1nullc(274),ro2nullc(274)
      double precision dkc0(10)
      double precision Y0, A, B, dM, dN, f1, f2, SS, Delta, a11, a22,s
      double precision dflux, factor
      integer nfix, ii

!      dflux=-1.0*vd
       dflux=0.0
      ii=101

      call SteadyState(ii,dflux,nfix,gamma0,ro0)
      write(6,*) 'nfix==',nfix

      do i=1,nfix
        write(6,*)  'gamma=',gamma0(i),'  ro=',ro0(i)
      enddo

!      do i=1,1450
      do i=1,Nx
             ro(i)=ro0(1)
            gamma(i)=gamma0(1)
      enddo

!      do i=1500, Nx
      do i=100,200
             ro(i)=ro0(1)
            gamma(i)=gamma0(1)+2
      enddo




      do i = 1,1
         gamma01=gamma0(i)

         ro01=ro0(i)

      s=s1*s2
      Y0=ro01*gamma01/(1.d0+gamma01)
      A=(dk*dL2-dL1)*dc/(1+dc*gamma01)**2
      B=(dk-1.d0)/(1+gamma01)**2

      write(6,*) '-----------------------------------------------------'
      write(6,*) 'A=',A,'B=',B

      dM=2*ro01*gamma01**2*(dlambda2-dlambda1)
     . /((1+gamma01)**2*(dlambda2+Y0**2)**2)
      dN=2*ro01**2*gamma01*(dlambda2-dlambda1)
     . /((1+gamma01)**3*(dlambda2+Y0**2)**2)
      write(6,*) '-----------------------------------------------------'
      write(6,*) 'M=',dM,'N=',dN,'s=',s
      f1=(1.d0+dk*gamma01)/(1.d0+gamma01)
      f2=(dL1+dk*dL2*dc*gamma01)/(1.d0+dc*gamma01)
      SS=-(f1+f2)+(dN*s-1)/depsilon
      Delta=((f1+f2)*(1.d0-dN*s)+s*dM*B*ro01
     . -s*dM*A*(1-ro01))/(depsilon)
      a11=(dN*s-1)/depsilon
      a22=-f1-f2
      write(6,*) '-----------------------------------------------------'
      write(6,*) 'a11=(N*s-1)/epsilon=',a11,'a22=-f1-f2=',a22
      write(6,*) '-----------------------------------------------------'
      write(6,*)'Stabilty of two-variable subsystem:S<0 and Delta>0'
      write(6,*) 'Delta (determinant)=',Delta
      write(6,*) 'S (trace)=',SS

      write(6,*) '-----------------------------------------------------'

      enddo



      return

      end

!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine SteadyState(ii,dflux,n,gamma0,ro0)

      implicit none
      integer ii, n, i
      double precision dflux, fixpoint, v, vnew, u

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma0(10),ro0(10)
      double precision gamma, ro
!     Apliying the constrain dtro=0 looks for when dtgamma changes
!     sign, when it does, calls zero gamma
      n = 0

      fixpoint = 0.d0

      do i = 1,90
      gamma=0.001+0.0001d0*i
      call function(ii,dflux,gamma,ro,vnew)
!     returns ro such that dtro=0 for that gamma
!     vnew=dtgamma

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then

            call zerogamma(ii,dflux,gamma-0.0001,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
            endif
      endif
      v=vnew
      enddo

      do i = 1,91
      gamma=0.01+0.001d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-0.001,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo

      do i = 1,91
      gamma=0.1+0.01d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-0.01,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo

      do i = 1,91
      gamma=1+0.1d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-0.1,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo

       do i = 1,91
      gamma=10+1.d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-1.d0,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo

       do i = 1,91
      gamma=100+10.d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-10.d0,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo

       do i = 1,91
      gamma=1000+100.d0*(i-1)
      call function(ii,dflux,gamma,ro,vnew)

      if ( i .gt. 1) then
            if (v*vnew .lt. 0.d0) then
            call zerogamma(ii,dflux,gamma-100.d0,gamma,fixpoint,ro,u)
            n=n+1
            gamma0(n)=fixpoint

            ro0(n)=ro
      endif
      endif
      v=vnew
      enddo



      return
      end
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine function(ii,dflux,gamma,ro,v)

      implicit none
      integer ii
      double precision dflux, v, x, s22, f1,f2, Y, Phi, gamma, ro
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      x=(1.d0*ii-50.d0)/10
      s22=s2*(1+dtanh(x))/2

      f1=(1.d0+dk*gamma)/(1.d0+gamma)
      f2=(dL1+dk*dL2*dc*gamma)/(1.d0+dc*gamma)
      ro=f2/(f1+f2)

      Y=ro*gamma/(1.d0+gamma)
      Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

      v=(s1*s22*Phi-(1+dtanh(x))/2*gamma)/depsilon+dflux


      return
      end

!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine zerogamma(ii,dflux,gamma1,gamma2,gamma,ro,v)

      implicit none
      integer ii
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision dflux, gamma, gamma1, gamma2, ro, ro1, ro2
      double precision v, v1, v2

      call function(ii,dflux,gamma1,ro1,v1)
      call function(ii,dflux,gamma2,ro2,v2)

 10   gamma=(gamma1+gamma2)/2
      call function(ii,dflux,gamma,ro,v)
      if(v1*v .lt. 0.d0) then
      gamma2=gamma

      ro2=ro
      v2=v
      else
      gamma1=gamma

      ro1=ro
      v1=v
      endif
      if (dabs(gamma2-gamma1) .lt. 1.d-12) then
        go to 25
      endif
      if(dabs(v) .gt. 1.d-6) then
         goto 10
      endif

 25   return
      end









