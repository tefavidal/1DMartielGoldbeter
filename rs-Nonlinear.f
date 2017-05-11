      subroutine rs(t,Nx,gamma,ro,gammaprime,roprime,vdx)
      
      implicit none

      integer i, Nx
      double precision t
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx),ro(Nx)
      double precision gammaprime(Nx),roprime(Nx)
      double precision gLaplace(Nx)
      double precision xgradeC(Nx)
      double precision vdx(Nx)
      double precision f1, f2, Y, aux, Phi,factor



      call functionLap(Nx,gamma,gLaplace,xgradeC)

!      do i=1,100
!        roprime(i)=0.0
!        gammaprime(i)=0.0
!      enddo

      do i=1,Nx
        if (i .gt. 500 .and. i .le. 513)then
        factor=0.0
        else
        factor=1.0
        endif

        aux=gamma(i)
        f1=(1.d0+dk*aux)/(1.d0+aux)
        f2=(dL1+dk*dL2*dc*aux)/(1.d0+dc*aux)
        Y=ro(i)*aux/(1.d0+aux)
        Phi=(dlambda1+Y**2)/(dlambda2+Y**2)

        roprime(i)=factor*(-f1*ro(i)+f2*(1.d0-ro(i)))
        gammaprime(i)=factor/depsilon*(s1*s2*Phi-aux)
     .                  +depsilon*gLaplace(i)-
     .            (vdx(i)*xgradeC(i))
      enddo

      return

      end
!     ***********************************************************
!      ***********************************************************

      subroutine functionLap(Nx,gamma,gLaplace,xgradeC)

      implicit none
      integer i, Nx

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx), gLaplace(Nx), xgradeC(Nx)
      double precision gammaim2, gammaim1,gammaip1
      double precision thetai, thetaim1, psii, psiim1

      do i=1,Nx
!       No-Flux boundary condition
!       if(i .eq. 1 .or. i .eq. 135) then
       if(i .eq. 1) then
        gammaim2=gamma(i+1)
        gammaim1=gamma(i+1)
        gammaip1=gamma(i+1)
!       elseif(i .eq. 2 .or. i .eq. 136) then
       elseif(i .eq. 2 ) then
        gammaim2=-gamma(i)+2*gamma(i-1)
        gammaim1=gamma(i-1)
        gammaip1=gamma(i+1)
       elseif(i .eq. Nx) then
!       elseif(i .eq. Nx .or. i .eq. 130) then
        gammaim2=gamma(i-2)
        gammaim1=gamma(i-1)
        gammaip1=gamma(i-1)
       else
        gammaim2=gamma(i-2)
        gammaim1=gamma(i-1)
        gammaip1=gamma(i+1)
       endif

!       Periodic boundary condition
!        if(i .eq. 1) then
!           gammaim2=gamma(Nx-1)
!           gammaim1=gamma(Nx)
!           gammaip1=gamma(i+1)
!        elseif(i .eq. 2) then
!            gammaim2=gamma(Nx)
!            gammaim1=gamma(1)
!            gammaip1=gamma(i+1)
!        elseif(i .eq. Nx) then
!           gammaim2=gamma(i-2)
!           gammaim1=gamma(i-1)
!           gammaip1=gamma(1)
!       else
!           gammaim2=gamma(i-2)
!           gammaim1=gamma(i-1)
!           gammaip1=gamma(i+1)
!       endif

        gLaplace(i)=(gammaip1+gammaim1-2*gamma(i))/(dx**2)




        if(gammaip1 .eq. gamma(i)) then
        thetai=1.d-10
        else
        thetai=(gamma(i)-gammaim1)/(gammaip1-gamma(i))+1.d-10
        endif

        if(gamma(i) .eq. gammaim1)then
        thetaim1=1.d-10
        else
        thetaim1=(gammaim1-gammaim2)/(gamma(i)-gammaim1)+1.d-10
        endif

        psii=max(0.0,min(1.0,1.0/3.0+thetai/6.0,thetai))
        psiim1=max(0.0,min(1.0,1.0/3.0+thetaim1/6.0,thetaim1))

       xgradeC(i)=(1.0-psiim1+psii/thetai)*(-gammaim1+gamma(i))/(dx)

      enddo

      return
      end

