      subroutine ODE(t,Nx,gamma,ro,vdx)
      
      implicit none
      double precision t
      integer Nx, i
      
      double precision gamma(Nx),ro(Nx)
      double precision gamma0(Nx),ro0(Nx)
      double precision gammak1(Nx),rok1(Nx)
      double precision gammak2(Nx),rok2(Nx)
      double precision gammak3(Nx),rok3(Nx)
      double precision gammak4(Nx),rok4(Nx)
      double precision gammak5(Nx),rok5(Nx)
      double precision vdx(Nx)
      double precision g1(Nx),r1(Nx)
      
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision tau,h,err
      integer iteration, index

      tau=0.d0
      h=dt

 13   do i=1,Nx
        gamma0(i)=gamma(i)
        ro0(i)=ro(i)
      enddo
      iteration=0

      call rs(t,Nx,gamma0,ro0,gammak1,rok1,vdx)
!     Runge-Kutta-Merson Method

 16   do i=1,Nx
        gamma(i)=gamma0(i)+h*gammak1(i)/3
        ro(i)=ro0(i)+h*rok1(i)/3
      enddo

      call rs(t+h/3,Nx,gamma,ro,gammak2,rok2,vdx)

      do i=1,Nx
        gamma(i)=gamma0(i)+h*(gammak1(i)+gammak2(i))/6
        ro(i)=ro0(i)+h*(rok1(i)+rok2(i))/6
      enddo

      call rs(t+h/3,Nx,gamma,ro,gammak3,rok3,vdx)

      do i=1,Nx
        gamma(i)=gamma0(i)+h*(gammak1(i)+3*gammak3(i))/8
        ro(i)=ro0(i)+h*(rok1(i)+3*rok3(i))/8
      enddo

      call rs(t+h/2,Nx,gamma,ro,gammak4,rok4,vdx)
       do i=1,Nx
        gamma(i)=gamma0(i)+h*(gammak1(i)-3*gammak3(i)
     .   +4*gammak4(i))/2
        ro(i)=ro0(i)+h*(rok1(i)-3*rok3(i)
     .   +4*rok4(i))/2
      enddo

      call rs(t+h,Nx,gamma,ro,gammak5,rok5,vdx)

      do i=1,Nx
        gamma(i)=gamma0(i)+h*(gammak1(i)+4*gammak4(i)
     .   +gammak5(i))/6
        ro(i)=ro0(i)+h*(rok1(i)+4*rok4(i)
     .   +rok5(i))/6
      enddo

      do i=1,Nx
!        g1(i)=gamma(i)-h*(gammak1(i)+gammak2(i)+gammak3(i)
!     .   +gammak4(i)+gammak5(i))/5-gamma0(i)
!        r1(i)=ro(i)-h*(rok1(i)+rok2(i)+rok2(i)+rok3(i)
!     .   +rok4(i)+rok5(i))/5-ro0(i)
         g1(i)=h*(2*gammak1(i)-9*gammak3(i)
     .   +8*gammak4(i)-gammak5(i))/30
         r1(i)=h*(2*rok1(i)-9*rok3(i)
     .   +8*rok4(i)-rok5(i))/30
      enddo

      err=0.d0
      index=0
      do i=1,Nx
        err=max(abs(g1(i)),abs(r1(i)))
        if (gamma(i) .lt. 0.d0 .or. ro(i) .lt. 0.d0)then
            index=1
            exit
        endif
      enddo

      if (err .gt. tol .or. index .eq. 1) then
        h=h/2
        iteration=iteration+1
        if (iteration .gt. 2) then
            write(6,*) 't =',t,' index =',index, 'iteration=',iteration
        endif
        if (iteration .gt. 40) then
            write(6,*) 'Emergency Exit'
            call EXIT(0)
        endif
        go to 16
      endif


      t=t+h
      tau=tau+h

      h=dt

      if (tau + h .le. tout+tol*dt) then


       go to 13
      elseif(tau .lt. tout-tol*dt)then
         h = tout - tau

         go to 13
      endif


      return
      end




