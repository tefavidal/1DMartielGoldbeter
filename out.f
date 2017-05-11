      subroutine out(t,Nx,gamma,ro)

      implicit none

      double precision t
      integer Nx,i

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx),ro(Nx)
      double precision dls



      dls=dk1/(dke0*Diffgamma)**0.5

      do i=1,Nx
        write(10,*) t/dk1,i*dx/dls,gamma(i),ro(i)
      enddo
      write(10,*)

      return
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine outFinal(t,Nx,gamma,ro)

      implicit none
      integer Nx, i
      double precision t


      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx),ro(Nx)


      do i=1,Nx
            write(42,*) t/dk1,gamma(i),ro(i)
      enddo
      write(42,*)
      return
      end
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine loadState(t,Nx,gamma,ro)

      implicit none

      integer Nx, i
      double precision t, aux

      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx),ro(Nx)



      do i=1,Nx
            read(7,*) aux,gamma(i),ro(i)
      enddo
      close(7)
      t=aux*dk1
      return
      end

