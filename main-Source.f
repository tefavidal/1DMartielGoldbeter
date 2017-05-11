      program main
      
      implicit none

      integer, parameter :: Nx=1200

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,tol,isf,itstart,pi,amplit,prob,tpulse

      double precision gamma01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision gamma(Nx),ro(Nx),vdx(Nx)
      double precision gammax(1000,Nx), rox(1000,Nx)
      double precision t
      double precision romin, romax
      integer ier,pgbeg, i, it, it1, printing
      character(len=30) ct1,ct2

      open(10,file ='OutputData2D/U15-Source-V1_5'
     . ,status = 'unknown',form = 'formatted')
!      ier=pgbeg(0,'OutputData2D/trajectory.ps/cps',1,1)
      t=0.d0

      call anfang(t,Nx,gamma,ro)
      call out(t,Nx,gamma,ro)

        if (tend/tout .ge. 1000)then
            printing= 0
        write(6,*) 'NO GRAPHIC OUTPUT'
        else
            printing = 1
        endif

!        if (printing .eq. 1) then
!            call vmap2(ier,1,Nx,gamma)
!
!
!!     ^^^^^^^^^^^^^First data Row^^^^^^^^^^^^^^^^^^^^^^^^^^^
!         do i=1,Nx
!          gammax(1,i)=gamma(i)
!          rox(1,i)=ro(i)
!         enddo
!!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!        endif


      call flow(t,Nx,vdx)
 5    continue

      it1=(t-tE)/tout
      call ODE(t,Nx,gamma,ro,vdx)
!      call noise(Nx,gamma)

      it=(t+0.000001-tE)/tout

      write(6,*) 'real t= '
      write(6,'(F6.2)') t/dk1



      if (it .gt. it1) then
      call out(t,Nx,gamma,ro)

         romin = ro(1)
         romax = ro(1)
         do i=1,Nx
!            if (printing .eq. 1) then
!                 gammax(it+1,i)=gamma(i)
!                    rox(it+1,i)=ro(i)
!            endif
          romin = min(romin,ro(i))
          romax = max(romax,ro(i))
         enddo
         write(6,*) 'ro=',romin,romax
!        if (printing .eq. 1) then
!            call vmap2(ier,it+1,Nx,gamma)
!        endif

      endif


      if (t+dt .lt. tend) then
         goto 5

      endif
      close(10)
      call pgend

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     WRITES FINAL STATE
      open(42,file ='OutputData2D/Final-State'
     . ,status = 'unknown',form = 'formatted')
      call outFinal(t,Nx,gamma,ro)
      close(42)

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!      SPACE-TIME PLOTS
!        if (printing .eq. 1) then
!      ier=pgbeg(0,'OutputData2D/Gamma-spacetime.ps/cps',1,1)
!
!      call vmap5(1,Nx,it+1,gammax)
!
!      call pgend
!
!      ier=pgbeg(0,'OutputData2D/Rho-spacetime.ps/cps',1,1)
!
!      call vmap5(1,Nx,it+1,rox)
!
!      call pgend
!        endif

      end
      

