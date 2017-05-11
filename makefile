all:SimGoldbeter
SimGoldbeter:	
	gfortran -o program5 main-Source.f anfang.f cg.f vmap2.f \
	vmap5.f rs-Nonlinear.f out.f ODE-Merson.f ic-Source.f flow.f \
	-L/usr/bmp/pgplot-5.2/ -lpgplot \
	-L/usr/bmp/slatec-4.1/lib -lslatec \
	-L/usr/bmp/lapack-3.4.0 \
	-lX11 \
	-lpng

	