all: zetaeq
obj = mod_params.o mod_zeta.o mod_popfun.o mod_ncio.o mod_derives.o
f90c = ifort
#flag = -I$(NETCDF_DIR)/include -lnetcdff -lnetcdf
flag = -L/Users/hewang/Library/netcdf/lib -lnetcdff  -I/Users/hewang/Library/netcdf/include

zetaeq: zetaeq.x
debug_zeta: debug_zeta.x

zetaeq.x: main.f90 $(obj) mod_io.o
	$(f90c) -o zetaeq.x main.f90 mod_io.o $(obj) $(flag)
debug_zeta.x: test_zeta.f90 $(obj) mod_io.o
	$(f90c) -o debug_zeta.x test_zeta.f90 mod_io.o $(obj) $(flag)


mod_io.o: mod_io.f90 mod_zeta.o mod_ncio.o mod_params.o
	$(f90c) -c mod_io.f90 $(flag)
mod_zeta.o: mod_zeta.f90 mod_ncio.o mod_derives.o mod_popfun.o mod_params.o
	$(f90c) -c mod_zeta.f90 $(flag)
mod_popfun.o: mod_popfun.f90 mod_derives.o mod_params.o
	$(f90c) -c mod_popfun.f90
mod_ncio.o: mod_ncio.f90 mod_params.o
	$(f90c) -c mod_ncio.f90 $(flag)
mod_derives.o: mod_derives.f90 mod_params.o
	$(f90c) -c mod_derives.f90
mod_params.o: mod_params.f90
	$(f90c) -c mod_params.f90 $(flag)
clean:
	rm *.o *.x *.mod
