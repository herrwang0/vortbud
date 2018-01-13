all: zetaeq meaneddy
obj = mod_params.o mod_zeta.o mod_popfun.o mod_popload.o mod_ncio.o mod_derives.o
f90c = ifort
# flag = -I$(NETCDF_DIR)/include -lnetcdff
flag = -L/Users/hewang/Library/netcdf/lib -lnetcdff  -I/Users/hewang/Library/netcdf/include

zetaeq: zetaeq.x
meaneddy: zetaeq_meaneddy.x

zetaeq.x: main.f90 $(obj) mod_control.o
	$(f90c) -o zetaeq.x main.f90 mod_control.o $(obj) $(flag)
zetaeq_meaneddy.x: main_eddy.f90 $(obj) mod_control_eddy.o
	$(f90c) -o zetaeq_meaneddy.x main_eddy.f90 mod_control_eddy.o $(obj) $(flag)

mod_control_eddy.o: mod_control_eddy.f90 mod_zeta.o mod_popload.o mod_ncio.o mod_params.o
	$(f90c) -c mod_control_eddy.f90 $(flag)
mod_control.o: mod_control.f90 mod_zeta.o mod_popload.o mod_ncio.o mod_params.o
	$(f90c) -c mod_control.f90 $(flag)
mod_popload.o: mod_popload.f90 mod_ncio.o mod_params.o
	$(f90c) -c mod_popload.f90 $(flag)
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
