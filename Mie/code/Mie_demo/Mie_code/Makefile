F90 = ifort #gfortran #f90
F77 = ifort #gfortran #ifort #f77
FLAGS= -r8 -O3
lib = MIEV0.o RDI1MACHmodule.o  ErrPack.o write_data.o
driver_lib = gamma_size.o REFWAT.o

#NetcdfHome = /Users/softwares/netcdf
NetcdfHome = /umbc/xfs1/zzbatmos/common/Exchange/MSCART_lib/netcdf/IFC
Netcdf_IncludeDir = $(NetcdfHome)/include
NetcdfLibs = -L$(NetcdfHome)/lib -lnetcdf  -lz -lm  -lcurl

all: single driver

single: mie_single_size.f90 ${lib}
	$(F90) $(FLAGS) -o mie_single_size.exe mie_single_size.f90 ${lib} -g -I$(Netcdf_IncludeDir) ${NetcdfLibs}

driver: driver.f90 ${driver_lib}
	$(F90) $(FLAGS) -o mie_single_size_driver.exe driver.f90 ${driver_lib}

gamma_size.o: gamma_size.f90
	$(F90) -c $(FLAGS) gamma_size.f90

REFWAT.o:REFWAT.f
	$(F90) -c $(FLAGS) REFWAT.f

MIEV0.o: MIEV0.f
	$(F77) -c $(FLAGS) MIEV0.f

RDI1MACHmodule.o: RDI1MACHmodule.F90
	$(F90) -c $(FLAGS) RDI1MACHmodule.F90

ErrPack.o: ErrPack.f
	$(F77) -c $(FLAGS) ErrPack.f

write_data.o : write_data.f90
	$(F90) -c $(FLAGS) write_data.f90

clean:
	rm *.o
	rm *.mod
	rm *.exe
