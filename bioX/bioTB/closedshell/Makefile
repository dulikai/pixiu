# make

bioparc: constant.o wavefunction.o sub_read.o sub_fock.o sub_site_parm.o sub_Lowdin.o sub_frontier.o main.o
	gfortran constant.o wavefunction.o sub_read.o sub_fock.o sub_site_parm.o sub_Lowdin.o sub_frontier.o main.o -o bioparc

constant.o: constant.f90
	gfortran -c constant.f90

wavefunction.o: wavefunction.f90
	gfortran -c wavefunction.f90

sub_read.o: sub_read.f90
	gfortran -c sub_read.f90

sub_fock.o: sub_fock.f90
	gfortran -c sub_fock.f90

sub_site_parm.o: sub_site_parm.f90
	gfortran -c sub_site_parm.f90

sub_Lowdin.o: sub_Lowdin.f90
	gfortran -c sub_Lowdin.f90

sub_frontier.o: sub_frontier.f90
	gfortran -c sub_frontier.f90

main.o: main.f90
	gfortran -c main.f90

clean:
	rm *.o

