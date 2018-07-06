# ============================================================================
# Name        : Makefile
# Author      : Roy Walters
# Version     :
# Copyright   : Copyright (C) 2012 under GNU General Public License
# Description : Makefile for RiCOMToolkit  (needs a lot of work yet)
# ============================================================================

.PHONY: all reorder rcmha refinex2 zoom_sigma decomp clean

# Change these lines if you are using different
# Fortran or C compilers and compiler options.
# Portland Group fortran
#FC = pgf90
#F_FLAGS = -r8 -O2 -fast
#DBG_FLAGS = -r8 -C -g -traceback
# intel fortran
#FC = ifort
#F_FLAGS = -r8
#DBG_FLAGS = -r8 -fpe0 -check bounds -check uninit -traceback
# gfortran
FC = gfortran
F_FLAGS = -O2
DBG_FLAGS = -g -O2 -fcheck=bounds -fbacktrace

libvtk0:
	$(FC) $(F_FLAGS) -c \
		src-vtk0/penf/penf_global_parameters_variables.F90 \
		src-vtk0/penf/penf_b_size.F90 \
		src-vtk0/penf/penf_stringify.F90
	$(FC) $(F_FLAGS) -c \
		src-vtk0/Lib_VTK_IO_Back_End.F90 \
		src-vtk0/Lib_VTK_IO_INI.f90 \
		src-vtk0/Lib_VTK_IO_GEO.f90 \
		src-vtk0/Lib_VTK_IO_CON.f90 \
		src-vtk0/Lib_VTK_IO_DAT_VAR.f90 \
		src-vtk0/Lib_VTK_IO_END.f90 \
		src-vtk0/Lib_VTK_IO.f90
	rm *.o
	rm penf_global_parameters_variables.mod
	rm penf_b_size.mod
	rm penf_stringify.mod
	rm lib_vtk_io_back_end.mod
	rm lib_vtk_io_ini.mod
	rm lib_vtk_io_geo.mod
	rm lib_vtk_io_con.mod
	rm lib_vtk_io_dat_var.mod
	rm lib_vtk_io_end.mod

libvtk0_dbg:
	$(FC) $(DBG_FLAGS) -c \
		src-vtk0/penf/penf_global_parameters_variables.F90 \
		src-vtk0/penf/penf_b_size.F90 \
		src-vtk0/penf/penf_stringify.F90
	$(FC) $(DBG_FLAGS) -c \
		src-vtk0/Lib_VTK_IO_Back_End.F90 \
		src-vtk0/Lib_VTK_IO_INI.f90 \
		src-vtk0/Lib_VTK_IO_GEO.f90 \
		src-vtk0/Lib_VTK_IO_CON.f90 \
		src-vtk0/Lib_VTK_IO_DAT_VAR.f90 \
		src-vtk0/Lib_VTK_IO_END.f90 \
		src-vtk0/Lib_VTK_IO.f90
	rm *.o
	rm penf_global_parameters_variables.mod
	rm penf_b_size.mod
	rm penf_stringify.mod
	rm lib_vtk_io_back_end.mod
	rm lib_vtk_io_ini.mod
	rm lib_vtk_io_geo.mod
	rm lib_vtk_io_con.mod
	rm lib_vtk_io_dat_var.mod
	rm lib_vtk_io_end.mod

rcm2vtk0:
	$(FC) $(F_FLAGS) -c \
		src-vtk0/penf/penf_global_parameters_variables.F90 \
		src-vtk0/penf/penf_b_size.F90 \
		src-vtk0/penf/penf_stringify.F90
	$(FC) $(F_FLAGS) -c \
		src-vtk0/Lib_VTK_IO_Back_End.F90 \
		src-vtk0/Lib_VTK_IO_INI.f90 \
		src-vtk0/Lib_VTK_IO_GEO.f90 \
		src-vtk0/Lib_VTK_IO_CON.f90 \
		src-vtk0/Lib_VTK_IO_DAT_VAR.f90 \
		src-vtk0/Lib_VTK_IO_END.f90 \
		src-vtk0/Lib_VTK_IO.f90
	$(FC) -c src/RCM2VTK.f90 
	mkdir -p bin
	$(FC) *.o -o bin/rcm2vtk0
	rm *.o
	rm *.mod

rcm2vtk0_dbg:
	$(FC) $(DBG_FLAGS) -c \
		src-vtk0/penf/penf_global_parameters_variables.F90 \
		src-vtk0/penf/penf_b_size.F90 \
		src-vtk0/penf/penf_stringify.F90
	$(FC) $(DBG_FLAGS) -c \
		src-vtk0/Lib_VTK_IO_Back_End.F90 \
		src-vtk0/Lib_VTK_IO_INI.f90 \
		src-vtk0/Lib_VTK_IO_GEO.f90 \
		src-vtk0/Lib_VTK_IO_CON.f90 \
		src-vtk0/Lib_VTK_IO_DAT_VAR.f90 \
		src-vtk0/Lib_VTK_IO_END.f90 \
		src-vtk0/Lib_VTK_IO.f90
	$(FC) $(DBG_FLAGS) -c src/RCM2VTK.f90 
	mkdir -p bin
	$(FC) *.o -o bin/rcm2vtk0_dbg
	rm *.o
	rm *.mod

rcm2vtk:
	$(FC) $(F_FLAGS) -c -cpp \
		src-vtk/IR_Precision.f90 \
		src-vtk/Lib_VTK_IO.f90
	$(FC) -c src/RCM2VTK.f90 
	mkdir -p bin
	$(FC) *.o -o bin/rcm2vtk
	rm *.o
	rm *.mod

rcm2vtk_dbg:
	$(FC) $(DBG_FLAGS) -c -cpp \
		src-vtk/IR_Precision.f90 \
		src-vtk/Lib_VTK_IO.f90
	$(FC) $(DBG_FLAGS) -c src/RCM2VTK.f90 
	mkdir -p bin
	$(FC) *.o -o bin/rcm2vtk_dbg
	rm *.o
	rm *.mod

clean:
	rm -f *.o *.mod
