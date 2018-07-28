#------------------------------------------------------------------------------#
# This is my first ever makefile
# Makefiles are cool because you don't have to re-compile the entire model
# every time you change a couple things in one file.
#------------------------------------------------------------------------------#

# First some variables
# := is an 'immediate set', i.e. expand ${} variables right now
# See: https://stackoverflow.com/a/448939/4970632
# Note don't expect the MKL utilities to change, so don't need to include
# them as dependencies
FC    := ifort # default is f77, ew
IOLIB := -L${NETCDFF_DIR}/lib
IOINC := -I${NETCDFF_DIR}/include
DFTI  := ${MKLROOT}/include/mkl_dfti.f90
TRIG  := ${MKLROOT}/include/mkl_trig_transforms.f90

# Main executable (depends on everything)
# Syntax is as follows (note literal tab character)
# target(s): dependency(ies)
# 	recipe
2layer.x: main.o initial.o cleanup.o forward.o prognostics.o \
          mkl_dfti.o mkl_trig_transforms.o \
          diagnostics.o transforms.o io.o global_variables.o
	$(FC) -O3 main.o initial.o cleanup.o forward.o prognostics.o \
            mkl_dfti.o mkl_trig_transforms.o \
            diagnostics.o transforms.o io.o global_variables.o \
            -mkl $(IOINC) $(IOLIB) -lnetcdff -o 2layer.x

# Dependencies
# * Note the 'use' statement in .f90 files causes compiler to search
#   for <module_name>.mod automatically; does not need to be specified in recipe
# * Note when we run <compiler> -c <file.f90>, the compiler will automatically
#   search for <module_name>.mod files in the same directory when it encounters
#   a 'use' statement.
# * Note the -mkl flag doesn't load libraries, just seems to link some mysterious
#   intel 'threading layer' that will make transforms run faster. Could still compile
#   without the flag it seems.
# Here's a good tutorial on designing fortran makefiles:
# http://genius2k.is-programmer.com/posts/40301.html
main.o: main.f90 initial.o cleanup.o forward.o prognostics.o \
        diagnostics.o transforms.o io.o global_variables.o mkl_dfti.o
	$(FC) -c -mkl $(IOINC) $(IOLIB) -lnetcdff main.f90

cleanup.o: cleanup.f90 global_variables.o mkl_dfti.o
	$(FC) -c -mkl $(IOINC) $(IOLIB) -lnetcdff cleanup.f90

transforms.o: transforms.f90 mkl_dfti.o mkl_trig_transforms.o
	$(FC) -c -mkl $(TRIG) transforms.f90

mkl_dfti.o: $(DFTI)
	$(FC) -c -mkl $(DFTI)

mkl_trig_transforms.o: $(TRIG)
	$(FC) -c -mkl $(TRIG)

global_variables.o: global_variables.f90 mkl_dfti.o
	$(FC) -c -mkl global_variables.f90

initial.o: initial.f90 transforms.o mkl_dfti.o global_variables.o
	$(FC) -c -mkl initial.f90

diagnostics.o: diagnostics.f90 transforms.o mkl_dfti.o global_variables.o
	$(FC) -c -mkl diagnostics.f90

prognostics.o: prognostics.f90 global_variables.o
	$(FC) -c prognostics.f90

forward.o: forward.f90 global_variables.o
	$(FC) -c forward.f90

io.o: io.f90 global_variables.o
	$(FC) -c $(IOINC) $(IOLIB) -lnetcdff io.f90

# Alternate version attempting to reference module files (.mod) directly
# * This one doesn't work (stuff re-compiles all the time) and sort of have
#   no idea why, oh well
# * The GFDL dry dynamical core model Makefile only references .o files as well, so even
#   though the .mod files are technically dependencies, guess it's ok to ignore them
#   as they are generated in concert with .o files
#
# main.o: main.f90 initial.mod cleanup.mod forward.mod prognostics.mod \
#         mkl_dfti.mod diagnostics.mod transforms.mod io.mod global_variables.mod
#   $(FC) -c $(IOINC) $(IOLIB) -lnetcdff main.f90
#
# mkl_dfti.o mkl_dfti.mod mkl_dft_type.mod: $(DFTI)
#   $(FC) -c -mkl $(DFTI)
#
# mkl_trig_transforms.o mkl_trig_transforms.mod mkl_tt_type.mod: $(TRIG)
#   $(FC) -c -mkl $(TRIG)
#
# transforms.o transforms.mod: transforms.f90 mkl_dfti.mod mkl_trig_transforms.mod
#   $(FC) -c -mkl transforms.f90
#
# global_variables.o global_variables.mod: global_variables.f90 mkl_dfti.mod
#   $(FC) -c -mkl global_variables.f90
#
# initial.o initial.mod: initial.f90 transforms.mod global_variables.mod mkl_dfti.mod
#   $(FC) -c -mkl initial.f90
#
# diagnostics.o diagnostics.mod: diagnostics.f90 transforms.mod global_variables.mod mkl_dfti.mod
#   $(FC) -c -mkl diagnostics.f90
#
# prognostics.o prognostics.mod: prognostics.f90 global_variables.mod
#   $(FC) -c prognostics.f90
#
# forward.o forward.mod: forward.f90 global_variables.mod
#   $(FC) -c forward.f90
#
# io.o io.mod: io.f90 global_variables.mod
#   $(FC) -c $(IOINC) $(IOLIB) -lnetcdff io.f90
#
# cleanup.o cleanup.mod: cleanup.f90 global_variables.mod mkl_dfti.mod
#   $(FC) -c -mkl $(IOINC) $(IOLIB) -lnetcdff cleanup.f90
