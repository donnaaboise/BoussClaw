
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw             # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:

LAPACK_DIR=/Users/calhoun/software/LAPACK/local

FFLAGS += -O0 -g -L${LAPACK_DIR}/lib -llapack -lblas

# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/1d
include $(AMRLIB)/Makefile.amr_1d

GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
# include $(GEOLIB)/Makefile.geoclaw

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = amr_module.f90 

EXCLUDE_SOURCES = b4step1.f90 rp1_geoclaw.f90


# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

BOUSSLIB= ../../../src/1d_new

MODULES =  $(BOUSSLIB)/amr_module.f90 \
	       $(BOUSSLIB)/bous_module.f90 \
	       $(BOUSSLIB)/qinit_module.f90 \
	       $(GEOLIB)/geoclaw_module.f90 \
	       $(GEOLIB)/utility_module.f90



SOURCES = \
  qinit.f90 \
  $(BOUSSLIB)/b4step1_maxrunup.f90 \
  $(BOUSSLIB)/src1.f90 \
  $(BOUSSLIB)/src1_bous.f90 \
  $(BOUSSLIB)/dlap.f90 \
  $(BOUSSLIB)/setprob.f90 \
  $(BOUSSLIB)/rp1_geoclaw.f90 \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f


#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

# Construct the topography data
.PHONY: topo all
topo:
	python maketopo.py

all: 
	$(MAKE) topo
	$(MAKE) .plots
	$(MAKE) .htmls

