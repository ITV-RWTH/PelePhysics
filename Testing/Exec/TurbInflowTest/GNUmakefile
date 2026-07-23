# AMReX
DIM          = 3
PRECISION    = DOUBLE
PROFILE      = FALSE
VERBOSE      = FALSE
DEBUG        = FALSE

# Compiler
COMP       = gnu
FCOMP      = gfortran
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_HIP    = FALSE

# PelePhysics
Eos_Model       = Fuego
Chemistry_Model = air
Transport_Model = Constant

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    ?= ../../..

Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
