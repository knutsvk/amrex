DEBUG = FALSE

USE_EB = TRUE

USE_MPI  = TRUE
USE_OMP  = FALSE

USE_HYPRE = FALSE
USE_PETSC = FALSE

COMP = gnu

DIM = 3

AMREX_HOME ?= ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs 	:= Base Boundary AmrCore 
Pdirs 	+= LinearSolvers/C_CellMG LinearSolvers/MLMG
Pdirs 	+= EB

ifeq ($(USE_HYPRE),TRUE)
    Pdirs += Extern/HYPRE
endif

ifeq ($(USE_PETSC),TRUE)
    Pdirs += Extern/PETSc
endif

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

