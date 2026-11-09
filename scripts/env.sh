# env.sh - environment for building and running everything in this repository.
# Load it in every shell first:   . scripts/env.sh      (from the repository root)
#
#   UGP_ROOT  this repository (found from this file's location)
#   UGP_DEPS  where the libraries are built and installed (default: deps/ inside the repository):
#               $UGP_DEPS/install     ADIOS2 and FFTW (bin, lib, include)
#               $UGP_DEPS/src/petsc   PETSc (PETSC_DIR), built with PETSC_ARCH below
#               $UGP_DEPS/src/slepc   SLEPc (SLEPC_DIR)
#             Set UGP_DEPS before loading this file to use libraries built somewhere else.
#   CC / CXX / FC and MPICH_*: the compilers MPICH's wrappers use (gcc-11 family by default).
export UGP_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export UGP_DEPS=${UGP_DEPS:-$UGP_ROOT/deps}
export CC=${CC:-gcc-11} CXX=${CXX:-g++-11} FC=${FC:-gfortran-11}
export MPICH_CC=$CC MPICH_CXX=$CXX MPICH_FC=$FC
export PETSC_DIR=$UGP_DEPS/src/petsc
export PETSC_ARCH=${PETSC_ARCH:-arch-complex-opt}
export SLEPC_DIR=$UGP_DEPS/src/slepc
export PATH=$UGP_DEPS/install/bin:$PATH
export LD_LIBRARY_PATH=$UGP_DEPS/install/lib:$UGP_DEPS/install/lib64:$PETSC_DIR/$PETSC_ARCH/lib:$SLEPC_DIR/$PETSC_ARCH/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
