# spod-dev environment.
# Uses the libraries in ~/ugp (ADIOS2, FFTW, PETSc/SLEPc, MPICH in /usr/local) read-only:
# the same ones the S4 job uses. Nothing here writes to ~/ugp.
export UGP=$HOME/ugp
export SPOD_DEV=$HOME/spod-dev
export PATH=$UGP/install/bin:/usr/local/bin:$PATH
export CC=gcc-11 CXX=g++-11 FC=gfortran-11
export MPICH_CC=gcc-11 MPICH_CXX=g++-11 MPICH_FC=gfortran-11
export PETSC_DIR=$UGP/src/petsc
export PETSC_ARCH=arch-complex-opt
export SLEPC_DIR=$UGP/src/slepc
export LD_LIBRARY_PATH=$UGP/install/lib:$UGP/install/lib64:$PETSC_DIR/$PETSC_ARCH/lib:$SLEPC_DIR/$PETSC_ARCH/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
