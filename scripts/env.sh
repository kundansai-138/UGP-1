export UGP=$HOME/ugp
export PATH=$UGP/install/bin:$PATH
export LD_LIBRARY_PATH=$UGP/install/lib:$UGP/install/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export CMAKE_PREFIX_PATH=$UGP/install
export CC=gcc-11 CXX=g++-11 FC=gfortran-11
export MPICH_CC=gcc-11 MPICH_CXX=g++-11 MPICH_FC=gfortran-11
export PETSC_DIR=$UGP/src/petsc
export PETSC_ARCH=arch-complex-opt
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
export SLEPC_DIR=$UGP/src/slepc
export LD_LIBRARY_PATH=$SLEPC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH

# --- Intel VTune (profiling only; does not affect normal runs) ---
VTUNE_VARS=/opt/intel/oneapi/vtune/2025.7/env/vars.sh
if [ -f "$VTUNE_VARS" ] && [ -z "$VTUNE_PROFILER_DIR" ]; then
    source "$VTUNE_VARS" > /dev/null 2>&1
fi
