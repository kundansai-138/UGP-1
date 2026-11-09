#!/bin/bash
# =====================================================================
# build.sh - build streaming SPOD with the old ADIOS2 reader.
#   source : spod/adios2-old/Streaming SPOD/streamingSPODModule.F90
#   result : spod/adios2-old/build/streamingSPOD
# Needs: scripts/env.sh (ADIOS2, FFTW, PETSc and SLEPc under UGP_DEPS).
# Usage: bash spod/adios2-old/build.sh
#        SRC="spod/adios2-old/Streaming SPOD/streamingSPODModule.F90.65grid" bash spod/adios2-old/build.sh
#        (the 65^3 variant; any other source file the same way)
# =====================================================================
set -euo pipefail
D=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)       # spod/adios2-old
[[ -n ${UGP_DEPS:-} ]] || . "$D/../../scripts/env.sh"
SRC=${SRC:-$D/Streaming SPOD/streamingSPODModule.F90}
B=$D/build
INST=$UGP_DEPS/install
A2MOD=$INST/include/adios2/fortran
[[ -f $SRC ]] || { echo "STOP: source not found: $SRC"; exit 1; }
[[ -f $A2MOD/adios2.mod ]] || { echo "STOP: adios2.mod not in $A2MOD (build ADIOS2 first)"; exit 1; }
[[ -f $INST/include/fftw3-mpi.f03 ]] || { echo "STOP: FFTW with MPI not in $INST (build FFTW first)"; exit 1; }
[[ -f $SLEPC_DIR/$PETSC_ARCH/include/slepcsvd.mod ]] || { echo "STOP: SLEPc not built in $SLEPC_DIR"; exit 1; }

rm -rf "$B"; mkdir -p "$B"
sed 's/\r$//' "$SRC" > "$B/streamingSPODModule.F90"      # the folder name has a space; build here
{
  echo "DEPS   = $INST"
  echo "A2_MOD = $A2MOD"
  echo 'FFLAGS = -O2 -g -fno-omit-frame-pointer -ffree-line-length-none -fallow-argument-mismatch -I$(DEPS)/include -I$(A2_MOD)'
  echo 'LDFLAGS = -L$(DEPS)/lib -Wl,-rpath,$(DEPS)/lib -lfftw3_mpi -lfftw3 -ladios2_fortran_mpi -ladios2_fortran -lm'
  echo 'streamingSPOD: streamingSPODModule.o'
  printf '\t$(FLINKER) $(FFLAGS) -o streamingSPOD streamingSPODModule.o $(SLEPC_SYS_LIB) $(LDFLAGS)\n'
  echo ''
  echo 'include ${SLEPC_DIR}/lib/slepc/conf/slepc_common'
} > "$B/Makefile"
cd "$B"
make > build.log 2>&1 || { tail -30 build.log; echo "STOP: build failed (full log: $B/build.log)"; exit 1; }
echo "built: $B/streamingSPOD"
