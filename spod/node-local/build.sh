#!/bin/bash
# =====================================================================
# build.sh - build both SPOD readers of this folder (re-run after every source edit).
#
#   build/spod_old/spod_old  from src/spod_old.F90 (copy of the old ADIOS2 reader)
#   build/spod_nl/spod_nl    from src/spod_nl.F90  (node-local reader)
# Both get the same environment-variable patch (tools/patch_env_overrides.py) and the same
# compiler flags and libraries as spod/adios2-old/build.sh. It first runs the stand-alone
# node-local mapping test (plain Fortran, no MPI).
# Needs: scripts/env.sh (ADIOS2, FFTW, PETSc and SLEPc under UGP_DEPS).
# Usage: bash spod/node-local/build.sh [map|old|nl]...   (default: map old nl)
# =====================================================================
set -euo pipefail
D=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)       # spod/node-local
[[ -n ${UGP_DEPS:-} ]] || . "$D/../../scripts/env.sh"
INST=$UGP_DEPS/install
A2MOD=$INST/include/adios2/fortran
[[ -f $A2MOD/adios2.mod ]] || { echo "STOP: adios2.mod not in $A2MOD (build ADIOS2 first)"; exit 1; }
[[ -f $INST/include/fftw3-mpi.f03 ]] || { echo "STOP: FFTW with MPI not in $INST (build FFTW first)"; exit 1; }
[[ -f $SLEPC_DIR/$PETSC_ARCH/include/slepcsvd.mod ]] || { echo "STOP: SLEPc not built in $SLEPC_DIR"; exit 1; }

build_one() {  # build_one <name> <source>
  local t=$1 s=$2 b=$D/build/$1
  rm -rf "$b"; mkdir -p "$b"
  cp "$s" "$b/$t.F90"
  python3 "$D/tools/patch_env_overrides.py" "$b/$t.F90"
  {
    echo "DEPS     = $INST"
    echo "A2_MOD   = $A2MOD"
    echo 'FFLAGS   = -O2 -g -fno-omit-frame-pointer -ffree-line-length-none -fallow-argument-mismatch -I$(DEPS)/include -I$(A2_MOD)'
    echo 'LDFLAGS  = -L$(DEPS)/lib -Wl,-rpath,$(DEPS)/lib -lfftw3_mpi -lfftw3 -ladios2_fortran_mpi -ladios2_fortran -lm'
    echo "$t: $t.o"
    printf '\t$(FLINKER) $(FFLAGS) -o %s %s.o $(SLEPC_SYS_LIB) $(LDFLAGS)\n' "$t" "$t"
    echo ''
    echo 'include ${SLEPC_DIR}/lib/slepc/conf/slepc_common'
  } > "$b/Makefile"
  if ! (cd "$b" && make > build.log 2>&1); then
    grep -n -E -i -B3 -A8 'error' "$b/build.log" | head -80
    echo "STOP: $t build failed (full log: $b/build.log)"; exit 1
  fi
  [[ -x $b/$t ]] || { echo "STOP: no $b/$t"; exit 1; }
  echo "built: $b/$t"
}

map_test() {
  local b=$D/build/maptest
  rm -rf "$b"; mkdir -p "$b/mod"
  # module node_local_map = first program unit of spod_nl.F90 (up to 'end module node_local_map')
  awk '{print} /^[ \t]*end[ \t]+module[ \t]+node_local_map/{exit}' "$D/src/spod_nl.F90" > "$b/node_local_map.F90"
  "$FC" -O2 -Wall -fcheck=all -J"$b/mod" -o "$b/test_node_local_map" \
      "$b/node_local_map.F90" "$D/src/test_node_local_map.f90" > "$b/build.log" 2>&1 \
    || { cat "$b/build.log"; echo "STOP: mapping test build failed"; exit 1; }
  "$b/test_node_local_map" | tail -5
}

targets=${*:-map old nl}
for t in $targets; do
  case $t in
    old) build_one spod_old "$D/src/spod_old.F90" ;;
    nl)  build_one spod_nl  "$D/src/spod_nl.F90" ;;
    map) map_test ;;
    *)   echo "unknown target $t"; exit 1 ;;
  esac
done
