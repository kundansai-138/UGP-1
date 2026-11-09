#!/bin/bash
set -euo pipefail
. ~/ugp/env.sh
U=$HOME/ugp
SRC="$U/src/SPOD/Streaming SPOD/streamingSPODModule.F90"
B=$U/build/spod

[ -f "$SRC" ] || { echo "STOP: SPOD source missing"; exit 1; }
[ -f $U/install/include/fftw3-mpi.f03 ] || { echo "STOP: FFTW-MPI missing"; exit 1; }
[ -f $SLEPC_DIR/$PETSC_ARCH/include/slepcsvd.mod ] || { echo "STOP: SLEPc modules missing"; exit 1; }

rm -rf $B; mkdir -p $B
sed 's/\r$//' "$SRC" > $B/streamingSPODModule.F90

{
echo 'FFTW_DIR ?= '"$U"'/install'
echo 'FOPT ?= -O2 -g -fno-omit-frame-pointer'
echo 'FFLAGS = $(FOPT) -ffree-line-length-none -fallow-argument-mismatch -I$(FFTW_DIR)/include'
echo 'LDFLAGS = -L$(FFTW_DIR)/lib -Wl,-rpath,$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 -lm'
echo 'TARGET = streamingSPOD'
echo 'SRC = streamingSPODModule.o'
echo '$(TARGET): $(SRC)'
printf '\t$(FLINKER) $(FFLAGS) -o $(TARGET) $(SRC) $(SLEPC_SYS_LIB) $(LDFLAGS)\n'
echo ''
echo 'include ${SLEPC_DIR}/lib/slepc/conf/slepc_common'
} > $B/Makefile

cd $B
make 2>&1 | tee build.log
[ -x streamingSPOD ] || { echo "STOP: no executable produced"; exit 1; }

echo "== linked libraries =="
ldd ./streamingSPOD | grep -E "petsc|slepc|fftw|libmpi" || true
echo "SPOD build finished: $B/streamingSPOD"
