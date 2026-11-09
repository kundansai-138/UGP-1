#!/bin/bash
set -euo pipefail
. ~/ugp/env.sh
U=$HOME/ugp
A=$U/src/SPOD-adios2
B=$U/build/spod-adios2
A2MOD=$(dirname "$(find $U/install -name adios2.mod | head -1)")
[ -f "$A2MOD/adios2.mod" ] || { echo "STOP: adios2.mod not found"; exit 1; }
rm -rf $B; mkdir -p $B
sed 's/\r$//' "$A/Streaming SPOD/streamingSPODModule.F90" > $B/streamingSPODModule.F90
{
echo "FFTW_DIR ?= $U/install"
echo "A2_DIR ?= $U/install"
echo "A2_MOD ?= $A2MOD"
echo 'FOPT ?= -O2 -g -fno-omit-frame-pointer'
echo 'FFLAGS = $(FOPT) -ffree-line-length-none -fallow-argument-mismatch -I$(FFTW_DIR)/include -I$(A2_MOD)'
echo 'LDFLAGS = -L$(FFTW_DIR)/lib -Wl,-rpath,$(FFTW_DIR)/lib -lfftw3_mpi -lfftw3 -L$(A2_DIR)/lib -Wl,-rpath,$(A2_DIR)/lib -ladios2_fortran_mpi -ladios2_fortran -lm'
echo 'TARGET = streamingSPOD'
echo 'SRC = streamingSPODModule.o'
echo '$(TARGET): $(SRC)'
printf '\t$(FLINKER) $(FFLAGS) -o $(TARGET) $(SRC) $(SLEPC_SYS_LIB) $(LDFLAGS)\n'
echo ''
echo 'include ${SLEPC_DIR}/lib/slepc/conf/slepc_common'
} > $B/Makefile
cd $B
make > build.log 2>&1 || { tail -30 build.log; echo "STOP: build failed"; exit 1; }
echo "built: $B/streamingSPOD"
