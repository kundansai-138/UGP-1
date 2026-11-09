#!/bin/bash
set -euo pipefail
. ~/ugp/env.sh
U=$HOME/ugp
PY=/usr/bin/python3
S=$U/src/SPOD
A=$U/src/SPOD-adios2
B=$U/build/spod-adios2
F="Streaming SPOD/streamingSPODModule.F90"
D=$U/data/synth_n20000_t96_d32

[ -f $U/tools/patch_userinput.py ] || { echo "STOP: run step5.sh first"; exit 1; }
[ -f $D/done ] || { echo "STOP: synthetic data missing (run step5.sh first)"; exit 1; }
A2MOD=$(dirname "$(find $U/install -name adios2.mod | head -1)")
[ -f "$A2MOD/adios2.mod" ] || { echo "STOP: adios2.mod not found under $U/install"; exit 1; }
echo "adios2 modules: $A2MOD"

# ---- 1. working copy on its own branch
[ -d $A ] || git clone --quiet $S $A
cd $A
git checkout -q -B adios2
$PY $U/tools/patch_userinput.py "$F"
git add -A
git -c user.name=ugp -c user.email=ugp@localhost commit -q -m "SPOD: environment-variable overrides in UserInput" || true
git tag -f env-baseline > /dev/null
git log --oneline -3

# ---- 2. reusable build script (re-run it after every source edit)
cat > $U/tools/build_spod_adios2.sh <<'BLD'
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
BLD
chmod +x $U/tools/build_spod_adios2.sh
$U/tools/build_spod_adios2.sh

# ---- 3. ADIOS2 Fortran link test with the same flags
cd $B
cat > a2test.f90 <<'F90'
program a2test
  use adios2
  implicit none
  include 'mpif.h'
  type(adios2_adios) :: adios
  integer :: ierr
  call MPI_Init(ierr)
  call adios2_init(adios, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) stop 1
  print *, 'adios2 init ok'
  call adios2_finalize(adios, ierr)
  call MPI_Finalize(ierr)
end program
F90
mpif90 a2test.f90 -I$A2MOD -L$U/install/lib -Wl,-rpath,$U/install/lib -ladios2_fortran_mpi -ladios2_fortran -o a2test
env -u DISPLAY -u XAUTHORITY mpiexec -n 2 ./a2test 2>&1 | grep -v Authorization

# ---- 4. regression: same smoke test on this build
R=$U/runs/adios2copy-baseline-$(date +%Y%m%d-%H%M%S); mkdir -p $R; cd $R
export SPOD_NPOINTS=20000 SPOD_NSNAP=96 SPOD_NDFT=32 SPOD_K=3 SPOD_OVERLAP=0.5 SPOD_NFREQSAVE=3
export SPOD_DATADIR=$D SPOD_WEIGHTS=none
env -u DISPLAY -u XAUTHORITY mpiexec -n 4 $B/streamingSPOD > spod.log 2>&1 \
  || { echo "RUN FAILED"; tail -30 spod.log; exit 1; }
got=$(grep "Eigen values at index" spod.log | awk '{print $5}' | tr '\n' ' ')
echo "expected: 4 10 14   got: $got"
if [ "$got" = "4 10 14 " ]; then echo "ADIOS2-COPY BASELINE: PASS"; else echo "DIFFERENT (see $R/spod.log)"; fi
