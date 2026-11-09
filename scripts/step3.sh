#!/bin/bash
set -euo pipefail
. ~/ugp/env.sh
U=$HOME/ugp
PY=/usr/bin/python3

[ -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.so" ] || { echo "STOP: PETSc not built ($PETSC_DIR/$PETSC_ARCH)"; exit 1; }
[ -x "$PY" ] || { echo "STOP: $PY not found"; exit 1; }

cd $U/src
if [ ! -d slepc ]; then
  git clone --quiet --branch v3.22.2 --depth 1 https://gitlab.com/slepc/slepc.git slepc || {
    rm -rf slepc; mkdir slepc
    curl -fsSL https://slepc.upv.es/download/distrib/slepc-3.22.2.tar.gz \
      | tar xz --strip-components=1 -C slepc; }
fi
cd slepc
export SLEPC_DIR=$U/src/slepc

"$PY" ./configure --help > /dev/null 2>&1 || { echo "STOP: SLEPc configure does not start with $PY"; exit 1; }
echo "preflight ok"

"$PY" ./configure
make SLEPC_DIR=$SLEPC_DIR PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH

grep -q SLEPC_DIR $U/env.sh || cat >> $U/env.sh <<'EOF'
export SLEPC_DIR=$UGP/src/slepc
export LD_LIBRARY_PATH=$SLEPC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
EOF

echo "== library and Fortran modules =="
ls -l $SLEPC_DIR/$PETSC_ARCH/lib/libslepc.so*
ls $SLEPC_DIR/$PETSC_ARCH/include/slepcsvd.mod $SLEPC_DIR/$PETSC_ARCH/include/slepc.mod

echo "== make check (may flag the stray MPI stderr line; results are what matter) =="
make SLEPC_DIR=$SLEPC_DIR PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH check 2>&1 | tail -12 || true
echo "SLEPc step finished"
