#!/bin/bash
set -euo pipefail
. ~/ugp/env.sh
U=$HOME/ugp
PY=/usr/bin/python3     # NOT anaconda: its python 3.13 segfaults in PETSc configure

[[ "$(mpicc -show)"  == *gcc-11*      ]] || { echo "STOP: mpicc is not gcc-11"; exit 1; }
[[ "$(mpif90 -show)" == *gfortran-11* ]] || { echo "STOP: mpif90 is not gfortran-11"; exit 1; }
[ -x "$PY" ] || { echo "STOP: $PY not found"; exit 1; }
"$PY" -c 'import sys; print("using python", sys.version.split()[0])'

cd $U/src
if [ ! -d petsc ]; then
  git clone --quiet --branch v3.22.2 --depth 1 https://gitlab.com/petsc/petsc.git petsc || {
    rm -rf petsc; mkdir petsc
    curl -fsSL https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.22.2.tar.gz \
      | tar xz --strip-components=1 -C petsc; }
fi
cd petsc

# preflight: does configure even start with this interpreter?
"$PY" ./configure --help > /dev/null 2>&1 || { echo "STOP: configure does not start with $PY"; exit 1; }
echo "preflight ok: configure starts"

export PETSC_DIR=$U/src/petsc
export PETSC_ARCH=arch-complex-opt
FL="-O3 -g -fno-omit-frame-pointer"

env -u CC -u CXX -u FC "$PY" ./configure PETSC_ARCH=$PETSC_ARCH \
    --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
    --with-scalar-type=complex --with-precision=double --with-debugging=0 \
    COPTFLAGS="$FL" CXXOPTFLAGS="$FL" FOPTFLAGS="$FL" \
    --download-fblaslapack --with-shared-libraries=1 --with-x=0

make MAKE_NP=8 all
make check

grep -q PETSC_DIR $U/env.sh || cat >> $U/env.sh <<'EOF'
export PETSC_DIR=$UGP/src/petsc
export PETSC_ARCH=arch-complex-opt
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
EOF
echo "PETSc step finished"
