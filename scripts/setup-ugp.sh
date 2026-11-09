#!/bin/bash
set -euo pipefail
U=$HOME/ugp
mkdir -p $U/{src,build,install,runs}

# --- guards -------------------------------------------------------------
for v in C_INCLUDE_PATH CPATH CPLUS_INCLUDE_PATH LIBRARY_PATH; do
  val="${!v:-}"
  case "$val" in
    *: | :* | *::*) echo "STOP: $v has an empty entry: [$val]"; exit 1;;
  esac
done
for t in gcc-11 g++-11 gfortran-11 mpicc mpicxx mpif90 cmake git; do
  command -v $t >/dev/null || { echo "MISSING: $t"; exit 1; }
done

# --- environment file ---------------------------------------------------
cat > $U/env.sh <<'ENVEOF'
export UGP=$HOME/ugp
export PATH=$UGP/install/bin:$PATH
export LD_LIBRARY_PATH=$UGP/install/lib:$UGP/install/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export CMAKE_PREFIX_PATH=$UGP/install
export CC=gcc-11 CXX=g++-11 FC=gfortran-11
export MPICH_CC=gcc-11 MPICH_CXX=g++-11 MPICH_FC=gfortran-11
ENVEOF
grep -q "ugp/env.sh" ~/.bashrc || echo '. $HOME/ugp/env.sh' >> ~/.bashrc
. $U/env.sh

mpicc -show | head -1; mpicxx -show | head -1; mpif90 -show | head -1
mpicc  -show | grep -q gcc-11      || { echo "STOP: mpicc does not use gcc-11"; exit 1; }
mpif90 -show | grep -q gfortran-11 || { echo "STOP: mpif90 does not use gfortran-11"; exit 1; }

# --- ADIOS2 2.10.2 ------------------------------------------------------
cd $U/src
[ -d ADIOS2 ] || git clone --quiet --branch v2.10.2 --depth 1 https://github.com/ornladios/ADIOS2.git
[ -z "$(git -C ADIOS2 status --porcelain)" ] || { echo "STOP: ADIOS2 source tree is modified"; exit 1; }
cmake -S ADIOS2 -B $U/build/adios2 -DCMAKE_INSTALL_PREFIX=$U/install \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$FC \
  -DADIOS2_USE_MPI=ON -DADIOS2_USE_Fortran=ON -DADIOS2_USE_SST=ON \
  -DADIOS2_USE_Python=OFF -DADIOS2_USE_HDF5=OFF -DADIOS2_USE_BZip2=OFF -DADIOS2_USE_Blosc2=OFF \
  -DADIOS2_USE_ZFP=OFF -DADIOS2_USE_SZ=OFF -DADIOS2_USE_MGARD=OFF -DADIOS2_USE_PNG=OFF \
  -DBUILD_TESTING=OFF
cmake --build $U/build/adios2 -j8 --target install
A2=$(dirname "$(find $U/install -name adios2-config.cmake | head -1)")
[ -n "$A2" ] && [ -f "$A2/adios2-config.cmake" ] || { echo "STOP: ADIOS2 cmake files not found"; exit 1; }
echo "ADIOS2 cmake dir: $A2"

# --- Incompact3d: clone once, patch the 2decomp template ----------------
cd $U/src
[ -d Incompact3d ] || git clone --quiet https://github.com/xcompact3d/Incompact3d.git
cd Incompact3d
F=cmake/decomp2d/downloadBuild2decomp.cmake.in
git checkout -- $F                      # always start from the original file
sed -i '/^ *GIT_TAG *"v2.0.4"/a\  CMAKE_ARGS        -DIO_BACKEND=@IO_BACKEND@ -DADIOS2_DIR=@ADIOS2_DIR@' $F
[ "$(grep -c '^ *CMAKE_ARGS' $F)" = 1 ] || { echo "STOP: patch did not apply cleanly"; sed -n 1,30p $F; exit 1; }
git diff > $U/adios2-decomp-forward.patch
echo "patch saved to $U/adios2-decomp-forward.patch"

FL="-g -fno-omit-frame-pointer"
cmake -S . -B $U/build/x3d_mpi -DIO_BACKEND=mpi -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="$FL"
cmake --build $U/build/x3d_mpi -j8

cmake -S . -B $U/build/x3d_adios -DIO_BACKEND=adios2 -DADIOS2_DIR=$A2 -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_Fortran_FLAGS="$FL"
cmake --build $U/build/x3d_adios -j8

# --- checks -------------------------------------------------------------
echo "adios2_init symbols (expect 0 for mpi, >0 for adios):"
nm $U/build/x3d_mpi/bin/xcompact3d   | grep -ci adios2_init || true
nm $U/build/x3d_adios/bin/xcompact3d | grep -ci adios2_init || true
echo DONE
