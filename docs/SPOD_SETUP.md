# Streaming SPOD setup

Two SPOD readers are in this repository. Both run the same streaming SPOD maths (PETSc/SLEPc,
complex scalars); they differ only in which part of each snapshot a rank reads.

| Folder | Reader | Server location |
|---|---|---|
| `spod/adios2-old/` | **old ADIOS2 reader**: each rank reads whole x-y planes | `~/ugp/src/SPOD-adios2` |
| `spod/node-local/` | **node-local reader**: each rank reads only data written on its own node | `~/spod-dev` |

Do [XCOMPACT3D_SETUP.md](XCOMPACT3D_SETUP.md) sections 1–3 first (toolchain, `env.sh`, ADIOS2).

## 1. Libraries

### FFTW 3.3.10 with MPI

The tarball is in `spod/adios2-old/Required libraries/fftw-3.3.10.tar.gz`.

```bash
mkdir -p ~/ugp/src/fftw ~/ugp/build/fftw
tar -xzf fftw-3.3.10.tar.gz -C ~/ugp/src/fftw
cd ~/ugp/build/fftw
~/ugp/src/fftw/fftw-3.3.10/configure --prefix=$HOME/ugp/install --enable-mpi --enable-shared MPICC=mpicc
make -j8 && make install
ls ~/ugp/install/include/fftw3-mpi.f03     # must exist
```

### PETSc 3.22.2, complex scalars (`scripts/step2.sh`)

```bash
. ~/ugp/env.sh
cd ~/ugp/src && git clone --branch v3.22.2 --depth 1 https://gitlab.com/petsc/petsc.git petsc
cd petsc
FL="-O3 -g -fno-omit-frame-pointer"
env -u CC -u CXX -u FC /usr/bin/python3 ./configure PETSC_ARCH=arch-complex-opt \
  --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
  --with-scalar-type=complex --with-precision=double --with-debugging=0 \
  COPTFLAGS="$FL" CXXOPTFLAGS="$FL" FOPTFLAGS="$FL" \
  --download-fblaslapack --with-shared-libraries=1 --with-x=0
make MAKE_NP=8 all && make check
```

Use `/usr/bin/python3`; anaconda's Python segfaults in PETSc's configure.

### SLEPc 3.22.2 (`scripts/step3.sh`)

```bash
cd ~/ugp/src && git clone --branch v3.22.2 --depth 1 https://gitlab.com/slepc/slepc.git slepc
cd slepc && export SLEPC_DIR=$PWD
/usr/bin/python3 ./configure
make SLEPC_DIR=$SLEPC_DIR PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH
make SLEPC_DIR=$SLEPC_DIR PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH check
```

`env.sh` already exports `PETSC_DIR`, `PETSC_ARCH` and `SLEPC_DIR`.

### Optional: the original file-based SPOD

`scripts/step4.sh` builds the upstream SPOD (reads `.bin` snapshot files) into `~/ugp/build/spod`.
`scripts/step5.sh` adds environment-variable settings to it (`ugp-tools/patch_userinput.py`)
and runs a smoke test on synthetic data with a known answer: the leading modes must come out
at frequency bins 4, 10 and 14 (`SMOKE TEST: PASS`).

## 2. Old ADIOS2 reader (`spod/adios2-old`)

`Streaming SPOD/streamingSPODModule.F90` reads the SST stream of Xcompact3d instead of files.

- The snapshot vector is `q = [ux; uy; uz; pp]`, 4 × 301³ = 109,083,604 rows.
- PETSc splits the rows into equal contiguous ranges. Each rank reads the whole x-y planes
  that cover its range, so every rank needs data from every Xcompact3d rank.
- It opens IO `spod-input`. If `spod_sst.xml` is not in the folder it runs in (the normal case),
  it uses a built-in SST engine.

**Its settings are written in `Subroutine UserInput`** (edit and rebuild to change them):

| Setting | Value |
|---|---|
| Grid | `gridNx = gridNy = gridNz = 301` |
| Variables | `ux, uy, uz, pp` |
| Stream | `streamName = 'data'` (run SPOD in the Xcompact3d run folder) |
| Snapshots | `numberOfSnapShot = 31` |
| SPOD | `nDFT = 8`, `percentageOfOverall = 0.5` (overlap 4), `k = 3` |
| Output | `nFreqToSave = 5` (5 mode files), equal weights |

`streamingSPODModule.F90.65grid` is the same reader set up for 65³ tests.

**Build** (`scripts/step6.sh` creates the working copy, `scripts/ugp-tools/build_spod_adios2.sh`
builds it):

```bash
bash ~/ugp/tools/build_spod_adios2.sh      # -> ~/ugp/build/spod-adios2/streamingSPOD
```

`spod/adios2-old/Makefile.ugp` is the Makefile it generates: `-O2 -g -ffree-line-length-none
-fallow-argument-mismatch`, FFTW and the ADIOS2 Fortran module on the include path, and
`-lfftw3_mpi -lfftw3 -ladios2_fortran_mpi -ladios2_fortran` plus SLEPc's `slepc_common`.

## 3. Node-local reader (`spod/node-local`)

### What it changes

The SPOD maths is unchanged. Only the layout of the rows changes:

1. On each node, the reader lists the pencils of the Xcompact3d ranks running on that node,
   line by line (one line = one run along x), and splits that list evenly over the SPOD ranks
   of the node. Every piece lies inside one Xcompact3d rank's pencil on the same node, so
   ADIOS2 serves it locally.
2. Every PETSc object with one row per grid point is created with this local row count instead of
   `PETSC_DECIDE`.
3. A `VecScatter` converts between this order and the natural order: once for the weights at
   start-up, and once per saved mode before writing. Mode files therefore have exactly the old
   layout.
4. At start-up it checks that the reordering covers every grid point exactly once (the run stops
   otherwise) and prints one line per node: host, SPOD ranks, and the Xcompact3d ranks it reads.

It assumes X-pencils, the same `p_row × p_col` as Xcompact3d, and both programs filling the same
nodes in the same order. If the last two are wrong, the results stay correct; only the reads stop
being local.

### Files

| File | Purpose |
|---|---|
| `src/spod_nl.F90` | node-local reader (module `node_local_map` + SPOD + program) |
| `src/spod_old.F90` | copy of the old reader, so both are built the same way |
| `src/test_node_local_map.f90` | stand-alone test of the mapping |
| `env.sh`, `setup.sh`, `build.sh` | environment, one-time set-up, build |
| `tools/patch_env_overrides.py` | lets both readers take their settings from environment variables |
| `tools/run_S4_nl.sh`, `tools/queue_S4.sh`, `run_compare.sh` | runs (see RUNNING.md) |
| `tools/compare_spod_outputs.py`, `tools/net_rates.py`, `tools/wait_analysis.py` | analysis |
| `design/` | design notes and the Python model used to check the mapping |

### Build

```bash
cp -r spod/node-local ~/spod-dev        # the scripts expect ~/spod-dev
bash ~/spod-dev/setup.sh                # copies the old reader, the Xcompact3d binary and the case from ~/ugp
bash ~/spod-dev/build.sh                # mapping test + both readers
```

`build.sh` first runs the mapping test (expect `node_local_map: 720 configurations PASS, 0 FAIL`),
then builds `build/spod_old/spod_old` and `build/spod_nl/spod_nl` with the same flags and
libraries as the old reader. Rebuild after any source change with `bash ~/spod-dev/build.sh nl`.

### Settings (environment variables)

Unset or empty variables keep the value in `UserInput`.

| Variable | Setting |
|---|---|
| `SPOD_NX`, `SPOD_NY`, `SPOD_NZ` | grid |
| `SPOD_STREAM` | stream name or path (contact file `<name>.sst`) |
| `SPOD_CONFIG` | ADIOS2 XML for IO `spod-input` (leave unset for the built-in SST engine) |
| `SPOD_NSNAP`, `SPOD_NDFT`, `SPOD_K` | snapshots, DFT length, modes kept |
| `SPOD_OVERLAP` | overlap as a **fraction** (0.5), not a snapshot count |
| `SPOD_NSAVE` | number of mode files written |
| `SPOD_WEIGHTS` | weight file, `none` = equal weights |
| `X3D_PROW`, `X3D_PCOL` | Xcompact3d's processor grid (node-local reader only) |
| `X3D_RANKS_PER_NODE` | 0 = ranks fill the same nodes in order (default) |

## 4. Output files

| File | Contents |
|---|---|
| `EnegySpectra.txt` | singular values per mode and frequency (PETSc ASCII matrix) |
| `SPODMODE_freq_<i>_sigma_<σ>.dat` | one mode: PETSc binary complex vector, natural order `[ux; uy; uz; pp]`, x fastest; 1.75 GB at 301³ |

Compare two runs with `python3 tools/compare_spod_outputs.py RUN_A RUN_B`. It allows for the
phase freedom of SPOD modes (a mode is only defined up to a factor e^iθ).
