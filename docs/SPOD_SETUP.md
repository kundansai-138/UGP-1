# Streaming SPOD setup

This guide builds the libraries SPOD needs and both SPOD readers. Do sections 1–3 of
[XCOMPACT3D_SETUP.md](XCOMPACT3D_SETUP.md) first (requirements, `scripts/env.sh`, ADIOS2).
All paths are relative to the cloned repository; load the environment first:

```bash
. scripts/env.sh
```

Both readers run the same streaming SPOD (PETSc/SLEPc with complex numbers). They differ only in
which part of each snapshot an MPI rank reads:

| Folder | Reader | What each SPOD rank reads |
|---|---|---|
| `spod/adios2-old/` | old ADIOS2 reader | whole x-y planes, from every Xcompact3d rank |
| `spod/node-local/` | node-local reader | only data written by Xcompact3d ranks on its own machine |

## 1. Libraries

### FFTW 3.3.10 with MPI

The source is included in `spod/adios2-old/Required libraries/`.

```bash
mkdir -p $UGP_DEPS/src $UGP_DEPS/build/fftw
tar -xzf "spod/adios2-old/Required libraries/fftw-3.3.10.tar.gz" -C $UGP_DEPS/src
cd $UGP_DEPS/build/fftw
$UGP_DEPS/src/fftw-3.3.10/configure --prefix=$UGP_DEPS/install --enable-mpi --enable-shared MPICC=mpicc
make -j8 && make install
cd $UGP_ROOT
ls $UGP_DEPS/install/include/fftw3-mpi.f03         # must exist
```

### PETSc 3.22.2 with complex numbers

```bash
cd $UGP_DEPS/src
git clone --branch v3.22.2 --depth 1 https://gitlab.com/petsc/petsc.git petsc
cd petsc
FL="-O3 -g -fno-omit-frame-pointer"
env -u CC -u CXX -u FC python3 ./configure PETSC_ARCH=$PETSC_ARCH \
  --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
  --with-scalar-type=complex --with-precision=double --with-debugging=0 \
  COPTFLAGS="$FL" CXXOPTFLAGS="$FL" FOPTFLAGS="$FL" \
  --download-fblaslapack --with-shared-libraries=1 --with-x=0
make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH check
cd $UGP_ROOT
```

`--download-fblaslapack` needs internet access during `configure`. Use a system Python 3;
an Anaconda Python crashed in PETSc's `configure` on the lab machines.

### SLEPc 3.22.2

```bash
cd $UGP_DEPS/src
git clone --branch v3.22.2 --depth 1 https://gitlab.com/slepc/slepc.git slepc
cd slepc
python3 ./configure
make SLEPC_DIR=$SLEPC_DIR PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH
make SLEPC_DIR=$SLEPC_DIR PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH check
cd $UGP_ROOT
```

## 2. Old ADIOS2 reader (`spod/adios2-old`)

The source is `spod/adios2-old/Streaming SPOD/streamingSPODModule.F90`.

How it reads: the snapshot vector is `q = [ux; uy; uz; pp]` (4 × 301³ = 109,083,604 rows).
PETSc gives every rank an equal, contiguous range of rows; each rank reads the whole x-y planes
that cover its range. It opens the ADIOS2 IO `spod-input`; when `spod_sst.xml` is not in the
folder it runs in (the normal case), it uses a built-in SST engine.

Its settings are written in `Subroutine UserInput`. To change them, edit that subroutine and
build again:

| Setting | Value |
|---|---|
| Grid | `gridNx = gridNy = gridNz = 301` |
| Variables | `ux, uy, uz, pp` |
| Stream | `streamName = 'data'`: SPOD runs in the Xcompact3d run folder |
| Snapshots | `numberOfSnapShot = 31` |
| SPOD | `nDFT = 8`, `percentageOfOverall = 0.5` (overlap 4 snapshots), `k = 3` |
| Output | `nFreqToSave = 5` mode files, equal weights |

`streamingSPODModule.F90.65grid` is the same reader set up for small tests: 65³ grid,
10 snapshots, `nDFT = 4`, 3 mode files.

Build:

```bash
bash spod/adios2-old/build.sh              # -> spod/adios2-old/build/streamingSPOD
```

`spod/adios2-old/Makefile.ugp` is the Makefile used on the lab machines, kept for reference;
`build.sh` writes the equivalent one for this repository.

For the 65³ variant (it replaces the build above):

```bash
SRC="spod/adios2-old/Streaming SPOD/streamingSPODModule.F90.65grid" bash spod/adios2-old/build.sh
```

## 3. Node-local reader (`spod/node-local`)

### What it changes

The SPOD maths is the same as in the old reader. Only the way rows are assigned to ranks changes:

1. On each machine, the reader lists the pencils of the Xcompact3d ranks running on that machine,
   line by line (a line is one run of points along x), and divides that list evenly among the
   SPOD ranks of the machine. Every piece lies inside one Xcompact3d rank's pencil on the same
   machine, so ADIOS2 delivers it without using the network.
2. Every PETSc object with one row per grid point is created with this number of local rows.
3. The weights are put into this order once at start-up, and every saved mode is put back into
   the normal order before it is written. Mode files therefore look exactly like the old
   reader's.
4. At start-up it checks that every grid point is read exactly once (the run stops otherwise) and
   prints one line per machine: machine name, SPOD ranks, and the Xcompact3d ranks it reads.

It assumes X-pencils, the same `p_row × p_col` as Xcompact3d, and both programs filling the same
machines in the same order. If the last two do not hold, the results are still correct; only the
reads are no longer local.

More detail: [spod/node-local/design/NODE_LOCAL_DESIGN.md](../spod/node-local/design/NODE_LOCAL_DESIGN.md).

### Files

| File | Purpose |
|---|---|
| `src/spod_nl.F90` | node-local reader: module `node_local_map`, the SPOD module and the program in one file |
| `src/spod_old.F90` | copy of the old reader, built the same way, for side-by-side tests |
| `src/test_node_local_map.f90` | stand-alone test of the mapping |
| `build.sh` | runs the mapping test and builds both readers |
| `run_compare.sh` | correctness test on one machine (see RUNNING.md) |
| `tools/patch_env_overrides.py` | lets both readers take their settings from environment variables |
| `tools/run_S4_nl.sh` | runs one case with the node-local reader (see RUNNING.md) |
| `tools/queue_S4.sh` | runs the five cases one after another (see RUNNING.md) |
| `tools/compare_spod_outputs.py` | compares the energies and modes of two runs |
| `tools/net_rates.py` | network MB/s per machine from a finished run |
| `tools/wait_analysis.py`, `tools/adios2_wait_shim.c` | wait analysis and its timing library (copies of the ones in `scripts/`) |
| `tools/rank_map.sh` | shows which MPI ranks run on which machine during a run |
| `design/` | design notes and the Python model used to check the mapping |
| `env.sh`, `setup.sh`, `origin.txt`, `hosts/coloc8`, `tools/status.sh`, `tools/watch_job.sh`, `tools/port_listen.sh`, `tools/serve_file.sh`, `tools/probe_nodes.sh`, `tools/run_S4.orig.sh` | from the original setup on the lab machines, kept for reference; the steps in these guides do not use them |

### Build

```bash
bash spod/node-local/build.sh
```

It first runs the mapping test, which must print
`node_local_map: 720 configurations PASS, 0 FAIL`. Then it builds
`spod/node-local/build/spod_old/spod_old` and `spod/node-local/build/spod_nl/spod_nl` with the
same compiler flags and libraries as the old reader. After changing a source file, rebuild one
reader with `bash spod/node-local/build.sh nl` (or `old`, or `map`).

### Settings (environment variables)

Both readers built by this script read these variables at start-up. A variable that is not set,
or is empty, keeps the value written in `UserInput`.

| Variable | Setting |
|---|---|
| `SPOD_NX`, `SPOD_NY`, `SPOD_NZ` | grid |
| `SPOD_STREAM` | stream name or path (the contact file is `<name>.sst`) |
| `SPOD_CONFIG` | ADIOS2 XML file for IO `spod-input`; leave unset for the built-in SST engine |
| `SPOD_NSNAP` | number of snapshots |
| `SPOD_NDFT` | snapshots per DFT block |
| `SPOD_K` | modes kept |
| `SPOD_OVERLAP` | overlap as a **fraction** of `SPOD_NDFT` (0.5), not a number of snapshots |
| `SPOD_NSAVE` | number of mode files written |
| `SPOD_WEIGHTS` | weight file; `none` = equal weights |
| `X3D_PROW`, `X3D_PCOL` | Xcompact3d's processor grid (node-local reader only) |
| `X3D_RANKS_PER_NODE` | 0 (default): Xcompact3d fills the same machines in the same order as SPOD |

The run scripts set all of these themselves.

## 4. Output files

Each SPOD run writes into the folder it runs in:

| File | Contents |
|---|---|
| `EnegySpectra.txt` | singular values per mode and frequency (PETSc ASCII matrix) |
| `SPODMODE_freq_<i>_sigma_<σ>.dat` | one mode: PETSc binary complex vector in the order `[ux; uy; uz; pp]`, x varying fastest; 1.75 GB each at 301³ |

To compare two runs (allowing for the arbitrary phase of each mode):

```bash
python3 spod/node-local/tools/compare_spod_outputs.py RUN_A RUN_B
```
