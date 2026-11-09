# Xcompact3d setup (ADIOS2 / SST build)

This builds the lab version of Xcompact3d (`xcompact3d/` in this repository) with ADIOS2, so
that it can stream snapshots over SST to streaming SPOD.

## 1. Machines and toolchain

| Item | Value on the lab machines |
|---|---|
| Machines | csews160–171 (and 172–199): i7-13700, 24 threads, 31 GB, Ubuntu 22.04 |
| Network | 1 GbE, interface `eno1` (also `docker0` and `virbr0`, which SST must not use) |
| Home | NFS, shared by all machines, 50 GB hard limit |
| Compilers | `gcc-11`, `g++-11`, `gfortran-11` |
| MPI | MPICH 4.3.1 in `/usr/local` (`mpicc`, `mpif90`, `mpirun`), device `ch4:ofi` |
| Other | `cmake`, `git`, `/usr/bin/python3` (not anaconda: its Python segfaults in PETSc configure) |

Check that the MPI wrappers use gcc-11:

```bash
mpicc -show     # must mention gcc-11
mpif90 -show    # must mention gfortran-11
```

## 2. Folder layout and environment

Everything lives under `~/ugp`:

```
~/ugp/src       sources (ADIOS2, PETSc, SLEPc, SPOD, x3d_lab_adios)
~/ugp/build     build folders
~/ugp/install   ADIOS2 and FFTW install prefix
~/ugp/runs      case folders and results
~/ugp/tools     run scripts (scripts/ugp-tools in this repository)
```

Copy `scripts/env.sh` to `~/ugp/env.sh` and load it from `~/.bashrc`:

```bash
echo '. $HOME/ugp/env.sh' >> ~/.bashrc
```

It puts `~/ugp/install/bin` on `PATH`, the ADIOS2, FFTW, PETSc and SLEPc libraries on
`LD_LIBRARY_PATH`, and sets `CC=gcc-11 FC=gfortran-11`, `PETSC_DIR`, `PETSC_ARCH=arch-complex-opt`
and `SLEPC_DIR`.

## 3. ADIOS2 2.10.2

`scripts/setup-ugp.sh` does this; the commands are:

```bash
cd ~/ugp/src
git clone --branch v2.10.2 --depth 1 https://github.com/ornladios/ADIOS2.git
cmake -S ADIOS2 -B ~/ugp/build/adios2 -DCMAKE_INSTALL_PREFIX=$HOME/ugp/install \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11 -DCMAKE_Fortran_COMPILER=gfortran-11 \
  -DADIOS2_USE_MPI=ON -DADIOS2_USE_Fortran=ON -DADIOS2_USE_SST=ON \
  -DADIOS2_USE_Python=OFF -DADIOS2_USE_HDF5=OFF -DADIOS2_USE_BZip2=OFF -DADIOS2_USE_Blosc2=OFF \
  -DADIOS2_USE_ZFP=OFF -DADIOS2_USE_SZ=OFF -DADIOS2_USE_MGARD=OFF -DADIOS2_USE_PNG=OFF \
  -DBUILD_TESTING=OFF
cmake --build ~/ugp/build/adios2 -j8 --target install
```

Check: `~/ugp/install/bin/adios2-config --version` and `~/ugp/install/bin/sst_conn_tool -h`.

## 4. Build Xcompact3d with ADIOS2

Put the `xcompact3d/` folder of this repository at `~/ugp/src/x3d_lab_adios`. Its sources already
contain four changes, made by `scripts/ugp-tools/setup_adios_case.sh` (the original files are
kept next to them):

| Change | File | Why |
|---|---|---|
| `int_to_str` fix | `src/visu.f90` (original: `visu.f90.orig`) | fixes the "End of record" crash |
| Snapshot-0 guard | `src/case.f90` (original: `case.f90.orig`) | snapshot at step 0 follows `ivisu`, as in the MPI build |
| BP5 and SST engine names | `decomp2d/io.f90` (original: `io.f90.pre_sst`) | the older 2decomp only knew BP4/HDF5; SST streams have no file extension |
| Synchronous snapshot writes | `src/visu.f90` (original: `visu.f90.pre_sync`) | deferred writes read the array only at `end_step` |

Build (serial `make`, about 5 minutes):

```bash
cd ~/ugp/src/x3d_lab_adios
make clean
make CMP=gcc FC=/usr/local/bin/mpif90 IO=adios2 ADIOS2DIR=$HOME/ugp/install > build.log 2>&1
```

`FFT` stays `generic` (2decomp's built-in FFT). Check the binary:

```bash
nm xcompact3d | grep -ci adios2_init      # must be > 0
ldd xcompact3d | grep -E 'adios|mpi|not found'
```

`setup_adios_case.sh` can redo all of this (patches, build, checks and the case folder):
`bash ~/ugp/tools/setup_adios_case.sh` (or `SKIP_BUILD=1` for the case folder only).

## 5. The test case

`case/testcase_adios/` holds the inputs; the scripts expect them in `~/ugp/runs/testcase_adios`.

| Setting | Value |
|---|---|
| Flow | `itype = 15`, nozzle jet with immersed boundary (`iibm = 3`) |
| Grid | `nx = ny = nz = 301`, box 20 × 20 × 20 |
| Physics | `re = 200`, `dt = 0.01`, RK3 (`itimescheme = 5`) |
| Length | `ilast = 150`, snapshot every `ioutput = 5` steps (31 snapshots with step 0) |
| Decomposition | `p_row = p_col = 0`: 2decomp chooses the grid |

### Processor grid

With `p_row = p_col = 0`, 2decomp's `best_2d_grid` lists the factors of the rank count in
ascending order and takes the middle one as `p_col`: 192 ranks give 12 × 16, 160 ranks give
10 × 16. Points that do not divide evenly go to the last ranks. Snapshots are written as
X-pencils: rank `r` holds all of x, the y-range of row `r / p_col` and the z-range of column
`r mod p_col`.

### ADIOS2 configurations

| File | Use |
|---|---|
| `adios2_config.sst.xml` | `solution-io` → SST stream `data`, everything else BP4 |
| `adios2_config.bp4.xml` | everything as BP4 files (not used: no snapshot files) |
| `adios2_config.xml` | the active one; the run scripts copy the SST config here |

SST settings for `solution-io`:

| Parameter | Value | Effect |
|---|---|---|
| `RendezvousReaderCount` | 1 | the solver waits at start until SPOD has connected |
| `QueueLimit` | 2 | at most 2 snapshots wait in memory |
| `QueueFullPolicy` | Block | when the queue is full, the solver's `end_step` waits for SPOD |
| `DataTransport` | WAN | TCP sockets over the normal network |
| `NetworkInterface` | eno1 | advertise and use `eno1`, not `docker0` |
| `OpenTimeoutSecs` | 3600 | a reader waits up to 1 h for the writer |

The stream's contact file `data.sst` appears in the run folder while the solver runs; SPOD
finds it through the shared home folder.
