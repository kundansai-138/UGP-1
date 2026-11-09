# Xcompact3d setup

This guide prepares the machines, builds ADIOS2, and builds the Xcompact3d version in
`xcompact3d/`, which can send its snapshots over an ADIOS2 SST stream.
All paths are relative to the cloned repository.

## 1. What you need

| Item | Notes |
|---|---|
| Linux machines | the results were measured on 24-thread machines with 31 GB of memory and 1 Gb/s Ethernet |
| Shared folder | every machine must see the repository at the same path (for example a shared home folder): the run scripts and the SST contact file rely on it |
| SSH | password-free SSH between all machines (the scripts log in to each machine) |
| Compilers | `gcc-11`, `g++-11`, `gfortran-11` by default; set `CC`, `CXX`, `FC` before loading `scripts/env.sh` to use others |
| MPI | MPICH with the Hydra launcher (`mpicc`, `mpif90`, `mpirun`); the results used MPICH 4.3.1 |
| Tools | `cmake` 3.12 or newer, `git`, `make`, `python3` (with `numpy` for the analysis tools) |

Check that the MPI compiler wrappers use the compilers you expect:

```bash
mpicc -show
mpif90 -show
```

## 2. Environment

Load the environment in every new shell, from the repository folder:

```bash
. scripts/env.sh
```

| Variable | Value |
|---|---|
| `UGP_ROOT` | the repository folder |
| `UGP_DEPS` | where the libraries are built and installed; default `deps/` inside the repository |
| `CC`, `CXX`, `FC` | compilers (default `gcc-11`, `g++-11`, `gfortran-11`); also passed to MPICH's wrappers |
| `PETSC_DIR`, `PETSC_ARCH` | `deps/src/petsc`, `arch-complex-opt` |
| `SLEPC_DIR` | `deps/src/slepc` |
| `PATH`, `LD_LIBRARY_PATH` | include `deps/install` and the PETSc and SLEPc libraries |

To use libraries that are already built somewhere else, set `UGP_DEPS` to that folder before
loading `env.sh`. It must have the same layout: `install/`, `src/petsc`, `src/slepc`.

## 3. ADIOS2 2.10.2

```bash
mkdir -p $UGP_DEPS/src && cd $UGP_DEPS/src
git clone --branch v2.10.2 --depth 1 https://github.com/ornladios/ADIOS2.git
cmake -S ADIOS2 -B $UGP_DEPS/build/adios2 -DCMAKE_INSTALL_PREFIX=$UGP_DEPS/install \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$FC \
  -DADIOS2_USE_MPI=ON -DADIOS2_USE_Fortran=ON -DADIOS2_USE_SST=ON \
  -DADIOS2_USE_Python=OFF -DADIOS2_USE_HDF5=OFF -DADIOS2_USE_BZip2=OFF -DADIOS2_USE_Blosc2=OFF \
  -DADIOS2_USE_ZFP=OFF -DADIOS2_USE_SZ=OFF -DADIOS2_USE_MGARD=OFF -DADIOS2_USE_PNG=OFF \
  -DBUILD_TESTING=OFF
cmake --build $UGP_DEPS/build/adios2 -j8 --target install
cd $UGP_ROOT
```

Check:

```bash
adios2-config --version
ls $UGP_DEPS/install/include/adios2/fortran/adios2.mod
```

## 4. Build Xcompact3d

```bash
cd xcompact3d
make CMP=gcc FC=mpif90 IO=adios2 ADIOS2DIR=$UGP_DEPS/install > build.log 2>&1
cd ..
```

This takes about 5 minutes (`make` without `-j`). The result is `xcompact3d/xcompact3d`.
`FFT` stays at its default `generic` (the FFT built into Xcompact3d's 2decomp library).

Check that ADIOS2 is linked in:

```bash
nm xcompact3d/xcompact3d | grep -ci adios2_init        # more than 0
ldd xcompact3d/xcompact3d | grep -E 'adios|mpi|not found'
```

### Changes already made in this source

Four changes were needed to stream snapshots. They were applied with
`scripts/original/setup_adios_case.sh` and are already in `xcompact3d/`; the original files
are kept next to the changed ones.

| Change | File (original kept as) | Why |
|---|---|---|
| `int_to_str` fix | `src/visu.f90` (`src/visu.f90.orig`) | fixes an "End of record" crash when writing snapshots |
| Start-up snapshot follows `ivisu` | `src/case.f90` (`src/case.f90.orig`) | `ivisu = 0` then gives no snapshots at all |
| BP5 and SST engine names | `decomp2d/io.f90` (`decomp2d/io.f90.pre_sst`) | the 2decomp I/O only knew BP4 and HDF5; an SST stream has no file extension |
| Snapshots written at once | `src/visu.f90` (`src/visu.f90.pre_sync`) | deferred writes read the array only at `end_step` |

## 5. The test case

`case/testcase_adios/` holds the inputs. The run scripts copy them into a new run folder for
every run.

| File | Contents |
|---|---|
| `input.i3d` | Xcompact3d input |
| `probes_points.prm` | probe positions |
| `adios2_config.sst.xml` | snapshots (`solution-io`) as SST stream `data`; everything else as BP4 files |
| `adios2_config.bp4.xml` | everything as BP4 files |
| `adios2_config.xml` | the configuration Xcompact3d reads; the run scripts copy the SST one here |
| `basecase.sh`, `hostfile`, `input.i3d.master`, `input.i3d.before_sst_test` | kept from the original case folder; the run scripts do not use them |

Main settings in `input.i3d`:

| Setting | Value |
|---|---|
| Flow | `itype = 15`: nozzle jet with immersed boundary (`iibm = 3`) |
| Grid | `nx = ny = nz = 301`, box 20 × 20 × 20 |
| Physics | `re = 200`, `dt = 0.01`, RK3 time stepping (`itimescheme = 5`) |
| Length | `ilast = 150` steps, a snapshot every `ioutput = 5` steps: 31 snapshots including step 0 |
| Processor grid | `p_row = p_col = 0`: chosen automatically (below) |

### Processor grid

With `p_row = p_col = 0`, 2decomp lists the factors of the number of ranks in ascending order
and takes the middle one as `p_col`, with `p_row = ranks / p_col`. 192 ranks give 12 × 16 and
160 ranks give 10 × 16. Grid points that do not divide evenly go to the last ranks.

Snapshots are written as X-pencils: rank `r` holds all of x, the y-range of row `r / p_col` and
the z-range of column `r mod p_col`.

### SST settings (`adios2_config.sst.xml`, IO `solution-io`)

| Parameter | Value | Effect |
|---|---|---|
| `RendezvousReaderCount` | 1 | Xcompact3d waits at start until SPOD has connected |
| `QueueLimit` | 2 | at most 2 snapshots wait in memory |
| `QueueFullPolicy` | Block | when the queue is full, Xcompact3d waits for SPOD |
| `DataTransport` | WAN | TCP sockets over the normal network |
| `NetworkInterface` | eno1 | network interface to use; the run scripts replace it with their `IFACE` setting |
| `OpenTimeoutSecs` | 3600 | a reader waits up to an hour for Xcompact3d to start |

While Xcompact3d runs, the contact file `data.sst` exists in the run folder. SPOD reads it
through the shared folder to find Xcompact3d.
