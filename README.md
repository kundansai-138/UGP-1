# UGP-1: streaming SPOD of Xcompact3d over ADIOS2 SST

Xcompact3d simulates a nozzle jet and, while it runs, sends every snapshot of the flow
(`ux, uy, uz, pp`) over an ADIOS2 **SST** stream straight into **streaming SPOD**
(PETSc/SLEPc). No snapshot is written to disk.

This repository contains the Xcompact3d version that can stream, two versions of the SPOD
reader, the scripts that build and run everything, and the measured results.

## Guides

Read them in this order:

| Guide | What it covers |
|---|---|
| [docs/XCOMPACT3D_SETUP.md](docs/XCOMPACT3D_SETUP.md) | Requirements, environment, ADIOS2, building Xcompact3d, the test case |
| [docs/SPOD_SETUP.md](docs/SPOD_SETUP.md) | FFTW, PETSc, SLEPc, building both SPOD readers |
| [docs/RUNNING.md](docs/RUNNING.md) | Correctness test, running the cases, the run queue, stopping, analysis |
| [results/README.md](results/README.md) | What each results folder holds and which runs to compare |

## Repository layout

```
xcompact3d/            Xcompact3d with ADIOS2/SST support (source)
spod/
  adios2-old/          streaming SPOD with the original ADIOS2 reader (source + build.sh)
  node-local/          streaming SPOD with the node-local ADIOS2 reader (source, build, test and run scripts)
scripts/
  env.sh               environment for every build and run (load it first)
  case_s4.sh           runs one case: Xcompact3d -> SST -> SPOD (old reader), with all measurements
  hostfile             machines for Xcompact3d (edit for your machines)
  hostfile_spod        machines for SPOD when it runs on separate machines (edit for your machines)
  wait_analysis.py     reader/writer wait analysis
  adios2_wait_shim.c   source of the timing library loaded into both programs
  original/            the scripts as first used on the lab machines, kept for reference only
case/testcase_adios/   inputs of the jet case and its ADIOS2 configurations
results/               measured results, split by SPOD reader version
docs/                  the guides
```

Two folders are created inside the repository when you build and run. Git ignores them:

```
deps/                  libraries you build: ADIOS2, FFTW, PETSc, SLEPc (see scripts/env.sh)
runs/                  results of your own runs
```

## Quick start

```bash
git clone https://github.com/kundansai-138/UGP-1.git
cd UGP-1
. scripts/env.sh                         # in every new shell
```

1. Build the libraries into `deps/`: ADIOS2 ([XCOMPACT3D_SETUP.md](docs/XCOMPACT3D_SETUP.md) section 3),
   then FFTW, PETSc and SLEPc ([SPOD_SETUP.md](docs/SPOD_SETUP.md) section 1).
2. Build the programs:

   ```bash
   (cd xcompact3d && make CMP=gcc FC=mpif90 IO=adios2 ADIOS2DIR=$UGP_DEPS/install > build.log 2>&1)
   bash spod/adios2-old/build.sh
   bash spod/node-local/build.sh
   ```

3. Check that everything works on one machine:

   ```bash
   bash spod/node-local/run_compare.sh      # ends with: RESULT: PASS
   ```

4. Run the cases: [docs/RUNNING.md](docs/RUNNING.md).

## Main results

301 × 301 × 301 grid, 31 snapshots, one run per case, 1 Gb/s Ethernet between machines.
All numbers: [results/README.md](results/README.md).

| Case | Xcompact3d ranks | SPOD ranks and reader | Time per step | SPOD read per snapshot | SPOD idle |
|---|---|---|---|---|---|
| Separate machines | 192 | 16, old reader | 35.3 s | 7.06 s | 94.6 % |
| Same machines | 160 | 32, old reader | 43.1 s | 3.27 s | 96.8 % |
| Same machines | 160 | 32, **node-local** | 42.8 s | **0.22 s** | 98.4 % |
| Separate machines | 160 | 32, old reader | 43.3 s | 7.11 s | 95.7 % |
| Separate machines | 160 | 32, **node-local** | 43.3 s | 6.98 s | 95.6 % |

- Xcompact3d is limited by the network: 93–95 % of its time is spent in `MPI_Alltoallv`,
  and the network cards reach their limit of 125 MB/s while it exchanges data.
- SPOD never slows the solver down: it waits for the next snapshot 95–98 % of the time.
- When SPOD runs on the same machines as Xcompact3d, the node-local reader makes SPOD's
  reads 15 times faster. Total time hardly changes, because the solver's own network
  traffic is far larger than the snapshots.

## Not in this repository

- Binaries and build files: build them as described in the guides.
- SPOD mode files: 1.75 GB each at 301³, above GitHub's 100 MB file limit.
- Restart files of the Xcompact3d examples (`checkpoint`, 870 MB).
