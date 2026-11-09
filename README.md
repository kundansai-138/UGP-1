# UGP-1: streaming SPOD of Xcompact3d over ADIOS2 SST

Xcompact3d computes a nozzle jet and, while it runs, sends every snapshot of the flow
(`ux, uy, uz, pp`) over an ADIOS2 **SST** stream straight into **streaming SPOD**
(PETSc/SLEPc). No snapshot is written to disk. This repository holds the solver build
used for that, both versions of the SPOD reader, the scripts that build and run
everything on the CSE lab machines (csews160–171), and the measured results.

## Guides

| Guide | What it covers |
|---|---|
| [docs/XCOMPACT3D_SETUP.md](docs/XCOMPACT3D_SETUP.md) | Toolchain, ADIOS2 2.10.2, the Xcompact3d ADIOS2 build and its patches, the test case and SST config |
| [docs/SPOD_SETUP.md](docs/SPOD_SETUP.md) | FFTW, PETSc/SLEPc (complex), the old ADIOS2 reader and the node-local reader |
| [docs/RUNNING.md](docs/RUNNING.md) | Running S0, S4 and the node-local cases, the run queue, correctness tests, monitoring, analysis |
| [results/README.md](results/README.md) | What each results folder holds and which runs to compare |

## Repository layout

```
xcompact3d/              Xcompact3d lab version, ADIOS2/SST build (source only, no binary)
spod/
  adios2-old/            streaming SPOD with the original ADIOS2 reader (the one S0/S4 used)
  node-local/            node-local ADIOS2 reader + build, test and run scripts (~/spod-dev)
scripts/                 environment, library build steps (step2-6.sh, setup-ugp.sh)
  ugp-tools/             S0/S4 run scripts, SST case setup, wait-timing library source, analysis
case/testcase_adios/     input.i3d and ADIOS2 configs of the jet case
results/                 S0 and S4 results, split by SPOD reader version
docs/                    the three guides
```

Not in the repository: binaries (`xcompact3d`, `streamingSPOD`, `spod_*`, `.so`), build files
(`.o`, `.mod`), Xcompact3d output data and SPOD mode files (1.75 GB each at 301³, over
GitHub's 100 MB limit). The plain MPI build of Xcompact3d used for S0 (`x3d_lab`) is also not
included.

## Main results (301³, 31 snapshots, one run per case)

| Case | Xcompact3d | SPOD | Step | SPOD read per snapshot | SPOD idle |
|---|---|---|---|---|---|
| S0, solver alone | 192 ranks | none | 41.2 s | – | – |
| S4 r4x4, separate | 192 ranks | 16, old reader | 35.3 s | 7.06 s | 94.6 % |
| Co-located | 160 ranks | 32, old reader | 43.1 s | 3.27 s | 96.8 % |
| Co-located | 160 ranks | 32, **node-local** | 42.8 s | **0.22 s** | 98.4 % |
| Separate | 160 ranks | 32, old reader | 43.3 s | 7.11 s | 95.7 % |
| Separate | 160 ranks | 32, **node-local** | 43.3 s | 6.98 s | 95.6 % |

- Xcompact3d is limited by the 1 GbE network: 93–95 % of its time is in `MPI_Alltoallv`,
  and the network cards reach their 125 MB/s ceiling during transposes.
- SPOD never limits the solver; it is idle 95–98 % of the time.
- The node-local reader makes SPOD's reads 15× faster when SPOD runs on the solver's nodes,
  but total time hardly changes because the solver's own traffic dominates the network.
- S0 was run on another day without APS's imbalance barrier, so it is not yet a fair
  baseline (see results/README.md).

## Paths

All scripts assume the layout used on the lab machines: libraries and builds under
`~/ugp` (`src`, `build`, `install`, `runs`, `tools`) and the node-local work under
`~/spod-dev`. Put the repository folders back in those places, or adjust the paths at the
top of each script.
