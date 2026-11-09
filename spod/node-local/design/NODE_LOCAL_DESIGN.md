# Node-local ADIOS2 input for streaming SPOD: design

The node-local reader (`spod/node-local/src/spod_nl.F90`) reads Xcompact3d's SST stream so that
**each SPOD rank only reads data written by Xcompact3d ranks on its own machine**. The PETSc/SLEPc
SPOD algorithm is unchanged.

| File | What it is |
|---|---|
| `../src/spod_nl.F90` | module `node_local_map` (which boxes each rank reads), the SPOD module with the new input layer, and the program |
| `../src/test_node_local_map.f90` | stand-alone test of the mapping (run by `../build.sh`) |
| `node_local_model.py` | Python model: mapping, SPOD equivalence, traffic numbers |
| `../tools/compare_spod_outputs.py` | compares two SPOD runs (energies and modes) |

## The problem with the old reader

Setup: `n` machines. Xcompact3d runs `k = p_row x p_col` ranks placed machine by machine
(ranks `0..k/n-1` on machine 0, ...). SPOD runs `l` ranks on the same machines.

* Xcompact3d (2decomp) gives rank `r` an X-pencil: all of `x`, the `y` range of row
  `r / p_col`, the `z` range of column `mod(r, p_col)`. Through ADIOS2 each rank writes
  exactly that box.
* The old SPOD reader keeps PETSc's default layout: rank `s` owns a contiguous range of the
  snapshot vector `q = [ux; uy; uz; pp]`. In 3-D that is a slab of whole `x-y` planes. A whole
  `x-y` plane is made of the pencils of `p_row` Xcompact3d ranks `c, c + p_col, c + 2 p_col, ...`,
  which are spread over every machine.

So almost every byte of every snapshot crosses the network: a fraction `1 - 1/n` of it.

```
 machine 0                      machine 1
 ┌─────────────────────┐        ┌─────────────────────┐
 │ X3D 0..k/2-1        │        │ X3D k/2..k-1        │
 │   writes pencils ───┼───┐ ┌──┼── writes pencils    │
 │                     │   ╳    │                     │   old reader:
 │ SPOD 0..l/2-1  <────┼───┘ └──┼──> SPOD l/2..l-1    │   each SPOD rank needs
 └─────────────────────┘        └─────────────────────┘   planes from all machines

 machine 0                      machine 1
 ┌─────────────────────┐        ┌─────────────────────┐
 │ X3D 0..k/2-1        │        │ X3D k/2..k-1        │   node-local reader:
 │      │              │        │      │              │   each SPOD rank reads
 │      v              │        │      v              │   only pencils written
 │ SPOD 0..l/2-1       │        │ SPOD l/2..l-1       │   on its own machine
 └─────────────────────┘        └─────────────────────┘
```

## Why reordering rows keeps the results the same

Every step of streaming SPOD is either

* **row by row**: running mean, windowed DFT accumulation, mean removal, scaling,
  `x = X_hat .* sqrt(W)`, `u_p = x - U (U^H x)`, `U_hat = [U u_new] * Umat`,
  `U_hat ./ sqrt(W)`; or
* **a sum over all rows**: `U^H x` and `||u_p||`, which feed the small `(k+1)x(k+1)`
  matrix that SLEPc decomposes.

Neither depends on the order of the rows, as long as the data, the weights and the modes all use
the same order. So the code can use any ordering of the grid points, and this version picks the
one that keeps the reads on each machine:

1. **Node-local ordering.** On each machine, take the Xcompact3d ranks of that machine in rank
   order and concatenate their pencils, line by line (a line is one run along `x`). Split that
   list evenly over the SPOD ranks of the machine. A rank's share is cut into at most three boxes
   per pencil: part of a first plane, whole planes, part of a last plane. Each box lies inside one
   pencil, so ADIOS2 serves it from that single writer, on the same machine.
2. **PETSc layout follows the data.** Every matrix or vector with one row per grid point
   (`x_new`, `mu`, `X_sum_BlockA/B`, `U_hat`, `U`, `Utemp`, `UtempMult`, `sqrtW`, `x`, `u_p`,
   `u_new`, `temp`, `X_hat_column`) is created with local row count `nLocalRows` instead of
   `PETSC_DECIDE`. Nothing else in the PETSc code changes. The small objects (`S_hat`, `Kmode`,
   `Umat`, `Ux`, `FourierRange`) keep `PETSC_DECIDE`.
3. **Natural order at the edges.** A `VecScatter` built from the natural index of every local row
   converts between the two orders:
   * weights: the weight file stays in natural order and is scattered once into `sqrtW`;
   * modes: before `VecView`, each saved mode is scattered back to natural order, so
     `SPODMODE_*.dat` has exactly the old layout (`[ux;uy;uz;pp]`, `i` fastest).

   Both happen once per run, not per snapshot.

At start-up the code checks that the reordering is an exact permutation (every natural row
receives exactly one local row; the run stops otherwise) and prints one line per machine with its
name, SPOD ranks and the Xcompact3d ranks it reads.

### What still crosses machines

* The snapshot data: **none**, if the placement assumption holds.
* SST metadata (block sizes per step), a few kB per step.
* The SPOD maths itself, as before: per block, `Allreduce`s of `k` numbers (`U^H x`, norms),
  broadcasts of one DFT row, the `(k+1)x(k+1)` SVD, `MatMatMult` with the small `Umat`. These are
  limited by latency, not bandwidth.

## Assumptions

1. Xcompact3d writes 3-D fields as X-pencils (`output2D = 0`, `nvisu = 1`).
2. The reader knows Xcompact3d's `p_row x p_col` (`X3D_PROW`, `X3D_PCOL`). The run scripts compute
   it with 2decomp's own rule and stop if Xcompact3d's log shows a different grid.
3. Xcompact3d ranks fill the same machines, in the same order, as the SPOD ranks
   (`X3D_RANKS_PER_NODE = 0`). If Xcompact3d uses a fixed number of ranks per machine, that
   number can be given in `X3D_RANKS_PER_NODE` instead.

If 2 or 3 is wrong, **the results stay correct**: the boxes still cover the grid exactly once and
ADIOS2 assembles them from whichever ranks hold them. Only the locality is lost.

SPOD ranks themselves can be placed in any way (block or round-robin); each rank finds its machine
with `MPI_Comm_split_type`.

Build and run: [docs/SPOD_SETUP.md](../../../docs/SPOD_SETUP.md) and [docs/RUNNING.md](../../../docs/RUNNING.md).

## How the correctness was checked

* `node_local_model.py` (Python, mirrors the Fortran routine by routine):
  * mapping: 724 configurations (grids, `p_row x p_col`, machines, SPOD ranks, all three pencils,
    both 2decomp split rules, block and round-robin SPOD placement, fixed `X3D_RANKS_PER_NODE`,
    SPOD ranks that end up with 0 rows). Every point is read exactly once, every box is inside one
    writer's pencil, and that writer is on the reader's machine;
  * maths: streaming SPOD run with natural row order and with node-local row order (random data,
    random weights). Energies agree to about 1e-16 and modes to about 1e-15. A mode can come out
    multiplied by `e^{i theta}` (often `-1`), the usual SVD phase freedom.
* `../src/test_node_local_map.f90`: the Fortran mapping on 720 configurations (`../build.sh map`).
* At run time: the permutation check described above.
* On real data: `../run_compare.sh` streams one Xcompact3d run to the old reader and to the
  node-local reader at the same time and compares the results (identical energies, modes equal to
  about 1e-14).

## Expected traffic (from `node_local_model.py`)

| Case | Old reader, across machines | Node-local, across machines |
|---|---|---|
| 301³ x 4 variables, X3D 8x16, 4 machines, SPOD 64 | 657 MiB per snapshot (75 %) | 0 |
| 301³ x 4 variables, X3D 16x16, 8 machines, SPOD 128 | 805 MiB per snapshot (87.5 %) | 0 |
| 256x128x256 x 4, X3D 4x16, 2 machines, SPOD 32 | 128 MiB per snapshot (50 %) | 0 |

The old reader also reads up to two extra partial planes per rank and variable (876 MiB instead
of 832 MiB per snapshot in the first case). Each node-local SPOD rank reads from at most 3
Xcompact3d ranks, instead of up to 32 with the old reader.
