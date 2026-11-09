# Streaming SPOD with node-local ADIOS2 input

This folder is the ADIOS2/SST version of streaming SPOD (`spod-adios2.txt`), with the
input rewritten so that **each SPOD rank only reads data written by Xcompact3d ranks
on its own node**. The PETSc/SLEPc SPOD algorithm is unchanged.

| File | What it is |
|---|---|
| `streamingSPODModule.F90` | SPOD program. New input layer, same SPOD maths. |
| `node_local_map.F90` | Decides which boxes each SPOD rank reads (plain Fortran, no MPI). |
| `test_node_local_map.f90` | Stand-alone test of the mapping (`make test_node_local_map`). |
| `Makefile` | Builds both. Set `FFTW_DIR` and `ADIOS2_DIR`. |
| `../tools/node_local_model.py` | Python model: mapping, SPOD equivalence, traffic numbers. |
| `../tools/compare_spod_outputs.py` | Compares two SPOD runs (energies + modes). |

## The problem with the current reader

Setup: `n` nodes. Xcompact3d (X3D) runs `k = p_row x p_col` ranks placed sequentially
(ranks `0..k/n-1` on node 0, ...). SPOD runs `l` ranks on the same nodes.

* X3D (2decomp-fft) gives rank `r` an X-pencil: all of `x`, the `y` range of row
  `r / p_col`, the `z` range of column `mod(r, p_col)`. Through ADIOS2 each rank
  writes exactly that box.
* The current SPOD reader keeps PETSc's default layout: rank `s` owns a contiguous
  range of the snapshot vector `q = [ux; uy; uz; pp]`. In 3-D that is a slab of whole
  `x-y` planes. A whole `x-y` plane is made of the pencils of `p_row` X3D ranks
  `c, c + p_col, c + 2 p_col, ...`, which are spread over every node.

So almost every byte of every snapshot crosses the network. With the numbers below
the fraction is `1 - 1/n`.

```
 node 0                         node 1
 ┌─────────────────────┐        ┌─────────────────────┐
 │ X3D 0..k/2-1        │        │ X3D k/2..k-1        │
 │   writes pencils ───┼───┐ ┌──┼── writes pencils    │
 │                     │   ╳    │                     │   current reader:
 │ SPOD 0..l/2-1  <────┼───┘ └──┼──> SPOD l/2..l-1    │   each SPOD rank needs
 └─────────────────────┘        └─────────────────────┘   planes from all nodes

 node 0                         node 1
 ┌─────────────────────┐        ┌─────────────────────┐
 │ X3D 0..k/2-1        │        │ X3D k/2..k-1        │   node-local reader:
 │      │              │        │      │              │   each SPOD rank reads
 │      v              │        │      v              │   only pencils written
 │ SPOD 0..l/2-1       │        │ SPOD l/2..l-1       │   on its own node
 └─────────────────────┘        └─────────────────────┘
```

## Why reordering rows keeps the results the same

Every step of streaming SPOD is either

* **row by row**: running mean, windowed DFT accumulation, mean removal, scaling,
  `x = X_hat .* sqrt(W)`, `u_p = x - U (U^H x)`, `U_hat = [U u_new] * Umat`,
  `U_hat ./ sqrt(W)`; or
* **a sum over all rows**: `U^H x` and `||u_p||`, which feed the small `(k+1)x(k+1)`
  matrix that SLEPc decomposes.

Neither depends on the order of the rows, as long as the data, the weights and the
modes all use the same order. So the code can use any ordering of the grid points,
and this version picks the one that keeps the reads on-node:

1. **Node-local ordering.** On each node, take the X3D ranks of that node in rank
   order and concatenate their pencils, line by line (a line is one run along `x`).
   Split that list evenly over the SPOD ranks of the node. A rank's share is cut into
   at most three boxes per pencil: part of a first plane, whole planes, part of a
   last plane. Each box lies inside one pencil, so ADIOS2 serves it from that single
   writer, on the same node.
2. **PETSc layout follows the data.** Every matrix or vector with one row per grid
   point (`x_new`, `mu`, `X_sum_BlockA/B`, `U_hat`, `U`, `Utemp`, `UtempMult`,
   `sqrtW`, `x`, `u_p`, `u_new`, `temp`, `X_hat_column`) is created with local row
   count `nLocalRows` instead of `PETSC_DECIDE`. Nothing else in the PETSc code
   changes. The small objects (`S_hat`, `Kmode`, `Umat`, `Ux`, `FourierRange`) keep
   `PETSC_DECIDE`.
3. **Natural order at the edges.** A `VecScatter` built from the natural index of
   every local row converts between the two orders:
   * weights: the weight file stays in natural order and is scattered once into
     `sqrtW`;
   * modes: before `VecView`, each saved mode is scattered back to natural order,
     so `SPODMODE_*.dat` has exactly the old layout (`[ux;uy;uz;pp]`, `i` fastest)
     and the MATLAB reader is unchanged.

   Both happen once per run, not per snapshot.

At start-up the code checks that the reordering is an exact permutation (every
natural row receives exactly one local row; the run aborts otherwise) and prints one
line per node with its host, SPOD ranks and the X3D ranks it reads.

### What still crosses nodes

* The snapshot data: **none**, if the placement assumption holds.
* SST metadata (block sizes per step), a few kB per step.
* The SPOD maths itself, as before: per block, `Allreduce`s of `k` numbers (`U^H x`,
  norms), broadcasts of one DFT row, the `(k+1)x(k+1)` SVD, `MatMatMult` with the
  small `Umat`. These are latency-bound, not bandwidth-bound.

## Assumptions

1. Xcompact3d writes 3-D fields: `output2D = 0`, `nvisu = 1` (full grid, X-pencils).
2. `p_row`, `p_col` are set explicitly in X3D `input.i3d` (not `0`) and copied into
   `x3dProw`, `x3dPcol` in `UserInput`.
3. X3D ranks are placed sequentially on the same nodes, in the same node order as
   the SPOD ranks (`x3dRanksPerNode = 0`). If X3D uses a fixed number of ranks per
   node, you can set `x3dRanksPerNode` to that number instead.

If 2 or 3 is wrong, **the results stay correct**. The boxes still cover the grid
exactly once and ADIOS2 assembles them from whichever ranks hold them. Only the
locality is lost.

SPOD ranks themselves can be placed in any way (block or round-robin); each rank
finds its node with `MPI_Comm_split_type`.

## Build

```bash
make                         # streamingSPOD
make test_node_local_map     # mapping test, only needs gfortran
./test_node_local_map        # expect: "... configurations PASS, 0 FAIL"
```

## Run (example: 4 nodes, X3D 8 x 16 = 128 ranks, SPOD 64 ranks)

The two launches must use the same nodes in the same order, with ranks filled node by
node.

MPICH (hydra), same host file for both:

```bash
mpiexec -f hosts -ppn 32 -n 128 ./xcompact3d &
mpiexec -f hosts -ppn 16 -n 64  ./streamingSPOD
```

Slurm, inside one allocation of 4 nodes (`-m block` keeps rank order node by node;
give each step its own cores):

```bash
srun -N 4 -n 128 -m block ./xcompact3d &
srun -N 4 -n 64  -m block ./streamingSPOD
```

Check the start-up table: the X3D ranks listed for a node must be the ones running on
that host.

## How the correctness was checked

* `tools/node_local_model.py` (Python, mirrors the Fortran routine by routine):
  * mapping: 724 configurations (grids, `p_row x p_col`, nodes, SPOD ranks, all
    three pencils, both 2decomp split rules, block and round-robin SPOD placement,
    explicit `x3dRanksPerNode`, SPOD ranks that end up with 0 rows).
    Every point is read exactly once, every box is inside one writer's pencil, and
    that writer is on the reader's node;
  * maths: streaming SPOD run with natural row order and with node-local row order
    (random data, random weights). Energies agree to ~1e-16 and modes to ~1e-15.
    A mode can come out multiplied by `e^{i theta}` (often `-1`). That is the usual
    SVD phase freedom, and the old code shows it too when you change the number of
    ranks.
* At run time: the permutation check described above.
* Recommended A/B test on the cluster: let X3D write a BP5 file instead of SST, run
  the old reader (`spod-adios2.txt`) and this version on that file with the BP5
  engine, then

  ```bash
  python tools/compare_spod_outputs.py run_old run_new
  ```

  It compares `EnegySpectra.txt` and every `SPODMODE_*.dat`, allowing for the
  phase freedom.

## Expected traffic (from `tools/node_local_model.py`)

| Case | Current reader, cross-node | Node-local, cross-node |
|---|---|---|
| 301³ x 4 vars, X3D 8x16, 4 nodes, SPOD 64 | 657 MiB / step (75%) | 0 |
| 301³ x 4 vars, X3D 16x16, 8 nodes, SPOD 128 | 805 MiB / step (87.5%) | 0 |
| 256x128x256 x 4, X3D 4x16, 2 nodes, SPOD 32 | 128 MiB / step (50%) | 0 |

The current reader also reads up to two extra partial planes per rank and variable
(876 MiB instead of 832 MiB per step in the first case). Each new SPOD rank talks to at
most 3 X3D ranks, instead of up to 32 before.
