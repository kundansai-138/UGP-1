# Results

All runs use the jet case of `case/testcase_adios`: 301 × 301 × 301 grid, 150 steps, a snapshot
every 5 steps (31 snapshots), SPOD with nDFT 8, overlap 4 and k 3. Machines: 24 threads, 31 GB,
1 Gb/s Ethernet; Xcompact3d on 8 machines, SPOD on the same 8 or on 4 other machines.

How to read a run folder: [docs/RUNNING.md](../docs/RUNNING.md), section 6. Mode files are not
included (1.75 GB each); runs that kept their energy spectrum have it in `run1/spod/EnegySpectra.txt`.

## Folders

`x3d192` = 192 Xcompact3d ranks (24 per machine), `x3d160` = 160 (20 per machine); `spodN` = N SPOD ranks.

| Folder | What ran | Runs | Note |
|---|---|---|---|
| `S0_solver_only/` | Xcompact3d alone, 192 ranks, no snapshots | 1 (run 2 failed at start) | baseline, see the note below |
| `spod-old-reader/separate_x3d192_spod48_A/` | old reader, 48 SPOD ranks on separate machines | 2 | first set, earlier version of the run script |
| `spod-old-reader/separate_x3d192_spod48_B/` | old reader, 48 SPOD ranks on separate machines | 2 | |
| `spod-old-reader/separate_x3d192_spod32/` | old reader, 32 SPOD ranks on separate machines | 2 | **run 2 overlapped another large run: its times are inflated** |
| `spod-old-reader/separate_x3d192_spod16/` | old reader, 16 SPOD ranks on separate machines | 1 | |
| `spod-old-reader/colocated_x3d160_spod32/` | old reader, 32 SPOD ranks on the same machines | 1 | |
| `spod-old-reader/separate_x3d160_spod32/` | old reader, 32 SPOD ranks on separate machines | 1 | |
| `spod-old-reader/spod_rank_sweep_console.log` | console output of the 48 / 32 / 16 SPOD-rank runs | | |
| `spod-node-local/colocated_x3d160_spod32/` | node-local reader, 32 SPOD ranks on the same machines | 1 | |
| `spod-node-local/colocated_x3d160_spod32_overlapped/` | node-local reader, same layout, on another set of 8 machines | 1 | **ran while another large run used the network: times inflated** |
| `spod-node-local/separate_x3d160_spod32/` | node-local reader, 32 SPOD ranks on separate machines | 1 | |
| `spod-node-local/quick65_colocated/` | node-local reader, 65³ grid, 15 steps, same machines | 1 | reports only |
| `spod-node-local/compare_301_nodelocal_vs_old.txt` | modes of `colocated_x3d160_spod32_overlapped` vs `separate_x3d192_spod32` run 2 | | see Correctness |
| `queue_S4.log` | log of the queue that ran the `x3d160` cases and `separate_x3d192_spod16` | | |

**S0 note.** `S0_solver_only` was produced with a plain MPI build of the same Xcompact3d, which is not
part of this repository, and its APS profiling ran without the imbalance barrier
(`APS_IMBALANCE_TYPE=2`) that all other runs used. Its time per step (41.2 s) is longer than the
streaming runs with the same 192 ranks (35 s), so it is not a fair baseline yet.

## Which runs to compare

Each pair ran one after the other on the same machines, with the same Xcompact3d grid (10 × 16):

| Question | Old reader | Node-local reader |
|---|---|---|
| SPOD on the same machines | `spod-old-reader/colocated_x3d160_spod32` | `spod-node-local/colocated_x3d160_spod32` |
| SPOD on separate machines | `spod-old-reader/separate_x3d160_spod32` | `spod-node-local/separate_x3d160_spod32` |

## Key numbers

| Run | Time per step | Extra time of a snapshot step | SPOD read per snapshot | SPOD idle | Total |
|---|---|---|---|---|---|
| S0, solver alone (192) | 41.2 s | – | – | – | 6,153 s |
| old, separate, 192 + 48 (B) | 35.1 s | 4.4 s | 6.6–6.8 s | 94.7 % | 5,516 s |
| old, separate, 192 + 16 | 35.3 s | 4.4 s | 7.06 s | 94.6 % | 5,530 s |
| old, same machines, 160 + 32 | 43.1 s | 5.7 s | 3.27 s | 96.8 % | 6,668 s |
| node-local, same machines, 160 + 32 | 42.8 s | 5.8 s | 0.22 s | 98.4 % | 6,655 s |
| old, separate, 160 + 32 | 43.3 s | 5.7 s | 7.11 s | 95.7 % | 6,719 s |
| node-local, separate, 160 + 32 | 43.3 s | 5.8 s | 6.98 s | 95.6 % | 6,694 s |

- Xcompact3d sets the pace: 192 ranks take about 35 s per step, 160 ranks about 43 s, whatever
  SPOD does.
- SPOD is idle 95–98 % of the time; the 2-snapshot SST queue absorbed every busy moment.
- With SPOD on the same machines, the node-local reader reads a snapshot in 0.22 s instead of
  3.27 s. With SPOD on separate machines every read crosses the network for both readers, and
  they take about the same time.

## Correctness

- **65³, one stream, three readers** (`spod/node-local/run_compare.sh`): the old reader and the
  node-local reader on 4 and 5 ranks give identical energies, and modes that agree to about 1e-14.
- **301³, two different runs** (`compare_301_nodelocal_vs_old.txt`): energies agree to 3e-7 (the
  precision of the file), modes point the same way (|cos| = 0.9999996), with local differences up
  to 3 % of the peak value. The two runs used different Xcompact3d grids (10 × 16 and 12 × 16), so
  their snapshots differ slightly; this comparison cannot separate reader effects from that.

## Why the network is the limit

- APS: 98.6–98.8 % of Xcompact3d's time is inside MPI, 93–95 % inside `MPI_Alltoallv`
  (in S0 and in every streaming run).
- `run1/node_<machine>.csv`: every Xcompact3d machine sends and receives 320–375 GB per run;
  the rate reaches 120–123.5 MB/s, the limit of 1 Gb/s Ethernet (125 MB/s)
  (`spod/node-local/tools/net_rates.py`).
- With 160 ranks, four machines carry about 372 GB and the other four 325 GB; a step takes 43 s
  instead of 35 s with 192 ranks, where every machine carries about 328 GB.
- When another large run used the network at the same time (on other machines),
  `separate_x3d192_spod32` run 2 slowed from 35.4 to 53.8 s per step.
