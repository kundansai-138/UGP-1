# Results

301³ nozzle jet, 150 steps, a snapshot every 5 steps (31 snapshots), SPOD nDFT 8, overlap 4,
k 3, on csews160–171 (i7-13700, 1 GbE). How to read a run folder: [docs/RUNNING.md](../docs/RUNNING.md)
section 6. SPOD mode files are not included (1.75 GB each); each run keeps `run1/spod/EnegySpectra.txt`
where it was saved.

## Folders

| Folder | Run | Date | Xcompact3d | SPOD | Notes |
|---|---|---|---|---|---|
| `S0/S0_20261104_1347` | S0, solver alone | 4 Nov | 192 (24/machine) | none | run 2 failed at start; APS without imbalance barrier |
| `spod-old-reader/S4/S4_20261104_1953` | S4 | 4 Nov | 192 | 48, separate | 2 runs |
| `spod-old-reader/S4/S4_r4x12_20261108_1809` | S4 r4x12 | 8 Nov | 192 | 48, separate | 2 runs |
| `spod-old-reader/S4/S4_r4x8_20261108_2120` | S4 r4x8 | 8 Nov | 192 | 32, separate | **run 2 overlapped another full-size run: its times are inflated** |
| `spod-old-reader/S4/S4_r4x4_20261109_0144` | S4 r4x4 | 9 Nov | 192 | 16, separate | |
| `spod-old-reader/S4_coloc/S4_r8x4_20261109_0320` | co-located | 9 Nov | 160 (20/machine) | 32 (4/machine), same machines | |
| `spod-old-reader/S4_sep160/S4_r4x8_20261109_0710` | separate | 9 Nov | 160 | 32 (8/machine) | |
| `spod-node-local/S4_coloc_nl/S4_r8x4_20261109_0516` | co-located | 9 Nov | 160 | 32, same machines | |
| `spod-node-local/S4_coloc_nl/S4_r8x4_20261108_2253` | co-located | 8 Nov | 160 on csews181–187, 190 | 32 | **overlapped the S4 r4x8 run 2: times inflated** |
| `spod-node-local/S4_sep_nl/S4_r4x8_20261109_0905` | separate | 9 Nov | 160 | 32 (8/machine) | |
| `spod-node-local/quick65_nl` | 65³ check, 15 steps | 8 Nov | 160 on csews181–187, 190 | 32 | reports only |
| `spod-node-local/compare_301_nl_vs_r4x8run2.txt` | | | | | modes of two different 301³ runs (different grids), see below |
| `queue_S4.log` | | 9 Nov | | | log of the 9 November queue |

## Which runs to compare

| Question | Old reader | Node-local reader |
|---|---|---|
| SPOD on the solver's machines | `S4_coloc/S4_r8x4_20261109_0320` | `S4_coloc_nl/S4_r8x4_20261109_0516` |
| SPOD on separate machines | `S4_sep160/S4_r4x8_20261109_0710` | `S4_sep_nl/S4_r4x8_20261109_0905` |

Each pair ran back to back on the same machines with the same Xcompact3d grid (10 × 16).

## Key numbers (one run each)

| Case | Step | Extra per snapshot step | SPOD read per snapshot | SPOD idle | Total |
|---|---|---|---|---|---|
| S0, solver alone | 41.2 s | – | – | – | 6,153 s |
| r4x12 (192 + 48) | 35.1 s | 4.4 s | 6.6–6.8 s | 94.7 % | 5,516 s |
| r4x4 (192 + 16) | 35.3 s | 4.4 s | 7.06 s | 94.6 % | 5,530 s |
| co-located, old | 43.1 s | 5.7 s | 3.27 s | 96.8 % | 6,668 s |
| co-located, node-local | 42.8 s | 5.8 s | 0.22 s | 98.4 % | 6,655 s |
| separate, old | 43.3 s | 5.7 s | 7.11 s | 95.7 % | 6,719 s |
| separate, node-local | 43.3 s | 5.8 s | 6.98 s | 95.6 % | 6,694 s |

## Correctness

- **65³, one stream, three readers** (8 Nov): old reader vs node-local reader on 4 and 5 ranks gave
  identical energies and modes equal to about 1e-14.
- **301³, two different runs** (`compare_301_nl_vs_r4x8run2.txt`): energies agree to 3e-7 (the
  printed precision), modes point the same way (|cos| = 0.9999996), local differences up to 3 % of
  the peak. The two runs used different Xcompact3d grids (10 × 16 and 12 × 16), so their snapshots
  differ slightly; this comparison cannot separate reader effects from that.

## Network evidence

- APS: 98.6–98.8 % of Xcompact3d's time is in MPI, 93–95 % in `MPI_Alltoallv` (S0 and all S4 runs).
- `node_<host>.csv`: each writer machine sends and receives 320–375 GB per run; peaks reach
  120–123.5 MB/s, the 125 MB/s limit of 1 GbE (`spod/node-local/tools/net_rates.py`).
- With 160 ranks, four machines carry about 372 GB and the others 325 GB; the step takes 43 s
  instead of 35 s.
- When two full-size runs overlapped on different machines (8–9 Nov), the S4 r4x8 run slowed from
  35.4 to 53.8 s per step.
