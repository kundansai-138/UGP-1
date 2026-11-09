# Running S0, S4 and the node-local cases

Set up first: [XCOMPACT3D_SETUP.md](XCOMPACT3D_SETUP.md) and [SPOD_SETUP.md](SPOD_SETUP.md).
Scripts live in `~/ugp/tools` (`scripts/ugp-tools` here) and `~/spod-dev` (`spod/node-local` here).

## 0. Before you start

- **Machines.** `~/ugp/tools/hostfile` lists csews160–167 (writer machines), `hostfile_spod`
  lists csews168–171 (SPOD machines). The scripts start ranks with MPICH's `mpirun -f <hostfile>`
  and log in to each machine with `ssh` (no password), so passwordless SSH between all machines
  must work. Check that nobody else is using them:
  `bash ~/spod-dev/tools/probe_nodes.sh 160 171` (shows CPU, memory, load and logged-in users).
- **One case at a time.** All machines share one 1 GbE network. When two full-size runs overlapped,
  both slowed down by about 50 %.
- **Disk.** Home is a shared 50 GB. Each SPOD mode file is 1.75 GB at 301³ and the old reader
  writes 5 per run. Use `case_s4.sh` (deletes them after each run) rather than `run_S4.sh`
  (keeps them).
- **Long runs.** Start them with `nohup setsid … &` so they keep running after you log out.

## 1. S0: the solver alone (baseline)

`basecase.sh` runs the case with zero snapshots under Intel APS: 192 ranks (24 per machine),
`ivisu = 0`, 2 runs. It uses the plain MPI build `~/ugp/src/x3d_lab/xcompact3d`, which is not in
this repository.

```bash
bash ~/ugp/tools/basecase.sh
```

The S0 of 4 November ran without APS's imbalance barrier, while all S4 runs used
`APS_IMBALANCE_TYPE=2`. Re-run S0 with the same APS settings before comparing step times.

## 2. S4: streaming with the old reader (`case_s4.sh`)

What one case does: fresh run folder from `~/ugp/runs/testcase_adios` (`ilast 150`, `ioutput 5`,
no restart files, one line per step) → `meta.txt` with the exact setup → for each repeat: start
the per-machine samplers, start Xcompact3d (under APS and the wait-timing library), wait for
`data.sst`, start SPOD, wait for both, write the APS reports, delete the mode files → summary →
wait analysis.

### Options (environment variables)

| Variable | Default | Meaning |
|---|---|---|
| `WHOSTS`, `RHOSTS` | `hostfile`, `hostfile_spod` | machines for Xcompact3d and SPOD |
| `WPPN`, `RPPN` | 24, 12 | ranks per machine for Xcompact3d and SPOD |
| `WPIN`, `RPIN` | none | CPU lists for `taskset -c`, e.g. `0-19` and `20-23` |
| `REPEATS` | 1 | runs per case |
| `ILAST`, `IOUTPUT` | 150, 5 | steps and snapshot interval |
| `APS`, `APS_SPOD` | 1, 1 | Intel APS on the solver and on SPOD |
| `WAIT`, `WAIT_WRITER` | 1, 1 | per-call timing of `begin_step`/`end_step` (LD_PRELOAD library) |
| `KEEP_SPOD`, `KEEP_MODES` | 0, 0 | keep the spectrum (`KEEP_SPOD=1`) and the first N modes of run 1 |
| `RUNROOT` | `~/ugp/runs/S4` | where the run folder is created |

### The cases that were run

```bash
# separate machines, 192 Xcompact3d ranks, 48 / 32 / 16 SPOD ranks
RPPN=12 bash ~/ugp/tools/case_s4.sh
RPPN=8  bash ~/ugp/tools/case_s4.sh
RPPN=4  bash ~/ugp/tools/case_s4.sh

# co-located: 20 Xcompact3d + 4 SPOD ranks on each of csews160-167
WPPN=20 RPPN=4 RHOSTS=~/ugp/tools/hostfile WPIN=0-19 RPIN=20-23 RUNROOT=~/ugp/runs/S4_coloc \
  bash ~/ugp/tools/case_s4.sh

# separate, 160 Xcompact3d ranks (20 per machine) and 32 SPOD ranks (8 per machine)
WPPN=20 RPPN=8 WPIN=0-19 RUNROOT=~/ugp/runs/S4_sep160 bash ~/ugp/tools/case_s4.sh
```

Add `KEEP_SPOD=1 KEEP_MODES=0` to keep each run's energy spectrum. The run folder is named
`S4_r<SPOD machines>x<SPOD ranks per machine>_<date>_<time>`.

## 3. S4 with the node-local reader (`run_S4_nl.sh`)

`~/spod-dev/tools/run_S4_nl.sh` is `case_s4.sh` with the node-local reader. It takes the same
options and adds three things: it gives the reader Xcompact3d's processor grid (computed with
2decomp's own rule) and stops the run if Xcompact3d's log shows a different grid; it passes the
reader's settings as environment variables (overlap as a fraction, 5 saved modes like the old
reader); and it limits each launch to `RUN_TIMEOUT` seconds (default 4 h).

```bash
# co-located, same layout as S4_coloc
WPPN=20 RPPN=4 WHOSTS=~/ugp/tools/hostfile RHOSTS=~/ugp/tools/hostfile WPIN=0-19 RPIN=20-23 \
  RUNROOT=~/spod-dev/runs/S4_coloc_nl KEEP_SPOD=1 KEEP_MODES=0 bash ~/spod-dev/tools/run_S4_nl.sh

# separate, same layout as S4_sep160
WPPN=20 RPPN=8 WPIN=0-19 WHOSTS=~/ugp/tools/hostfile RHOSTS=~/ugp/tools/hostfile_spod \
  RUNROOT=~/spod-dev/runs/S4_sep_nl KEEP_SPOD=1 KEEP_MODES=0 bash ~/spod-dev/tools/run_S4_nl.sh

# quick 65^3 check of the whole chain (about 3 minutes)
GRID=65 ILAST=15 IOUTPUT=1 WPPN=20 RPPN=4 WHOSTS=~/ugp/tools/hostfile RHOSTS=~/ugp/tools/hostfile \
  WPIN=0-19 RPIN=20-23 RUNROOT=~/spod-dev/runs/quick65_nl bash ~/spod-dev/tools/run_S4_nl.sh
```

In `run1/spod_ts.log` the reader prints which Xcompact3d ranks each machine reads and
`Node-local layout: permutation check passed`.

## 4. Everything in one queue (`queue_S4.sh`)

`~/spod-dev/tools/queue_S4.sh` runs the cases one after another, each with `REPEATS=1` and
`KEEP_SPOD=1 KEEP_MODES=0`:

| # | Case | Script |
|---|---|---|
| 0 | compare saved 301³ modes, then delete them (specific to 9 Nov; skip with `ONLY`) | – |
| 1 | `r4x4`, old reader | `case_s4.sh RPPN=4` |
| 2 | co-located, old reader | `case_s4.sh` |
| 3 | co-located, node-local | `run_S4_nl.sh` |
| 4 | separate 160, old reader | `case_s4.sh` |
| 5 | separate 160, node-local | `run_S4_nl.sh` |

```bash
# start on csews160 (keeps running after you log out)
ONLY="1 2 3 4 5" nohup setsid bash ~/spod-dev/tools/queue_S4.sh > /dev/null 2>&1 < /dev/null &

# progress: one line per case start/end; full output per case in ~/spod-dev/runs/queue_S4/
tail -f ~/spod-dev/runs/queue_S4.log
```

For 2 repeats, change `REPEATS=1` to `REPEATS=2` in `COMMON=(…)` at the top of the script
(not while it is running). The 9 November queue took about 9.5 hours.

### Stopping a queue or a case

```bash
kill $(cat ~/spod-dev/runs/queue_S4.pid)                 # the queue (no new case starts)
pkill -u $USER -f '[r]un_S4|[c]ase_s4'                   # the running case script
for h in csews{160..171}; do ssh -n $h "pkill -u \$USER -f '[x]compact3d|[s]treamingSPOD|[s]pod_nl|[h]ydra_pmi_proxy'"; done
touch <run folder>/run1/.stop                            # ends the per-machine samplers
for h in csews{160..171}; do echo "$h: $(ssh -n $h "pgrep -u \$USER -c -f '[x]compact3d|[s]treamingSPOD|[s]pod_nl|[h]ydra_pmi'")"; done   # all 0
```

The `[x]compact3d` form matches the programs but not the `pkill` command itself.

## 5. Correctness tests

**One machine, three readers on one stream** (`~/spod-dev/run_compare.sh`): Xcompact3d at 65³
(2 × 4 ranks) streams to the old reader (4 ranks) and the node-local reader (4 and 5 ranks) at the
same time, then the outputs are compared. Options: `N`, `PROW`, `PCOL`, `ILAST`, `IOUTPUT`, `NDFT`,
`OVERLAP`, `K`, `NSAVE`, `NP_OLD`, `NP_NL`. A size guard stops runs over 1 GB of output or 8 GB of
memory unless `ALLOW_BIG=1`.

```bash
bash ~/spod-dev/run_compare.sh      # expect: SAME RESULTS ... RESULT: PASS
```

Result on 8 November: energies identical, modes equal to about 1e-14 after removing the arbitrary
phase, for both rank counts.

## 6. What a run folder contains

| File | Contents |
|---|---|
| `summary.txt` / `.csv` | total time, start-up, median step, extra time per snapshot step, peak memory |
| `wait_summary.txt`, `wait_steps.csv`, `bursts.csv` | SPOD wait / read / compute per snapshot, writer hand-off, queue check |
| `meta.txt` | binaries with md5, launch commands, settings, machines |
| `aps_report_*.html` | full APS report |
| `run1/x3d_ts.log`, `run1/spod_ts.log` | time-stamped output of both programs |
| `run1/aps_report.txt`, `run1/aps_spod_report.txt`, `run1/aps_*mpi_*.txt` | APS summaries, MPI per function and per rank |
| `run1/node_<host>.csv` | CPU busy, memory, network bytes, once per second |
| `run1/wait/*.csv` | every `begin_step` / `end_step` call of every rank |
| `run1/spod/EnegySpectra.txt` | the SPOD energy spectrum (with `KEEP_SPOD=1`) |

How the timings are measured: step times from the time stamps between `Time step = N` lines;
SPOD's *wait* = time inside `adios2_begin_step`, *read* = until `adios2_end_step` returns (the data
moves here), *compute* = until the next `begin_step`.

## 7. Analysis tools

```bash
# reader/writer wait for several cases side by side
python3 ~/ugp/tools/wait_analysis.py ~/ugp/runs/S4_coloc/S4_* ~/spod-dev/runs/S4_coloc_nl/S4_*

# network MB/s per machine and second (1 GbE = 125 MB/s)
python3 ~/spod-dev/tools/net_rates.py ~/ugp/runs/S4_coloc/S4_r8x4_*

# compare two runs' spectra and modes (needs KEEP_SPOD=1, KEEP_MODES>=1)
python3 ~/spod-dev/tools/compare_spod_outputs.py RUN_A/run1/spod RUN_B/run1/spod
```

Watching a run: `tail -f <run folder>/run1/x3d_ts.log`, and
`bash ~/spod-dev/tools/rank_map.sh ~/ugp/tools/hostfile` to see which ranks run on which machine.

## 8. Things learned the hard way

- `run_S4.sh` keeps all SPOD outputs; a queue of them would have filled the 50 GB home.
- Two full-size runs at the same time slow each other by about 50 % through the shared network,
  even on different machines.
- The stop commands match programs by name on the listed machines; don't run other programs
  called `xcompact3d` or `streamingSPOD` there at the same time.
- The old reader's settings are compiled in (always 5 modes); the node-local reader reads them
  from the environment.
- HPC2013 compute nodes have no route to the CSE lab network, so SST between csews and HPC compute
  nodes needs a route from the HPC admins or a relay through an HPC login node.
