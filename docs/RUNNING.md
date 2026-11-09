# Running

Build everything first: [XCOMPACT3D_SETUP.md](XCOMPACT3D_SETUP.md) and [SPOD_SETUP.md](SPOD_SETUP.md).
Run every command from the repository folder, in a shell where the environment is loaded:

```bash
. scripts/env.sh
```

Your own results are written under `runs/` in the repository.

## 0. Before you start

**Machines.** Two host files list the machines, one `name:slots` per line:

| File | Used for |
|---|---|
| `scripts/hostfile` | Xcompact3d; also SPOD when it runs on the same machines |
| `scripts/hostfile_spod` | SPOD when it runs on separate machines |

They contain the machines the results were measured on (8 in `hostfile`, 4 in `hostfile_spod`);
replace them with yours. Every machine must see the repository at the same path and accept
password-free SSH from the others.

**Network interface.** MPI and SST must use the same network interface on every machine. The
scripts use `eno1`; if your interface has another name (see `ip -o -4 addr`), add `IFACE=<name>`
to every command below.

**Profiling.** The run scripts profile both programs with Intel APS (part of Intel oneAPI VTune),
loaded from `/opt/intel/oneapi/vtune/2025.7/env/vars.sh`. If it is installed elsewhere, set
`APS_ENV=<path>/env/vars.sh`. Without it, add `APS=0` to the commands.

**One case at a time.** All machines share the network, and Xcompact3d uses it heavily. Two large
cases at the same time slow each other down by about 50 %.

**Disk space.** At 301³ every SPOD mode file is 1.75 GB and each run writes 5 of them. The run
scripts delete them at the end of each run unless told otherwise (`KEEP_SPOD`, `KEEP_MODES`).

**Long runs.** A 301³ case takes about 1.5–2 hours. Start long runs with `nohup setsid … &` so they
keep running after you log out (see section 4).

## 1. Correctness test on one machine

`spod/node-local/run_compare.sh` runs Xcompact3d on a 65³ grid (2 × 4 ranks) on the machine you are
logged in to. It streams the snapshots to three SPOD runs at the same time: the old reader on 4 ranks,
and the node-local reader on 4 and on 5 ranks. All three receive the same snapshots, so their
results must agree.

```bash
bash spod/node-local/run_compare.sh
```

Expected end of the output: `SAME RESULTS` twice and `RESULT: PASS`. The energies are identical and
the modes agree to about 1e-14 (each mode is compared after removing its arbitrary phase). It takes
about a minute and writes about 100 MB to `runs/compare/`.

Options: `N` (grid, default 65), `PROW`, `PCOL` (Xcompact3d grid, default 2 × 4), `ILAST`, `IOUTPUT`,
`NDFT`, `OVERLAP` (fraction), `K`, `NSAVE`, `NP_OLD`, `NP_NL` (SPOD ranks, default `"4 5"`). Runs
that would write more than 1 GB or need more than 8 GB of memory stop before starting, unless you
add `ALLOW_BIG=1`.

## 2. Small test on two machines

Before running the full cases, test the multi-machine runs on a 65³ grid. Make a host file with
two of your machines (example names):

```bash
printf 'machineA:24\nmachineB:24\n' > runs/test_hosts
```

**Node-local reader** (about 1 minute): 8 Xcompact3d ranks and 2 SPOD ranks on each machine, 15
steps with a snapshot every step:

```bash
GRID=65 ILAST=15 IOUTPUT=1 APS=0 WPPN=8 RPPN=2 WPIN=0-7 RPIN=20-21 \
  WHOSTS=runs/test_hosts RHOSTS=runs/test_hosts KEEP_SPOD=1 \
  bash spod/node-local/tools/run_S4_nl.sh
```

In `run1/spod_ts.log` of the new folder under `runs/S4_nl/`, the reader lists, for each machine,
the Xcompact3d ranks it reads, followed by `Node-local layout: permutation check passed`.

**Old reader** (about 1 minute): the old reader has its grid built in, so this test uses its 65³
variant and a copy of the case with a 65³ grid:

```bash
SRC="spod/adios2-old/Streaming SPOD/streamingSPODModule.F90.65grid" bash spod/adios2-old/build.sh
cp -r case/testcase_adios runs/case65
sed -i -E 's/^(n[xyz]\s*=\s*)[0-9]+/\165/' runs/case65/input.i3d
CASE=$PWD/runs/case65 ILAST=45 IOUTPUT=5 SPOD_NDFT=4 SPOD_OVERLAP=2 APS=0 WPPN=8 RPPN=2 \
  WPIN=0-7 RPIN=20-21 WHOSTS=runs/test_hosts RHOSTS=runs/test_hosts KEEP_SPOD=1 \
  bash scripts/case_s4.sh
bash spod/adios2-old/build.sh              # back to the 301^3 reader for the full cases
```

`ILAST=45 IOUTPUT=5` gives the 10 snapshots this variant expects; `SPOD_NDFT=4 SPOD_OVERLAP=2`
tells the script its DFT settings.

Both tests should end with `run1: writer rc 0, SPOD rc 0` in their output.

## 3. Full cases (301³)

### What one case does

`scripts/case_s4.sh` (old reader) and `spod/node-local/tools/run_S4_nl.sh` (node-local reader) work
the same way:

1. Check the settings, the programs, the host files and the CPU lists.
2. Make a new run folder `runs/<RUNROOT>/S4_r<SPOD machines>x<SPOD ranks per machine>_<time>`, copy
   the case into it and set `ilast`, `ioutput`, one output line per step and no restart files.
3. Write `meta.txt`: programs with checksums, launch commands, settings, every machine's hardware
   and load.
4. For each repeat: start a sampler on every machine (CPU, memory, network once per second), start
   Xcompact3d, wait for the SST contact file `data.sst`, start SPOD, wait for both to finish, write
   the APS reports, then delete the mode files (unless `KEEP_SPOD`/`KEEP_MODES` say otherwise).
5. Summarise all repeats (`summary.txt`).
6. Analyse how long SPOD waited, read and computed per snapshot (`wait_summary.txt`).

`run_S4_nl.sh` also gives the reader Xcompact3d's processor grid and stops the run if Xcompact3d's
log shows a different grid.

### Options

Set them in front of the command, for example `RPPN=4 bash scripts/case_s4.sh`.

| Variable | Default | Meaning |
|---|---|---|
| `WHOSTS`, `RHOSTS` | `scripts/hostfile`, `scripts/hostfile_spod` | machines for Xcompact3d and for SPOD |
| `WPPN`, `RPPN` | 24, 12 | ranks per machine for Xcompact3d and for SPOD |
| `WPIN`, `RPIN` | none | CPUs each rank may use (`taskset -c`), for example `0-19` and `20-23` |
| `REPEATS` | 1 | runs of the same case |
| `ILAST`, `IOUTPUT` | 150, 5 | number of steps and steps between snapshots |
| `APS`, `APS_SPOD` | 1, 1 | Intel APS on Xcompact3d and on SPOD |
| `WAIT`, `WAIT_WRITER` | 1, 1 | time every snapshot hand-over (needed for `wait_summary.txt`) |
| `KEEP_SPOD`, `KEEP_MODES` | 0, 0 | `KEEP_SPOD=1` keeps the energy spectrum; `KEEP_MODES=N` also keeps the first N modes |
| `RUNROOT` | `runs/S4` (old), `runs/S4_nl` (node-local) | where the run folder is created |
| `IFACE` | `eno1` | network interface for MPI and SST |
| `GRID` | none | node-local script only: `nx = ny = nz = GRID` for tests |
| `RUN_TIMEOUT` | 14400 | node-local script only: seconds before each launch is stopped |

### Old reader

```bash
# Xcompact3d 24 per machine on hostfile (192 ranks); SPOD on hostfile_spod with 12, 8 or 4 per machine
RPPN=12 KEEP_SPOD=1 bash scripts/case_s4.sh
RPPN=8  KEEP_SPOD=1 bash scripts/case_s4.sh
RPPN=4  KEEP_SPOD=1 bash scripts/case_s4.sh

# same machines: 20 Xcompact3d ranks on CPUs 0-19 and 4 SPOD ranks on CPUs 20-23 of each hostfile machine
WPPN=20 RPPN=4 RHOSTS=scripts/hostfile WPIN=0-19 RPIN=20-23 RUNROOT=runs/S4_coloc KEEP_SPOD=1 \
  bash scripts/case_s4.sh

# separate machines: 20 Xcompact3d ranks per hostfile machine (160), 8 SPOD ranks per hostfile_spod machine (32)
WPPN=20 RPPN=8 WPIN=0-19 RUNROOT=runs/S4_sep160 KEEP_SPOD=1 bash scripts/case_s4.sh
```

### Node-local reader

```bash
# same machines, same layout as S4_coloc above
WPPN=20 RPPN=4 WHOSTS=scripts/hostfile RHOSTS=scripts/hostfile WPIN=0-19 RPIN=20-23 \
  RUNROOT=runs/S4_coloc_nl KEEP_SPOD=1 bash spod/node-local/tools/run_S4_nl.sh

# separate machines, same layout as S4_sep160 above
WPPN=20 RPPN=8 WPIN=0-19 WHOSTS=scripts/hostfile RHOSTS=scripts/hostfile_spod \
  RUNROOT=runs/S4_sep160_nl KEEP_SPOD=1 bash spod/node-local/tools/run_S4_nl.sh
```

The CPU lists `0-19` and `20-23` fit 24-thread machines; check yours with `lscpu -e`.

## 4. All five cases in one queue

`spod/node-local/tools/queue_S4.sh` runs these cases one after another, with `REPEATS=1` and
`KEEP_SPOD=1 KEEP_MODES=0`:

| # | Case | Script | Results |
|---|---|---|---|
| 1 | separate machines, 192 + 16 ranks, old reader | `scripts/case_s4.sh` | `runs/S4_sep192` |
| 2 | same machines, 160 + 32 ranks, old reader | `scripts/case_s4.sh` | `runs/S4_coloc` |
| 3 | same machines, 160 + 32 ranks, node-local reader | `run_S4_nl.sh` | `runs/S4_coloc_nl` |
| 4 | separate machines, 160 + 32 ranks, old reader | `scripts/case_s4.sh` | `runs/S4_sep160` |
| 5 | separate machines, 160 + 32 ranks, node-local reader | `run_S4_nl.sh` | `runs/S4_sep160_nl` |

All five take about 9.5 hours. Start the queue so that it keeps running after you log out:

```bash
nohup setsid bash spod/node-local/tools/queue_S4.sh > /dev/null 2>&1 < /dev/null &
```

Only some cases, for example 3 and 5:

```bash
ONLY="3 5" nohup setsid bash spod/node-local/tools/queue_S4.sh > /dev/null 2>&1 < /dev/null &
```

Progress, one line per case start and end:

```bash
tail -f runs/queue_S4.log
```

The full output of each case is in `runs/queue_S4/<number>_<case>.log`. For 2 repeats per case,
change `REPEATS=1` to `REPEATS=2` in the `COMMON=(…)` line of the script before starting it.

## 5. Watching and stopping a run

Follow a running case:

```bash
tail -f runs/S4_coloc/S4_*/run1/x3d_ts.log
```

See which MPI ranks run on which machine:

```bash
bash spod/node-local/tools/rank_map.sh scripts/hostfile
```

Stop everything (queue, case script, programs on all machines, samplers):

```bash
kill $(cat runs/queue_S4.pid)                     # the queue, if you started one
pkill -u $USER -f '[c]ase_s4.sh|[r]un_S4_nl.sh'   # the case script
for h in $(cut -d: -f1 scripts/hostfile scripts/hostfile_spod); do
  ssh -n $h "pkill -u \$USER -f '[x]compact3d|[s]treamingSPOD|[s]pod_nl|[h]ydra_pmi_proxy'"
done
touch <run folder>/run1/.stop                     # ends the per-machine samplers
```

The `[x]compact3d` form matches the programs but not the `pkill` command itself. Check that nothing
is left (every machine should show 0):

```bash
for h in $(cut -d: -f1 scripts/hostfile scripts/hostfile_spod); do
  echo "$h: $(ssh -n $h "pgrep -u \$USER -c -f '[x]compact3d|[s]treamingSPOD|[s]pod_nl|[h]ydra_pmi'")"
done
```

These commands stop programs by name on the listed machines. Do not run other programs called
`xcompact3d` or `streamingSPOD` there at the same time.

## 6. What a run folder contains

| File | Contents |
|---|---|
| `summary.txt`, `summary.csv` | total time, start-up time, median time per step, extra time of a snapshot step, slowest step, peak memory |
| `wait_summary.txt`, `wait_steps.csv`, `bursts.csv` | SPOD wait / read / compute per snapshot, Xcompact3d hand-over time, whether the 2-snapshot queue was enough |
| `meta.txt` | programs with checksums, launch commands, settings, machines |
| `aps_report_*.html` | full APS report, to open in a browser |
| `run1/x3d_ts.log`, `run1/spod_ts.log` | output of both programs, each line with its time |
| `run1/aps_report.txt`, `run1/aps_spod_report.txt` | APS summaries of Xcompact3d and SPOD |
| `run1/aps_mpi_functions.txt`, `run1/aps_mpi_per_rank.txt` (and `aps_spod_…`) | time in each MPI function and per rank |
| `run1/node_<machine>.csv` | CPU busy, memory used, network bytes received and sent, once per second |
| `run1/wait/*.csv` | every snapshot hand-over of every rank |
| `run1/spod/EnegySpectra.txt` | SPOD energy spectrum (with `KEEP_SPOD=1`) |

How the times are measured:

- **Time per step:** the time between consecutive `Time step = N` lines in `x3d_ts.log`, median over
  steps without a snapshot.
- **SPOD wait:** time spent inside `adios2_begin_step`, waiting for the next snapshot.
- **SPOD read:** from `adios2_begin_step` returning until `adios2_end_step` returns; the data moves here.
- **SPOD compute:** from `adios2_end_step` returning until the next `adios2_begin_step`.

## 7. Analysis tools

```bash
# wait / read / compute of several runs side by side
python3 scripts/wait_analysis.py runs/S4_coloc/S4_* runs/S4_coloc_nl/S4_*

# network MB/s per machine and second (1 Gb/s Ethernet carries at most 125 MB/s)
python3 spod/node-local/tools/net_rates.py runs/S4_coloc/S4_*

# energies and modes of two runs (needs KEEP_SPOD=1 and KEEP_MODES of 1 or more)
python3 spod/node-local/tools/compare_spod_outputs.py RUN_A/run1/spod RUN_B/run1/spod
```
