#!/usr/bin/env bash
# =====================================================================
# basecase.sh - S0 base case, end to end, in one script.
#
# Runs the Xcompact3d case with zero snapshots (ivisu = 0, ioutput = ilast + 1),
# every run under Intel APS (MPI profiling).
#
# Steps
#   1. settings and checks
#   2. fresh run folder (inputs only); input.i3d: ioutput = ilast+1, ivisu = 0 (zero snapshots),
#      icheckpoint = ilast+1 (no restart writes), ilist = 1 (one "Time step" line per step)
#   3. machine information (meta.txt)
#   4. REPEATS runs, each under APS, each with its own APS report
#   5. summary table over all runs (summary.csv / summary.txt)
#
# Your case: nozzle jet, 301^3, ilast = 150 (taken from the master input), 8 machines x 24 ranks,
#            2 runs. The master input keeps ioutput = 5 for S1-S4; only the S0 copy is changed.
#
# Usage:   bash basecase.sh
# Options: CASE=~/ugp/runs/testcase_mpi  BIN=~/ugp/src/x3d_lab/xcompact3d
#          HOSTFILE=~/ugp/tools/hostfile  PPN=24  REPEATS=2
#          NOCKPT=1 (no restart writes)  ILIST=1 (step line every step)
#          NX= ILAST=   (only for a small test grid; empty = your case unchanged)
# =====================================================================
set -euo pipefail

# ---------------------------------------------------------------- 1. settings
CASE=${CASE:-$HOME/ugp/runs/testcase_mpi}            # folder with the master input.i3d
BIN=${BIN:-$HOME/ugp/src/x3d_lab/xcompact3d}     # binary (MPI-IO build)
HOSTFILE=${HOSTFILE:-$HOME/ugp/tools/hostfile}       # all 8 machines csews160..167 (192 ranks)
PPN=${PPN:-24}                                       # MPI ranks per machine
REPEATS=${REPEATS:-2}                                # measured runs (~1 h 40 min each)
NX=${NX:-}                                           # test only: nx=ny=nz
ILAST=${ILAST:-}                                     # test only: number of steps
NOCKPT=${NOCKPT:-1}                                  # 1 = no restart file during the run (pure compute)
ILIST=${ILIST:-1}                                    # print step info every ILIST steps (1 = every step, needed for step times)
MPIRUN=${MPIRUN:-/usr/local/bin/mpirun}              # MPICH on csews
APS_ENV=${APS_ENV:-/opt/intel/oneapi/vtune/2025.7/env/vars.sh}
RUNROOT=${RUNROOT:-$HOME/ugp/runs/S0}

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -f $CASE/input.i3d ]] || die "no input.i3d in $CASE"
[[ -x $BIN ]]            || die "binary not found: $BIN"
[[ -f $HOSTFILE ]]       || die "hostfile not found: $HOSTFILE"
[[ -x $MPIRUN ]]         || die "mpirun not found: $MPIRUN"
command -v perl >/dev/null    || die "perl is needed for timestamps"
command -v python3 >/dev/null || die "python3 is needed for the summary"
[[ -f $APS_ENV ]]             || die "APS environment not found: $APS_ENV"
# shellcheck disable=SC1090
source "$APS_ENV"
command -v aps >/dev/null     || die "aps not found after sourcing $APS_ENV"

mapfile -t HOSTS < <(grep -v '^\s*#' "$HOSTFILE" | grep -v '^\s*$' | cut -d: -f1)
NP=$(( ${#HOSTS[@]} * PPN ))
OUTPUT_DIRS="data probes statistics streaklines out"

# read / change one number in input.i3d (keeps comments; works with "p_row=0, p_col=0")
getv() { grep -ioE "(^|[ ,])$1\s*=\s*[-0-9]+" "$2" | head -1 | sed -E 's/.*=\s*//'; }
setv() {
  grep -qiE "(^|[ ,])$1\s*=" "$3" || die "$1 not found in $3"
  sed -i -E "s/(^|[ ,])($1\s*=\s*)[-0-9]+/\1\2$2/I" "$3"
  [[ $(getv "$1" "$3") == "$2" ]] || die "could not set $1=$2 in $3"
}

# ---------------------------------------------------------------- 2. run folder
RUN=$RUNROOT/S0${NX:+_n$NX}_$(date +%Y%m%d_%H%M)
mkdir -p "$RUN"
# copy the case folder WITHOUT old outputs; the original folder is not changed
tar -C "$CASE" --exclude='./data' --exclude='./probes' --exclude='./statistics' \
    --exclude='./streaklines' --exclude='./out' --exclude='./prof_*' --exclude='*.log' \
    --exclude='*.log1' --exclude='*.bin' --exclude='*.bp4' --exclude='*.bp5' --exclude='*.sst' \
    --exclude='restart*' --exclude='checkpoint*' -cf - . | tar -C "$RUN" -xf -
cd "$RUN"
for d in $OUTPUT_DIRS; do mkdir -p "$d"; done
cp input.i3d input.i3d.orig

[[ -n $NX ]]    && { setv nx "$NX" input.i3d; setv ny "$NX" input.i3d; setv nz "$NX" input.i3d; }
[[ -n $ILAST ]] && setv ilast "$ILAST" input.i3d
IL=$(getv ilast input.i3d); [[ -n $IL ]] || die "ilast not found in input.i3d"
setv ioutput $((IL + 1)) input.i3d                   # no snapshots in the time loop
# snapshot 0 (written at start-up) is skipped only if the visu switch is off
# (needs the patched case.f90: snapshot 0 in preprocessing is guarded by ivisu)
if grep -qiE '(^|[ ,])ivisu\s*=' input.i3d; then
  setv ivisu 0 input.i3d                             # zero snapshots, also at start-up
elif grep -qiE '^\s*&BasicParam' input.i3d; then
  # ivisu belongs to the &BasicParam namelist; add it right after the group header
  sed -i -E '0,/^\s*&BasicParam/I{s/^(\s*&BasicParam.*)$/\1\nivisu=0   ! added by basecase.sh (zero snapshots)/I}' input.i3d
  [[ $(getv ivisu input.i3d) == 0 ]] || die "could not add ivisu=0 to &BasicParam"
else
  die "no &BasicParam group found in input.i3d to put ivisu=0 into"
fi
# no restart file inside the run (it is I/O too: ~40 s per write at 301^3 over NFS)
[[ $NOCKPT == 1 ]] && setv icheckpoint $((IL + 1)) input.i3d
# Xcompact3d prints "Time step = n/N" only every ilist steps; 1 gives a timestamp per step
[[ -n $ILIST ]] && setv ilist "$ILIST" input.i3d
echo "run folder : $RUN"
echo "launch     : $MPIRUN -f $HOSTFILE -np $NP -ppn $PPN aps -r <run>/aps $BIN input.i3d"
echo "changes to input.i3d:"; diff input.i3d.orig input.i3d || true

# ---------------------------------------------------------------- 3. machine information
{
  echo "date    : $(date -Is)"
  echo "binary  : $BIN (md5 $(md5sum "$BIN" | cut -d' ' -f1))"
  echo "launch  : $MPIRUN -f $HOSTFILE -np $NP -ppn $PPN"
  echo "case fs : $(findmnt -n -o FSTYPE,SOURCE -T "$RUN" 2>/dev/null || echo unknown)"
  echo "--- libraries"; ldd "$BIN" | grep -iE 'mpi|adios|fftw|not found' || true
  echo "--- input (S0)"
  grep -iE '(^|[ ,])(nx|ny|nz|p_row|p_col|ilast|ilist|ioutput|ivisu|icheckpoint|iibm|nraf)\s*=' input.i3d || true
  echo "--- machines (name, cpu, logical cpus, RAM GB, network card speed Mb/s, load average)"
  for h in "${HOSTS[@]}"; do
    ssh -o BatchMode=yes "$h" 'IF=$(ip -o route show default | awk "{print \$5; exit}");
      echo "$(hostname) | $(lscpu | sed -n "s/^Model name:\s*//p") | $(nproc) |" \
           "$(awk "/MemTotal/{printf \"%.1f\", \$2/1048576}" /proc/meminfo) | $IF $(cat /sys/class/net/$IF/speed 2>/dev/null) |" \
           "$(cut -d" " -f1-3 /proc/loadavg)"' 2>/dev/null || echo "$h | ssh failed"
  done
} > meta.txt
cat meta.txt

# ---------------------------------------------------------------- helpers for one run
# memory + CPU + network sampler, 1 line per second on one machine, stops when $2 exists
SAMPLER='
STOP=$1
IF=$(ip -o route show default 2>/dev/null | awk "{print \$5; exit}")
echo "t,cpu_busy_pct,mem_used_gb,net_rx_bytes,net_tx_bytes"
read -r _ a b c d e f g h _ < /proc/stat; pt=$((a+b+c+d+e+f+g+h)); pi=$((d+e))
end=$(( $(date +%s) + 86400 ))
while [ ! -e "$STOP" ] && [ "$(date +%s)" -lt "$end" ]; do
  sleep 1
  read -r _ a b c d e f g h _ < /proc/stat; tt=$((a+b+c+d+e+f+g+h)); ti=$((d+e))
  busy=$(awk -v t=$((tt-pt)) -v i=$((ti-pi)) "BEGIN{ if (t>0) printf \"%.1f\", 100*(t-i)/t; else print 0 }")
  pt=$tt; pi=$ti
  mem=$(awk "/MemTotal/{t=\$2} /MemAvailable/{a=\$2} END{printf \"%.3f\", (t-a)/1048576}" /proc/meminfo)
  net=$(awk -v IF="$IF" "NR>2{sub(/:/,\" \"); if (\$1!=\"lo\" && (IF==\"\" || \$1==IF)) {r+=\$2; x+=\$10}} END{print r+0\",\"x+0}" /proc/net/dev)
  echo "$(date +%s.%N | cut -c1-14),$busy,$mem,$net"
done'

clean_out() { for d in $OUTPUT_DIRS; do rm -rf "${d:?}"/*; done; rm -f restart.info checkpoint*; }

# run_one <name>  -> folder <name>/ with x3d_ts.log, node_<host>.csv, APS result + report
run_one() {
  local name=$1 dir=$RUN/$1 h pids=() rc
  mkdir -p "$dir"; clean_out
  for h in "${HOSTS[@]}"; do              # start one sampler per machine
    ssh -o BatchMode=yes "$h" "bash -s -- '$dir/.stop'" <<< "$SAMPLER" > "$dir/node_$h.csv" 2>/dev/null &
    pids+=($!)
  done
  # aps sits between mpirun and the program: every rank is profiled
  local cmd=("$MPIRUN" -f "$HOSTFILE" -np "$NP" -ppn "$PPN" aps -r "$dir/aps" "$BIN" input.i3d)
  echo "$(date +%s.%N) # START $name" > "$dir/x3d_ts.log"
  set +e
  # every output line gets a wall-clock timestamp; gfortran must not buffer stdout
  GFORTRAN_UNBUFFERED_PRECONNECTED=y "${cmd[@]}" 2>&1 \
    | perl -MTime::HiRes=time -ne '$|=1; printf "%.3f %s", time, $_' >> "$dir/x3d_ts.log"
  rc=${PIPESTATUS[0]}
  set -e
  echo "$(date +%s.%N) # END rc=$rc" >> "$dir/x3d_ts.log"
  touch "$dir/.stop"; wait "${pids[@]}" 2>/dev/null || true
  # APS report (summary: MPI time share, top MPI calls); syntax differs between versions
  aps --report "$dir/aps" > "$dir/aps_report.txt" 2>&1 \
    || aps-report "$dir/aps" > "$dir/aps_report.txt" 2>&1 \
    || echo "WARNING: aps report failed for $name - try: aps --report $dir/aps"
  echo "$name: exit code $rc, $(awk '/# START/{s=$1} /# END/{e=$1} END{printf "%.1f", e-s}' "$dir/x3d_ts.log") s"
  return "$rc"
}

# ---------------------------------------------------------------- 4. runs (all under APS)
RUNS=()
for ((r = 1; r <= REPEATS; r++)); do
  echo; echo "=== run $r of $REPEATS (APS) ==="
  run_one "run$r" || echo "WARNING: run$r failed - see $RUN/run$r/x3d_ts.log"
  RUNS+=("run$r")
done

# ---------------------------------------------------------------- 5. summary
echo; echo "=== summary ==="
python3 - "$RUN" "${RUNS[@]}" <<'PY'
import csv, glob, os, re, statistics as st, sys
run_dir, names = sys.argv[1], sys.argv[2:]
rows = []
for n in names:
    d = os.path.join(run_dir, n)
    log = os.path.join(d, "x3d_ts.log")
    if not os.path.exists(log):
        continue
    start = end = rc = None; steps = []; snaps = 0; rst = []; rst_t0 = None
    for line in open(log, errors="replace"):
        p = line.split(" ", 1)
        try: t = float(p[0])
        except ValueError: continue
        txt = p[1] if len(p) > 1 else ""
        if "# START" in txt: start = t
        if "# END" in txt:
            end = t; m = re.search(r"rc=(\d+)", txt); rc = int(m.group(1)) if m else None
        m = re.search(r"Time step\s*=\s*(\d+)\s*/", txt)
        if m: steps.append((int(m.group(1)), t))
        if "Writing snapshots" in txt: snaps += 1
        if "Writing restart point" in txt: rst_t0 = t
        if "saved successfully" in txt and rst_t0 is not None:
            rst.append((rst_t0, t)); rst_t0 = None
    # time per step: gap between consecutive printed steps / number of steps between them
    # (the last printed step has no next line, so it is not used; restart-write time inside a gap is subtracted)
    dur = []
    for (n1, t1), (n2, t2) in zip(steps, steps[1:]):
        if n2 > n1:
            io = sum(b - a for a, b in rst if t1 <= a <= t2)
            dur.append((t2 - t1 - io) / (n2 - n1))
    last_t = steps[-1][1] if steps else None
    mem, cpu = {}, []
    for f in glob.glob(os.path.join(d, "node_*.csv")):
        h = os.path.basename(f)[5:-4]
        rr = [r for r in csv.DictReader(open(f)) if r.get("mem_used_gb")]
        if rr:
            mem[h] = max(float(r["mem_used_gb"]) for r in rr)
            cpu += [float(r["cpu_busy_pct"]) for r in rr if steps and steps[0][1] <= float(r["t"]) <= steps[-1][1]]
    first_t = steps[0][1] if steps else None
    rows.append({
        "run": n, "exit_code": rc,
        "wall_s": round(end - start, 2) if start and end else "",
        "init_s": round(first_t - start, 2) if first_t and start else "",
        "loop_s": round(last_t - first_t, 2) if steps else "",
        "end_s": round(end - last_t, 2) if steps and end else "",
        "steps": steps[-1][0] if steps else 0, "step_lines": len(steps), "snapshots": snaps,
        "restarts": len(rst), "restart_s": round(sum(b - a for a, b in rst), 2),
        "step_median_s": round(st.median(dur), 3) if dur else "",
        "step_min_s": round(min(dur), 3) if dur else "",
        "step_max_s": round(max(dur), 3) if dur else "",
        "cpu_busy_mean_pct": round(sum(cpu) / len(cpu), 1) if cpu else "",
        "mem_peak_max_gb": round(max(mem.values()), 2) if mem else "",
        "mem_peak_host": max(mem, key=mem.get) if mem else "",
    })
keys = list(rows[0].keys()) if rows else []
with open(os.path.join(run_dir, "summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, keys); w.writeheader(); w.writerows(rows)
lines = ["  ".join(f"{k:>16s}" for k in keys)] + ["  ".join(f"{str(r[k]):>16s}" for k in keys) for r in rows]
meas = [r for r in rows if r["run"].startswith("run") and r["step_median_s"] != ""]
if meas:
    lines += ["", f"S0 result (median of {len(meas)} runs): "
              f"wall {st.median([r['wall_s'] for r in meas]):.1f} s, "
              f"step {st.median([r['step_median_s'] for r in meas]):.4f} s, "
              f"peak memory {max([r['mem_peak_max_gb'] for r in meas if r['mem_peak_max_gb'] != ''] or ['?'])} GB"]
if any(r["step_lines"] and r["step_lines"] < r["steps"] for r in rows):
    lines += ["NOTE: fewer step lines than steps (ilist > 1) - step times are averages over the gaps"]
if any(r["snapshots"] for r in rows):
    lines += ["WARNING: snapshots were written (see the 'snapshots' column) - this is not a zero-snapshot run"]
txt = "\n".join(lines); print(txt)
open(os.path.join(run_dir, "summary.txt"), "w").write(txt + "\n")
PY
echo
echo "APS reports: $RUN/run*/aps_report.txt"
echo "All results: $RUN"
