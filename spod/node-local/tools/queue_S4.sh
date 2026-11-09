#!/bin/bash
# =====================================================================
# queue_S4.sh - the five S4 cases one after another, 1 run each.
#
#   1. sep192_old  old reader   Xcompact3d on hostfile 24/machine (192), SPOD on hostfile_spod 4/machine (16)
#   2. coloc_old   old reader   Xcompact3d 20/machine on CPUs 0-19 + SPOD 4/machine on CPUs 20-23, hostfile
#   3. coloc_nl    node-local   same as 2
#   4. sep160_old  old reader   Xcompact3d on hostfile 20/machine (160), SPOD on hostfile_spod 8/machine (32)
#   5. sep160_nl   node-local   same as 4
#
# Every case: scripts/case_s4.sh measurements (APS both sides, reader/writer wait timing, per-machine
# CSVs, summary), REPEATS=1, KEEP_SPOD=1 KEEP_MODES=0 -> mode files deleted after each run, energy
# spectrum kept. Cases 3 and 5 use spod/node-local/tools/run_S4_nl.sh. A failed case is logged and the
# queue goes on. Results: runs/S4_sep192, runs/S4_coloc, runs/S4_coloc_nl, runs/S4_sep160, runs/S4_sep160_nl.
#
# Start (keeps running after you log out):
#     nohup setsid bash spod/node-local/tools/queue_S4.sh > /dev/null 2>&1 < /dev/null &
# Subset:   ONLY="3 5" nohup setsid bash spod/node-local/tools/queue_S4.sh ...
# Progress: tail -f runs/queue_S4.log      (full output of each case in runs/queue_S4/)
# Stop:     kill $(cat runs/queue_S4.pid), then stop the running case (docs/RUNNING.md)
# =====================================================================
set -u
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)   # this repository
[[ -n ${UGP_DEPS:-} ]] || . "$ROOT/scripts/env.sh"
R=$ROOT/runs; Q=$R/queue_S4; mkdir -p "$Q"
LOG=$R/queue_S4.log
ONLY=${ONLY:-1 2 3 4 5}
COMMON=(REPEATS=1 KEEP_SPOD=1 KEEP_MODES=0)
H=$ROOT/scripts/hostfile; HS=$ROOT/scripts/hostfile_spod
OLD=$ROOT/scripts/case_s4.sh; NL=$ROOT/spod/node-local/tools/run_S4_nl.sh

if [[ -f $R/queue_S4.pid ]] && kill -0 "$(cat "$R/queue_S4.pid")" 2>/dev/null; then
  echo "queue already running (pid $(cat "$R/queue_S4.pid"))"; exit 1
fi
echo $$ > "$R/queue_S4.pid"
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
want() { [[ " $ONLY " == *" $1 "* ]]; }

case_run() {  # case_run <n> <name> <script> VAR=value ...
  local n=$1 name=$2 script=$3 rc t0; shift 3
  want "$n" || return 0
  local clog=$Q/${n}_${name}.log
  say "=== $n $name: start  $*"
  t0=$(date +%s)
  env "${COMMON[@]}" "$@" bash "$script" > "$clog" 2>&1; rc=$?
  say "=== $n $name: end rc=$rc after $(( ($(date +%s) - t0) / 60 )) min; $(grep -E '^All results' "$clog" | tail -1)"
  grep -E '^run[0-9]+: (writer rc|Xcompact3d grid|removed)|WARNING|ERROR' "$clog" | sed 's/^/    /' | tee -a "$LOG"
}

case_run 1 sep192_old "$OLD" RPPN=4 RUNROOT="$R/S4_sep192"
case_run 2 coloc_old  "$OLD" WPPN=20 RPPN=4 RHOSTS="$H" WPIN=0-19 RPIN=20-23 RUNROOT="$R/S4_coloc"
case_run 3 coloc_nl   "$NL"  WPPN=20 RPPN=4 WHOSTS="$H" RHOSTS="$H" WPIN=0-19 RPIN=20-23 RUNROOT="$R/S4_coloc_nl"
case_run 4 sep160_old "$OLD" WPPN=20 RPPN=8 WPIN=0-19 RUNROOT="$R/S4_sep160"
case_run 5 sep160_nl  "$NL"  WPPN=20 RPPN=8 WPIN=0-19 WHOSTS="$H" RHOSTS="$HS" RUNROOT="$R/S4_sep160_nl"

say "=== queue finished"
rm -f "$R/queue_S4.pid"
