#!/bin/bash
# =====================================================================
# queue_S4.sh - the S4 cases one after another, 1 run each (replaces the stopped job of 2026-11-08).
#
#   0. full-size check: 301^3 SPOD modes 1-2 of the node-local reader (spod-dev run of 2026-11-08)
#      vs the old reader (S4_r4x8 run2), then delete those mode files (~11.7 GB, spectra kept)
#      and the cut-off r4x4 folder S4_r4x4_20261109_0123
#   1. r4x4        old reader   X3D csews160-167 24/node (192), SPOD csews168-171 4/node (16)   ~/ugp/runs/S4
#   2. coloc_old   old reader   X3D 20/node CPU 0-19 + SPOD 4/node CPU 20-23, csews160-167    ~/ugp/runs/S4_coloc
#   3. coloc_nl    node-local   same as 2                                                   ~/spod-dev/runs/S4_coloc_nl
#   4. sep_old     old reader   X3D csews160-167 20/node CPU 0-19 (160), SPOD csews168-171 8/node (32)
#                                                                                           ~/ugp/runs/S4_sep160
#   5. sep_nl      node-local   same as 4                                                   ~/spod-dev/runs/S4_sep_nl
#
# Every case: case_s4.sh measurements (APS both sides, reader/writer wait timing, node CSVs, summary),
# REPEATS=1, KEEP_SPOD=1 KEEP_MODES=0 -> mode files deleted after each run, EnegySpectra.txt kept.
# Cases 3 and 5 use run_S4_nl.sh = case_s4.sh with the node-local reader (see its header).
# A failed case is logged and the queue goes on.
#
# Start on csews160:  nohup setsid bash ~/spod-dev/tools/queue_S4.sh > /dev/null 2>&1 < /dev/null &
# Subset:             ONLY="3 5" nohup setsid bash ~/spod-dev/tools/queue_S4.sh ...
# Progress:           tail ~/spod-dev/runs/queue_S4.log      (one full log per case in ~/spod-dev/runs/queue_S4/)
# Stop:               kill $(cat ~/spod-dev/runs/queue_S4.pid), then stop the running case
# =====================================================================
set -u
U=$HOME/ugp; D=$HOME/spod-dev
Q=$D/runs/queue_S4; mkdir -p "$Q"
LOG=$D/runs/queue_S4.log
ONLY=${ONLY:-0 1 2 3 4 5}
COMMON=(REPEATS=1 KEEP_SPOD=1 KEEP_MODES=0)

[[ ":${LD_LIBRARY_PATH:-}:" == *"/ugp/install/lib:"* ]] || . "$U/env.sh"     # same environment as the old job
if [[ -f $D/runs/queue_S4.pid ]] && kill -0 "$(cat "$D/runs/queue_S4.pid")" 2>/dev/null; then
  echo "queue already running (pid $(cat "$D/runs/queue_S4.pid"))"; exit 1
fi
echo $$ > "$D/runs/queue_S4.pid"
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
want() { [[ " $ONLY " == *" $1 "* ]]; }

# ---------------------------------------------------------------- 0. full-size check, then clean up
if want 0; then
  A=$D/runs/S4_coloc_nl/S4_r8x4_20261108_2253/run1/spod        # node-local, 10 x 16, modes 1-2
  B=$U/runs/S4/S4_r4x8_20261108_2120/run2/spod                 # old reader, 12 x 16, modes 1-5
  OUT=$D/runs/compare_301_nl_vs_r4x8run2.txt
  say "=== 0 compare 301^3: node-local $A vs old reader $B"
  if ls "$A"/SPODMODE_freq_1_* "$B"/SPODMODE_freq_1_* > /dev/null 2>&1; then
    python3 "$D/tools/compare_spod_outputs.py" "$A" "$B" > "$OUT" 2>&1
    { echo "(different Xcompact3d grids 10x16 vs 12x16: expect rounding-level, not bitwise, agreement;"
      echo " modes 3-5 exist only in the old run)"; } >> "$OUT"
    sed 's/^/    /' "$OUT" | tee -a "$LOG"
  fi
  if grep -q '^mode 2:' "$OUT" 2>/dev/null; then
    rm -f "$A"/SPODMODE_freq_* "$B"/SPODMODE_freq_*
    say "deleted the 301^3 mode files in both folders (EnegySpectra.txt kept)"
  else
    say "WARNING: comparison did not finish - mode files NOT deleted"
  fi
  [[ -d $U/runs/S4/S4_r4x4_20261109_0123 ]] && rm -rf "$U/runs/S4/S4_r4x4_20261109_0123" \
    && say "deleted the cut-off r4x4 folder S4_r4x4_20261109_0123"
fi

# ---------------------------------------------------------------- 1-5. the cases
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

case_run 1 r4x4      "$U/tools/case_s4.sh" RPPN=4
case_run 2 coloc_old "$U/tools/case_s4.sh" WPPN=20 RPPN=4 RHOSTS="$U/tools/hostfile" WPIN=0-19 RPIN=20-23 \
                     RUNROOT="$U/runs/S4_coloc"
case_run 3 coloc_nl  "$D/tools/run_S4_nl.sh" WPPN=20 RPPN=4 WHOSTS="$U/tools/hostfile" RHOSTS="$U/tools/hostfile" \
                     WPIN=0-19 RPIN=20-23 RUNROOT="$D/runs/S4_coloc_nl"
case_run 4 sep_old   "$U/tools/case_s4.sh" WPPN=20 RPPN=8 WPIN=0-19 RUNROOT="$U/runs/S4_sep160"
case_run 5 sep_nl    "$D/tools/run_S4_nl.sh" WPPN=20 RPPN=8 WPIN=0-19 WHOSTS="$U/tools/hostfile" \
                     RHOSTS="$U/tools/hostfile_spod" RUNROOT="$D/runs/S4_sep_nl"

say "=== queue finished"
rm -f "$D/runs/queue_S4.pid"
