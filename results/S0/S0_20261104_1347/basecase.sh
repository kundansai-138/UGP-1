#!/usr/bin/env bash
# =====================================================================
# basecase.sh : the same Xcompact3d case with No snapshot output
# during the time loop (ioutput > ilast)
# runs the MPI based base case implementation with aps based profiling

set -euo pipefail

CASE=${CASE:-$HOME/ugp/runs/testcase_mpi}
BIN=${BIN:-$HOME/ugp/src/x3d_lab_mpi/xcompact3d}
HOSTFILE=${HOSTFILE:-$HOME/ugp/tools/hostfile7}
PPN_LIST=${PPN_LIST:-24}
REPEATS=${REPEATS:-3}
WARMUP=${WARMUP:-1}
NX=${NX:-}
ILAST=${ILAST:-}
NOCKPT=${NOCKPT:-0}
FIX_PROW=${FIX_PROW:-1}
RUNROOT=${RUNROOT:-$HOME/ugp/runs/S0}
TOOLS=$(cd "$(dirname "$0")" && pwd)
PROFILE=$TOOLS/x3d_profile.sh
ANALYZE=$TOOLS/analyze_x3d_profile.py

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -f $CASE/input.i3d ]] || die "no input.i3d in $CASE"
[[ -x $BIN ]]            || die "binary not found: $BIN"
[[ -f $HOSTFILE ]]       || die "hostfile not found: $HOSTFILE"
[[ -f $PROFILE && -f $ANALYZE ]] || die "x3d_profile.sh / analyze_x3d_profile.py must be next to this script ($TOOLS)"

OUTPUT_DIRS="data probes statistics streaklines out"
NHOSTS=$(grep -v '^\s*#' "$HOSTFILE" | grep -cv '^\s*$')

# ---- namelist helpers (input.i3d lines like "ioutput=500   ! comment") ----
getv() { grep -ioE "(^|[ ,])$1\s*=\s*[-0-9]+" "$2" | head -1 | sed -E 's/.*=\s*//'; }
setv() {  # setv name value file   (only replaces an existing line)
  grep -qiE "(^|[ ,])$1\s*=" "$3" || die "$1 not found in $3"
  sed -i -E "s/(^|[ ,])($1\s*=\s*)[-0-9]+/\1\2$2/I" "$3"
  [[ $(getv "$1" "$3") == "$2" ]] || die "could not set $1=$2 in $3"
}

# ---- 1. fresh run folder with inputs only ----
GRIDTAG=${NX:+_n$NX}
STAMP=$(date +%Y%m%d_%H%M)
RUN=$RUNROOT/S0${GRIDTAG}_$STAMP
mkdir -p "$RUN"
tar -C "$CASE" --exclude='./data' --exclude='./probes' --exclude='./statistics' --exclude='./streaklines' \
    --exclude='./out' --exclude='./prof_*' --exclude='*.log' --exclude='*.log1' --exclude='*.bin' \
    --exclude='*.bp4' --exclude='*.bp5' --exclude='*.sst' --exclude='restart*' --exclude='checkpoint*' \
    -cf - . | tar -C "$RUN" -xf -
for d in $OUTPUT_DIRS; do mkdir -p "$RUN/$d"; done
cd "$RUN"
cp input.i3d input.i3d.orig

# ---- 2. S0 edits ----
[[ -n $NX ]]    && { setv nx "$NX" input.i3d; setv ny "$NX" input.i3d; setv nz "$NX" input.i3d; }
[[ -n $ILAST ]] && setv ilast "$ILAST" input.i3d
IL=$(getv ilast input.i3d); [[ -n $IL ]] || die "ilast not found"
setv ioutput $((IL + 1)) input.i3d
if [[ $NOCKPT == 1 ]]; then setv icheckpoint $((IL + 1)) input.i3d; fi
echo "S0 run folder: $RUN"
echo "changes to input.i3d:"; diff input.i3d.orig input.i3d || true

clean_out() { for d in $OUTPUT_DIRS; do rm -rf "${d:?}"/*; done; rm -f restart.info checkpoint* *.sst; }

# ---- 3-5. runs ----
TAGS=()
for PPN in $PPN_LIST; do
  NP=$((NHOSTS * PPN))
  echo; echo "=== PPN=$PPN  NP=$NP on $NHOSTS hosts ==="
  for ((w = 1; w <= WARMUP; w++)); do
    clean_out
    TAG=ppn${PPN}_warmup$w
    STAGE=S0 BIN=$BIN HOSTFILE=$HOSTFILE NP=$NP PPN=$PPN OUTPUT_DIRS="$OUTPUT_DIRS" TAG=$TAG \
      bash "$PROFILE" > "warmup_$TAG.out" 2>&1 || die "warm-up failed, see $RUN/warmup_$TAG.out and prof_S0_$TAG/x3d_ts.log"
    tail -8 "prof_S0_$TAG/summary.txt" | head -6
  done
  if [[ $FIX_PROW == 1 && $(getv p_row input.i3d) == 0 ]]; then
    PR=$(python3 -c "import json;print(json.load(open('prof_S0_ppn${PPN}_warmup1/summary.json'))['p_row_x_p_col'] or '')" 2>/dev/null | tr -d '[]() ' || true)
    if [[ -n $PR ]]; then
      setv p_row "${PR%%,*}" input.i3d; setv p_col "${PR##*,}" input.i3d
      echo "fixed p_row x p_col = ${PR%%,*} x ${PR##*,} in input.i3d"
      # same case for all later stages: the ORIGINAL input (same ioutput, ilast, grid) + this fixed split
      cp input.i3d.orig "input_stages_ppn$PPN.i3d"
      [[ -n $NX ]]    && { setv nx "$NX" "input_stages_ppn$PPN.i3d"; setv ny "$NX" "input_stages_ppn$PPN.i3d"; setv nz "$NX" "input_stages_ppn$PPN.i3d"; }
      [[ -n $ILAST ]] && setv ilast "$ILAST" "input_stages_ppn$PPN.i3d"
      setv p_row "${PR%%,*}" "input_stages_ppn$PPN.i3d"; setv p_col "${PR##*,}" "input_stages_ppn$PPN.i3d"
      echo "input for S1-S4 (identical case, only ioutput differs from S0): $RUN/input_stages_ppn$PPN.i3d"
    else
      echo "WARNING: could not read p_row x p_col from the warm-up log; left at auto (0 0)"
    fi
  fi
  for ((r = 1; r <= REPEATS; r++)); do
    clean_out
    TAG=ppn${PPN}_run$r
    STAGE=S0 BIN=$BIN HOSTFILE=$HOSTFILE NP=$NP PPN=$PPN OUTPUT_DIRS="$OUTPUT_DIRS" TAG=$TAG \
      bash "$PROFILE" > "log_$TAG.out" 2>&1 || echo "WARNING: run $TAG exited non-zero, see prof_S0_$TAG/"
    TAGS+=("prof_S0_$TAG")
    echo "--- $TAG"; sed -n '2,5p' "prof_S0_$TAG/summary.txt" 2>/dev/null || true
  done
  # p_row/p_col depend on NP: reset to auto for the next PPN
  if [[ $FIX_PROW == 1 ]]; then cp input.i3d "input.i3d.ppn$PPN"; setv p_row 0 input.i3d; setv p_col 0 input.i3d; fi
done

echo; echo "=== S0 comparison ==="
python3 "$ANALYZE" --compare "${TAGS[@]}" | tee S0_compare.txt
echo
echo "Results in $RUN  (input used per PPN: input.i3d.ppn<PPN>)"
