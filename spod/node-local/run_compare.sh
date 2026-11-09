#!/bin/bash
# =====================================================================
# run_compare.sh - correctness check of the node-local SPOD reader, SST only.
#
# One Xcompact3d run streams its snapshots over SST to several SPOD readers at
# the same time (RendezvousReaderCount = number of readers). Every reader gets
# exactly the same snapshots, so the SPOD outputs must agree:
#   old     build/spod_old/spod_old on NP_OLD ranks (the reader the S4 job uses)
#   nl<N>   build/spod_nl/spod_nl   on N ranks      (node-local reader)
# tools/compare_spod_outputs.py then compares old against every nl run.
#
# Everything runs on this machine only (mpirun -launcher fork). No BP4/BP5 data
# files: snapshots only travel over SST. Disk use = logs + SPOD outputs.
#
# Defaults = the jet case (case/testcase_adios, ilast 150, ioutput 5 -> 31 snapshots,
# SPOD nDFT 8, overlap 4 of 8, k 3, 2 saved frequencies, ux uy uz pp) on a smaller
# grid, so writer and all readers fit on one machine.
# Size guard: a run whose estimated output (SPOD modes) is over 1 GB, or whose estimated
# memory is over 8 GB, stops before starting unless ALLOW_BIG=1. 65^3 writes ~100 MB
# and needs ~3.5 GB. Results go to runs/compare/ in this repository.
# Needs: xcompact3d/xcompact3d built, spod/node-local/build.sh run, scripts/env.sh loaded.
# Usage:   bash spod/node-local/run_compare.sh                 (65^3)
# Options: N=65 PROW=2 PCOL=4 ILAST=150 IOUTPUT=5   Xcompact3d: N^3 grid, PROW x PCOL ranks
#          NDFT=8 OVERLAP=0.5 K=3 NSAVE=2            SPOD
#          NP_OLD=4 NP_NL="4 5"                       SPOD ranks per reader
#          TLIMIT=7200                                s, everything is stopped after this
#          IFACE=eno1                                 network interface for MPI and SST
#          ALLOW_BIG=1                                skip the size guard
# =====================================================================
set -euo pipefail
D=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)       # spod/node-local
ROOT=$(cd "$D/../.." && pwd)                           # this repository
[[ -n ${UGP_DEPS:-} ]] || . "$ROOT/scripts/env.sh"
N=${N:-65}; PROW=${PROW:-2}; PCOL=${PCOL:-4}; ILAST=${ILAST:-150}; IOUTPUT=${IOUTPUT:-5}
NDFT=${NDFT:-8}; OVERLAP=${OVERLAP:-0.5}; K=${K:-3}; NSAVE=${NSAVE:-2}
NP_OLD=${NP_OLD:-4}; NP_NL=${NP_NL:-4 5}
TLIMIT=${TLIMIT:-7200}
MPIRUN=${MPIRUN:-$(command -v mpirun)}
IFACE=${IFACE:-eno1}
X3D=$ROOT/xcompact3d/xcompact3d; OLD=$D/build/spod_old/spod_old; NL=$D/build/spod_nl/spod_nl
NET=(-launcher fork -iface "$IFACE" -genv FI_PROVIDER tcp -genv FI_TCP_IFACE "$IFACE" -genv UCX_LOG_LEVEL error)

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -n $MPIRUN && -x $MPIRUN ]] || die "mpirun not found (load scripts/env.sh or set MPIRUN)"
for b in "$X3D" "$OLD" "$NL"; do [[ -x $b ]] || die "missing $b (build xcompact3d, then run spod/node-local/build.sh)"; done

getv() { grep -ioE "(^|[ ,])$1\s*=\s*[-0-9]+" "$2" | head -1 | sed -E 's/.*=\s*//'; }
setv() {
  grep -qiE "(^|[ ,])$1\s*=" "$3" || die "$1 not found in $3"
  sed -i -E "s/(^|[ ,])($1\s*=\s*)[-0-9]+/\1\2$2/I" "$3"
  [[ $(getv "$1" "$3") == "$2" ]] || die "could not set $1=$2 in $3"
}

# readers: name:binary:ranks
READERS=("old:$OLD:$NP_OLD")
for n in $NP_NL; do READERS+=("nl$n:$NL:$n"); done
NR=${#READERS[@]}
NW=$((PROW * PCOL))
NTOT=$NW; for r in "${READERS[@]}"; do NTOT=$((NTOT + ${r##*:})); done
(( NTOT <= $(nproc) )) || echo "NOTE: $NTOT ranks on $(nproc) logical CPUs (oversubscribed, only slower)"
NSNAP=$(( ILAST / IOUTPUT + 1 ))

# size guard: modes = readers x NSAVE x 4 N^3 complex; memory measured at 65^3:
# ~1.1 GB per reader + 0.75 GB writer, scaled with N^3
OUT_MB=$(( NR * NSAVE * 4 * N * N * N * 16 / 1048576 ))
MEM_MB=$(( (NR * 1100 + 750) * N * N * N / (65 * 65 * 65) ))
echo "estimate: SPOD output ${OUT_MB} MB, memory ${MEM_MB} MB"
if [[ ${ALLOW_BIG:-0} != 1 ]] && (( OUT_MB > 1024 || MEM_MB > 8192 )); then
  die "too big while the S4 job shares this home (limit 1 GB output, 8 GB memory); set ALLOW_BIG=1 to override"
fi

# ---------------------------------------------------------------- run folder
RUN=$ROOT/runs/compare/cmp_${N}_$(date +%Y%m%d_%H%M%S)
mkdir -p "$RUN"
cp -r "$ROOT/case/testcase_adios/." "$RUN/"
cp "$RUN/adios2_config.sst.xml" "$RUN/adios2_config.xml"
sed -i -E "s/(key=\"NetworkInterface\" value=\")[^\"]*/\1$IFACE/" "$RUN/adios2_config.xml"
cd "$RUN"
for d in data probes statistics streaklines out; do mkdir -p "$d"; done
cp input.i3d input.i3d.orig
for v in nx ny nz; do setv $v "$N" input.i3d; done
setv p_row "$PROW" input.i3d; setv p_col "$PCOL" input.i3d
setv ilast "$ILAST" input.i3d; setv ioutput "$IOUTPUT" input.i3d
setv icheckpoint $((ILAST + 1)) input.i3d     # no restart writes
setv ilist 10 input.i3d
grep -qiE '(^|[ ,])ivisu\s*=' input.i3d && setv ivisu 1 input.i3d
# SST writer waits for all readers; every reader gets every step
sed -i -E 's/(key="RendezvousReaderCount" value=")[0-9]+/\1'"$NR"'/' adios2_config.xml
grep -q "key=\"RendezvousReaderCount\" value=\"$NR\"" adios2_config.xml || die "could not set RendezvousReaderCount"
grep -A1 '<io name="solution-io">' adios2_config.xml | grep -q 'type="SST"' || die "solution-io is not SST"

# SPOD settings (same for every reader). No SPOD config file: the readers use
# their hard-coded SST engine, as in the S4 job.
export SPOD_NX=$N SPOD_NY=$N SPOD_NZ=$N SPOD_STREAM=$RUN/data
export SPOD_NSNAP=$NSNAP SPOD_NDFT=$NDFT SPOD_K=$K SPOD_OVERLAP=$OVERLAP SPOD_NSAVE=$NSAVE
export SPOD_WEIGHTS=none X3D_PROW=$PROW X3D_PCOL=$PCOL
unset SPOD_CONFIG

{
  echo "date    : $(date -Is) on $(hostname)"
  echo "x3d     : $X3D, ${N}^3, p_row $PROW p_col $PCOL ($NW ranks), ilast $ILAST ioutput $IOUTPUT -> $NSNAP snapshots"
  echo "readers : ${READERS[*]}"
  echo "SPOD    : nDFT $NDFT overlap $OVERLAP k $K nSave $NSAVE"
  echo "md5     : x3d $(md5sum "$X3D" | cut -c1-12) old $(md5sum "$OLD" | cut -c1-12) nl $(md5sum "$NL" | cut -c1-12)"
  echo "changes to input.i3d:"; diff input.i3d.orig input.i3d || true
} | tee meta.txt

# ---------------------------------------------------------------- memory sampler
# every 2 s: machine 'used' MB and resident MB per program (x3d = writer, else the reader name)
sample_mem() {
  set +e
  local run_base; run_base=$(basename "$RUN")
  while [[ ! -e $RUN/.stop ]]; do
    {
      printf '%s used=%s' "$(date +%s)" "$(free -m | awk '/^Mem:/{print $3}')"
      for p in $(pgrep -u "$USER" -x 'xcompact3d|spod_old|spod_nl'); do
        echo "$(basename "$(readlink /proc/$p/cwd 2>/dev/null)") $(awk '/VmRSS/{print $2}' /proc/$p/status 2>/dev/null)"
      done | awk -v rb="$run_base" '{k = ($1 == rb) ? "x3d" : $1; s[k] += $2} END {for (k in s) printf " %s=%d", k, s[k] / 1024}'
      echo
    } >> "$RUN/mem.log"
    sleep 2
  done
}
sample_mem &
spid=$!

# ---------------------------------------------------------------- writer
echo; echo "== Xcompact3d (writer, $NW ranks)"
( env -u DISPLAY -u XAUTHORITY timeout "$TLIMIT" "$MPIRUN" "${NET[@]}" -np "$NW" "$X3D" input.i3d > x3d.log 2>&1
  echo $? > x3d.rc ) &
wpid=$!
for ((i = 0; i < 300; i += 2)); do
  [[ -s data.sst ]] && break; kill -0 $wpid 2>/dev/null || break; sleep 2
done
if [[ ! -s data.sst ]]; then
  tail -20 x3d.log; pkill -u "$USER" -f "$X3D" || true; die "writer made no data.sst"
fi
echo "data.sst ready after ${i}s"

# ---------------------------------------------------------------- readers
pids=()
for r in "${READERS[@]}"; do
  IFS=: read -r name bin np <<< "$r"
  mkdir -p "$RUN/$name"
  echo "== $name ($np ranks): $bin"
  ( cd "$RUN/$name" && env -u DISPLAY -u XAUTHORITY timeout "$TLIMIT" "$MPIRUN" "${NET[@]}" -np "$np" "$bin" \
      > spod.log 2>&1; echo $? > rc ) &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || true; done
for r in "${READERS[@]}"; do
  name=${r%%:*}; rc=$(cat "$RUN/$name/rc" 2>/dev/null || echo "?")
  [[ $rc == 0 ]] || { fail=1; echo "$name: SPOD rc $rc"; tail -15 "$RUN/$name/spod.log"; }
done
if (( fail )); then pkill -u "$USER" -f "$X3D" || true; fi
wait $wpid || true
touch "$RUN/.stop"; wait $spid || true
wrc=$(cat x3d.rc 2>/dev/null || echo "?")
echo "Xcompact3d rc $wrc, $(grep -c 'Writing snapshots' x3d.log || true) snapshot writes," \
     "$(grep -E 'Averaged time per step' x3d.log | awk '{print $NF}') s/step"
[[ $wrc == 0 ]] || tail -15 x3d.log
echo "peak memory (MB):" $(awk '{for (i = 2; i <= NF; i++) {split($i, a, "="); if (NR == 1 && a[1] == "used") b = a[2];
  if (a[2] + 0 > m[a[1]]) m[a[1]] = a[2] + 0}} END {for (k in m) printf "%s=%d ", k, m[k]; printf "(used before start %d)", b}' "$RUN/mem.log")

# ---------------------------------------------------------------- compare
echo; echo "== compare (old vs node-local)"
status=0
for r in "${READERS[@]:1}"; do
  name=${r%%:*}
  echo "--- old vs $name"
  python3 "$D/tools/compare_spod_outputs.py" "$RUN/old" "$RUN/$name" | tee "compare_old_vs_$name.txt" || status=1
done
leftover=$(find "$RUN" -name '*.bp4' -o -name '*.bp5' | head -3)
[[ -z $leftover ]] || echo "NOTE: ADIOS2 files found: $leftover"
echo; echo "run folder: $RUN ($(du -sh "$RUN" | cut -f1))"
(( fail == 0 && status == 0 )) && [[ $wrc == 0 ]] && echo "RESULT: PASS" || { echo "RESULT: FAIL"; exit 1; }
